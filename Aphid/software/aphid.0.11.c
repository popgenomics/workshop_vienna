#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXLLINE 1000  //max line length in option and taxon files
#define MAXLNAME 500 //max taxon name length 
#define MAXNSP 500 //max nb of species in input trees
#define MAX_PARAMETER_NAME_LG 50

#define QuasiZero 0.000001 

#define NTOPOL 4 //distinct topologies: (C(A,B)) ; (B(A,C)) ; (A(B,C)) ; (A,B,C)
#define NBRANCH 4 //relevant ranches: a, b, c, d
#define NTIME 3 //time-dimensioned parameter: tau1, tau2, theta

#define MAX_PHGT 0.666  //maximal value for pab+pac+pbc 


#define DEFAULT_POSTERIOR_THRESHOLD 0.95
#define DEFAULT_nb_ILS_times 1
#define DEFAULT_nb_HGT_times 2
#define DEFAULT_HGT_TIME_1 1
#define DEFAULT_HGT_TIME_2 0.5
#define DEFAULT_STRICT 0
#define DEFAULT_REQUIRE_TO_MONOPH 1
#define DEFAULT_REQUIRE_OUTG_MONOPH 1
#define DEFAULT_MAX_ALPHAI 10
#define DEFAULT_MAX_BL 200
#define DEFAULT_MAX_CLOCK_RATIO 2
#define DEFAULT_TRIPLET_CLOCK_THRESH 0.0
#define DEFAULT_MIN_D 0.5
#define DEFAULT_MAX_BACKWARD_MOVES 10
#define DEFAULT_MAX_ITERATIONS 1000
#define DEFAULT_CONVERGENCE 0.000001
#define DEFAULT_DO_CI1 0
#define DEFAULT_DO_CI2 0
#define DEFAULT_CI_THR 1.92

#define MAX_nb_HGT_times 100
#define NB_FILTERING_CASES 11

int debug;

/*
int debug_deriv;
double delta_t=0.00000001;
double delta_p=0.000001;
double debL, debLp, debdLt1, debdLpab, dlup, dldown, debLup, debLdo, debd2L;  
*/

double* log_facto;

int** topol_scenario;
double *p, *pr;
double *dpdt1, *dpdt2, *dpdth; 
double *dprdt1, *dprdt2, *dprdth;
double *d2pdt1, *d2pdt2, *d2pdth; 
double *d2prdt1, *d2prdt2, *d2prdth;
double **exp_lg;
double ***exp_lg_coeff; 
//exp_lg_coeff[i][j][0]*t1+exp_lg_coeff[i][j][1]*t2+exp_lg_coeff[i][j][2]*th is the expected length of branch j under scenario i
//branch in order a b c d l
//scenarios from 1 to nb_scenario_plus1 (0 not used)
//detailed coeffs are needed for calculating derivatives

#define NB_THETA_START 4
double theta_start[NB_THETA_START]={0.0001, 0.001, 0.005, 0.01};

int nb_scenario_plus1;
char** scenario_name;
//int** scenario_type; //1=basic; 2=ILS; 3=HGT
int it;

int nb_parameters;
char** parameter_name;


struct option{
  int nb_ILS_times; //number of times for the second (oldest) ILS time; uniform distribution over these times assumed, with mean equal to t2+theta/6+theta/2
  double* discr_expo_mean1;  //means of equal-weights bins of a discretized exponential distribution of mean 1
  int nb_HGT_times;
  double* HGT_time_coeff;  //the ith possible HGT time is HGT_time_coeff[i]*t1
  int strict_no_event;  //if true, under the no-event scenario, the younger coalescence time must be <t2; free otherwise (I didn't find this useful on simul)

  int require_triplet_other_monophyly;
  int require_outgroup_monophyly;
  double max_alphai; 
  double max_bl;
  double max_clock_ratio;
  double triplet_clock_thresh;
  
  int max_backward_moves;
  int max_iterations;
  double convergence;
    
  double posterior_threshold;  //a posterior above this threshold => a scenario or scenario type is considered well supported, for a specific locus
  
  double min_d; //if d*lgseq is below this threshold triplet topology will be condidered unresolved
  
  int optimize_tau1, optimize_tau2, optimize_theta;
  int optimize_pab, optimize_pac, optimize_pbc, optimize_pone;
  double fixed_tau1, fixed_tau2, fixed_theta;
  double fixed_pab, fixed_pac, fixed_pbc, fixed_pone;
  
  int do_ILS_GF_CI, do_param_CI; //confidence interval 
  double CI_max_lnL_diff; //max difference in log likelihood between optimum and any theta value in confidence interval  

  int verbose; //1=normal output; 0=one-line output
};


struct gene_data{
  char* name;
  int topology; //0=canonical=(C,(B,A)) ; 1=anomalous1=(B,(C,A)); 2=anomalous2=(A,(B,C)); 3=unresolved=(A,B,C)
  double a, b, c, d, e, l;
  double triplet_l, non_triplet_l;  //for clock-likeness check
  int lgseq;
  int to_ignore;  //0=ok ; 1=missing >=1 triplet sp ; 2=missing all outgroup sp ; 3=non-monophyletic triplet ; 5=non-monophyletic triplet + other ; 6=non-monophyletic outgroup ; 7=long branch;  8=extreme alpha_i ; 9=slow-fast triplet ; 10=non-clock-like triplet
};


struct relevant_taxa{
  char** triplet; //A, B, C, in this order; must be monophyletic
  int nb_outg;
  char** outg; //outgroup; must be monophyletic and sister group to triplet
  int nb_oth;
  char** oth;  //other taxa to be used, together with outgroup, when calculating reference length
};


struct param{
  double tau1, tau2;
  double theta;
  double pab, pac, pbc;
  double pone;
  double* alpha;
  int nbg;
  double l_ref;
};


typedef struct noeud{
	struct noeud *v1, *v2, *v3;	/* neighbors */
	double l1, l2, l3;		/* branch lengths */
	int alive1, alive2, alive3;
	char nom[MAXLNAME+1];		/* name */	
	struct noeud *remain;		/* unchanged neighbor */
	double lremain;			/* unchanged length */
	int alremain;			/* unchanged status */
	int to_remove;
	double av_l_to_tip;
	int n_to_tip;
} *noeud;


typedef struct branche{
	noeud bout1;
	noeud bout2;
	double length;
	double bootstrap;
} *branche;


struct log_likelihood{
  double lnL;
  double dlnLdt1, dlnLdt2, dlnLdth;
  double dlnLdpab, dlnLdpac, dlnLdpbc, dlnLdpone;
  double d2lnLdt1, d2lnLdt2, d2lnLdth;  
  double d2lnLdpab, d2lnLdpac, d2lnLdpbc, d2lnLdpone;
};


struct gene_likelihood{
  double L;
  double *scenar_p, *scenar_l;
  double *posterior_scenar, posterior_ILS, posterior_HGT;
  int predicted_scenario, predicted_scenario_type;
    
  double dLdt1, dLdt2, dLdth;
  double dLdpab, dLdpac, dLdpbc, dLdpone;
  double d2Ldt1, d2Ldt2, d2Ldth;
  double d2Ldpab, d2Ldpac, d2Ldpbc, d2Ldpone;  
};


struct scenario_likelihood{
  double L;
  double dLdt1, dLdt2, dLdth;  
  double d2Ldt1, d2Ldt2, d2Ldth;  
};


struct posterior{
  double* scenar;
  double basic;
  double noconflict_ILS, conflict_ILS; 
  double noconflict_HGT, conflict_HGT;
  double asymmetry_ILS, asymmetry_HGT;
  int dominant_ILS, dominant_HGT; //1: ((AC)B) > ((BC)A) ; 2: ((AC)B) < ((BC)A) 
};


struct confidence_interval{
  double max_lnL_diff; //max difference in lnL between the optimum and any point in the CI
  char master_param[10]; //either "tau1", "tau2", "theta", "pab", "pac", "pbc", or "pone"; refers to the parameter that will be varied around its optimum 
  double master_param_min, master_param_max; //min and max value of the master parameter (boundaries of the CI)
  struct param par_min, par_max; //optimized parameters at min and max
  struct posterior post_min, post_max; //posterior annotations at min and max
};

noeud root;


/* usage */
void usage(){
  printf("usage: aphid tree_file taxon_file option_file outfile \n");
  exit(0);
}


/* fin */
void fin(char* s){
  printf("%s\n", s);
  exit(0);
}


/* check_alloc */
/* Return a pointer with the required allocated memory or leave program */
/* if not enough memory. */

void *check_alloc(int nbrelt, int sizelt)
{
void *retval;

if( (retval=(char*)calloc(nbrelt,sizelt)) != NULL ) { return retval; }
printf("Not enough memory\n");
exit(0);
}


/*******************************/
/* tree manipulation functions */
/*******************************/







/***********************************/
/* end tree manipulation functions */
/***********************************/




/* allocate_global_var */
void allocate_global_var(){

  int i, j;
  
  topol_scenario=(int**)check_alloc(NTOPOL, sizeof(int*));
  for(i=0;i<NTOPOL;i++) topol_scenario[i]=(int*)check_alloc(nb_scenario_plus1, sizeof(int));
  
  p=(double*)check_alloc(nb_scenario_plus1, sizeof(double));
  pr=(double*)check_alloc(nb_scenario_plus1, sizeof(double));
  dpdt1=(double*)check_alloc(nb_scenario_plus1, sizeof(double));
  dpdt2=(double*)check_alloc(nb_scenario_plus1, sizeof(double));
  dpdth=(double*)check_alloc(nb_scenario_plus1, sizeof(double));
  dprdt1=(double*)check_alloc(nb_scenario_plus1, sizeof(double));
  dprdt2=(double*)check_alloc(nb_scenario_plus1, sizeof(double));
  dprdth=(double*)check_alloc(nb_scenario_plus1, sizeof(double));
  d2pdt1=(double*)check_alloc(nb_scenario_plus1, sizeof(double));
  d2pdt2=(double*)check_alloc(nb_scenario_plus1, sizeof(double));
  d2pdth=(double*)check_alloc(nb_scenario_plus1, sizeof(double));
  d2prdt1=(double*)check_alloc(nb_scenario_plus1, sizeof(double));
  d2prdt2=(double*)check_alloc(nb_scenario_plus1, sizeof(double));
  d2prdth=(double*)check_alloc(nb_scenario_plus1, sizeof(double));
  
  exp_lg=(double**)check_alloc(nb_scenario_plus1, sizeof(double*));
  for(i=0;i<nb_scenario_plus1;i++) exp_lg[i]=(double*)check_alloc(NBRANCH, sizeof(double));
  exp_lg_coeff=(double***)check_alloc(nb_scenario_plus1, sizeof(double**));
  for(i=0;i<nb_scenario_plus1;i++){
    exp_lg_coeff[i]=(double**)check_alloc(NBRANCH, sizeof(double*));
    for(j=0;j<NBRANCH;j++) exp_lg_coeff[i][j]=(double*)check_alloc(NTIME, sizeof(double));
  }
  
}




/* max_rank */

void max_rank(double* tab, int n, double* max, int* max_r){

  int i;
  
  *max=tab[0]; *max_r=0;
  for(i=1;i<n;i++){
    if(tab[i]>*max){
      *max=tab[i];
      *max_r=i;
    }
  }
}



/* discretize_exponential */
/* returns an pointer towards an array of nb_bin positive values, which are the means of a discretization of the exponential distribution of mean mean in nb_bin equally weight bins */
double* discretize_exponential(double mean, int nb_bin){
  int i;
  double lambda;
  double* boundaries;
  double* bin_means;
 
  lambda=1/mean;
  boundaries=(double*)check_alloc(nb_bin, sizeof(double));
  bin_means=(double*)check_alloc(nb_bin, sizeof(double));  
  
  //calculate boundaries
  boundaries[0]=0;
  for(i=1;i<nb_bin;i++)
    boundaries[i]=-log(exp(-lambda*boundaries[i-1])-1./nb_bin)/lambda;
    
  //calculate bin means
  for(i=0;i<nb_bin-1;i++)
    bin_means[i]=nb_bin*(boundaries[i]*exp(-lambda*boundaries[i])-boundaries[i+1]*exp(-lambda*boundaries[i+1])+(exp(-lambda*boundaries[i])-exp(-lambda*boundaries[i+1]))/lambda);
  bin_means[nb_bin-1]=nb_bin*(boundaries[nb_bin-1]*exp(-lambda*boundaries[nb_bin-1])+exp(-lambda*boundaries[i])/lambda);

  return bin_means;    
}



/* has_long_branches */
int has_long_branches(struct gene_data* gd, struct option *opt){

  if(gd->a*gd->lgseq>opt->max_bl) return 1;
  if(gd->b*gd->lgseq>opt->max_bl) return 1;
  if(gd->c*gd->lgseq>opt->max_bl) return 1;
  if(gd->d*gd->lgseq>opt->max_bl) return 1;
  if(gd->e*gd->lgseq>opt->max_bl) return 1;
  return 0;
}


/* has_unequal_branches */
int has_unequal_branches(struct gene_data* gd, struct option* opt){

  int i, b1, b2, b3;
  int maxb;
  double meanb, pval;
  
  if(gd->topology==0){
    b1=(int)((gd->a+gd->d)*gd->lgseq);
    b2=(int)((gd->b+gd->d)*gd->lgseq);
    b3=(int)(gd->c*gd->lgseq);
  }
  else if(gd->topology==1){
    b1=(int)((gd->a+gd->d)*gd->lgseq);
    b2=(int)((gd->c+gd->d)*gd->lgseq);
    b3=(int)(gd->b*gd->lgseq);
  }
  else if(gd->topology==2){
    b1=(int)((gd->b+gd->d)*gd->lgseq);
    b2=(int)((gd->c+gd->d)*gd->lgseq);
    b3=(int)(gd->a*gd->lgseq);
  }
  else if(gd->topology==3){
    b1=(int)(gd->a*gd->lgseq);
    b2=(int)(gd->b*gd->lgseq);
    b3=(int)(gd->c*gd->lgseq);
  }
 
  maxb=b1; if(b2>maxb) maxb=b2; if(b3>maxb) maxb=b3;
  meanb=(b1+b2+b3)/3.;
  
  if(maxb==0) return 0;
  
  pval=0;
  
  for(i=maxb;i<opt->max_bl;i++){
    pval+=exp(i*log(meanb)-meanb-log_facto[i]);
  }
  
  if(pval<opt->triplet_clock_thresh)
    return 1;

  return 0;
}


/* copy_data */
void copy_data(struct gene_data* from, struct gene_data* to){

  to->a=from->a;
  to->b=from->b;
  to->c=from->c;
  to->d=from->d;
  to->e=from->e;
  to->l=from->l;
  to->triplet_l=from->triplet_l;
  to->non_triplet_l=from->non_triplet_l;
  to->topology=from->topology;
  to->lgseq=from->lgseq;
  if(strlen(to->name)<strlen(from->name)){
    to->name=(char*)check_alloc((int)strlen(from->name)+1, sizeof(char));
  }
  sprintf(to->name, "%s", from->name);  
}




/* find_dubious_genes */

void find_dubious_genes(struct gene_data* gd, int nb_gene, struct option* opt){

  int i, j;
  double meanl, meantl, meanntl, alpha, ratio_i, ratio_mean;


  //genes with exceedingly long branches a, b, c or d
  for(i=0;i<nb_gene;i++){
    if(gd[i].to_ignore==0 && has_long_branches(gd+i, opt))
      gd[i].to_ignore=7;
  }

  //genes with anomalously long/short reference l
  meanl=0.; for(i=0;i<nb_gene;i++) meanl+=gd[i].l; meanl/=nb_gene;
  for(i=0;i<nb_gene;i++){
    if(gd[i].to_ignore==0 && (gd[i].l/meanl>opt->max_alphai || meanl/gd[i].l>opt->max_alphai))
      gd[i].to_ignore=8;
  }

  //clock departure 1: slow/fast triplet
  meantl=0; for(i=0;i<nb_gene;i++) meantl+=gd[i].triplet_l; meantl/=nb_gene;
  meanntl=0; for(i=0;i<nb_gene;i++) meanntl+=gd[i].non_triplet_l; meanntl/=nb_gene;
  ratio_mean=meantl/meanntl;
  for(i=0;i<nb_gene;i++){
    ratio_i=gd[i].triplet_l/gd[i].non_triplet_l;
//if(gd[i].to_ignore==0) printf("ratio_i %f, ratio_mean %f, max %f\n", ratio_i, ratio_mean, opt->max_clock_ratio);
    if(gd[i].to_ignore==0 && (ratio_i/ratio_mean>opt->max_clock_ratio || ratio_mean/ratio_i>opt->max_clock_ratio))
      gd[i].to_ignore=9;
//if(gd[i].to_ignore==0 && (ratio_i>opt->max_clock_ratio || 1./ratio_i>opt->max_clock_ratio)){
//  gd[i].to_ignore=9;
//  printf("slow fast\n");
//}

  }
  
  //clock departure 2: non-clock triplet

  for(i=0;i<nb_gene;i++){
    if(gd[i].to_ignore==0 && has_unequal_branches(gd+i, opt))
      gd[i].to_ignore=10;
  }   


}



/* print_ignored */

void print_ignored(struct gene_data* gd, int nb_gene){

  int i, nbign_tot, nbign[NB_FILTERING_CASES];
  
  for(i=0;i<NB_FILTERING_CASES;i++) nbign[i]=0;
  for(i=0;i<nb_gene;i++) nbign[gd[i].to_ignore]++;
  nbign_tot=0;
  for(i=1;i<NB_FILTERING_CASES;i++) nbign_tot+=nbign[i]; 
  
  printf("ignoring %d genes: ", nbign_tot);
  printf("%d missing triplet; %d no outgroup; %d unexpected topology; %d unexpected branch lengths\n", nbign[1], nbign[2], nbign[3]+nbign[4]+nbign[5]+nbign[6], nbign[7]+nbign[8]+ nbign[9]+nbign[10]);
  if(nbign[3]+nbign[4]+nbign[5]+nbign[6]) printf("[unexpected topology: %d triplet paraphyly; %d wrong root; %d triplet+other paraphyly; %d outgroup paraphyly]\n", nbign[3], nbign[4], nbign[5], nbign[6]);
  if(nbign[10]>0) 
    printf("[unexpected branch lengths: %d long branch; %d extreme alpha_i; %d slow/fast triplet; %d non-clock triplet]\n", nbign[7], nbign[8], nbign[9], nbign[10]); //user is expert and has activated triplet-clock filter (not mentioned in the doc)
  else
    if(nbign[7]+nbign[8]+nbign[9]) printf("[unexpected branch lengths: %d long branch; %d extreme alpha_i; %d slow/fast triplet]\n", nbign[7], nbign[8], nbign[9]); //normal case, no triplet-clock filter
  printf("analyzing %d gene trees\n", nb_gene-nbign_tot);
   
}



/* do_ignore */

void do_ignore(struct gene_data* gd, int* nb_gene){

  int i, j, nbg;
  nbg=*nb_gene;
  
  for(i=nbg-1;i>=0;i--){
    if(gd[i].to_ignore){
      for(j=i;j<nbg-1;j++){
        copy_data(gd+j+1, gd+j);
      }
      *nb_gene=*nb_gene-1;
    }
  }
}



/* set_unresolved_topology */
int set_unresolved_topology(struct gene_data* gd, int nb_gene, double min){

  int i;
  
  for(i=0;i<nb_gene;i++)
    if(gd[i].d*gd[i].lgseq<min)
      gd[i].topology=3;

}


/* print_options */
void print_options(struct option * opt){
  int i;
  int optimize_all;


  printf("HGT times (relative to tau1):");
  for(i=0;i<opt->nb_HGT_times;i++) printf(" %f", opt->HGT_time_coeff[i]);
  printf("\n");
  if(opt->nb_ILS_times>1){
    printf("ILS coalescence times (relative to theta/2):");
    for(i=0;i<opt->nb_ILS_times;i++) printf(" %f", opt->discr_expo_mean1[i]);
    printf("\n");
  }
  printf("topology requirements:");
  if(opt->require_triplet_other_monophyly) printf(" triplet+other monophyly");
  if(opt->require_outgroup_monophyly){
    if(opt->require_triplet_other_monophyly) printf(";");
    printf(" outgroup monophyly");
  }
  if(!opt->require_triplet_other_monophyly && !opt->require_outgroup_monophyly) printf(" none");
  printf("\n");
  printf("min_d: %f\n", opt->min_d);  
    
  printf("filter thresholds: max_alphai=%.3f, max_bl=%.1f, max_clock_ratio=%.3f\n", opt->max_alphai, opt->max_bl, opt->max_clock_ratio);

  //printf("optimization: max_iterations=%d, max_backward_moves=%d, convergence=%f\n", opt->max_iterations, opt->max_backward_moves, opt->convergence);

  optimize_all=(opt->optimize_tau1 || opt->optimize_tau2 || opt->optimize_theta || opt->optimize_pab || opt->optimize_pac || opt->optimize_pbc || opt->optimize_pone);
  if(!optimize_all){
    printf("do not optimize ");
    if(!opt->optimize_tau1) printf("tau1 [%f] ", opt->fixed_tau1);
    if(!opt->optimize_tau2) printf("tau2 [%f] ", opt->fixed_tau2);
    if(!opt->optimize_theta) printf("theta [%f] ", opt->fixed_theta);
    if(!opt->optimize_pab) printf("pab [%f] ", opt->fixed_pab);
    if(!opt->optimize_pac) printf("pac [%f] ", opt->fixed_pac);
    if(!opt->optimize_pbc) printf("pbc [%f] ", opt->fixed_pbc);
    if(!opt->optimize_pone) printf("pone [%f] ", opt->fixed_pone);
  }
  
  printf("Confidence intervals:" );
  if(opt->do_ILS_GF_CI) printf("ILS/GF ");
  if(opt->do_param_CI) printf("parameters ");
  if(!opt->do_ILS_GF_CI && !opt->do_param_CI) printf("none");
  else printf("[%f]", opt->CI_max_lnL_diff);
  printf("\n");
  

  //printf("posterior_threshold: %.3f\n", opt->posterior_threshold);  
}



/* set_default_options */
void set_default_options(struct option * opt){

  opt->verbose=1;
  opt->posterior_threshold=DEFAULT_POSTERIOR_THRESHOLD;
  opt->nb_ILS_times=DEFAULT_nb_ILS_times;
  opt->strict_no_event=DEFAULT_STRICT;
  opt->require_triplet_other_monophyly=DEFAULT_REQUIRE_TO_MONOPH;
  opt->require_outgroup_monophyly=DEFAULT_REQUIRE_OUTG_MONOPH;
  opt->max_alphai=DEFAULT_MAX_ALPHAI;
  opt->max_bl=DEFAULT_MAX_BL;
  opt->max_clock_ratio=DEFAULT_MAX_CLOCK_RATIO;
  opt->triplet_clock_thresh=DEFAULT_TRIPLET_CLOCK_THRESH;
  opt->max_backward_moves=DEFAULT_MAX_BACKWARD_MOVES;
  opt->max_iterations=DEFAULT_MAX_ITERATIONS;
  opt->convergence=DEFAULT_CONVERGENCE;
  opt->nb_HGT_times=DEFAULT_nb_HGT_times;
  opt->HGT_time_coeff=(double*)check_alloc(MAX_nb_HGT_times, sizeof(double));
  opt->HGT_time_coeff[0]=DEFAULT_HGT_TIME_1;
  opt->HGT_time_coeff[1]=DEFAULT_HGT_TIME_2;
  opt->optimize_tau1=opt->optimize_tau2=opt->optimize_theta=1;
  opt->optimize_pab=opt->optimize_pac=opt->optimize_pbc=opt->optimize_pone=1;
  opt->fixed_tau1=opt->fixed_tau2=opt->fixed_theta=-1.;
  opt->fixed_pab=opt->fixed_pac=opt->fixed_pbc=opt->fixed_pone=-1.;
  opt->min_d=DEFAULT_MIN_D;
  opt->do_ILS_GF_CI=DEFAULT_DO_CI1;
  opt->do_param_CI=DEFAULT_DO_CI2;
  opt->CI_max_lnL_diff=DEFAULT_CI_THR;
  
}




/* set_scenario_names */
void set_scenario_names(struct option* opt){
  int i;
  
  scenario_name=(char**)check_alloc(nb_scenario_plus1, sizeof(char*));
  for(i=1;i<=nb_scenario_plus1;i++) scenario_name[i]=(char*)check_alloc(101, sizeof(char));
  sprintf(scenario_name[1], "no_event");
  
  if(opt->nb_ILS_times>1){
    for(i=0;i<opt->nb_ILS_times;i++){
      sprintf(scenario_name[2+i*3],"ILS_%.3f_((AB)C)", opt->discr_expo_mean1[i]);
      sprintf(scenario_name[3+i*3],"ILS_%.3f_((AC)B)", opt->discr_expo_mean1[i]);
      sprintf(scenario_name[4+i*3],"ILS_%.3f_((BC)A)", opt->discr_expo_mean1[i]);
    }
  }
  else{
    sprintf(scenario_name[2],"ILS_((AB)C)");
    sprintf(scenario_name[3],"ILS_((AC)B)");
    sprintf(scenario_name[4],"ILS_((BC)A)");  
  }
  
  if(opt->nb_HGT_times>1){
    for(i=0;i<opt->nb_HGT_times;i++){
      sprintf(scenario_name[2+3*opt->nb_ILS_times+i*5],"GF_%.3f_A<->B", opt->HGT_time_coeff[i]);
      sprintf(scenario_name[3+3*opt->nb_ILS_times+i*5],"GF_%.3f_A->C", opt->HGT_time_coeff[i]);
      sprintf(scenario_name[4+3*opt->nb_ILS_times+i*5],"GF_%.3f_C->A", opt->HGT_time_coeff[i]);
      sprintf(scenario_name[5+3*opt->nb_ILS_times+i*5],"GF_%.3f_B->C", opt->HGT_time_coeff[i]);
      sprintf(scenario_name[6+3*opt->nb_ILS_times+i*5],"GF_%.3f_C->B", opt->HGT_time_coeff[i]);
    }
  }
  else{
    sprintf(scenario_name[2+3*opt->nb_ILS_times],"GF_A<->B");
    sprintf(scenario_name[3+3*opt->nb_ILS_times],"GF_A->C");
    sprintf(scenario_name[4+3*opt->nb_ILS_times],"GF_C->A");  
    sprintf(scenario_name[5+3*opt->nb_ILS_times],"GF_B->C");
    sprintf(scenario_name[6+3*opt->nb_ILS_times],"GF_C->B");  
  }
    
}


/* read_1_param */

double read_1_param(char* line, char* deb){

  char* prov;
  int i, reti;
  double par;

  i=0;
  
  while(line[i]==deb[i]) i++;
  prov=line+i;
  while(*prov==' ' || *prov==':' || *prov=='=') prov++;
  reti=sscanf(prov, "%le", &par);
  if(reti==1) return par;
  return -1.; 

}

/* read_n_param */

double* read_n_param(char* line, char* deb, int nb){

  char* prov;
  int i, j, reti;
  double* par;

  par=(double*)check_alloc(nb, sizeof(double));

  i=0;
  while(line[i]==deb[i]) i++;
  
  prov=strtok(line+i, " :=,;\t");
  
  for(j=0;j<nb;j++){
    reti=sscanf(prov, "%le", par+j);
    if(reti!=1) return NULL;
    prov=strtok(NULL, " :=,;\t");
  }
  
  return par; 
}



/* read_option_file */
void read_option_file(FILE* optf, struct option * opt){
  char line[MAXLLINE+1];

  set_default_options(opt);
  
  while(fgets(line, MAXLLINE, optf)){
    if(strncmp(line, "verbose", 7)==0) opt->verbose=read_1_param(line, "verbose");
    if(strncmp(line, "posterior_threshold", 13)==0) opt->posterior_threshold=read_1_param(line, "posterior_threshold");
    if(strncmp(line, "nb_ILS_times", 11)==0) opt->nb_ILS_times=read_1_param(line, "nb_ILS_times");
    if(strncmp(line, "nb_GF_times", 11)==0) opt->nb_HGT_times=read_1_param(line, "nb_GF_times");
    if(strncmp(line, "GF_time_coeff", 11)==0) opt->HGT_time_coeff=read_n_param(line, "GF_time_coeff", opt->nb_HGT_times);
    if(strncmp(line, "strict_no_event", 15)==0) opt->strict_no_event=read_1_param(line, "strict_no_event");
    if(strncmp(line, "require_triplet_other_monophyly", 30)==0) opt->require_triplet_other_monophyly=read_1_param(line, "require_triplet_other_monophyly");
    if(strncmp(line, "require_outgroup_monophyly", 25)==0) opt->require_outgroup_monophyly=read_1_param(line, "require_outgroup_monophyly");
    if(strncmp(line, "max_alphai", 9)==0) opt->max_alphai=read_1_param(line, "max_alphai");
    if(strncmp(line, "max_bl", 6)==0) opt->max_bl=read_1_param(line, "max_bl");
    if(strncmp(line, "max_clock_ratio", 15)==0) opt->max_clock_ratio=read_1_param(line, "max_clock_ratio");
    if(strncmp(line, "triplet_clock_thresh", 20)==0) opt->triplet_clock_thresh=read_1_param(line, "triplet_clock_thresh");
    if(strncmp(line, "max_backward_moves", 15)==0) opt->max_backward_moves=read_1_param(line, "max_backward_moves");
    if(strncmp(line, "max_iterations", 12)==0) opt->max_iterations=read_1_param(line, "max_iterations");
    if(strncmp(line, "convergence", 10)==0) opt->convergence=read_1_param(line, "convergence");
    if(strncmp(line, "fixed_tau1", 10)==0) {opt->fixed_tau1=read_1_param(line, "fixed_tau1"); if(opt->fixed_tau1>0.) opt->optimize_tau1=0;}
    if(strncmp(line, "fixed_tau2", 10)==0) {opt->fixed_tau2=read_1_param(line, "fixed_tau2"); if(opt->fixed_tau2>0.) opt->optimize_tau2=0;}
    if(strncmp(line, "fixed_theta", 11)==0) {opt->fixed_theta=read_1_param(line, "fixed_theta"); if(opt->fixed_theta>0.) opt->optimize_theta=0;}
    if(strncmp(line, "fixed_pab", 9)==0) {opt->fixed_pab=read_1_param(line, "fixed_pab"); if(opt->fixed_pab>0.) opt->optimize_pab=0;}
    if(strncmp(line, "fixed_pac", 9)==0) {opt->fixed_pac=read_1_param(line, "fixed_pac"); if(opt->fixed_pac>0.) opt->optimize_pac=0;}
    if(strncmp(line, "fixed_pbc", 9)==0) {opt->fixed_pbc=read_1_param(line, "fixed_pbc"); if(opt->fixed_pbc>0.) opt->optimize_pbc=0;}
    if(strncmp(line, "fixed_panc", 10)==0) {opt->fixed_pone=read_1_param(line, "fixed_panc"); if(opt->fixed_pone>0.) opt->optimize_pone=0;}
    if(strncmp(line, "min_d", 5)==0) opt->min_d=read_1_param(line, "min_d");
    if(strncmp(line, "CI_ILS_GF", 9)==0) opt->do_ILS_GF_CI=read_1_param(line, "CI_ILS_GF");
    if(strncmp(line, "CI_param", 8)==0) opt->do_param_CI=read_1_param(line, "CI_param");
  }  
  
}



/* set_parameter_names */
void set_parameter_names(){

  int i;
  
  nb_parameters=7;
  
  parameter_name=(char**)check_alloc(nb_parameters, sizeof(char*));
  for(i=0;i<nb_parameters;i++)
    parameter_name[i]=(char*)check_alloc(MAX_PARAMETER_NAME_LG+1, sizeof(char));

  sprintf(parameter_name[0], "tau1");
  sprintf(parameter_name[1], "tau2");
  sprintf(parameter_name[2], "theta");
  sprintf(parameter_name[3], "pab");
  sprintf(parameter_name[4], "pac");
  sprintf(parameter_name[5], "pbc");
  sprintf(parameter_name[6], "pone");

}



/* count tokens */
int count_tokens(char* s, char* sep){

  char* s2;
  int n;
  
  s2=check_alloc(strlen(s)+1, sizeof(char));
  sprintf(s2, "%s", s);
  
  if(strtok(s2, sep)==NULL) return 0;
  n=1;
  while(strtok(NULL, sep)) n++;
  return n;
}



/* read_taxon_file */
void read_taxon_file(FILE* taxf, struct relevant_taxa* tax){

  char line[MAXLLINE+1];
  char* prov, *prov2;
  int i, triplet_ok;
  
  triplet_ok=tax->nb_outg=tax->nb_oth=0;
  
  while(fgets(line, MAXLLINE, taxf)){
    prov=line;
    while(*prov==' ' || *prov=='\t') prov++;
    if(strncmp(prov, "triplet", 7)==0){
      tax->triplet=(char**)check_alloc(3, sizeof(char*));
      strtok(prov, " :=,;\t\n");
      for(i=0;i<3;i++){
        prov2=strtok(NULL, " :=,;\t\n"); if(prov2==NULL) fin("unexpected triplet line in taxon file");
        tax->triplet[i]=check_alloc(strlen(prov2)+1, sizeof(char)); sprintf(tax->triplet[i], "%s", prov2);
      }
      triplet_ok=1;
    }
    else if(strncmp(prov, "outgroup", 8)==0){
      tax->nb_outg=count_tokens(prov, " :=,;\t\n")-1;
      tax->outg=check_alloc(tax->nb_outg, sizeof(char*));
      strtok(prov, " :=,;\t");
      for(i=0;i<tax->nb_outg;i++){
        prov2=strtok(NULL, " :=,;\t\n");
        tax->outg[i]=check_alloc(strlen(prov2)+1, sizeof(char));
        sprintf(tax->outg[i], "%s", prov2);
      }
    }
    else if(strncmp(prov, "other", 5)==0){
      tax->nb_oth=count_tokens(prov, " :=,;\t\n")-1;
      tax->oth=check_alloc(tax->nb_oth, sizeof(char*));
      strtok(prov, " :=,;\t");
      for(i=0;i<tax->nb_oth;i++){
        prov2=strtok(NULL, " :=,;\t\n");
        tax->oth[i]=check_alloc(strlen(prov2)+1, sizeof(char));
        sprintf(tax->oth[i], "%s", prov2);
      }
    }
  }
  
  if(triplet_ok==0) fin("missing triplet line in taxon file");
  if(tax->nb_outg==0) fin("missing outgroup line in taxon file");
  

}


/* print_taxa */
void print_taxa(struct relevant_taxa *tax){

  int i;

  printf("found triplet: %s %s %s\n", tax->triplet[0], tax->triplet[1], tax->triplet[2]);

  printf("found %d outgroup taxa", tax->nb_outg);
  if(tax->nb_outg<6){
    printf(": ");
    for(i=0;i<tax->nb_outg-1;i++) printf("%s, ", tax->outg[i]);
    printf("%s\n", tax->outg[tax->nb_outg-1]);
  }
  else printf("\n");
  
  if(tax->nb_oth){
    printf("found %d other relevant taxa", tax->nb_oth);
    if(tax->nb_oth<6){
      printf(": ");
      for(i=0;i<tax->nb_oth-1;i++) printf("%s, ", tax->oth[i]);
      printf("%s\n", tax->oth[tax->nb_oth-1]);
    }
    else printf("\n");
  }
}



int retder(int *liste){
  int i=0, j;
  while (liste[i] != 0) i++;
  j = *(liste + i - 1);
  *(liste + i - 1) = 0;
  return j;
}

void aj(int *liste, int nb) {
  int  i=0;
  while (liste[i] != 0) i++;
  *(liste + i) = nb;
  return;
}


int is_in_list(char* s, char** sel, int nbsel){

  int i;
  
  for(i=0;i<nbsel;i++)
    if(strcmp(s, sel[i])==0) return 1;
  return 0;

}




/* samename */
/* TRUE if names name1 and name2 are equal, FALSE otherwise */

int samename(char* name1, char* name2){

  int i=0;
  char c1, c2;

  while(i<20){
    c1=name1[i]; c2=name2[i];
    if((c1==' ' || c1=='\n' || c1=='\t' || c1==0) && (c2==' ' || c2=='\n' || c2=='\t' || c2==0)) return 1;
    if(c1!=c2) return 0;
    i++;
  }
  return 1;
}



/* equallist */
/* Check the list of taxa down a branch. */

int equallist(noeud from, noeud nd, int nb, char** list){
  noeud n1, n2;
  int i, flag;


  if(from==nd->v1) { n1=nd->v2; n2=nd->v3; }
  else if(from==nd->v2) { n1=nd->v1; n2=nd->v3; }
  else if(from==nd->v3) { n1=nd->v1; n2=nd->v2; }
  else { printf("Erreur\n"); }

  if(!n1){
    for(i=0;i<nb;i++){
      flag=0;
      if(samename(nd->nom, list[i])) {flag=1; break; }
    }
    return flag;
  }

  return(equallist(nd, n1, nb, list) && equallist(nd, n2, nb, list));
}








/* isbranch */
/* Check wether lists of taxa name list1 and list2 (respectively l1 */
/* and l2 taxa) exactly match the lists of taxa laying left-side  */
/* and right-side branch br. */

int isbranch(branche br, int l1, char** list1, int l2, char** list2){

  if( equallist(br->bout1, br->bout2, l1, list1) && equallist(br->bout2, br->bout1, l2, list2)) return 2;
  if( equallist(br->bout2, br->bout1, l1, list1) && equallist(br->bout1, br->bout2, l2, list2)) return 1;
  return 0;
}


/* organize_tree */

void organize_tree(noeud from, noeud nd){
  noeud prov;
  double lprov;
  int alprov;
  static int numappel=0;

  if(!nd->v1) return;
  if(nd!=root){
    numappel++;
    sprintf(nd->nom, "int%d", numappel);
  }
  if(from==nd->v1) {
    prov=nd->v3; nd->v3=nd->v1; nd->v1=prov;
    lprov=nd->l3; nd->l3=nd->l1; nd->l1=lprov;
    alprov=nd->alive3; nd->alive3=nd->alive1; nd->alive1=alprov;
  }
  if(from==nd->v2) {
    prov=nd->v3; nd->v3=nd->v2; nd->v2=prov;
    lprov=nd->l3; nd->l3=nd->l2; nd->l2=lprov;
    alprov=nd->alive3; nd->alive3=nd->alive2; nd->alive2=alprov;
  }
  organize_tree(nd, nd->v1);
  organize_tree(nd, nd->v2);

  if(nd==root) numappel=0;
}




/* taxa_count */
/* Return the number of taxa in c_tree ctree. */

int taxa_count(char* ctree){
  char* prov;
  int nbpo=0, nbpf=0, nbv=0;

  prov=ctree;
  while(*prov!='[' && *prov!='(' && *prov) prov++;
  if(*prov=='\0'){
    printf("Trees should begin with '[' or '('\n");
    return -1;
  }
  if(*prov=='['){
    while(*prov!=']') prov++;
    while(*prov!='(') prov++;
  }

/* debut arbre */

  while(*prov && *prov!=';'){
    if(*prov=='(') nbpo++;
    if(*prov==')') nbpf++;
    if(nbpo-nbpf==1 && *prov==',') nbv++;
    prov++;
  }

  if(nbv==1){ /* rooted */
    return (nbpo+1);
  }
  if(nbv==2){ /* unrooted */
    return (nbpo+2);
  }
  printf("Bad number of commas ',' in tree\n");
  return -1;
}





/* ctot */
/* ChaineTOTable. Read c_tree input (string) and write t_tree arbre (int**), */
/* branch length lgbi (internal) and lgbp (terminal), bootstrap values, */
/* species names (nom) and rooted/unrooted (racine-> r (rooted) or n (not). */

int ctot(char *input, int *arbre[], double *lgbi, double *lgbp, double* bootstrap, char *nom[], char *racine, int* nbbi){


  int   i=0, j, k, fdf=0, fin=0, nbpo=0,
	nbpf = 0, cptv1 = 0, nbotu, br_ouverte, *listecour,
	otu = -1, pomoinspf, cpttree=0, t=1, nbnom=0;
  char  c, cc, cas, dejalu = '\0';
  double  f;

  sscanf(input+(cpttree++), "%c", &c);
  if (c == '[') {
    while ((c != ']') && (fdf != EOF)) fdf = sscanf(input+(cpttree++), "%c", &c);
    if (c != ']'){
      printf("Unmatched '[' ']'\n");
      return -1;
    }
  }
  else
    if (c == '(') cpttree=0;
    else{
      printf("Tree file 1st character must be '(' or '['\n");
      return -1;
    }
  while ((c != ';') && (fdf != EOF)) {
    fdf = sscanf(input+(cpttree++), "%c", &c);
      if (c == '(') nbpo++;
      if (c == ')') nbpf++;
      if ((nbpo == nbpf + 1) && (c == ',')) cptv1++;
  }

  if (c != ';'){
    printf("';' missing at end of tree file\n");
    return -1;
  }

  if (nbpo != nbpf){
    printf("Unmatched parenthesis\n");
    return -1;
  }

  if (cptv1 == 1) cas = 'c';
  if (cptv1 == 2) cas = 'a';

  if ((cptv1!=1) && (cptv1!=2)){
    printf("Bad number of ',' in tree file\n");
    return -1;
  }

  nbotu = nbpo + 2;


  listecour=(int*) check_alloc(nbotu, sizeof(int));

  if (cas == 'a') *racine='n'; else *racine='r';

  if (cas=='c'){
    if (lgbp) lgbp[0] = 0.;
    strcpy(nom[0],"ROOT");
    for(i=0;i<nbotu-3;i++) arbre[0][i]=0;
    otu++;
  }


  cpttree=0;
  sscanf(input+(cpttree++), "%c", &c);
  if (c == '['){
    while (c != ']') sscanf(input+(cpttree++), "%c", &c);
    while((c==']') || (c==' ') || (c=='\n') || (c=='\t')) sscanf(input+(cpttree++), "%c", &c);
    if (c!='(') return -1;
  }
  else
    while(c!='(') sscanf(input+(cpttree++), "%c", &c);

  pomoinspf=1;
  for (i = 0; i < nbotu; i++) listecour[i] = 0;
  br_ouverte = 0;


  for (k = 0; t==1; k++) {
    if (dejalu == '\0') sscanf(input+(cpttree++), "%c", &c);
    else {
      c = dejalu;
      dejalu = '\0';
    }
    switch (c) {
    case ';': fin = 1; break;
    case ',': case '\n': case '\t': case '\'': case ' ': break;
    case '(':
      pomoinspf ++;
      br_ouverte++;
      aj(listecour, br_ouverte); break;
    case ')':
      pomoinspf--;
      sscanf(input+(cpttree++), "%c", &cc);
      if (cc == ';' || pomoinspf==0) {
	fin = 1; break;
      }
      j = retder(listecour);
      while (cc=='\n' || cc==' ' || cc=='\t') sscanf(input+(cpttree++),"%c",&cc);
      if (strpbrk(input+cpttree-1, "-0123456789.Ee")==input+cpttree-1){
        sscanf(input+cpttree-1, "%le", bootstrap+j-1);
        cpttree+=strspn(input+cpttree-1,"-.0123456789Ee");
        cc=*(input+cpttree-1);
        while (cc=='\n' || cc==' ' || cc=='\t') sscanf(input+(cpttree++),"%c",&cc);
      }
      if (cc == ':') {
	while(input[cpttree]==' ') cpttree++;
	sscanf(input+cpttree, "%le", &f);
	cpttree+=strspn(input+cpttree,"-0123456789.Ee");
	lgbi[j - 1] = f;
      }
      else dejalu=cc;
      break;
    default:
      nbnom++;
      otu++; cc = c; i = 0;
      while ((cc != ':') && (cc != ',') && (cc != ')') && (cc != '\n') && (cc!=' ')) {
	if (cc != '\'') { nom[otu][i] = cc; i++; nom[otu][i]='\0'; }
	sscanf(input+cpttree, "%c", &cc);
        cpttree++;
      }
      while ((cc != ':') && (cc != ',') && (cc != ')') && (cc != '\n')){
	sscanf(input+cpttree, "%c", &cc);
        cpttree++;
      }
      if (cc == ':') {
	while(input[cpttree]==' ') cpttree++;
	sscanf(input+(cpttree), "%le", &f);
	cpttree+=strspn(input+cpttree,"-0123456789.Ee");
	lgbp[otu] = f;
      }
      else dejalu = cc;
      for (i = 0; i < nbotu - 3; i++) arbre[otu][i] = 0;
      for (i = 0; i < nbotu; i++)
	if (listecour[i] != 0)
	  arbre[otu][listecour[i] - 1] = 1;
    }
    if (fin == 1) break;
  }

  if(nbbi) *nbbi=br_ouverte;

  free(listecour);

  return nbnom;

/*  if (cas=='a') return nbotu; else return (nbotu-1);*/
}




/* unroot */
/* Unroot t_tree treer. All parameters are modified (branch lengths, names, ...). */
/* Root info may be kept in list1-list2 l1-l2 fracroot1 if l1 is non NULL and list1 and list2  */
/* are allocated as nbtaxa-sized char* arrays. */

int unroot(int** treer, int notu, double* lgbi, double* lgbp, double* bootvals, char** nom, char** list1, char** list2, int* l1, int* l2, double* fracroot1){

int i, j, j1, j2, br_a_virer=-1, bp_a_modifier=-1, bi_a_modifier=-1, somme1, somme2, drapeau;
  double l;

  for(j=0;j<notu-2;j++) if(treer[0][j]!=0) return 1;
  for(j1=0;j1<notu-3;j1++){
    somme1=0;
    for(i=0;i<notu+1;i++) somme1 += treer[i][j1];
    if (somme1==notu-1){
      br_a_virer=j1;
      if(lgbp){
        l=lgbi[j1];
        for(i=1;i<notu+1;i++) if(treer[i][j1]==0) {bp_a_modifier=i; break;}
      }
      break;
    }
    for(j2=j1+1;j2<notu-2;j2++){
      somme2=0;
      for(i=0;i<notu+1;i++) somme2+=treer[i][j2];
      if (somme1==notu-somme2){
	drapeau=0;
	for(i=1;i<notu;i++) if(abs(treer[i][j1]-treer[i][j2])!=1) drapeau=1;
	if (drapeau==0){
          br_a_virer=j1;
          if(lgbi) { bi_a_modifier=j2; l=lgbi[j1]; }
	  break;
	}
      }

      if (br_a_virer!=-1) break;
    }
    if (br_a_virer!=-1) break;
  }


  if(l1){
    j1=j2=0;
    for(i=1;i<=notu;i++){
      if(treer[i][br_a_virer]==0) {list1[j1]=nom[i]; j1++;}
      else {list2[j2]=nom[i]; j2++;}
    }
    if(j1<1 || j2<1) printf("probleme\n");
    *l1=j1; *l2=j2;
    if (bp_a_modifier!=-1) *fracroot1=lgbp[bp_a_modifier]/(lgbp[bp_a_modifier]+l);
    if (bi_a_modifier!=-1) *fracroot1=lgbi[bi_a_modifier]/(lgbi[bi_a_modifier]+l);
  }


  if (bp_a_modifier!=-1) lgbp[bp_a_modifier]+=l;
  if (bi_a_modifier!=-1) lgbi[bi_a_modifier]+=l;
  for(i=br_a_virer;i<notu-3;i++){
    if (lgbi) lgbi[i]=lgbi[i+1];
    if (bootvals) bootvals[i]=bootvals[i+1];
  }
  if(lgbp) for(i=0;i<notu;i++) lgbp[i]=lgbp[i+1];
  if(nom){
    free(nom[0]);
    for(i=0;i<notu;i++) nom[i]=nom[i+1];
    nom[notu]=check_alloc(MAXLNAME+1, sizeof(char));
  }

  for(i=1;i<notu+1;i++)
    for(j=0;j<notu-2;j++)
      treer[i-1][j]=treer[i][j];

  for(j=br_a_virer+1;j<notu-2;j++)
    for(i=0;i<notu;i++)
      treer[i][j-1]=treer[i][j];

  return 0;
}




void makelistbr_unrooted(noeud from, noeud nd, branche* br, int* alivebr,int nbbranch){
  noeud to1, to2;
  double lto1, lto2;
  int i=0, alive1, alive2;



  if (nd==root){
    if(from==root->v1) makelistbr_unrooted(root, nd->v2, br, alivebr, nbbranch);
    if(from==root->v2) makelistbr_unrooted(root, nd->v1, br, alivebr, nbbranch);
    return;
  }


  if (from==nd->v1) {to1=nd->v2; to2=nd->v3; lto1=nd->l2; lto2=nd->l3; alive1=nd->alive2; alive2=nd->alive3;}
  else if (from==nd->v2) {to1=nd->v1; to2=nd->v3; lto1=nd->l1; lto2=nd->l3; alive1=nd->alive1; alive2=nd->alive3;}
  else if (from==nd->v3) {to1=nd->v1; to2=nd->v2; lto1=nd->l1; lto2=nd->l2; alive1=nd->alive1; alive2=nd->alive2;}
  else {
    printf("erreur making branch list.\n");
  }


  while(i<nbbranch && br[i]->bout1) i++;

  if (to1 && to1==root) {
    if(nd==root->v1){
      lto1+=root->l2;
      br[i]->bout1=nd; br[i]->bout2=root->v2; br[i]->length=lto1; i++;
    }
    else{
      lto1+=root->l1;
      br[i]->bout1=nd; br[i]->bout2=root->v1; br[i]->length=lto1; i++;
    }
    if(alivebr) alivebr[i-1]=1;
  }


  if (to2 && to2==root) {
    if(nd==root->v1){
      lto2+=root->l2;
      br[i]->bout1=nd; br[i]->bout2=root->v2; br[i]->length=lto2; i++;
    }
    else{
      lto2+=root->l1;
      br[i]->bout1=nd; br[i]->bout2=root->v1; br[i]->length=lto2; i++;
    }
    if(alivebr) alivebr[i-1]=1;
  }



  if (to1 && to1!=root) {
    br[i]->bout1=nd;
    br[i]->bout2=to1;
    br[i]->length=lto1;
    if(alivebr) alivebr[i]=alive1;
    i++;
  }

  if (to2 && to2!=root) {
    br[i]->bout1=nd;
    br[i]->bout2=to2;
    br[i]->length=lto2;
    if(alivebr) alivebr[i]=alive2;
  }

  if (to1) makelistbr_unrooted(nd, to1, br, alivebr, nbbranch);
  if (to2) makelistbr_unrooted(nd, to2, br, alivebr, nbbranch);
}




/* create_node */
/* Return a node whose "parent", "childs" and values are */
/* v1, v2, v3, l1, l2, l3, nom, order. */

noeud create_node(noeud v1, noeud v2, noeud v3, double l1, double l2, double l3, double b1, double b2, double b3, char* nom) {

  noeud nd;

  nd = (struct noeud*) check_alloc(1, sizeof(struct noeud));
  nd->v1 = v1; nd->v2 = v2; nd->v3 = v3;
  nd->l1 = l1; nd->l2 = l2; nd->l3 = l3;
  nd->alive1 = (int)b1; nd->alive2 = (int)b2; nd->alive3 = (int)b3;

  if (nom!=NULL) {strncpy(nd->nom, nom, MAXLNAME); nd->nom[MAXLNAME]='\0'; }
  return nd;
}



/* bottomnode */
/* Return the root node of rooted s_tree including node nd. */
/*    !!BEWARE!!    */
/* If the s_tree including node nd is unrooted, program will bug */
/* (infinite loop). */

noeud bottomnode(noeud nd){
  if(nd->v3 == NULL) return nd;
  return (bottomnode(nd->v3));
}


/* ttos */
/* TableTOStruct. Create s_tree arbre_s from : t_tree arbre_t, leaves number notu, */
/* branch lengths lgbi (internal branches) , lgbp (terminal branches). */
/* If tree has no branch length, arguments 3 and 4 must be NULL. */


int ttos(int** arbre_t, int notu, double* lgbi, double* lgbp, double* alive, char** nom, noeud* arbre_s){

  noeud   	nd, p1, p2, p3;
  int           notuv, i, j, k, sommebi, tax1, tax2, n1, n2, n3, cpt1=-1;
  int           *kill_tax, *kill_bi;
  double	arg6, arg9;

  kill_tax = (int *)check_alloc(notu, sizeof(int));
  kill_bi = (int *)check_alloc(notu, sizeof(int));

  notuv = notu;
  for (i = 0; i <notu; i++){
    kill_tax[i] = 0;
    kill_bi[i] = 0;
  }

	/* terminal nodes */

  for (i = 0; i < notu; i++){
    if(lgbp==NULL) arg6=-1.0; else arg6=lgbp[i];
    arbre_s[i] = create_node(NULL, NULL, NULL, 0., 0., arg6, -1., -1., 1., nom[i]);
  }


	/* internal nodes */

  for (i = 0; i < notu - 3; i++) {
	/* determination de la bi a creer */
    for (j = 0; j < notu - 3 ; j++) {
      if (kill_bi[j] == 0) {
	sommebi = 0;
	for (k = 0; k < notu; k++)
	  if (kill_tax[k] == 0) {
	    sommebi += arbre_t[k][j];
	  }
	if (sommebi == 2 || notuv - sommebi == 2)
	  break;
      }
    }

	/* determination des 2 otus/noeuds fils */
    if (sommebi == 2) {
      for (k = 0; k < notu; k++)
	if (arbre_t[k][j] == 1 && kill_tax[k] == 0) {
	  tax1 = k;
	  break;
	}
      for (k = tax1 + 1; k < notu; k++)
	if (arbre_t[k][j] == 1 && kill_tax[k] == 0) {
	  tax2 = k;
	  break;
	}
    } else {
      for (k = 0; k < notu; k++)
	if (arbre_t[k][j] == 0 && kill_tax[k] == 0) {
	  tax1 = k;
	  break;
	}
      for (k = tax1 + 1; k < notu; k++)
	if (arbre_t[k][j] == 0 && kill_tax[k] == 0) {
	  tax2 = k;
	  break;
	}
    }

    p1 = bottomnode(arbre_s[tax1]);
    p2 = bottomnode(arbre_s[tax2]);

    if (lgbp) arg6=lgbi[j]; else arg6=0.;
    if (alive) arg9=alive[j]; else arg9=-1.;

    arbre_s[notu + (++cpt1)] = create_node(p1, p2, NULL, p1->l3, p2->l3, arg6, (double)p1->alive3, (double)p2->alive3, arg9, NULL);


    p1->v3 = arbre_s[notu + cpt1];
    p2->v3 = arbre_s[notu + cpt1];

    kill_tax[tax1] = kill_bi[j] = 1;
    notuv--;
  }


	/* last node */


  for (i = 0; i < 2 * notu - 3; i++)
    if (arbre_s[i]->v3 == NULL) {
      n1 = i;
      break;
    }

  for (i = n1 + 1; i < 2 * notu - 3; i++)
    if (arbre_s[i]->v3 == NULL) {
      n2 = i;
      break;
    }

  for (i = n2 + 1; i < 2 * notu - 3; i++)
    if (arbre_s[i]->v3 == NULL) {
      n3 = i;
      break;
    }


  arbre_s[2 * notu - 3] = create_node(arbre_s[n1], arbre_s[n2], arbre_s[n3], arbre_s[n1]->l3, arbre_s[n2]->l3, arbre_s[n3]->l3, arbre_s[n1]->alive3, arbre_s[n2]->alive3, arbre_s[n3]->alive3, NULL);
  arbre_s[n1]->v3 = arbre_s[2 * notu - 3];
  arbre_s[n2]->v3 = arbre_s[2 * notu - 3];
  arbre_s[n3]->v3 = arbre_s[2 * notu - 3];


  free(kill_tax);
  free(kill_bi);
  return 0;
}



/* root_s */

void root_s(branche br, double frac1){
  noeud nd1, nd2;
  double l;

  nd1=br->bout1;
  nd2=br->bout2;
  l=br->length;

  root->v1=nd1; root->v2=nd2; root->l1=l*frac1; root->l2=l*(1-frac1);
  root->alive1=root->alive2=1; root->alive3=-1;
  nd1->remain=nd2; nd2->remain=nd1;
  nd1->lremain=l; nd2->lremain=l;
  if (nd1->v1==nd2) { nd1->v1=root; nd1->l1=l*frac1; nd1->alremain=nd1->alive1; nd1->alive1=1;}
  else if (nd1->v2==nd2) { nd1->v2=root; nd1->l2=l*frac1; nd1->alremain=nd1->alive2; nd1->alive2=1; }
  else if (nd1->v3==nd2) { nd1->v3=root; nd1->l3=l*frac1; nd1->alremain=nd1->alive3; nd1->alive3=1; }
  else printf("Error rooting.\n");
  if (nd2->v1==nd1) { nd2->v1=root; nd2->l1=l*(1-frac1); nd2->alremain=nd2->alive1; nd2->alive1=1; }
  else if (nd2->v2==nd1) { nd2->v2=root; nd2->l2=l*(1-frac1); nd2->alremain=nd2->alive2; nd2->alive2=1; }
  else if (nd2->v3==nd1) { nd2->v3=root; nd2->l3=l*(1-frac1); nd2->alremain=nd2->alive3; nd2->alive3=1; }
  else printf("Error rooting.\n");
}




/* ctos */

noeud* ctos(char* ctree, int* notu, int *tree_is_rooted){

  int i, nb, **ttree, l1, l2, nbbi;
  double *lgbp, *lgbi, *bootvals, fracroot1;
  char **nom, racine, *list1[MAXNSP], *list2[MAXNSP];
  noeud* stree;
  branche* br_list;

  nb=taxa_count(ctree);

  lgbp=(double*)check_alloc(nb+1, sizeof(double));
  lgbi=(double*)check_alloc(nb+1, sizeof(double));
  bootvals=(double*)check_alloc(nb+1, sizeof(double));
  /*list1=check_alloc(2*nb+2, sizeof(char*));
  list2=check_alloc(2*nb+2, sizeof(char*));*/

  nom=(char**)check_alloc(nb+1, sizeof(char*));
  ttree=(int**)check_alloc(nb+1, sizeof(int*));
  for(i=0;i<nb+1;i++){
    nom[i]=(char*)check_alloc(MAXLNAME+1, sizeof(char));
    ttree[i]=(int*)check_alloc(nb, sizeof(int));
  }
  stree=(noeud*)check_alloc(2*nb, sizeof(noeud));
  br_list=(branche*)check_alloc(2*nb, sizeof(branche));
  for(i=0;i<2*nb-2;i++){
    br_list[i]=(struct branche*)check_alloc(1, sizeof(struct branche));
  }

  nb=ctot(ctree, ttree, lgbi, lgbp, bootvals, nom, &racine, &nbbi);

  if(racine=='r'){
    unroot(ttree, nb, lgbi, lgbp, bootvals, nom, list1, list2, &l1, &l2, &fracroot1);
  }
  ttos(ttree, nb, lgbi, lgbp, bootvals, nom, stree);
  makelistbr_unrooted(NULL, stree[0], br_list, NULL, 2*nb-3);
  root=create_node(NULL, NULL, NULL, -1., -1., -1., -1., -1., -1., "ROOT");

  if(racine=='r'){
    *tree_is_rooted=1;
    for(i=0;i<2*nb-3;i++){
      if(isbranch(br_list[i], l1, list1, l2, list2)==1){
 	root_s(br_list[i], fracroot1);
	break;
      }
      if(isbranch(br_list[i], l1, list1, l2, list2)==2){
 	root_s(br_list[i], 1.-fracroot1);
	break;
      }
    }
    if(i==2*nb-3) printf("Erreur racinage\n");
  }
  else{
    *tree_is_rooted=0;
    root_s(br_list[0], 0.5);
  }
  organize_tree(NULL, root);
  *notu=nb;

  free(lgbp);
  free(lgbi);
  free(bootvals);

  for(i=0;i<nb;i++){
    free(ttree[i]); free(nom[i]);
  }
  free(ttree); free(nom);
  /* free(list1); free(list2); */

  for(i=0;i<2*nb-2;i++)
    free(br_list[i]);
  free(br_list);

  return stree;

}


/* stree_retriction */

int stree_restriction(noeud nd, char** sel, int nbsel){

  int ret1, ret2;
  noeud keep, *me_in_parent;
  double newlg;
  double *my_lg_in_parent;

  if(nd->v1==NULL){
    if(nd->to_remove) return 1;
    if(is_in_list(nd->nom, sel, nbsel)) return 0;
    return 1;
  }
  
  ret1=stree_restriction(nd->v1, sel, nbsel);
  ret2=stree_restriction(nd->v2, sel, nbsel);
  
  if(ret1==0 && ret2==0) return 0;  //to keep
  if(ret1==1 && ret2==1) return 1;  //to remove
  
  if(ret1==0) keep=nd->v1; else keep=nd->v2;
  if(nd==root) {root=keep; root->v3=NULL; root->l3=-1; return -1;}
  
  newlg=nd->l3+keep->l3;
  keep->v3=nd->v3;
  keep->l3=newlg;
  if(nd->v3->v1==nd) {me_in_parent=&(nd->v3->v1); my_lg_in_parent=&(nd->v3->l1);}
  else  {me_in_parent=&(nd->v3->v2); my_lg_in_parent=&(nd->v3->l2);}
  *me_in_parent=keep;
  *my_lg_in_parent=newlg;
  return 0;
}



/* stoc */
char* stoc(noeud nd){
  char* ctree, *c1, *c2;
  
  if(nd->v1==0) return nd->nom;
  c1=stoc(nd->v1);
  c2=stoc(nd->v2);
  ctree=(char*)check_alloc(strlen(c1)+strlen(c2)+5+50, sizeof(char));
  sprintf(ctree, "(%s:%f,%s:%f)", c1, nd->l1, c2, nd->l2);
  if(nd->v3==NULL){
    char* prov;
    prov=ctree;
    while(*prov) prov++;
    *prov=';';
    *(prov+1)=0;
  }
  return ctree;
}



/* is_in_tree */

int is_in_tree(noeud nd, char* tax){

  int is1, is2;

  if(nd->v1==NULL)
    if(strcmp(nd->nom, tax)==0) return 1; else return 0;
    
  is1=is_in_tree(nd->v1, tax);
  if(is1) return 1;
  is2=is_in_tree(nd->v2, tax);
  if(is2) return 1;
  
  return 0;

}


/* is_monophyletic_rooted */
int is_monophyletic_rooted(int** ttree, int nb, int nbbi, char** names, char** taxa, int nbtax){

  int i, j, k, n, in;
  
  if(nbtax==1) return 1;

/*
for(i=0;i<nbtax;i++) printf("%s\n", taxa[i]);
for(j=0;j<=nb;j++){
  for(i=0;i<nbbi;i++){
    printf("%d", ttree[j][i]);
  }
  printf(" %s\n", names[j]);
}
*/

  for(i=0;i<nbbi;i++){
    n=0;
    for(j=1;j<nb+1;j++) n+=ttree[j][i];
    if(n!=nbtax) continue;
    for(j=1;j<nb+1;j++){
      in=is_in_list(names[j], taxa, nbtax);
      if(in==1 && ttree[j][i]==0) break;
      if(in==0 && ttree[j][i]==1) break;
    }
/*
printf("branch %d:%d\n", i, j);
getchar();    
*/
    if(j==nb+1) return 1;
  }
  
  return 0;
}



/* is_root_mrca_rooted */
/* return 1 if the common ancestor to the two groups of taxa taxa1 and taxa2 is the tre root, 0 if not */
/* inefficient version in o(n2) but called once per run so that's ok*/
int is_root_mrca_rooted(int** ttree, int nb, int nbbi, char** names, char** taxa1, int nbtax1, char** taxa2, int nbtax2){

  int i, j, k, n, n1;
  
  n=nbtax1+nbtax2;
  
  for(i=0;i<nbbi;i++){
    n1=0;
    for(j=0;j<nbtax1;j++){
      for(k=1;k<=nb;k++){
        if(strcmp(names[k],taxa1[j])==0){
          n1+=ttree[k][i];
          break;
        }
      }
    }
    for(j=0;j<nbtax2;j++){
      for(k=1;k<=nb;k++){
        if(strcmp(names[k],taxa2[j])==0){
          n1+=ttree[k][i];
          break;
        }
      }
    }
    if(n1==n) return 0;
  }

  return 1;
}



/* get_triplet_topol */
int get_triplet_topol(int** ttree, int nb, int nbbi, char** names, char** triplet){

  int i, j, n;
  int rka, rkb, rkc;
  
  //find triplet species pos
  for(j=1;j<nb+1;j++){
    if(strcmp(names[j], triplet[0])==0) rka=j;
    if(strcmp(names[j], triplet[1])==0) rkb=j;
    if(strcmp(names[j], triplet[2])==0) rkc=j;
  }
    
  //find relevant branch and get topol
  for(i=0;i<nbbi;i++){
    n=0;
    for(j=1;j<nb+1;j++) n+=ttree[j][i];
    if(n!=2) continue;
    if(ttree[rka][i]==1 && ttree[rkb][i]==1) return 0;
    if(ttree[rka][i]==1 && ttree[rkc][i]==1) return 1; 
    if(ttree[rkb][i]==1 && ttree[rkc][i]==1) return 2;         
  }
  
  return -1;
}


/* get_triplet_bl */
void get_triplet_bl(int** ttree, int nb, int nbbi, char** names, double* lgbi, double* lgbp, char** triplet, int topol, struct gene_data* gd){


  int i, j, n;
  int rka, rkb, rkc;
  
  //find triplet species pos
  for(j=1;j<nb+1;j++){
    if(strcmp(names[j], triplet[0])==0) rka=j;
    if(strcmp(names[j], triplet[1])==0) rkb=j;
    if(strcmp(names[j], triplet[2])==0) rkc=j;
  }

  //terminal bl
  gd->a=lgbp[rka];
  gd->b=lgbp[rkb];
  gd->c=lgbp[rkc];

  //d branch
  for(i=0;i<nbbi;i++){
    n=0;
    for(j=1;j<nb+1;j++) n+=ttree[j][i];
    if(n!=2) continue;
    if(ttree[rka][i]==1 && ttree[rkb][i]==1 && topol==0) {gd->d=lgbi[i]; break;}
    if(ttree[rka][i]==1 && ttree[rkc][i]==1 && topol==1) {gd->d=lgbi[i]; break;}
    if(ttree[rkb][i]==1 && ttree[rkc][i]==1 && topol==2) {gd->d=lgbi[i]; break;}
  }

  //e branch
  for(i=0;i<nbbi;i++){
    n=0;
    for(j=1;j<nb+1;j++) n+=ttree[j][i];
    if(n!=3) continue;
    if(ttree[rka][i]==1 && ttree[rkb][i]==1 && ttree[rkc][i]==1) {gd->e=lgbi[i]; break;}
  }
  
}



/* get_triplet_topol_branch_lengths */
void get_triplet_topol_branch_lengths(char* ctree, struct relevant_taxa* tax, struct option* opt, struct gene_data* gd){

  int i, nb, nbtax, is_rooted, nbsp, nbbi;
  int noutg, noth;
  int **ttree;
  double* lgbi, *lgbp, *bootvals;
  char** all_tax, **triplet_other, **outgroup, **names, rac;
  noeud* stree;
  
  gd->to_ignore=0;

  //count taxa in input tree
  nb=taxa_count(ctree);

  lgbp=(double*)check_alloc(nb+1, sizeof(double));
  lgbi=(double*)check_alloc(nb+1, sizeof(double));
  bootvals=(double*)check_alloc(nb+1, sizeof(double));

  names=(char**)check_alloc(nb+1, sizeof(char*));
  ttree=(int**)check_alloc(nb+1, sizeof(int*));
  for(i=0;i<nb+1;i++){
    names[i]=(char*)check_alloc(MAXLNAME+1, sizeof(char));
    ttree[i]=(int*)check_alloc(nb, sizeof(int));
  }

  //convert to struct
  stree=ctos(ctree, &nbsp, &is_rooted);
  
  //check rooted
  if(!is_rooted) {printf("rooted tree expected"); exit(0);}    
  
  //list all relevant taxa
  nbtax=3+tax->nb_outg+tax->nb_oth;
  all_tax=(char**)check_alloc(nbtax, sizeof(char*));
  for(i=0;i<3;i++) all_tax[i]=tax->triplet[i];
  for(i=0;i<tax->nb_outg;i++) all_tax[3+i]=tax->outg[i];
  for(i=0;i<tax->nb_oth;i++) all_tax[3+tax->nb_outg+i]=tax->oth[i];
  triplet_other=(char**)check_alloc(3+tax->nb_oth, sizeof(char*));
  for(i=0;i<3;i++) triplet_other[i]=tax->triplet[i];
  outgroup=(char**)check_alloc(tax->nb_outg, sizeof(char*));
  
  //restrict tree
  /* restrict */  
  stree_restriction(root, all_tax, nbtax);

  //check all triplet species
  for(i=0;i<3;i++){
    if(is_in_tree(root, tax->triplet[i])==0){
      gd->to_ignore=1;
      return;
    }
  } 

  //count+identify other
  noth=0;
  for(i=0;i<tax->nb_oth;i++){
    if(is_in_tree(root, tax->oth[i])){
      triplet_other[3+noth]=tax->oth[i];
      noth++;
    }
  }
  
  //count+identify outgroup  
  noutg=0;
  for(i=0;i<tax->nb_outg;i++){
    if(is_in_tree(root, tax->outg[i])){
      outgroup[noutg]=tax->outg[i];
      noutg++;
    }
  }
  
  //check >=1 outgroup species
  if(noutg==0){
    gd->to_ignore=2; 
    return;
  }
    
  //convert to table
  nb=ctot(stoc(root), ttree, lgbi, lgbp, bootvals, names, &rac, &nbbi);  
  
  //check triplet monophyly  
  if(is_monophyletic_rooted(ttree, nb, nbbi, names, tax->triplet, 3)==0){
    gd->to_ignore=3;
    return;    
  }

  //check MRCA(triplet, outg)==root
  if(is_root_mrca_rooted(ttree, nb, nbbi, names, tax->triplet, 3, outgroup, noutg)==0){
    gd->to_ignore=4;
    return;      
  }

  //check triplet+other monophyly (optional)
  if(opt->require_triplet_other_monophyly && is_monophyletic_rooted(ttree, nb, nbbi, names, triplet_other, 3+noth)==0){
    gd->to_ignore=5;
    return;    
  }

  //check outgroup monophyly (optional)  
  if(opt->require_outgroup_monophyly && is_monophyletic_rooted(ttree, nb, nbbi, names, outgroup, noutg)==0){
    gd->to_ignore=6;
    return;    
  }    

  //get triplet topology
  gd->topology=get_triplet_topol(ttree, nb, nbbi, names, tax->triplet);
  
  //get triplet branch lengths
  get_triplet_bl(ttree, nb, nbbi, names, lgbi, lgbp, tax->triplet, gd->topology, gd);

}



/* calculate_root_to_tip_taxa */
void calculate_root_to_tip_taxa(noeud nd, char** tax_list, int ntax){

  int n1, n2;
  double avl1, avl2;
  
  if(nd->v1==NULL) { 
    if(is_in_list(nd->nom, tax_list, ntax)){ nd->n_to_tip=1; nd->av_l_to_tip=0.;}
    else {nd->n_to_tip=0; nd->av_l_to_tip=-1.; }
    return;
  }
  
  calculate_root_to_tip_taxa(nd->v1, tax_list, ntax);
  calculate_root_to_tip_taxa(nd->v2, tax_list, ntax);
  
  n1=nd->v1->n_to_tip;
  n2=nd->v2->n_to_tip;
  avl1=nd->v1->av_l_to_tip;
  avl2=nd->v2->av_l_to_tip;
      
  if(n1+n2) nd->av_l_to_tip=((nd->l1+avl1)*n1+(nd->l2+avl2)*n2)/(n1+n2);
  else nd->av_l_to_tip=-1.;
  nd->n_to_tip=n1+n2;
}





/* calculate_root_to_tip */
void calculate_root_to_tip(noeud nd){

  int n1, n2;
  double avl1, avl2;
  
  if(nd->v1==NULL) {nd->av_l_to_tip=0.; nd->n_to_tip=1; return;}
  
  calculate_root_to_tip(nd->v1);
  calculate_root_to_tip(nd->v2);
  
  n1=nd->v1->n_to_tip;
  n2=nd->v2->n_to_tip;
  avl1=nd->v1->av_l_to_tip;
  avl2=nd->v2->av_l_to_tip;
      
  nd->n_to_tip=n1+n2;
  nd->av_l_to_tip=((nd->l1+avl1)*n1+(nd->l2+avl2)*n2)/(n1+n2);
}



/* read_data_file */

struct gene_data* read_data_file(FILE* inf, struct option* opt, struct relevant_taxa* tax, int* nbg){ 

  struct gene_data* gdata;
  int i, nb_gene, iret, lline, maxlline, lgt;
  char c, *line, *ret, *tree, *lgs, *name, **non_triplet;
  double all_rtt, triplet_rtt, non_triplet_rtt;

  //first pass: max line length
  maxlline=lline=0;
  while(1){
    c=fgetc(inf);
    if(c==EOF) break;
    if(c=='\n'){
      if(lline>maxlline) maxlline=lline;
      lline=0;
    }
    else lline++;
  }
  if(lline>maxlline) maxlline=lline;


  //allocate 1
  line=check_alloc(maxlline+2, sizeof(char));
  
  //second pass: count gene trees
  rewind(inf);
  nb_gene=0;
  while(1){
    ret=fgets(line, maxlline+1, inf);
    if(ret==NULL) break;
    if(line[0]=='\n' || line[0]==' ' || line[0]=='\t' || line[0]=='#')
      continue;
    nb_gene++;
  }  
  
  //allocate 2
  gdata=(struct gene_data*)check_alloc(nb_gene, sizeof(struct gene_data));
    
  //third pass: read gene trees + seq lg + name, check, store relevant branch lengths
  rewind(inf);
  non_triplet=(char**)check_alloc(tax->nb_outg+tax->nb_oth, sizeof(char*));
  for(i=0;i<tax->nb_outg;i++) non_triplet[i]=tax->outg[i];
  for(i=0;i<tax->nb_oth;i++) non_triplet[tax->nb_outg+i]=tax->oth[i];
  i=0;
  while(1){
    ret=fgets(line, maxlline+1, inf);
    if(ret==NULL) break;
    if(line[0]=='\n' || line[0]==' ' || line[0]=='\t' || line[0]=='#')
      continue;
    tree=strtok(line, "\t\n");
    
    lgs=strtok(NULL, "\t\n");
    if(lgs==NULL) fin("unexpected line in data file (1)");
    iret=sscanf(lgs, "%d", &(gdata[i].lgseq));
    if(iret!=1) fin("unexpected line in data file (2)");
    
    name=strtok(NULL, "\t\n");
    if(name==NULL) fin("unexpected line in data file (3)"); 
    gdata[i].name=(char*)check_alloc(strlen(name)+1, sizeof(char)); 
    sprintf(gdata[i].name, "%s", name);
      
    get_triplet_topol_branch_lengths(tree, tax, opt, gdata+i);

    calculate_root_to_tip(root);
    all_rtt=root->av_l_to_tip;

    calculate_root_to_tip_taxa(root, tax->triplet, 3);
    triplet_rtt=root->av_l_to_tip;

    calculate_root_to_tip_taxa(root, non_triplet, tax->nb_outg+tax->nb_oth);
    non_triplet_rtt=root->av_l_to_tip;
        
    gdata[i].l=all_rtt;
    gdata[i].triplet_l=triplet_rtt;
    gdata[i].non_triplet_l=non_triplet_rtt;
        
    i++;
  }    

  *nbg=nb_gene;  
  return gdata;
}





/* set_log_facto_data */
void set_log_facto_data(struct gene_data* gdata, int nbg){
  int max, i;

  max=0;
  for(i=0;i<nbg;i++)
    if(gdata[i].lgseq>max) max=gdata[i].lgseq;
  max++; 

  log_facto=(double*)check_alloc(max+1, sizeof(double));
    
  log_facto[0]=log_facto[1]=0;
  for(i=2;i<=max;i++)
    log_facto[i]=log_facto[i-1]+log(i);
}  



/* set_log_facto */
void set_log_facto(int max){
  int i;

  log_facto=(double*)check_alloc(max+1, sizeof(double));
    
  log_facto[0]=log_facto[1]=0;
  for(i=2;i<=max;i++)
    log_facto[i]=log_facto[i-1]+log(i);
}  




/* set_topol_scenario */

void  set_topol_scenario(struct option opt){

  int i, j;
  
  for(i=0;i<3;i++) for(j=0;j<nb_scenario_plus1;j++) topol_scenario[i][j]=0;
  
  //basic
  topol_scenario[0][1]=1;

  //ILS
  for(i=0;i<opt.nb_ILS_times;i++){
    topol_scenario[0][2+3*i]=1; 
    topol_scenario[1][3+3*i]=1;
    topol_scenario[2][4+3*i]=1; 
  }
  
  //HGT
  for(i=0;i<opt.nb_HGT_times;i++){
    topol_scenario[0][3*opt.nb_ILS_times+2+5*i]=1;
    topol_scenario[1][3*opt.nb_ILS_times+3+5*i]=1;    
    topol_scenario[1][3*opt.nb_ILS_times+4+5*i]=1;    
    topol_scenario[2][3*opt.nb_ILS_times+5+5*i]=1;    
    topol_scenario[2][3*opt.nb_ILS_times+6+5*i]=1;   
  }
  
  //unresolved topology
  for(j=1;j<nb_scenario_plus1;j++) topol_scenario[3][j]=1;  
}


/* print_data */

void print_data(struct gene_data* gd){

  printf("topol:%d\n", gd->topology);
  printf("bl: %f %f %f %f %f\n", gd->a, gd->b, gd->c, gd->d, gd->l);
  printf("lgseq:%d\n", gd->lgseq);
  printf("\n");
}

/* print_par */

void print_par(struct param* par, struct confidence_interval* ci, int nbci){

  if(nbci==0){
    printf("tau1: %f, tau2: %f, theta: %f\n", par->tau1, par->tau2, par->theta);
    printf("Pac: %f, Pbc: %f, Pab:%f\n", par->pac, par->pbc, par->pab);
    printf("Pa: %f\n", par->pone);
  }
  
  if(nbci==1){
    printf("tau1: %f, tau2: %f, theta: %f [%f-%f]\n", par->tau1, par->tau2, par->theta, ci[0].master_param_min, ci[0].master_param_max);
    printf("Pac: %f, Pbc: %f, Pab:%f\n", par->pac, par->pbc, par->pab);
    printf("Pa: %f\n", par->pone);  
  }

  if(nbci==nb_parameters){
    printf("tau1: %f [%f-%f], tau2: %f [%f-%f], theta: %f [%f-%f]\n", par->tau1, ci[0].master_param_min, ci[0].master_param_max, par->tau2, ci[1].master_param_min, ci[1].master_param_max, par->theta, ci[2].master_param_min, ci[2].master_param_max);
    printf("Pab: %f [%f-%f], Pac: %f [%f-%f], Pbc:%f [%f-%f]\n", par->pab, ci[3].master_param_min, ci[3].master_param_max, par->pac, ci[4].master_param_min, ci[4].master_param_max, par->pbc, ci[5].master_param_min, ci[5].master_param_max);
    printf("Pa: %f [%f-%f]\n", par->pone, ci[6].master_param_min, ci[6].master_param_max);  
  }

}


/* print_lh */

void print_lh(struct log_likelihood* lh){

  printf("lnL:%f\n", lh->lnL);
  printf("1st order derivatives: [t1]%f  [t2]%f  [th]%f  [pab]%f  [pac]%f  [pbc]%f  [py]%f\n", lh->dlnLdt1, lh->dlnLdt2, lh->dlnLdth, lh->dlnLdpab, lh->dlnLdpac, lh->dlnLdpbc, lh->dlnLdpone);
  printf("2nd order derivatives: [t1]%f  [t2]%f  [th]%f  [pab]%f  [pac]%f  [pbc]%f  [py]%f\n", lh->d2lnLdt1, lh->d2lnLdt2, lh->d2lnLdth, lh->d2lnLdpab, lh->d2lnLdpac, lh->d2lnLdpbc, lh->d2lnLdpone);
  printf("\n");
}



/* apply_constraints */
void apply_constraints(struct param* par){
  double rat;
  
  //times
  if(par->tau1<QuasiZero) par->tau1=QuasiZero;
  if(par->tau2<par->tau1) par->tau2=par->tau1+QuasiZero;
  if(par->theta<QuasiZero) par->theta=QuasiZero;
      
  //PHGT
  if(par->pab<QuasiZero) par->pab=QuasiZero;  
  if(par->pac<QuasiZero) par->pac=QuasiZero;  
  if(par->pbc<QuasiZero) par->pbc=QuasiZero;      
  if(par->pab+par->pac+par->pbc>MAX_PHGT){
    rat=(par->pab+par->pac+par->pbc)/MAX_PHGT;
    par->pab/=rat; par->pac/=rat; par->pbc/=rat; 
  }
  if(par->pone<QuasiZero) par->pone=QuasiZero; 
  if(par->pone>1-QuasiZero) par->pone=1-QuasiZero;

}



/* set_starting */
void set_starting(struct gene_data* gdata, int nbg, struct option opt, double theta, struct param* par){

  int i;
  double mean_a, mean_b, mean_c, mean_d, mean_l;
  double pils;
  double prop_topol[NTOPOL];

  // theta: argument
  par->theta=theta;
  
  //constants
  par->nbg=nbg;
  par->l_ref=0.;
  for(i=0;i<nbg;i++) par->l_ref+=gdata[i].l;
  par->l_ref/=nbg;

  //calculate topo proportions, not counting unresolved
  for(i=0;i<NTOPOL;i++) prop_topol[i]=0;
  for(i=0;i<nbg;i++) prop_topol[gdata[i].topology]++;
  for(i=0;i<3;i++) prop_topol[i]/=(nbg-prop_topol[3]);
    
  //branch lengths: based on mean observed 
  mean_a=mean_b=mean_c=mean_d=0.;
  for(i=0;i<nbg;i++){
    mean_a+=gdata[i].a; mean_b+=gdata[i].b; mean_c+=gdata[i].c; mean_d+=gdata[i].d; 
  }
  mean_a/=nbg; mean_b/=nbg; mean_c/=nbg; mean_d/=nbg; 
  par->tau1=(mean_a+mean_b)/2.-theta/2.;
  if(opt.optimize_tau1==0) par->tau1=opt.fixed_tau1;
  par->tau2=(mean_c+par->tau1+mean_d)/2. - theta/2.; 
  if(opt.optimize_tau1==0) par->tau2=opt.fixed_tau2;
  if(par->tau2<par->tau1) par->tau2=par->tau1+QuasiZero;

  //alpha[i]: observed l[i] / mean observed l
  for(i=0;i<nbg;i++) {par->alpha[i]=gdata[i].l/par->l_ref; if(par->alpha[i]==0) par->alpha[i]=QuasiZero;}

  //starting pac pbc : anomalous topology freq time, having substracted probable ILS ones
  pils=exp(-2*(par->tau2-par->tau1)/par->theta);
  par->pac=prop_topol[1]-pils/3.; if(par->pac<QuasiZero) par->pac=QuasiZero; 
  par->pbc=prop_topol[2]-pils/3.; if(par->pbc<QuasiZero) par->pbc=QuasiZero; 
  //par->pac=prop_topol[1]/2.;
  if(opt.optimize_pac==0) par->pac=opt.fixed_pac;
  //par->pbc=prop_topol[2]/2.;  
  if(opt.optimize_pbc==0) par->pbc=opt.fixed_pbc;
  
  //starting pab: mean of pac and pbc
  par->pab=(par->pac+par->pbc)/2.;
  if(opt.optimize_pab==0) par->pab=opt.fixed_pab;
  
  //starting pone: 1/nb HGT times
  par->pone=1./opt.nb_HGT_times;
  if(opt.optimize_pone==0) par->pone=opt.fixed_pone;

}



/* calculate_scenario_likelihood */
void calculate_scenario_likelihood(struct gene_data* gd, struct param* par, int gene, int scenar, struct scenario_likelihood* slh){

  int i, obse[NBRANCH];
  double expe[NBRANCH];
  double logL, dlogLdt1, dlogLdt2, dlogLdth, d2logLdt1, d2logLdt2, d2logLdth;
  double lal;
//struct scenario_likelihood slh;

  obse[0]=round(gd->lgseq*gd->a); 
  obse[1]=round(gd->lgseq*gd->b); 
  obse[2]=round(gd->lgseq*gd->c); 
  obse[3]=round(gd->lgseq*gd->d); 
  //obse[4]=round(gd->lgseq*gd->l);
  
  lal=gd->lgseq*par->alpha[gene];
  for(i=0;i<NBRANCH;i++) expe[i]=lal*exp_lg[scenar][i];

//for(i=0;i<NBRANCH;i++) printf("%d obs=%d exp=%f\n", scenar, obse[i], expe[i]);
//getchar();

  logL=dlogLdt1=dlogLdt2=dlogLdth=d2logLdt1=d2logLdt2=d2logLdth=0.;
  for(i=0;i<NBRANCH;i++){
    logL+=(obse[i]*log(expe[i])-expe[i]-log_facto[obse[i]]);    
    dlogLdt1+=(obse[i]*lal*exp_lg_coeff[scenar][i][0]/expe[i] - lal*exp_lg_coeff[scenar][i][0]);
    dlogLdt2+=(obse[i]*lal*exp_lg_coeff[scenar][i][1]/expe[i] - lal*exp_lg_coeff[scenar][i][1]);
    dlogLdth+=(obse[i]*lal*exp_lg_coeff[scenar][i][2]/expe[i] - lal*exp_lg_coeff[scenar][i][2]);
    d2logLdt1-=obse[i]*lal*exp_lg_coeff[scenar][i][0]*lal*exp_lg_coeff[scenar][i][0]/(expe[i]*expe[i]);
    d2logLdt2-=obse[i]*lal*exp_lg_coeff[scenar][i][1]*lal*exp_lg_coeff[scenar][i][1]/(expe[i]*expe[i]);
    d2logLdth-=obse[i]*lal*exp_lg_coeff[scenar][i][2]*lal*exp_lg_coeff[scenar][i][2]/(expe[i]*expe[i]);
  }
  
  slh->L=exp(logL);
  slh->dLdt1=dlogLdt1*slh->L;
  slh->dLdt2=dlogLdt2*slh->L;
  slh->dLdth=dlogLdth*slh->L;
  slh->d2Ldt1=(d2logLdt1*slh->L*slh->L+slh->dLdt1*slh->dLdt1)/slh->L;
  slh->d2Ldt2=(d2logLdt2*slh->L*slh->L+slh->dLdt2*slh->dLdt2)/slh->L;
  slh->d2Ldth=(d2logLdth*slh->L*slh->L+slh->dLdth*slh->dLdth)/slh->L; 
  if(slh->L==0.) slh->d2Ldt1=slh->d2Ldt2=slh->d2Ldth=0.;
//to calculte the derivatives of L we first calculate the derivatives of logL then use dL=dlogL*L
//and similarly for 2nd derivatives: d2L=(d2logL*L*L+dL*dL)/L=d2logL*L+(dL*dL)/L

  return;
}



/* free_gene_lkh_content */
void free_gene_lkh_content(struct gene_likelihood* glh, int post_too){

  if(post_too) free(glh->posterior_scenar);
  free(glh->scenar_p);
  free(glh->scenar_l);
}



/*
struct gene_likelihood{
  double L;
  double *scenar_p, *scenar_l;
  double *posterior_scenar, posterior_ILS, posterior_HGT;
  int predicted_scenario, predicted_scenario_type;
    
  double dLdt1, dLdt2, dLdth;
  double dLdpab, dLdpac, dLdpbc, dLdpone;
  double d2Ldt1, d2Ldt2, d2Ldth;
  double d2Ldpab, d2Ldpac, d2Ldpbc, d2Ldpone;  
};
*/

/* calculate_gene_likelihood */
/* !! uses global variables p, pr, exp_lg !! */

void calculate_gene_likelihood(struct gene_data* gd, struct param* par, struct option* opt, int gene, struct gene_likelihood* glh){

  //struct gene_likelihood glh;
  struct scenario_likelihood slh;
  double L;
  double pils;
  double dpilsdt1, dpilsdt2, dpilsdth, d2pilsdt1, d2pilsdt2, d2pilsdth;
  double t1, t2, th;
  double pab, pac, pbc, pone;
  int i;

    
  t1=par->tau1; t2=par->tau2; th=par->theta; 
  pab=par->pab; pac=par->pac; pbc=par->pbc; pone=par->pone;
  pils=exp(-2*(t2-t1)/th);
  
//printf("pils=%f\n", pils);
  
  //below we note L=sum_i p[i] pr[i] , where p[i] is the probability of a scenario, pr[i] is the likelihood of a scenario, 1<=i<=nb scenario
  
  //scenario probabilities
  p[1]=(1-pab-pac-pbc)*(1.-pils);
  
  for(i=2;i<3*opt->nb_ILS_times+2;i++)
    p[i]=(1-pab-pac-pbc)*pils/(3*opt->nb_ILS_times);
  
  p[3*opt->nb_ILS_times+2]=pab*pone; p[3*opt->nb_ILS_times+3]=p[3*opt->nb_ILS_times+4]=pac*pone/2.; p[3*opt->nb_ILS_times+5]=p[3*opt->nb_ILS_times+6]=pbc*pone/2.;

  for(i=1;i<opt->nb_HGT_times;i++){
    p[3*opt->nb_ILS_times+2+5*i]=pab*(1-pone)/(opt->nb_HGT_times-1);    
    p[3*opt->nb_ILS_times+3+5*i]=p[3*opt->nb_ILS_times+4+5*i]=pac*(1-pone)/(2.*(opt->nb_HGT_times-1));
    p[3*opt->nb_ILS_times+5+5*i]=p[3*opt->nb_ILS_times+6+5*i]=pbc*(1-pone)/(2.*(opt->nb_HGT_times-1));  
  }

//for(i=1;i<nb_scenario_plus1;i++){
//  printf("scenario proba: %f\n", p[i]);
//}

  //scenario likelihoods and scenario likelihood derivatives wrt t1, t2, th
  for(i=1;i<nb_scenario_plus1;i++){
    pr[i]=dprdt1[i]=dprdt2[i]=dprdth[i]=d2prdt1[i]=d2prdt2[i]=d2prdth[i]=0.;
    if(topol_scenario[gd->topology][i]){
      calculate_scenario_likelihood(gd, par, gene, i, &slh);
      pr[i]=slh.L;
      dprdt1[i]=slh.dLdt1; dprdt2[i]=slh.dLdt2; dprdth[i]=slh.dLdth;
      d2prdt1[i]=slh.d2Ldt1; d2prdt2[i]=slh.d2Ldt2; d2prdth[i]=slh.d2Ldth;
    }
  }

  //store scenario proba and likelihood
//  glh->scenar_p=(double*)check_alloc(nb_scenario_plus1, sizeof(double));
//  glh->scenar_l=(double*)check_alloc(nb_scenario_plus1, sizeof(double));
  for(i=1;i<nb_scenario_plus1;i++){
    glh->scenar_p[i]=p[i];
    glh->scenar_l[i]=pr[i];
  }

  //gene likelihood
  glh->L=0.;
  for(i=1;i<nb_scenario_plus1;i++){
    glh->L+=p[i]*pr[i];
  }

  //gene likelihood derivatives wrt pab pac pbc pone
  glh->dLdpab=-pr[1]*(1.-pils);
  for(i=0;i<3*opt->nb_ILS_times;i++) glh->dLdpab-=(pr[2+i]*pils/(3*opt->nb_ILS_times));
  glh->dLdpab+=pr[2+3*opt->nb_ILS_times]*pone;
  for(i=1;i<opt->nb_HGT_times;i++) glh->dLdpab+=pr[2+3*opt->nb_ILS_times+5*i]*(1.-pone)/(opt->nb_HGT_times-1);
  glh->dLdpac=-pr[1]*(1.-pils);
  for(i=0;i<3*opt->nb_ILS_times;i++) glh->dLdpac-=(pr[2+i]*pils/(3*opt->nb_ILS_times));
  glh->dLdpac+=(pr[3+3*opt->nb_ILS_times]+pr[4+3*opt->nb_ILS_times])*pone/2.; 
  for(i=1;i<opt->nb_HGT_times;i++) glh->dLdpac+=(pr[3+3*opt->nb_ILS_times+5*i]+pr[4+3*opt->nb_ILS_times+5*i])*(1.-pone)/(2*(opt->nb_HGT_times-1));  
  glh->dLdpbc=-pr[1]*(1.-pils);
  for(i=0;i<3*opt->nb_ILS_times;i++) glh->dLdpbc-=(pr[2+i]*pils/(3*opt->nb_ILS_times));
  glh->dLdpbc+=(pr[5+3*opt->nb_ILS_times]+pr[6+3*opt->nb_ILS_times])*pone/2.;  
  for(i=1;i<opt->nb_HGT_times;i++) glh->dLdpbc+=(pr[5+3*opt->nb_ILS_times+5*i]+pr[6+3*opt->nb_ILS_times+5*i])*(1.-pone)/(2*(opt->nb_HGT_times-1));  

  glh->dLdpone=pr[2+3*opt->nb_ILS_times]*pab+(pr[3+3*opt->nb_ILS_times]+pr[4+3*opt->nb_ILS_times])*pac/2.+(pr[5+3*opt->nb_ILS_times]+pr[6+3*opt->nb_ILS_times])*pbc/2.;
  for(i=1;i<opt->nb_HGT_times;i++) glh->dLdpone-=((pr[2+3*opt->nb_ILS_times+5*i]*pab)+(pr[3+3*opt->nb_ILS_times+5*i]+pr[4+3*opt->nb_ILS_times+5*i])*pac/2.+(pr[5+3*opt->nb_ILS_times+5*i]+pr[6+3*opt->nb_ILS_times+5*i])*pbc/2.)/(opt->nb_HGT_times-1);
  
  glh->d2Ldpab=glh->d2Ldpac=glh->d2Ldpbc=glh->d2Ldpone=0;
  
  //gene likelihood derivatives wrt t1, t2, th (harder)

    //pils derivatives
  dpilsdt1=2*pils/th;
  dpilsdt2=-2*pils/th;
  dpilsdth=2*(t2-t1)*pils/(th*th);
  d2pilsdt1=4*pils/(th*th);
  d2pilsdt2=4*pils/(th*th);
  d2pilsdth=2*(t2-t1)*(dpilsdth*th*th-2*th*pils)/(th*th*th*th);


    //scenario proba derivatives
  dpdt1[1]=-(1-pab-pac-pbc)*dpilsdt1;
  dpdt2[1]=-(1-pab-pac-pbc)*dpilsdt2;  
  dpdth[1]=-(1-pab-pac-pbc)*dpilsdth;
  for(i=0;i<opt->nb_ILS_times;i++){
    dpdt1[2+3*i]=dpdt1[3+3*i]=dpdt1[4+3*i]=(1-pab-pac-pbc)*dpilsdt1/(3*opt->nb_ILS_times);  
    dpdt2[2+3*i]=dpdt2[3+3*i]=dpdt2[4+3*i]=(1-pab-pac-pbc)*dpilsdt2/(3*opt->nb_ILS_times);  
    dpdth[2+3*i]=dpdth[3+3*i]=dpdth[4+3*i]=(1-pab-pac-pbc)*dpilsdth/(3*opt->nb_ILS_times);
  }
  for(i=3*opt->nb_ILS_times+2;i<nb_scenario_plus1;i++) dpdt1[i]=dpdt2[i]=dpdth[i]=0.;

  d2pdt1[1]=-(1-pab-pac-pbc)*d2pilsdt1;
  d2pdt2[1]=-(1-pab-pac-pbc)*d2pilsdt2;  
  d2pdth[1]=-(1-pab-pac-pbc)*d2pilsdth;
  for(i=0;i<opt->nb_ILS_times;i++){
    d2pdt1[2+3*i]=d2pdt1[3+3*i]=d2pdt1[4+3*i]=(1-pab-pac-pbc)*d2pilsdt1/(3*opt->nb_ILS_times);  
    d2pdt2[2+3*i]=d2pdt2[3+3*i]=d2pdt2[4+3*i]=(1-pab-pac-pbc)*d2pilsdt2/(3*opt->nb_ILS_times);  
    d2pdth[2+3*i]=d2pdth[3+3*i]=d2pdth[4+3*i]=(1-pab-pac-pbc)*d2pilsdth/(3*opt->nb_ILS_times);
  }
  for(i=3*opt->nb_ILS_times+2;i<nb_scenario_plus1;i++) d2pdt1[i]=d2pdt2[i]=d2pdth[i]=0.;  
        
    //t1
  glh->dLdt1=glh->d2Ldt1=0.;
  for(i=1;i<nb_scenario_plus1;i++){
    glh->dLdt1+=p[i]*dprdt1[i]+pr[i]*dpdt1[i];
    glh->d2Ldt1+=p[i]*d2prdt1[i]+pr[i]*d2pdt1[i]+2*dpdt1[i]*dprdt1[i];
  }
    
    //t2
  glh->dLdt2=glh->d2Ldt2=0.;
  for(i=1;i<nb_scenario_plus1;i++){
    glh->dLdt2+=p[i]*dprdt2[i]+pr[i]*dpdt2[i];
    glh->d2Ldt2+=p[i]*d2prdt2[i]+pr[i]*d2pdt2[i]+2*dpdt2[i]*dprdt2[i];
  }
    
    //th  
  glh->dLdth=glh->d2Ldth=0.;
  for(i=1;i<nb_scenario_plus1;i++){
    glh->dLdth+=p[i]*dprdth[i]+pr[i]*dpdth[i];
    glh->d2Ldth+=p[i]*d2prdth[i]+pr[i]*d2pdth[i]+2*dpdth[i]*dprdth[i];
  }
  
      
  //return glh;  

}


/* set_expected_length */
/* !! uses global variables exp_lg !! */
/* sets expected branch lengths a b c and d for each scenaro given a parameter set */
/* also sets the linear coefficients relating expected lengths to parameters (needed for derivative calculation) */

void set_expected_length(struct param* par, struct option* opt){
  double t1, t2, tg, th, ths2, ths6;
  double coeff;
  int i, k;
  
  t1=par->tau1; t2=par->tau2; th=par->theta; ths2=th/2; ths6=th/6;
  
  //scenario 1: basic
  exp_lg[1][0]=t1+ths2;  //a   
  exp_lg_coeff[1][0][0]=1; exp_lg_coeff[1][0][1]=0; exp_lg_coeff[1][0][2]=0.5; 
  exp_lg[1][1]=t1+ths2;  //b
  exp_lg_coeff[1][1][0]=1; exp_lg_coeff[1][1][1]=0; exp_lg_coeff[1][1][2]=0.5;   
  exp_lg[1][2]=t2+ths2;  //c
  exp_lg_coeff[1][2][0]=0; exp_lg_coeff[1][2][1]=1; exp_lg_coeff[1][2][2]=0.5;   
  exp_lg[1][3]=t2-t1;  //d
  exp_lg_coeff[1][3][0]=-1; exp_lg_coeff[1][3][1]=1; exp_lg_coeff[1][3][2]=0;   
  
  if(opt->strict_no_event && t1+ths2>t2){
    exp_lg[1][0]=t2;  //a   
    exp_lg_coeff[1][0][0]=0; exp_lg_coeff[1][0][1]=1; exp_lg_coeff[1][0][2]=0; 
    exp_lg[1][1]=t2;  //b
    exp_lg_coeff[1][1][0]=0; exp_lg_coeff[1][1][1]=1; exp_lg_coeff[1][1][2]=0;   
    exp_lg[1][2]=t2+ths2;  //c
    exp_lg_coeff[1][2][0]=0; exp_lg_coeff[1][2][1]=1; exp_lg_coeff[1][2][2]=0.5;   
    exp_lg[1][3]=ths2;  //d
    exp_lg_coeff[1][3][0]=0; exp_lg_coeff[1][3][1]=0; exp_lg_coeff[1][3][2]=0.5;  
  }


  for(i=0;i<opt->nb_ILS_times;i++){ 
    //scenario 2 modulo 3: ILS topology c
    k=2+3*i;
    exp_lg[k][0]=t2+ths6;  //a
    exp_lg_coeff[k][0][0]=0; exp_lg_coeff[k][0][1]=1; exp_lg_coeff[k][0][2]=1./6.;   
    exp_lg[k][1]=t2+ths6;  //b
    exp_lg_coeff[k][1][0]=0; exp_lg_coeff[k][1][1]=1; exp_lg_coeff[k][1][2]=1./6.;   
    exp_lg[k][2]=t2+ths6+ths2*opt->discr_expo_mean1[i];  //c
    exp_lg_coeff[k][2][0]=0; exp_lg_coeff[k][2][1]=1; exp_lg_coeff[k][2][2]=1./6.+0.5*opt->discr_expo_mean1[i];   
    exp_lg[k][3]=ths2*opt->discr_expo_mean1[i];  //d
    exp_lg_coeff[k][3][0]=0; exp_lg_coeff[k][3][1]=0; exp_lg_coeff[k][3][2]=0.5*opt->discr_expo_mean1[i];   
  
    //scenario 3 modulo 3: ILS topology b
    k=3+3*i;
    exp_lg[k][0]=t2+ths6;  //a
    exp_lg_coeff[k][0][0]=0; exp_lg_coeff[k][0][1]=1; exp_lg_coeff[k][0][2]=1./6.;     
    exp_lg[k][1]=t2+ths6+ths2*opt->discr_expo_mean1[i];  //b
    exp_lg_coeff[k][1][0]=0; exp_lg_coeff[k][1][1]=1; exp_lg_coeff[k][1][2]=1./6.+0.5*opt->discr_expo_mean1[i];     
    exp_lg[k][2]=t2+ths6;  //c
    exp_lg_coeff[k][2][0]=0; exp_lg_coeff[k][2][1]=1; exp_lg_coeff[k][2][2]=1./6.;       
    exp_lg[k][3]=ths2*opt->discr_expo_mean1[i];  //d
    exp_lg_coeff[k][3][0]=0; exp_lg_coeff[k][3][1]=0; exp_lg_coeff[k][3][2]=0.5*opt->discr_expo_mean1[i];     

    //scenario 4 modulo 3: ILS topology a
    k=4+3*i;
    exp_lg[k][0]=t2+ths6+ths2*opt->discr_expo_mean1[i];  //a
    exp_lg_coeff[k][0][0]=0; exp_lg_coeff[k][0][1]=1; exp_lg_coeff[k][0][2]=1./6.+0.5*opt->discr_expo_mean1[i];     
    exp_lg[k][1]=t2+ths6;  //b
    exp_lg_coeff[k][1][0]=0; exp_lg_coeff[k][1][1]=1; exp_lg_coeff[k][1][2]=1./6.;     
    exp_lg[k][2]=t2+ths6;  //c
    exp_lg_coeff[k][2][0]=0; exp_lg_coeff[k][2][1]=1; exp_lg_coeff[k][2][2]=1./6.;     
    exp_lg[k][3]=ths2*opt->discr_expo_mean1[i];  //d
    exp_lg_coeff[k][3][0]=0; exp_lg_coeff[k][3][1]=0; exp_lg_coeff[k][3][2]=0.5*opt->discr_expo_mean1[i];     
  }


  for(i=0;i<opt->nb_HGT_times;i++){

    coeff=opt->HGT_time_coeff[i];
    tg=par->tau1*coeff; 
        
    //scenario 5 modulo 5: HGT A<->B
    k=2+3*opt->nb_ILS_times+5*i;
    exp_lg[k][0]=tg;  //a
    exp_lg_coeff[k][0][0]=coeff; exp_lg_coeff[k][0][1]=0; exp_lg_coeff[k][0][2]=0;     
    exp_lg[k][1]=tg;  //b
    exp_lg_coeff[k][1][0]=coeff; exp_lg_coeff[k][1][1]=0; exp_lg_coeff[k][1][2]=0;     
    exp_lg[k][2]=t2+ths2;  //c
    exp_lg_coeff[k][2][0]=0; exp_lg_coeff[k][2][1]=1; exp_lg_coeff[k][2][2]=0.5;     
    exp_lg[k][3]=t2+ths2-tg;  //d
    exp_lg_coeff[k][3][0]=-coeff; exp_lg_coeff[k][3][1]=1; exp_lg_coeff[k][3][2]=0.5;     

  //scenario 6 modulo 5: HGT A->C
    k=3+3*opt->nb_ILS_times+5*i;
    exp_lg[k][0]=tg;  //a
    exp_lg_coeff[k][0][0]=coeff; exp_lg_coeff[k][0][1]=0; exp_lg_coeff[k][0][2]=0;     
    exp_lg[k][1]=t1+ths2;  //b
    exp_lg_coeff[k][1][0]=1; exp_lg_coeff[k][1][1]=0; exp_lg_coeff[k][1][2]=0.5;     
    exp_lg[k][2]=tg;  //c
    exp_lg_coeff[k][2][0]=coeff; exp_lg_coeff[k][2][1]=0; exp_lg_coeff[k][2][2]=0;     
    exp_lg[k][3]=t1+ths2-tg;  //d
    exp_lg_coeff[k][3][0]=1-coeff; exp_lg_coeff[k][3][1]=0; exp_lg_coeff[k][3][2]=0.5;     
  
  //scenario 7 modulo 5: HGT C->A
    k=4+3*opt->nb_ILS_times+5*i;
    exp_lg[k][0]=tg;  //a
    exp_lg_coeff[k][0][0]=coeff; exp_lg_coeff[k][0][1]=0; exp_lg_coeff[k][0][2]=0;     
    exp_lg[k][1]=t2+ths2;  //b
    exp_lg_coeff[k][1][0]=0; exp_lg_coeff[k][1][1]=1; exp_lg_coeff[k][1][2]=0.5;     
    exp_lg[k][2]=tg;  //c
    exp_lg_coeff[k][2][0]=coeff; exp_lg_coeff[k][2][1]=0; exp_lg_coeff[k][2][2]=0;     
    exp_lg[k][3]=t2+ths2-tg;  //d
    exp_lg_coeff[k][3][0]=-coeff; exp_lg_coeff[k][3][1]=1; exp_lg_coeff[k][3][2]=0.5;     

  //scenario 8 modulo 5: HGT B->C
    k=5+3*opt->nb_ILS_times+5*i;
    exp_lg[k][0]=t1+ths2;  //a
    exp_lg_coeff[k][0][0]=1; exp_lg_coeff[k][0][1]=0; exp_lg_coeff[k][0][2]=0.5;     
    exp_lg[k][1]=tg;  //b
    exp_lg_coeff[k][1][0]=coeff; exp_lg_coeff[k][1][1]=0; exp_lg_coeff[k][1][2]=0;     
    exp_lg[k][2]=tg;  //c
    exp_lg_coeff[k][2][0]=coeff; exp_lg_coeff[k][2][1]=0; exp_lg_coeff[k][2][2]=0;     
    exp_lg[k][3]=t1+ths2-tg;  //d
    exp_lg_coeff[k][3][0]=1-coeff; exp_lg_coeff[k][3][1]=0; exp_lg_coeff[k][3][2]=0.5;     

  //scenario 9 modulo 5: HGT C->B
    k=6+3*opt->nb_ILS_times+5*i;
    exp_lg[k][0]=t2+ths2;  //a
    exp_lg_coeff[k][0][0]=0; exp_lg_coeff[k][0][1]=1; exp_lg_coeff[k][0][2]=0.5;     
    exp_lg[k][1]=tg;  //b
    exp_lg_coeff[k][1][0]=coeff; exp_lg_coeff[k][1][1]=0; exp_lg_coeff[k][1][2]=0;     
    exp_lg[k][2]=tg;  //c
    exp_lg_coeff[k][2][0]=coeff; exp_lg_coeff[k][2][1]=0; exp_lg_coeff[k][2][2]=0;     
    exp_lg[k][3]=t2+ths2-tg;  //d
    exp_lg_coeff[k][3][0]=-coeff; exp_lg_coeff[k][3][1]=1; exp_lg_coeff[k][3][2]=0.5;   
  }  
}



/* calculate_log_likelihood */

struct log_likelihood calculate_log_likelihood(struct gene_likelihood* glh, int nbg){

  struct log_likelihood lh;
  int i;
  double gL, gL2;
  
  //lnL
  lh.lnL=0.;
  for(i=0;i<nbg;i++) lh.lnL+=log(glh[i].L);
  
  //1st order derivatives
  lh.dlnLdt1=lh.dlnLdt2=lh.dlnLdth=lh.dlnLdpab=lh.dlnLdpac=lh.dlnLdpbc=lh.dlnLdpone=0;

  for(i=0;i<nbg;i++){
    gL=glh[i].L; 
    lh.dlnLdt1+=glh[i].dLdt1/gL;
    lh.dlnLdt2+=glh[i].dLdt2/gL;
    lh.dlnLdth+=glh[i].dLdth/gL;
    lh.dlnLdpab+=glh[i].dLdpab/gL;
    lh.dlnLdpac+=glh[i].dLdpac/gL;
    lh.dlnLdpbc+=glh[i].dLdpbc/gL;
    lh.dlnLdpone+=glh[i].dLdpone/gL;
    if(gL==0.) lh.dlnLdt1=lh.dlnLdt2=lh.dlnLdth=lh.dlnLdpab=lh.dlnLdpac=lh.dlnLdpbc=lh.dlnLdpone=0.;
  }
    
  //2nd order derivatives
  lh.d2lnLdt1=lh.d2lnLdt2=lh.d2lnLdth=lh.d2lnLdpab=lh.d2lnLdpac=lh.d2lnLdpbc=lh.d2lnLdpone=0;    

  for(i=0;i<nbg;i++){
    gL2=glh[i].L*glh[i].L;
    lh.d2lnLdt1+=(glh[i].d2Ldt1*glh[i].L-glh[i].dLdt1*glh[i].dLdt1)/gL2; 
    lh.d2lnLdt2+=(glh[i].d2Ldt2*glh[i].L-glh[i].dLdt2*glh[i].dLdt2)/gL2; 
    lh.d2lnLdth+=(glh[i].d2Ldth*glh[i].L-glh[i].dLdth*glh[i].dLdth)/gL2;     
    lh.d2lnLdpab+=-(glh[i].dLdpab*glh[i].dLdpab)/gL2; 
    lh.d2lnLdpac+=-(glh[i].dLdpac*glh[i].dLdpac)/gL2;
    lh.d2lnLdpbc+=-(glh[i].dLdpbc*glh[i].dLdpbc)/gL2;    
    lh.d2lnLdpone+=-(glh[i].dLdpone*glh[i].dLdpone)/gL2;    
    if(gL2==0.) lh.d2lnLdt1=lh.d2lnLdt2=lh.d2lnLdth=lh.d2lnLdpab=lh.d2lnLdpac=lh.d2lnLdpbc=lh.d2lnLdpone=0.;
  }        

  return lh;
}



/* copy_param */
void copy_param(struct param* from, struct param* to){
  int i;
  to->tau1=from->tau1;
  to->tau2=from->tau2;
  to->theta=from->theta;
  to->pab=from->pab;
  to->pac=from->pac;
  to->pbc=from->pbc;  
  to->pone=from->pone;  
  to->nbg=from->nbg;
  to->l_ref=from->l_ref;
  to->alpha=from->alpha;
}



/* optimize_lnL_starting */
struct gene_likelihood* optimize_lnL_starting(struct gene_data* gdata, int nbg, struct option opt, struct param starting_par, struct param* arg_par, double* maxlnL){

  int i, nbdec;
  struct log_likelihood lh, lhp, lhm;
  struct gene_likelihood* glh;
  struct param par, mem_par;
//struct param true_par;
  double gL, gL2;
  double lnL, new_lnL;
  double dt1, dt2, dth, dpab, dpac, dpbc, dpone;

  //initialize
  //set_log_facto_data(gdata, nbg);  
  copy_param(&starting_par, &par);
  //par.alpha=starting_par.alpha;
  //par.nbg=starting_par.nbg;
  //par.l_ref=starting_par.l_ref;
  


  //allocate
  glh=(struct gene_likelihood*)check_alloc(nbg, sizeof(struct gene_likelihood));
  for(i=0;i<nbg;i++){
    glh[i].scenar_p=(double*)check_alloc(nb_scenario_plus1, sizeof(double));
    glh[i].scenar_l=(double*)check_alloc(nb_scenario_plus1, sizeof(double));  
    glh[i].posterior_scenar=(double*)check_alloc(nb_scenario_plus1, sizeof(double));
  }


//printf("Starting params:\n");  
//print_par(&par, NULL, 0);

  // fill in gene likelihood and derivatives   
  set_expected_length(&par, &opt);


  for(i=0;i<nbg;i++){
    //glh[i]=calculate_gene_likelihood(gdata+i, &par, &opt, i);    
    calculate_gene_likelihood(gdata+i, &par, &opt, i, glh+i);    
  }


  // calculate starting L + derivatives 
  lh=calculate_log_likelihood(glh, nbg);

  if(opt.verbose) {printf("starting lnL: %f -> ", lh.lnL); fflush(stdout);}  
    

  // Newton Raphson loop 
  it=0;
  nbdec=0;
  
  while(it<opt.max_iterations){
    if(nbdec==0) lnL=lh.lnL;
    if(nbdec==opt.max_backward_moves){
      if(opt.verbose) printf("Optimization problem (too many backward moves), maximum lnL:%f\n", new_lnL);
      break; 
    }

    copy_param(&par, &mem_par);

//print_par(&par, NULL, 0);

    dt1=dt2=dth=dpab=dpac=dpbc=dpone=0.;
    if(opt.optimize_tau1) dt1=-(lh.dlnLdt1/lh.d2lnLdt1)/pow(2, nbdec); if(lh.d2lnLdt1==0.) dt1=0.;
    if(opt.optimize_tau2) dt2=-(lh.dlnLdt2/lh.d2lnLdt2)/pow(2, nbdec); if(lh.d2lnLdt2==0.) dt2=0.;
    if(opt.optimize_theta) dth=-(lh.dlnLdth/lh.d2lnLdth)/pow(2, nbdec); if(lh.d2lnLdth==0.) dth=0.;
    if(opt.optimize_pab) dpab=-(lh.dlnLdpab/lh.d2lnLdpab)/pow(2, nbdec); if(lh.d2lnLdpab==0.) dpab=0.;
    if(opt.optimize_pac) dpac=-(lh.dlnLdpac/lh.d2lnLdpac)/pow(2, nbdec); if(lh.d2lnLdpac==0.) dpac=0.;
    if(opt.optimize_pbc) dpbc=-(lh.dlnLdpbc/lh.d2lnLdpbc)/pow(2, nbdec); if(lh.d2lnLdpbc==0.) dpbc=0.;
    if(opt.optimize_pone) dpone=-(lh.dlnLdpone/lh.d2lnLdpone)/pow(2, nbdec); if(lh.d2lnLdpone==0.) dpone=0.;
    
    par.tau1+=dt1; par.tau2+=dt2; par.theta+=dth; 
    par.pab+=dpab; par.pac+=dpac; par.pbc+=dpbc; 
    par.pone+=dpone; 

    apply_constraints(&par);

    set_expected_length(&par, &opt);
    for(i=0;i<nbg;i++){
      calculate_gene_likelihood(gdata+i, &par, &opt, i, glh+i);
    }
    lh=calculate_log_likelihood(glh, nbg);    
    new_lnL=lh.lnL;
    
    if(new_lnL<lnL) { //decreasing -lnL: parameters move back
      copy_param(&mem_par, &par);
      nbdec++; 
      continue;
    } else nbdec=0;

    it++;
        
    if(new_lnL-lnL<opt.convergence){
      if(opt.verbose) printf("Maximum lnL:%f\n", new_lnL);
      break;
    }
    
    if(it>opt.max_iterations){
      if(opt.verbose) printf("Optimization problem (too many iterationss), maximum lnL:%f\n", new_lnL);
      break;    
    }

//if(it%100==0) printf("it %d\n", it);

  }


  //return;
  
  *maxlnL=new_lnL;
  copy_param(&par, arg_par);
  arg_par->alpha=(double*)check_alloc(nbg, sizeof(double));
  for(i=0;i<nbg;i++) arg_par->alpha[i]=par.alpha[i]; 
  return glh;  

}


/* optimize_lnL */
struct gene_likelihood* optimize_lnL(struct gene_data* gdata, int nbg, struct option opt, struct param* final_par, double* max_maxlnL){
  struct param start_par, arg_par;
  int i, j, best;
  double maxlnL;
  struct gene_likelihood** glh;
  
  glh=(struct gene_likelihood**)check_alloc(NB_THETA_START, sizeof(struct gene_likelihood*));
  start_par.alpha=(double*)check_alloc(nbg, sizeof(double));

  if(opt.optimize_theta){ //loop across starting values
    for(i=0;i<NB_THETA_START;i++){  
      set_starting(gdata, nbg, opt, theta_start[i], &start_par);
      apply_constraints(&start_par);  
      glh[i]=optimize_lnL_starting(gdata, nbg, opt, start_par, &arg_par, &maxlnL);
      if(i==0 || maxlnL>*max_maxlnL){
        *max_maxlnL=maxlnL;
        copy_param(&arg_par, final_par);
        best=i;
      }
    }
    
    for(i=0;i<NB_THETA_START;i++){  
      if(i!=best){
        for(j=0;j<nbg;j++)
          free_gene_lkh_content(glh[i]+j, 0);
          free(glh[i]);
      }
    }
  }
  else{ //single starting values
    set_starting(gdata, nbg, opt, opt.fixed_theta, &start_par);
    apply_constraints(&start_par);
    glh[0]=optimize_lnL_starting(gdata, nbg, opt, start_par, &arg_par, &maxlnL);
    *max_maxlnL=maxlnL;
    copy_param(&arg_par, final_par);
    best=0;
  }
  

  
  return glh[best];
}



/* annotate_genes */
void annotate_genes(struct gene_likelihood* glh, int nb_gene, struct option opt){

  int i, j, bests;
  double sum_post, maxps;

  for(i=0;i<nb_gene;i++){
    sum_post=0;
    for(j=1;j<nb_scenario_plus1;j++) {
      glh[i].posterior_scenar[j]=glh[i].scenar_p[j]*glh[i].scenar_l[j];
      sum_post+=glh[i].posterior_scenar[j];
    }
    for(j=1;j<nb_scenario_plus1;j++) glh[i].posterior_scenar[j]/=sum_post;
    glh[i].posterior_ILS=0.;
    for(j=0;j<opt.nb_ILS_times;j++) glh[i].posterior_ILS+=glh[i].posterior_scenar[2+3*j]+glh[i].posterior_scenar[3+3*j]+glh[i].posterior_scenar[4+3*j];
    glh[i].posterior_HGT=1.-glh[i].posterior_ILS-glh[i].posterior_scenar[1];
    
    glh[i].predicted_scenario=glh[i].predicted_scenario_type=-1;
    max_rank(glh[i].posterior_scenar+1, nb_scenario_plus1-1, &maxps, &bests);
    bests++;
    if(maxps>opt.posterior_threshold) glh[i].predicted_scenario=bests;
    if(glh[i].posterior_scenar[1]>opt.posterior_threshold) glh[i].predicted_scenario_type=0;
    if(glh[i].posterior_ILS>opt.posterior_threshold) glh[i].predicted_scenario_type=1;
    if(glh[i].posterior_HGT>opt.posterior_threshold) glh[i].predicted_scenario_type=2;
  }
}


/* annotate_global */
void annotate_global(struct gene_likelihood* glh, int nb_gene, struct posterior* post, struct option opt){

  int i, j;

  post->scenar=(double*)check_alloc(nb_scenario_plus1, sizeof(double));
  post->basic=post->noconflict_ILS=post->noconflict_HGT=post->conflict_ILS=post->conflict_HGT=0;

  for(j=1;j<nb_scenario_plus1;j++)
    post->scenar[j]=0.;
  for(i=0;i<nb_gene;i++){
    for(j=1;j<nb_scenario_plus1;j++)
      post->scenar[j]+=glh[i].posterior_scenar[j];
  }
  for(j=1;j<nb_scenario_plus1;j++)
    post->scenar[j]/=nb_gene;
       
  post->basic=post->scenar[1]; 
  post->noconflict_ILS=0.;
  for(j=0;j<opt.nb_ILS_times;j++) post->noconflict_ILS+=post->scenar[2+3*j]; 
  post->conflict_ILS=0;
  for(j=0;j<opt.nb_ILS_times;j++) post->conflict_ILS+=post->scenar[3+3*j]+post->scenar[4+3*j];
  post->noconflict_HGT=0; 
  for(j=0;j<opt.nb_HGT_times;j++) post->noconflict_HGT+=post->scenar[2+3*opt.nb_ILS_times+5*j]; 
  post->conflict_HGT=0;  
  for(j=0;j<opt.nb_HGT_times;j++) post->conflict_HGT+=post->scenar[3+3*opt.nb_ILS_times+5*j]+post->scenar[4+3*opt.nb_ILS_times+5*j]+post->scenar[5+3*opt.nb_ILS_times+5*j]+post->scenar[6+3*opt.nb_ILS_times+5*j];
  post->asymmetry_ILS=0;
  for(j=0;j<opt.nb_ILS_times;j++) post->asymmetry_ILS+=post->scenar[3];
  post->asymmetry_ILS/=post->conflict_ILS;
  if(post->asymmetry_ILS<0.5) {post->asymmetry_ILS=1.-post->asymmetry_ILS; post->dominant_ILS=2;}
  else post->dominant_ILS=1;
  post->asymmetry_HGT=0;
  for(j=0;j<opt.nb_HGT_times;j++) post->asymmetry_HGT+=post->scenar[3+3*opt.nb_ILS_times+5*j]+post->scenar[4+3*opt.nb_ILS_times+5*j];
  post->asymmetry_HGT/=post->conflict_HGT;
  if(post->asymmetry_HGT<0.5) {post->asymmetry_HGT=1.-post->asymmetry_HGT; post->dominant_HGT=2;}
  else post->dominant_HGT=1;

}


/* calculate_logLn_degenerate */
double calculate_logLn_degenerate(struct gene_data* gdata, int nbg){

  double lnLd, Ld_gene, lnLd_scenar;
  int i, j;
  int oa, ob, oc, od;
  double ea, eb, ec, ed;
  double oalogea, oblogeb, oclogec, odloged;
  
  lnLd=0;
  
  for(i=0;i<nbg;i++){ //for ech gene tree

    oa=(round)(gdata[i].lgseq*gdata[i].a);
    ob=(round)(gdata[i].lgseq*gdata[i].b);
    oc=(round)(gdata[i].lgseq*gdata[i].c);
    od=(round)(gdata[i].lgseq*gdata[i].d);
    Ld_gene=0.;
    
    for(j=0;j<nbg;j++){ //consider every (other or not) gene tree as a scenario, with expected = observed
      lnLd_scenar=0.;
      ea=gdata[j].lgseq*gdata[j].a;
      eb=gdata[j].lgseq*gdata[j].b;
      ec=gdata[j].lgseq*gdata[j].c;
      ed=gdata[j].lgseq*gdata[j].d;
      
      if(ea>0) oalogea=oa*log(ea); else oalogea=0;
      if(eb>0) oblogeb=ob*log(eb); else oblogeb=0;
      if(ec>0) oclogec=oc*log(ec); else oclogec=0;
      if(ed>0) odloged=od*log(ed); else odloged=0;
    
      lnLd_scenar+=oalogea-ea-log_facto[oa];
      lnLd_scenar+=oblogeb-eb-log_facto[ob];
      lnLd_scenar+=oclogec-ec-log_facto[oc];
      lnLd_scenar+=odloged-ed-log_facto[od];

      Ld_gene+=exp(lnLd_scenar)/nbg;    //each scenario has probability 1./nbg
    }        


    lnLd+=log(Ld_gene);  
  }
  
   
  return lnLd;

}


/* calculate_ci */
int calculate_ci(struct gene_data* gdata, int nb_gene, struct option opt, struct param optimal_par, double max_lnL, char* master_s, struct confidence_interval* ci){

  int i, j, verb;
  double* master; 
  double param_down, param_up, param_med, param_opt, param_min, param_max, lnL, diff;
  double tol=0.1;
  struct gene_likelihood* glh;
  struct param par, start_par;
  struct posterior post;


  verb=opt.verbose;
  
  //identify master parameter
  if(strncmp(master_s, "tau1", 4)==0){param_opt=optimal_par.tau1; param_min=QuasiZero; param_max=2*optimal_par.tau1; master=&(start_par.tau1); opt.optimize_tau1=0;}
  else if(strncmp(master_s, "tau2", 4)==0){param_opt=optimal_par.tau2; param_min=QuasiZero; param_max=2*optimal_par.tau2; master=&(start_par.tau2); opt.optimize_tau2=0;}
  else if(strncmp(master_s, "theta", 5)==0){param_opt=optimal_par.theta; param_min=QuasiZero; param_max=10*optimal_par.theta; master=&(start_par.theta); opt.optimize_theta=0;}
  else if(strncmp(master_s, "pab", 3)==0){param_opt=optimal_par.pab; param_min=QuasiZero; param_max=0.5; master=&(start_par.pab); opt.optimize_pab=0;}
  else if(strncmp(master_s, "pac", 3)==0){param_opt=optimal_par.pac; param_min=QuasiZero; param_max=0.33; master=&(start_par.pac); opt.optimize_pac=0;}
  else if(strncmp(master_s, "pbc", 3)==0){param_opt=optimal_par.pbc; param_min=QuasiZero; param_max=0.33; master=&(start_par.pbc); opt.optimize_pbc=0;}
  else if(strncmp(master_s, "pone", 4)==0){param_opt=optimal_par.pone; param_min=QuasiZero; param_max=1.-QuasiZero; master=&(start_par.pone); opt.optimize_pone=0;}
  else {if(verb) printf("unexpected parameter name in CI calculation\n"); return 0;}
  
  //initialize
  ci->max_lnL_diff=opt.CI_max_lnL_diff;  
  opt.verbose=0;
  copy_param(&optimal_par, &start_par);
  
  //1.lower bound
//printf("lower\n");
  param_down=param_min; param_up=param_opt; param_med=param_down;
  *master=param_med;
  apply_constraints(&start_par);
  glh=optimize_lnL_starting(gdata, nb_gene, opt, start_par, &par, &lnL);
  diff=max_lnL-lnL;
  if(diff>opt.CI_max_lnL_diff+tol){
    for(j=0;j<nb_gene;j++) free_gene_lkh_content(glh+j, 0);
    i=0;
    while(1){    
      param_med=(param_up+param_down)/2.;
      *master=param_med;
      apply_constraints(&start_par);
      glh=optimize_lnL_starting(gdata, nb_gene, opt, start_par, &par, &lnL);
      diff=max_lnL-lnL;
//printf("%f %f\n", *master, diff);
      if(diff>opt.CI_max_lnL_diff+tol) param_down=param_med; 
      else if(diff<opt.CI_max_lnL_diff-tol) param_up=param_med; 
      else break;   
      for(j=0;j<nb_gene;j++) free_gene_lkh_content(glh+j, 0);
      if(i>100) {if(verb) printf("problem caculating confidence interval\n"); return 0;}   
      i++;
    }    
  }
  ci->master_param_min=*master;
  copy_param(&par, &(ci->par_min));
  annotate_genes(glh, nb_gene, opt);
  annotate_global(glh, nb_gene, &(ci->post_min), opt);

  //2.upper bound
//printf("upper\n");
  param_down=param_opt; param_up=param_max; param_med=param_up;
  *master=param_med;
  apply_constraints(&start_par);
  glh=optimize_lnL_starting(gdata, nb_gene, opt, start_par, &par, &lnL);
  diff=max_lnL-lnL;
  if(diff>opt.CI_max_lnL_diff+tol){
    for(j=0;j<nb_gene;j++) free_gene_lkh_content(glh+j, 0);
    i=0;
    while(1){    
      param_med=(param_up+param_down)/2.;
      *master=param_med;
      apply_constraints(&start_par);
      glh=optimize_lnL_starting(gdata, nb_gene, opt, start_par, &par, &lnL);
      diff=max_lnL-lnL;
//printf("%f %f\n", *master, diff);
      if(diff>opt.CI_max_lnL_diff+tol) param_up=param_med;
      else if(diff<opt.CI_max_lnL_diff-tol) param_down=param_med; 
      else break;   
      for(j=0;j<nb_gene;j++) free_gene_lkh_content(glh+j, 0);
      if(i>100) {if(verb) printf("problem caculating confidence interval\n"); return 0;}   
      i++;
    } 
  }
  ci->master_param_max=*master;
  copy_param(&par, &(ci->par_max));
  annotate_genes(glh, nb_gene, opt);
  annotate_global(glh, nb_gene, &(ci->post_max), opt);
  
  return 1;
}



/* count_topol */
void  count_topol(struct gene_data* gdata, int nb_gene, int* ntopo0, int* ntopo1, int* ntopo2, int* ntopo3){

  int i;
 
  *ntopo1=*ntopo2=*ntopo3=*ntopo0=0;
  for(i=0;i<nb_gene;i++){
    if(gdata[i].topology==0) (*ntopo0)++;
    else if(gdata[i].topology==1) (*ntopo1)++;
    else if(gdata[i].topology==2) (*ntopo2)++;
    else if(gdata[i].topology==3) (*ntopo3)++;  
  }
}
  
  
/* av_branch_lengths */
void  av_branch_lengths(struct gene_data* gdata, int nb_gene, double* ava, double* avb, double* avc, double* avd){

  int i;
  
  *ava=*avb=*avc=*avd=0;  
  for(i=0;i<nb_gene;i++){
    *ava+=gdata[i].a;
    *avb+=gdata[i].b;
    *avc+=gdata[i].c;
    *avd+=gdata[i].d;
  }
  (*ava)/=nb_gene;
  (*avb)/=nb_gene;
  (*avc)/=nb_gene;
  (*avd)/=nb_gene;
}


/* ctopol */
char* ctopol(int to){

  if(to==0) return "((AB)C)";
  if(to==1) return "((AC)B)"; 
  if(to==2) return "((BC)A)"; 
  return "";
}


/* main */
int main(int argc, char** argv){

  struct option opt;
  struct relevant_taxa tax;
  struct gene_data* gdata;
  struct param optimal_par;
  struct posterior post;
  struct confidence_interval* ci;
  int nb_gene, nb_ignored, nb_sp, reti, nbok, nbci, i, j;
  int ntopo1, ntopo2, ntopo3, ntopo0;
  int gof_ddl;
  double max_lnL, lnL_degn;
  double opt_pab, opt_pac, opt_pbc, opt_pils;
  int npred_ils, npred_hgt, nils_a, nils_b, nhgt_a, nhgt_b;
  double av_lg, av_a, av_b, av_c, av_d, sum_post, thresh;
  struct gene_likelihood* glh;
  FILE* inf, *optf, *taxf, *outf;


  if(argc!=5) usage();

  //read option file

  optf=fopen(argv[3], "r");
  if(optf==NULL) {printf("Cannot open file: %s\n", argv[3]); exit(0);}  
  read_option_file(optf, &opt);  
  opt.discr_expo_mean1=discretize_exponential(1, opt.nb_ILS_times);
  set_log_facto((int)opt.max_bl+1);

  if(opt.verbose){  
    printf("Reading option file...\n");    
    print_options(&opt);
    printf("\n");
  }
  
  //read taxon file
  if(opt.verbose) printf("Reading taxon file...\n");
  taxf=fopen(argv[2], "r");
  if(taxf==NULL) {printf("Cannot open file: %s\n", argv[2]); exit(0);}  
  read_taxon_file(taxf, &tax); 
  if(opt.verbose){ 
    print_taxa(&tax); 
    printf("\n");
  }

  // read data file
  inf=fopen(argv[1], "r");
  if(inf==NULL) {printf("Cannot open file: %s\n", argv[1]); exit(0);}
  if(opt.verbose) printf("Reading data file...\n");
  gdata=read_data_file(inf, &opt, &tax, &nb_gene);
  if(opt.verbose) printf("found %d gene trees\n", nb_gene);
  set_unresolved_topology(gdata, nb_gene, opt.min_d);
  find_dubious_genes(gdata, nb_gene, &opt);
  if(opt.verbose) print_ignored(gdata, nb_gene);
  do_ignore(gdata, &nb_gene);
  
  // data summary
  count_topol(gdata, nb_gene, &ntopo0, &ntopo1, &ntopo2, &ntopo3);
  if(opt.verbose) printf("((AB)C):%d [%.1f%%]  ((AC)B):%d [%.1f%%]  ((BC)A):%d [%.1f%%]  (ABC):%d [%.1f%%]\n", ntopo0, 100.*(double)ntopo0/nb_gene, ntopo1, 100.*(double)ntopo1/nb_gene, ntopo2, 100.*(double)ntopo2/nb_gene, ntopo3, 100.*(double)ntopo3/nb_gene);
  av_lg=0.;
  for(i=0;i<nb_gene;i++) av_lg+=gdata[i].lgseq;
  av_lg/=nb_gene;
  if(opt.verbose) printf("average sequence length: %.1f\n", av_lg);
  av_branch_lengths(gdata, nb_gene, &av_a, &av_b, &av_c, &av_d);
  if(opt.verbose) printf("average branch lengths: a=%f, b=%f, c=%f, d=%f\n", av_a, av_b, av_c, av_d);
  if(opt.verbose) printf("\n");
  
  //initialize
  set_parameter_names();
  nb_scenario_plus1=1+3*opt.nb_ILS_times+5*opt.nb_HGT_times+1;
  allocate_global_var();
  set_scenario_names(&opt);
  set_topol_scenario(opt);
  
  // optimize lnL
  if(opt.verbose) printf("Optimizing likelihood...\n");  
  glh=optimize_lnL(gdata, nb_gene, opt, &optimal_par, &max_lnL);
  if(opt.verbose) printf("\n");

  // annotate
  if(opt.verbose) printf("Posterior annotation...\n");    
  
  annotate_genes(glh, nb_gene, opt);
  annotate_global(glh, nb_gene, &post, opt);

  if(opt.verbose) printf("\n");
  
  // gof (not nested models so null ditrib to be approached via simulations)
/*  
  lnL_degn=calculate_logLn_degenerate(gdata, nb_gene);
  gof_ddl=5*nb_gene-1-nb_parameters;
  if(opt.verbose) printf("maxlnL:%f, deglnL:%f, 2delta_lnL:%f, ddl:%d, normal_deviate:%f\n", max_lnL, lnL_degn, 2*(lnL_degn-max_lnL), gof_ddl, (2*(lnL_degn-max_lnL)-gof_ddl)/sqrt(2*gof_ddl));
  //getchar();
*/

  // confidence intervals
  ci=check_alloc(nb_parameters, sizeof(struct confidence_interval));
  for(i=0;i<nb_parameters;i++){
    ci[i].master_param_min=ci[i].master_param_max=-1.;
    ci[i].post_min.conflict_ILS=ci[i].post_max.conflict_ILS=ci[i].post_min.conflict_HGT=ci[i].post_max.conflict_HGT=-1.;
  }
  if(opt.do_param_CI){
    if(opt.verbose) printf("Confidence intervals:\n");     
    nbci=nb_parameters;
    for(i=0;i<nb_parameters;i++){
      if(opt.verbose) printf("%s...\n", parameter_name[i]);
      reti=calculate_ci(gdata, nb_gene, opt, optimal_par, max_lnL, parameter_name[i], ci+i);
      //if(reti==0) {nbci=0; break;}
    }
  }
  else if(opt.do_ILS_GF_CI){
    if(opt.verbose) printf("Confidence intervals...\n");
    nbci=1;
    reti=calculate_ci(gdata, nb_gene, opt, optimal_par, max_lnL, "theta", ci+2);
    if(reti==0) nbci=0;
  }
  else nbci=0;
  
  
  // output detail to file    
  if(opt.verbose) printf("Writing to file...\n");    
  outf=fopen(argv[4], "w");
  if(outf==NULL) {printf("Cannot write file\n"); exit(0);}
  
  fprintf(outf, "gene,topology,a,b,c,d,l,lgseq,");
  for(j=1;j<nb_scenario_plus1;j++) fprintf(outf, "%s,", scenario_name[j]);
  fprintf(outf, "p_ILS,p_HGT\n");
  for(i=0;i<nb_gene;i++){
    fprintf(outf, "%s,%s,%f,%f,%f,%f,%f,%d,", gdata[i].name, ctopol(gdata[i].topology), gdata[i].a, gdata[i].b, gdata[i].c, gdata[i].d, gdata[i].l, gdata[i].lgseq);
    for(j=1;j<nb_scenario_plus1;j++) fprintf(outf, "%f,", glh[i].posterior_scenar[j]);
    fprintf(outf, "%f,%f\n", glh[i].posterior_ILS, glh[i].posterior_HGT);
  }
  if(opt.verbose) printf("\n");

  //summary stats
  opt_pab=optimal_par.pab; opt_pac=optimal_par.pac; opt_pbc=optimal_par.pbc; 
  opt_pils=exp(-2*(optimal_par.tau2-optimal_par.tau1)/optimal_par.theta);

  npred_ils=npred_hgt=0;
  for(i=0;i<nb_gene;i++){
    if(glh[i].predicted_scenario_type==1){ npred_ils++; }
    if(glh[i].predicted_scenario_type==2){ npred_hgt++; }
  }

  nils_a=nils_b=nhgt_a=nhgt_b=0;
  for(i=0;i<nb_gene;i++){
    if(glh[i].predicted_scenario_type==1){
      if(gdata[i].topology==1) nils_b++;
      if(gdata[i].topology==2) nils_a++;
    }
    if(glh[i].predicted_scenario_type==2){
      if(gdata[i].topology==1) nhgt_b++;
      if(gdata[i].topology==2) nhgt_a++;
    }
  }  
  
  //output summary
  if(opt.verbose){
    printf("*** Summary results: ***\n\n");
    printf("maximum likelihood: %.3f\n", max_lnL);
    printf("optimized parameters:\n");  
    print_par(&optimal_par, ci, nbci);
    printf("\n");

    printf("posterior no-event: %.3f\n", post.basic);
    printf("posterior no-conflict ILS: %.3f\n", post.noconflict_ILS);
    printf("posterior no-conflict GF : %.3f\n\n", post.noconflict_HGT);
        
    printf("posterior conflict ILS: %.3f", post.conflict_ILS);
    if(nbci!=0) printf(" [%.3f-%.3f]", ci[2].post_min.conflict_ILS, ci[2].post_max.conflict_ILS);
    printf("\n");
    printf("posterior conflict GF : %.3f", post.conflict_HGT);
    if(nbci!=0) printf(" [%.3f-%.3f]", ci[2].post_max.conflict_HGT, ci[2].post_min.conflict_HGT);
    printf("\n\n");   
    printf("posterior imbalance ILS: %.3f     [%s dominant]\n", post.asymmetry_ILS, ctopol(post.dominant_ILS));
    printf("posterior imbalance GF : %.3f     [%s dominant]\n", post.asymmetry_HGT, ctopol(post.dominant_HGT));
    printf("\n");  
  }
  
  //one-line output
  if(!opt.verbose){ 
    //printf("%d,%d,%d,%d,%d,%.1f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%s,%f,%s,%f\n", nb_gene, ntopo0, ntopo1, ntopo2, ntopo3, av_lg, optimal_par.tau1, optimal_par.tau2, optimal_par.theta, optimal_par.pab, optimal_par.pac, optimal_par.pbc, optimal_par.pone, post.basic, post.noconflict_ILS, post.noconflict_HGT, post.conflict_ILS, post.conflict_HGT, post.asymmetry_ILS, ctopol(post.dominant_ILS), post.asymmetry_HGT, ctopol(post.dominant_HGT), max_lnL);
    printf("%d,%d,%d,%d,%d,%.1f,", nb_gene, ntopo0, ntopo1, ntopo2, ntopo3, av_lg);
    printf("%f,%f,%f,", optimal_par.tau1, ci[0].master_param_min, ci[0].master_param_max);
    printf("%f,%f,%f,", optimal_par.tau2, ci[1].master_param_min, ci[1].master_param_max);
    printf("%f,%f,%f,", optimal_par.theta, ci[2].master_param_min, ci[2].master_param_max);
    printf("%f,%f,%f,", optimal_par.pab, ci[3].master_param_min, ci[3].master_param_max); 
    printf("%f,%f,%f,", optimal_par.pac, ci[4].master_param_min, ci[4].master_param_max);
    printf("%f,%f,%f,", optimal_par.pbc, ci[5].master_param_min, ci[5].master_param_max);
    printf("%f,%f,%f,", optimal_par.pone, ci[6].master_param_min, ci[6].master_param_max);
    printf("%f,%f,%f,", post.basic, post.noconflict_ILS, post.noconflict_HGT);
    printf("%f,%f,%f,", post.conflict_ILS, ci[2].post_min.conflict_ILS, ci[2].post_max.conflict_ILS);
    printf("%f,%f,%f,", post.conflict_HGT, ci[2].post_max.conflict_HGT, ci[2].post_min.conflict_HGT);
    printf("%f,%s,%f,%s,%f\n", post.asymmetry_ILS, ctopol(post.dominant_ILS), post.asymmetry_HGT, ctopol(post.dominant_HGT), max_lnL);
  }
}




