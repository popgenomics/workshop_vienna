import numpy as np
import pandas as pd
import os
from itertools import groupby
import argparse

def all_equal(iterable):
    """
    Test if all the elements in a list are the same or not

    Parameters
    ----------
    iterable : list
        List of element

    Returns
    -------
    bool
        Return True if all the elements are the same and False if not.
    """
    g = groupby(iterable)
    return next(g, True) and not next(g, False)

def get_tgf_estimated(geneflow_df, topology, times):
    """
    Calculate the mean probability to observe gene flow at each possible timing and return the timing where the probability
    is higher.

    Parameters
    ----------
    gene_output : pandas.DataFrame
        Dataframe with the probability of different geneflow for each gene.
    topology : int
        Number corresponding at the current topology, 0: ((AB)C), 1: ((AC)B) and 2: ((BC)A)

    Returns
    -------
    max_timing : float
        Return the timing where the probabilty of significaive gene flow is maximum.
        posible value between 0.1 and 1.0 with 0.1 incrementation.
    """
    df = geneflow_df
    mean_proba_t = []
    timing = times
    #timing = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9, 1.0] ## ajouter une façon de le récupérer depuis le fichier de config
    if topology == 0: #((AB)C)
        topo_df = df[df.columns[pd.Series(df.columns).str.endswith('A<->B')]]
        for i in range(0, len(timing)):
            crt_df = topo_df.iloc[:,i]
            mean_proba_t.append(crt_df.mean())

    if topology == 1: #((AC)B)
        topo_df = df[df.columns[(pd.Series(df.columns).str.endswith('A->C') | pd.Series(df.columns).str.endswith('C->A'))]]
        CA = df[df.columns[(pd.Series(df.columns).str.endswith('A->C'))]]  ## The name of the variable is express backward and aphid express direction forward
        AC = df[df.columns[(pd.Series(df.columns).str.endswith('C->A'))]]  ## The name of the variable is express backward and aphid express direction forward
        for i in range(0, len(timing)*2, 2): # vérifier si bien len(timing)*2
            crt_df = topo_df.iloc[:, i:i+2]
            crt_df_sum = crt_df.sum(axis=1)
            mean_proba_t.append(crt_df_sum.mean())

    if topology == 2: #((BC)A)
        topo_df = df[df.columns[(pd.Series(df.columns).str.endswith('B->C') | pd.Series(df.columns).str.endswith('C->B'))]]
        CB = df[df.columns[(pd.Series(df.columns).str.endswith('B->C'))]]  ## The name of the variable is express backward and aphid express direction forward
        BC = df[df.columns[(pd.Series(df.columns).str.endswith('C->B'))]]  ## The name of the variable is express backward and aphid express direction forward
        for i in range(0, len(timing)*2, 2): # vérifier si bien len(timing)*2
            crt_df = topo_df.iloc[:, i:i+2]
            crt_df_sum = crt_df.sum(axis=1)
            mean_proba_t.append(crt_df_sum.mean())

    sum_pGF = topo_df.sum(axis=1).sum()
    if all_equal(mean_proba_t) and (mean_proba_t[0] == 0 or mean_proba_t[0] is np.nan):
        return None, None, sum_pGF, None
    else:
        ## get weighted tgf ##
        proba_timing = topo_df.groupby(lambda x: x.rsplit('_', 1)[0], axis=1).sum()
        gf_timing =  []
        proba_gf_per_gene = []
        for i in range(proba_timing.shape[0]):
            crt_line = proba_timing.iloc[i,::]
            ## ensure that crt_line and timing are in the same order ##
            labels = [f"GF_{t:.3f}" for t in timing]
            crt_ordered = crt_line[labels]
            sum_proba = crt_line.sum()
            proba_gf_per_gene.append(sum_proba)
            if sum(crt_line) != 0:
                gf_timing.append(np.round(np.average(a=timing, weights=crt_ordered),1))
            else:
                gf_timing.append(0)
        if sum(proba_gf_per_gene) != 0:
            weighted_t = round(np.average(a=gf_timing, weights=proba_gf_per_gene),1) ## V2 with averaging per proba GF
        else:
            weighted_t = None

        if topology == 0:
            return weighted_t, 'A<->B', sum_pGF, None
        else:
            ## determine direction by taking into account all the timing
            if topology == 1:
                proba_CA = CA.sum(axis=1)
                proba_AC = AC.sum(axis=1)
                mean_CA = proba_CA.mean()
                mean_AC = proba_AC.mean()
                AC_sup = sum(proba_AC > proba_CA)
                CA_sup = sum(proba_CA > proba_AC)
                dir_total = np.argmax([AC_sup, CA_sup])
                if np.min([mean_CA, mean_AC]) != 0:
                    bayes_f = np.max([mean_CA, mean_AC]) / np.min([mean_CA, mean_AC])
                else:
                    bayes_f = None
            elif topology == 2:
                proba_CB = CB.sum(axis=1)
                proba_BC = BC.sum(axis=1)
                mean_CB = proba_CB.mean()
                mean_BC = proba_BC.mean()
                BC_sup = sum(proba_BC > proba_CB)
                CB_sup = sum(proba_CB > proba_BC)
                dir_total = np.argmax([BC_sup, CB_sup])
                if np.min([mean_CB, mean_BC]) != 0:
                    bayes_f = np.max([mean_CB, mean_BC]) / np.min([mean_CB, mean_BC])
                else:
                    bayes_f = None
            ## prepare return
            dir_string_global = {1: ['A->C', 'C->A'], 2: ['B->C',
                                                          'C->B']}  ## prepare string for direction when determine with all timing
            dir_global = dir_string_global[topology][dir_total]
            return weighted_t, dir_global, sum_pGF, bayes_f

def get_proba(taxon, aphid_standard, processed_aphid, aphid_output, times):
    ## prepare data
    df = pd.read_csv(aphid_standard, header=0)
    df = df[df['nb_gene'] != 'missing triplet line in taxon file']
    df['nb_gene'] = df['nb_gene'].astype(float)
    df = df.reset_index(drop=True)
    n_row = df.shape[0]
    summary = {
        'taxon': [None] * (n_row*3),
        'triplet': [None] * (n_row*3),
        'species_A': [None] * (n_row*3),
        'species_B': [None] * (n_row*3),
        'topology': [None] * (n_row*3),
        '((AB)C)_majority': [None] * (n_row*3),
        'nb_gene': [None] * (n_row*3),
        'contribution_noevent' : [None] * (n_row*3),
        'contribution_ILS': [None] * (n_row*3),
        'imbalance_ILS': [None] * (n_row*3),
        'contribution_GF': [None] * (n_row*3),
        'imbalance_GF': [None] * (n_row*3),
        'tgf_estimated': [None] * (n_row*3),
        'direction_GF': [None] * (n_row*3),
        'tau1': [None] * (n_row*3),
        'tau2': [None] * (n_row*3),
        'theta': [None] * (n_row*3),
        'p_noevent_topology': [None] * (n_row*3),
        'pILS_topology': [None] * (n_row * 3),
        'pGF_topology': [None] * (n_row*3),
        'bayes_factor': [None] * (n_row*3)
    }
    
    last_index = 0
    for i in range(df.shape[0]):
        ntopo3 = df['ntopo3'][i] # get number of polytomies
        t1 = df['tau1'][i] # get tau1 inferred by aphid, tau1 = t1 * µ
        t2 = df['tau2'][i] # get tau2 inferred by aphid, tau1 = t2 * µ
        theta = df['theta'][i] # get theta inferred by aphid, theta = 4Ne * µ
        file = df['dataset'][i] # get the path to the taxa file
        name = os.path.basename(file) # change path to the file name
        name = name.replace('.tax', '.csv') # change to have output file name

        with open(file) as filein:
            triplet = filein.readline().rstrip()
            triplet_string = triplet.replace('triplet: ', '')
            triplet = triplet_string.split(', ')
            triplet_string = triplet_string.replace(', ', '_')
        #aphid_output = os.path.join(aphid_output, f'out_{name}')
        gene_df = pd.read_csv(aphid_output, header=0)
        index_triplet = [(0, 1, 2), (0, 2, 1), (1, 2, 0)]
        topology = ['((AB)C)', '((AC)B)', '((BC)A)']
        ## Count the number of concordant trees and the total number of trees
        counter = gene_df[gene_df['topology'] == topology[0]].shape[0]
        total_gene = gene_df.shape[0]
        ## calculate contribution of each topology
        df_ils = gene_df[gene_df.columns[pd.Series(gene_df.columns).str.startswith('ILS_')]]
        df_gf = gene_df[gene_df.columns[pd.Series(gene_df.columns).str.startswith('GF_')]]
        df_noevent = gene_df[gene_df.columns[pd.Series(gene_df.columns).str.startswith('no_event')]]
        for it in range(len(index_triplet)):
            summary['triplet'][last_index] = triplet_string
            summary['topology'][last_index] = topology[it]
            summary['tau1'][last_index] = t1
            summary['tau2'][last_index] = t2
            summary['theta'][last_index] = theta
            ## Add if the current triplet have a majority of concordant trees
            if counter / (total_gene - ntopo3) >= 0.5:
                summary['((AB)C)_majority'][last_index] = 1
            else:
                summary['((AB)C)_majority'][last_index] = 0
            ## Get each species 
            spe_A = triplet[index_triplet[it][0]]
            spe_B = triplet[index_triplet[it][1]]
            spe_C = triplet[index_triplet[it][2]]
            ## 
            df_topo_ils = df_ils.loc[gene_df[gene_df['topology'] == topology[it]].index]
            df_topo_gf = df_gf.loc[gene_df[gene_df['topology'] == topology[it]].index]
            df_topo_noevent = df_noevent.loc[gene_df[gene_df['topology'] == topology[it]].index]
            summary['p_noevent_topology'][last_index] = df_topo_noevent.mean().values[0] ### À rajouter dans summary pour avoir les contribution au sein d'une topo
            summary['pILS_topology'][last_index] = df_topo_ils.sum(axis=1).mean()
            summary['pGF_topology'][last_index] = df_topo_gf.sum(axis=1).mean()
            n_gene = df_topo_gf.shape[0]
            if n_gene:
                ## calculate timing and direction where the probability of gene flow is the higher
                weighted_tgf, dir_gf_global, sum_pgf, bayes_f = get_tgf_estimated(geneflow_df=df_topo_gf, topology=it, times=times)
                ## calculate number of gene with the majority of gene flow
                # signi_gf = gene_df[(gene_df['topology'] == topology[it]) & (gene_df['p_HGT'] > (gene_df['p_ILS'] + gene_df['no_event']))]
                # signi_ils = gene_df[(gene_df['topology'] == topology[it]) & (gene_df['p_ILS'] > (gene_df['p_HGT'] + gene_df['no_event']))]
                # signi_noevent = gene_df[(gene_df['topology'] == topology[it]) & (gene_df['no_event'] > (gene_df['p_ILS'] + gene_df['p_HGT']))]
                # summary['nb_gene_GF'][last_index] = signi_gf.shape[0]
                # summary['nb_gene_ILS'][last_index] = signi_ils.shape[0]
                # summary['nb_gene_noevent'][last_index] = signi_noevent.shape[0]
                ## calculate contribution / probability of each event for each topology
                contrib_noevent = (df_topo_noevent.sum()) / (df_noevent.shape[0]-ntopo3) # n_gene less gene with undefined topology
                contrib_noevent = contrib_noevent.values[0]
                contrib_ils = (df_topo_ils.iloc[:, it].sum()) / (df_ils.shape[0]-ntopo3)  # n_gene less gene with undefined topology
                contrib_gf = sum_pgf / (df_gf.shape[0]-ntopo3)  # n_gene less gene with undefined topology
                if it == 0 :
                    if weighted_tgf is not None:
                        pass
                    else:
                        weighted_tgf = -1
                        dir_gf_global = -1
                elif it == 1:
                    if weighted_tgf is not None:
                        pass
                    else:
                        weighted_tgf = -1
                        dir_gf_global = -1
                elif it == 2:
                    if weighted_tgf is not None:
                        pass
                    else:
                        weighted_tgf = -1
                        dir_gf_global = -1

                ## add probability for each topology
                summary['contribution_noevent'][last_index] = contrib_noevent
                summary['contribution_ILS'][last_index] = contrib_ils
                summary['contribution_GF'][last_index] = contrib_gf
                summary['tgf_estimated'][last_index] = weighted_tgf
                summary['direction_GF'][last_index] = dir_gf_global
            ## add general information
            summary['taxon'][last_index] = taxon
            summary['species_A'][last_index] = spe_A
            summary['species_B'][last_index] = spe_B
            summary['nb_gene'][last_index] = n_gene
            ## add imbalance
            summary['imbalance_ILS'][last_index] = df['imbalance_ILS'][i]
            summary['imbalance_GF'][last_index] = df['imbalance_GF'][i]
            ## add bayes factor
            summary['bayes_factor'][last_index] = bayes_f
            last_index += 1
    output_df = pd.DataFrame.from_dict(summary)
    output_df.to_csv(processed_aphid, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--taxon', '-t', type=str)
    parser.add_argument('--aphid_standard', '-f', type=str)
    parser.add_argument('--output', '-o', type=str)
    parser.add_argument('--aphid_output', '-g', type=str)
    parser.add_argument('--times', nargs='+', type=float)
    #parser.add_argument('--workingdir', '-w', type=str)
    args = parser.parse_args()

    #os.chdir(args.workingdir)
    get_proba(taxon=args.taxon, aphid_standard=args.aphid_standard, processed_aphid=args.output, aphid_output=args.aphid_output, times=args.times)
