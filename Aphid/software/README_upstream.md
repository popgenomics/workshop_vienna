# Aphid

This is a C implementation of the "Aphid" method,
following its specifications within Nicolas Galtier's draft paper
*An approximate likelihood method reveals ancient gene flow
 between human, chimpanzee and gorilla* (2023 *in prep.*),
and the original C code.

Sources available at: https://gitlab.mbb.cnrs.fr/ibonnici/aphid/.

## Build instructions

### With Apptainer/Singularity

Aphid can be run from an apptainer
provided you have [Apptainer] (aka. Singularity) installed.

First, build the container image as root:

```sh
# Get definition file.
curl -L \
  https://gitlab.mbb.univ-montp2.fr/ibonnici/aphid/-/jobs/artifacts/main/raw/singularity.def?job=container-recipes \
  > aphid.def # (or download by hand)

# Build image.
sudo apptainer build aphid.sif aphid.def
```

Then run the container as regular user:

```sh
./aphid.sif <my_input_files>
```

[Apptainer]: https://apptainer.org/

### With Docker

Aphid can be run from a docker container,
provided you have [Docker] installed.

First, build the container image:
```sh
docker buildx build -t aphid \
  https://gitlab.mbb.univ-montp2.fr/ibonnici/aphid/-/jobs/artifacts/main/raw/Dockerfile?job=container-recipes
```

Then run the container:
```sh
docker run --rm -it -v ${PWD}:/home/aphid aphid <my_input_files>
```

[Docker]: https://www.docker.com/

### Manual compilation

The program is still fairly easy to compile yet,
provided you have [git] and [gcc] installed.

First clone and compile the project:

```sh
git clone --recursive https://gitlab.mbb.univ-montp2.fr/ibonnici/aphid/
cd aphid
gcc -lm aphid.*.c -o aphid
```

Then run the binary obtained:

```sh
./aphid <my_input_files>
```

[git]: https://git-scm.com/
[gcc]: https://gcc.gnu.org/

## Running Aphid

## Input

Aphid reads data from three distinct files (see examples below):
- A __gene-tree__ file.
- A __taxon__ file.
- A __control__ file.

### The Gene-Tree File

The gene-tree file includes one line per gene tree.
Each line contains three tab-separated entries:
- Gene tree.
- Locus length.
- Locus name.

The *gene tree* must be in Newick format, rooted,
with branch lengths expressed in per site number of mutation.  
The *locus length* corresponds to the number of sites of the alignment
from which the gene tree was reconstructed,
such that the product of branch length by locus length
is the total number of mutations having occurred in the considered branch
at the considered locus.  
The *locus name* is an unconstrained string supposed to identify the tree.

### The Taxon File

The taxon file includes:
- A line defining the focal `triplet` (mandatory).
- A line defining the `outgroup` (mandatory).
- An optional line defining `other` relevant species
  according to the following format:

```
triplet: species_A, species_B, species_C
outgroup: <comma-separated species list>
other: <comma-separated species list>
```

The `triplet` line defines the three focal species of the analysis.
Order matters: the third species of the list is the earliest branching one,
*i.e.* it corresponds to species `C` of the assumed `((A,B),C)` species tree.  
The `outgroup` and `other` lines
define the species used to estimate locus-specific mutation rates.
Gene trees not including the three species of the `triplet`
and at least one `outgroup` species are not analyzed.
Species not listed in the taxon file, if any, are not considered.

The `other` species should branch somewhere
beween the `outgroup` and the focal `triplet` in the species tree.
Their inclusion is intended to increase the accuracy
of relative mutation rate estimation,
especially when the `outgroup` includes few species.
Consider for instance the following species tree:

```
(out,(other1,(other2,(other3,(other4,(C,(B,A)))))));
```

If we only specify the triplet `A, B, C`
and the outgroup `out` in the taxon file,
the relative mutation rate of a locus will be estimated
based on these four species.
If, in the taxon file we add the following line:
```
other = other1, other2, other3, other4
```

Then the relative mutation rate of a locus will be estimated
based on up to eight species,
depending on how many of the `other` species
are present in the considered gene tree.

### The Control File

Options are passed *via* the control file
using the `<option_name> = <option_value>` syntax (see example below):

- `require_triplet_other_monophyly`:
  If set to `1` (default),
  only gene trees in which the `triplet` species and `other` species (if any)
  form a monophyletic group will be analyzed.

- `require_outgroup_monophyly`:
  If set to `1` (default),
  only gene trees in which the `outgroup` species form a monophyletic group
  will be analyzed.

- `max_alphai`:
  Gene trees with an estimated relative mutation rate higher than this value,
  or lower than 1 / this value, will not be analyzed (default = `10`).

- `max_clock_ratio`:
  Gene trees in which the ratio of `triplet` estimated mutation rate to
  `outgroup` + `other` estimated mutation rate exceeds this value,
  or is below 1 / this value, will not be analyzed (default = `2`).

- `nb_GF_times`:
  Number of assumed gene flow times (default = `2`).

- `GF_time_coeff`:
  Comma-separated list of `nb_GF_times` coefficients
  specifying the gene flow times,
  expressed as fractions of the `A|B` speciation time.
  Each time must be between `0` and `1`.
  (default = `1.0, 0.5`).

- `min_d`:
  If internal branch length times locus length is below this value,
  then the gene tree will be considered topologically unresolved
  (default = `0.5`).

- `CI_param`:
  If set to `1`,
  confidence intervals around parameter estimates will be calculated
  (default = `0`).

- `CI_ILS_GF`:
  If set to `1`,
  confidence intervals around the estimated contribution of ILS and GF
  to the conflict will be calculated (default = `0`).

- `verbose`:
  If set to `1` (default),
  results are output in a human-readable way.
  If set to `0` results are output as a single line (see below).

- `fixed_tau1`:
  If set to a positive value,
  parameter `τ_1` will not be optimized but rather fixed to that value
  (default = `-1`).

- `fixed_tau`, `fixed_theta`, `fixed_pab`, `fixed_pac`, `fixed_pbc`, `fixed_panc`:
  Same as above for parameters `τ_2`, `θ`, `p_ab`, `p_ac`, `p_bc`, and `p_a`.
  If `nb_GF_times` is set to `1`,
  then `fixed_panc` is automatically set to `1.0`.

Lines starting with `#` are ignored
and can be used for inserting comments in the control file.

## Output

### Standard output

Parameter estimates,
ILS/GF estimated contributions,
discordant topology imbalance
and various descriptors of the data and analysis
are written to the standard output (terminal).
If `verbose` is set to 1, the output is detailed and human-readable.
If `verbose` is set to 0, the output is less detailed
and written as a single line in the following format:

```  
_nb_gene,          \
ntopo0,            \
ntopo1,            \
ntopo2,            \
ntopo3,            \
av_lg,             \
tau1,              \
tau1_low,          \
tau1_high,         \
tau2,              \
tau2_low,          \
tau2_high,         \
theta,             \
theta_low,         \
theta_high,        \
pab,               \
pab_low,           \
pab_high,          \
pac,               \
pac_low,           \
pac_high,          \
pbc,               \
pbc_low,           \
pbc_high,          \
pa,                \
pa_low,            \
pa_high,           \
no_event,          \
noconflict_ILS,    \
noconflict_GF,     \
conflict_ILS,      \
conflict_ILS_low,  \
conflict_ILS_high, \
conflict_GF,       \
conflict_GF_low,   \
conflict_GF_high,  \
imbalance_ILS,     \
dominant_ILS,      \
imbalance_GF,      \
dominant_GF,       \
max_lnL_
```

where:

- `nb_gene`:
  Number of analyzed gene trees, after filtering
  (detailed output gives mor einformation on filtering steps).

- `ntopo, ntopo1, ntopo2, ntopo3`:
  Number of analyzed gene trees
  with an `((A,B),C)`, `((A,C),B)`, `((B,C),A)`
  and `(A,B,C)` (unresolved) topology, respectively.

- `av_lg`:
  Mean locus length of analyzed gene trees.

- `tau1, tau2, theta, pab, pac, pbc, pa`:
  Parameters estimates,
  with the `_low` and `_high` suffix indicating confidence interval boundaries.

- `no_event`:
  Estimated contribution of the `no_event` scenario

- `no_conflict_ILS`:
  Estimated contribution of the ILS scenario
  with an `((A,B),C)` topology.

- `no_conflict_GF`:
  Estimated contribution of GF scenarios
  with an `((A,B),C)` topology.

- `conflict_ILS`:
  Estimated contribution of ILS scenarios
  with a topology different from `((A,B),C)`.

- `conflict_GF`:
  Estimated contribution of GF scenarios
  with a topology different from `((A,B),C)`.

- `imbalance_ILS`:
  Estimated discordant topology imbalance associated with ILS.

- `dominant_ILS`:
  The most prevalent discordant topology associated with ILS.

- `imbalance_GF`:
  Estimated discordant topology imbalance associated with GF.

- `dominant_ILS`:
  The most prevalent discordant topology associated with GF.

- `max_lnL`:
  Maximum likelihood.

The `no_event`, `no_conflict_ILS` and `no_conflict_GF` posterior estimates
are provided.
But it is unclear how reliable these estimates generally are.
Their sum correspond to the overall estimated prevalence
of the canonical `((A,B),C)` topology among gene trees.
The discordant topology imbalance is a number
between 0.5 (no imbalance) and 1 (maximal imbalance)
indicating whether one of the two discordant topologies
`((A,C),B)` and `((B,C),A)` dominates over the other one.
Which one dominates is specified by `dominant_ILS` and `dominant_GF`.

### Output File

In addition to the information displayed on standard output,
Aphid generates a `.csv` file including one line per analyzed gene tree.
This file provides gene tree-specific descriptors
(name, topology, branch lengths, locus length)
and Aphid annotations
— specifically, the posterior probabilities
of each scenario for each gene tree.

## Example Input Files

### Gene-Tree File

`example.in`
```
((out1:0.07,out2:0.07):0.03,(oth1:0.06,((A:0.02,B:0.02):0.01,C:0.03):0.03):0.04);	1000	canonical
((out1:0.068,out2:0.066):0.028,(oth1:0.059,((A:0.042,C:0.044):0.002,B:0.045):0.017):0.036);	500	discordant_tall
((out1:0.068,out2:0.066):0.028,(oth1:0.059,((B:0.042,C:0.044):0.002,A:0.045):0.017):0.036);	500	other_discordant_tall
((out1:0.064,out2:0.075):0.027,(oth1:0.058,((B:0.01,C:0.008):0.01,A:0.02):0.027):0.036);	500	discordant_short
((out1:0.064,out2:0.075):0.027,(oth1:0.058,((B:0.019,C:0.022):0.01,A:0.028):0.037):0.036);	555	ancient_GF
((out1:0.07,out2:0.07):0.03,(oth1:0.06,((A:0.035,B:0.037):0.003,C:0.041):0.02):0.04);	1212	canonical_ILS
((out1:0.04,out2:0.03):0.02,(oth1:0.04,((A:0.035,B:0.038):0.026,C:0.062):0.052):0.05);	666	non_clock_like
((out1:0.07,out2:0.07):0.03,(oth1:0.06,((A:0.008,B:0.009):0.02,C:0.03):0.03):0.04);	421	canonical_GF
((out1:0.07,out2:0.07):0.03,((A:0.02,B:0.02):0.01,C:0.03):0.07);	1500	missing_other=OK
((out1:0.07,out2:0.07):0.03,(oth1:0.06,(A:0.03,C:0.03):0.03):0.04);	1001	incomplete_triplet=notOK
(out2:0.1,(oth1:0.06,((A:0.02,B:0.02):0.01,C:0.03):0.03):0.04);	1234	missing_1outgroup=OK
(out1:0.1,((A:0.02,B:0.02):0.01,C:0.03):0.07);	100	minimal_sampling
(out2:0.15,((A:0.03,B:0.02):0.015,C:0.025):0.06);	1200	other_minimal_sampling
(oth1:0.06,((A:0.02,B:0.02):0.01,C:0.03):0.03);	600	missing_all_outgroup=notOK
((out1:0.07,out2:0.07):0.03,((oth1:0.04,additional1:0.04):0.02,((A:0.02,B:0.02):0.01,C:0.03):0.03):0.04);	321	additional_species=OK
((additional2:0.04,additional3:0.04):0.06,((out1:0.07,out2:0.07):0.03,((oth1:0.04,additional1:0.04):0.02,((A:0.02,B:0.02):0.01,C:0.03):0.03):0.04):0.03);	1111	more_additional_species=stillOK
((additional2:0.04,additional3:0.04):0.06,((out1:0.07,(out2:0.06,additional4:0.06):0.01):0.03,((oth1:0.04,additional1:0.04):0.02,((A:0.02,B:0.02):0.01,(C:0.02,additional5:0.02):0.01):0.03):0.04):0.03);	200	even_more_additional_species=stillOK
(additional1:0.1,(oth1:0.06,((A:0.02,B:0.02):0.01,C:0.03):0.03):0.04);	2000	missing_all_outgroup=notOK_even_if_additional
(out1:0.2,(out2:0.14,(oth1:0.12,((A:0.04,B:0.04):0.02,C:0.06):0.06):0.08):0.06);	500	paraphyletic_outgroup=OK_or_not_depending_options
(other1:0.1,((out1:0.07,out2:0.07):0.03,((A:0.02,B:0.02):0.01,C:0.03):0.04):0.01);	500	paraphyletic_triplet+other=OK_or_not_depending_options
```

### Taxon file

`example.tax`
```
triplet: A, B, C
outgroup: out1, out2
other: oth1
```

### Control File

`example.opt`
```
# Gene tree filtering.
require_triplet_other_monophyly = 1
require_outgroup_monophyly = 1
max_alphai = 10
max_clock_ratio = 2

# Likelihood calculation.
nb_GF_times = 2
GF_time_coeff = 1.0, 0.5
min_d = 0.5
# fixed_pac = 0
# (uncomment the above line to calculate the likelihood
#  assuming no GF between A and C)

# Output.
CI_param = 0
CI_ILS_GF = 0
verbose = 1
```
