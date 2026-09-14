# DILS practical: from dice to demographic histories

Which model could have produced our observations? We start with dice to explore how simulations help us compare competing explanations. We then apply the same approach to two mussel species, *Mytilus edulis* and *Mytilus galloprovincialis*: did they exchange genes, when did this happen, and did all parts of their genomes respond in the same way?

The practical combines short explanations, R exercises and questions to discuss. We will compare models, estimate parameters, check whether the selected model reproduces the observations, and identify candidate barrier loci.

Download this folder, keeping `workshop.html` and `data/` together. Open `workshop.html` in your browser and run the R code from the `DILS/` directory. Download the HTML file before opening it: the GitHub file viewer does not display the practical as a webpage.

## R packages

```r
install.packages(c(
  "tidyverse",
  "abcrf",
  "ggdensity",
  "ggthemes",
  "knitr",
  "scales"
))
```

- **tidyverse**: reading and manipulating tables, and plotting.
- **abcrf**: model comparison and parameter estimation using random forests.
- **ggdensity**, **ggthemes** and **scales**: additional plotting tools.
- **knitr**: formatting the parameter estimates as a table.

The demographic simulations are already provided. No simulator or Quarto installation is needed to follow the practical.

## Contents

Introduction, followed by a first exercise: **Before Mytilus: which die was drawn?**

1. Explore the data
2. Is there ongoing gene flow?
3. When did gene flow occur?
4. Does effective population size vary across the genome?
5. Does migration vary across the genome?
6. Putting the model choices together
7. Estimating effective population size, followed by estimates for all parameters
8. Can the selected model reproduce the data?
9. Looking at the joint pattern
10. Which loci resist gene flow?
