
**Code and data from Telenský, T., Storch, D., Klvaňa, P., Reif, J. Extension of Pradel capture-recapture survival-recruitment model accounting for transients. Methods in Ecology and Evolution.**

Implementations of basic as well extended model in BUGS language for both Nimble and JAGS, along with simple examples.

 - [code](/code) - R scripts to run the models using Nimble or JAGS
 - [data](/data) - the data from the Case study (Czech Constant Effort Sites scheme)
 - [model](/model) - basic and extended model in BUGS language

Install instructions

1. Install either Nimble (https://r-nimble.org/) or [JAGS](https://sourceforge.net/projects/mcmc-jags) (along with the R package 'runjags')
2. Install R packages:

```r
install.packages(c("coda", "reshape2", "stringr"))
```

3. Install the talonplot R package - either from a binary, if available, or from the source repository:

```r
require(devtools)
install_github("https://github.com/telenskyt/talonplot/")
```

This work is available under [CC-BY license](https://creativecommons.org/licenses/by/4.0/deed.en). If you reuse these scripts for your work, please cite our paper.



