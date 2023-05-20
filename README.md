# modtcov

## Overview
This software implements moderated t-covariates to improve large-scale inferences under the limma EB test. A machine learning model is available to fit local false discoveries rate to rank differentially expressed genes.

## Installation         

```
# source("https://bioconductor.org/biocLite.R")
# biocLite("limma")
# install.packages("devtools")
devtools::install_github("ninhtran02/modtcov")
```

### An Example

```
# Load package
library("modtcov")

mod.t.cov.df <- create.t.stat.cov.df(exampleY, exampleX, c(1,-1,0,0,0))
fit <- fit.t.model(mod.t.cov.df)
mlfdr <- mlfdr_fn_tmm(mod.t.stat = mod.t.cov.df$mod.t.stat, prob_mat = fit$prob_mat, 
                      param = fit$param, df = mod.t.cov.df$df, alt.df = mod.t.cov.df$df)	
```
