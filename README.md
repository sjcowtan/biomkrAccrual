# biomkrAccrual
## Simulating recruitment at time of randomisation to adaptive trials with arm eligibility determined by biomarker status.

<!-- badges: start -->
[![R-CMD-check](https://github.com/sjcowtan/biomkrAccrual/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sjcowtan/biomkrAccrual/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/github/sjcowtan/biomkrAccrual/graph/badge.svg?token=1LIAWRVBU3)](https://codecov.io/github/sjcowtan/biomkrAccrual)
<!-- badges: end -->

The `{biomkrAccrual}` package uses a Poisson-Gamma-Dirichlet model to simulate
trial recruitment for multi-site, multi-region, multi-arm trials. Recruitment per 
site is modelled with the Poisson-Gamma model (Anisimov and Federov, 2007).
A hierarchical Dirichlet model is used to model biomarker proportions for sites
within regions. Recruitment to a given site in a given week is then randomised
to biomarker status using the prevalences drawn from the Dirichlet model for that
site.

## Running a single simulation

`biomkrAccrual()`

The default settings will use the configuration files in the `extdata` directory, and will
keep the resulting data files and recruitment plots. The location can be specified with 
`output_path` and `figs_path`.

## Running a set of simulations

`biomkrAccrualSim(n = 250)`

The datafiles and recruitment plots from the individual runs will not be kept (this
can be changed with `quietly = FALSE`) but will preserve the summary datafiles and
distribution plots.

## Practical notes

There are a very large number of arguments to both commands, and three configuation files,
one of which (the relationship of treatment arms to biomarker recruitment arms) is a JSON.  
This is because flexibility is required, and they are intended to be driven by a dashboard 
in future.

This package will not pass R CMD Check because it is written using the new object orientation 
system for R, `{S7}`, and CMD Check does not yet understand the syntax.


Anisimov, V.V., Fedorov, V.V., 2007. Modelling, prediction and adaptive adjustment of 
recruitment in multicentre trials. Statistics in Medicine 26, 4958–4975. 
https://doi.org/10.1002/sim.2956

