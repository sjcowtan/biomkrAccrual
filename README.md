# biomkrAccrual <img src="man/figures/logo.png" align="right" height="139" alt="" />
## Simulating recruitment at time of randomisation to adaptive trials with arm eligibility determined by biomarker status.

<!-- badges: start -->
[![R-CMD-check](https://github.com/sjcowtan/biomkrAccrual/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sjcowtan/biomkrAccrual/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/github/sjcowtan/biomkrAccrual/graph/badge.svg?token=1LIAWRVBU3)](https://codecov.io/github/sjcowtan/biomkrAccrual)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

# Overview

The `{biomkrAccrual}` package uses a Poisson-Gamma-Dirichlet model to simulate
trial recruitment for multi-site, multi-region, multi-arm trials. Recruitment per 
site is modelled with the Poisson-Gamma model (Anisimov and Federov, 2007).
A hierarchical Dirichlet model is used to model biomarker proportions for sites
within regions. Recruitment to a given site in a given week is then randomised
to biomarker status using the prevalences drawn from the Dirichlet model for that
site.

# Installation

```
devtools::install_github("sjcowtan/biomkrAccrual")
```

# Contacts

# Use cases

# Usage

## Running a single simulation

```
biomkrAccrual(precision = 10, var_lambda = 0.25)
```

The default settings will use the example configuration files in the `extdata` directory, and will
keep the resulting data files and recruitment plots. The location can be specified with 
`output_path` and `figs_path`.

## Running a set of simulations

```
biomkrAccrualSim(n = 250, precision = 10, var_lambda = 0.25)`
```

The datafiles and recruitment plots from the individual runs will not be kept (this
can be changed with `quietly = FALSE`) but will preserve the summary datafiles and
distribution plots.

## Configuration files

### Trial structure: how biomarkers map onto the trial arms{#arms-json}

The example `arms.json` file in `biomkrAccrual/inst/extdata/` looks like this:

```
{
  "T1":[1,2,3,4],
  "T2":[2,4,5,7],
  "T3":[3,4,6,7],
  "T4":[2,3,4,5,6,7],
  "T5":[2,3,4,5,6,7]
}
```

This describes a trial with five experimental arms and seven [biomarkers](#biomarkers).  Multiple biomarkers can be eligible for the same arm, and biomarkers can be eligible for multiple arms.  In this case the first four biomarkers are eligible for trial arm `T1`.

A file like this can be created by hand, or you can create one from a list in R using the `{jsonlite}` package:

```
r$> arms_ls <- list(
    T1 = 1:4, 
    T2 = c(2, 4, 6, 7), 
    T3 = c(3, 4, 6, 7), 
    T4 = 2:7, 
    T5 = 2:7
  )
r$> jsonlite::write_json(arms_ls, "arms.json")
```

This file will look a bit uglier:

```
{"T1":[1,2,3,4],"T2":[2,4,6,7],"T3":[3,4,6,7],"T4":[2,3,4,5,6,7],"T5":[2,3,4,5,6,7]}
```

Make sure you stay consistent with which biomarker is which in the other configuration files!

### Recruitment centres (sites)

Expectations about each planned recruitment centre
are encapsulated in `centres.csv`.  The example in `biomkrAccrual/inst/extdata/centres.csv` is rather long.  A shorter example with five sites might look like:

```
"site","start_month","mean_rate","region","site_cap"
1,1,3,1,15
2,1,5,1,10
3,1,2,1,9
4,3,1.6,1,8
5,4,1.3,1,10
```
- The `site` column is just an identifier for each site.
- Counting `1` as the month trial recruitment opens, `start_month` contains the month each centre opens to recruitment.
- The expected number of patients per month for each site is in the `mean_rate` column.
- If trial centres are likely to draw from very different populations, you might want to group those as `regions`, e.g. a trial that has some sites in the UK and some in India. This allows for different expectations of biomarker prevalences for the different regions.  In this example they're all in one region, region 1.  You do need this column, even if every entry is the same.
- The `site_cap` column is optional. In practice, trial centres often close to recruitment once they have recruited their agreed minimum number of patients.  Using this column allows you to take this into account when predicting recruitment.

### Biomarker prevalences{#biomarkers}

This example from `biomkrAccrual/inst/extdata/centres.csv` uses the same seven biomarkers as in the [trial structure](#arms-json) section.

```
"category","region_1","region_2"
"B1+, B2",0.47,0.41
"B1+, B3",0.1,0.21
"B1+, B4",0.05,0.04
"B1+, B3, B4",0.13,0.14
"B1-, B3",0.08,0.11
"B1-, B4",0.04,0.03
"B1-, B3, B4",0.13,0.06
```

- The `category` column is just some text which will tell you which biomarker is which. This example was actually using combinations of two biomarker statuses for each "biomarker", `B1` which was either positive or negative, in combination with `B2`, `B3` or `B4`.
- The `region_1` column contains the **expected prevalence** of each biomarker in the recruited population of the first region.
- Each region referenced in the [`centres.csv`](#biomarkers) file needs to have a column here with a column name in the same format, e.g. `region_2`.

### Recruitment targets

The recruitment targets for each [*experimental arm*](#arms-json) are defined in the `targets.csv` file.  This is the example from `biomkrAccrual/inst/extdata/targets.csv`:

```
"arm","interim","final"
"T1",30,60
"T2",30,60
"T3",30,60
"T4",30,60
"T5",30,60
```

- The names of the arms in the `arm` column correspond to those used in the [trial structure file](#arms-json).
- Recruitment targets can be specified for any number of interim monitoring points; here there is one called `interim`, but for a trial with monitoring at 6 and 24 months you could use `interim_6` and `interim_24`.
- The final recruitment target for each arm is in a column called `final`.
- The numbers in the columns represent the number of patients *in the experimental arms only*.  The number of patients in the control arms are specified by the run-time argument `control_ratio`.

## Run-time arguments

### Practical notes

There are a very large number of arguments to both commands, and four configuation files,
one of which (the relationship of treatment arms to biomarker recruitment arms) is a JSON.  
This is because flexibility is required, and they are intended to be driven by a dashboard 
in future.







# References

Anisimov, V.V., Fedorov, V.V., 2007. Modelling, prediction and adaptive adjustment of 
recruitment in multicentre trials. Statistics in Medicine 26, 4958–4975. 
https://doi.org/10.1002/sim.2956

