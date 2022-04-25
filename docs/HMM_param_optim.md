# HMM parameter optimisation

Using the code from `rule compare_params` in https://github.com/brettellebi/pilot_paper/blob/master/workflow/rules/05_param_optim.smk.

## Read in files


```r
CONC = list.files("/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/",full.names = T, recursive = T)
KW = "/hps/nobackup/birney/users/ian/pilot/kruskal_wallis/out.rds"
```

## Process



