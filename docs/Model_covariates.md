# Model covariates


```r
library(tidyverse)
library(DT)
```
## Read in data


```r
IN = "/hps/nobackup/birney/users/ian/pilot/hmm_out/0.05/dist_angle/18.csv"
```


```r
# Read 

df = readr::read_csv(IN)
#> Rows: 14651580 Columns: 15
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr (6): assay, ref_fish, test_fish, tank_side, quadrant...
#> dbl (9): date, time, frame, seconds, x, y, distance, ang...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

# Controls (iCab-iCab pairings only)

df_control = df %>% 
  # Filter for only iCabs paired with iCabs
  dplyr::filter(test_fish == "icab") %>% 
  # Get individual
  tidyr::unite(date, time, quadrant, fish,
               col = "indiv", 
               remove = F) %>% 
  # Group by assay and individual to get mean speed
  dplyr::group_by(assay, date, time, quadrant, indiv) %>% 
  # Calculate mean speed
  dplyr::summarise(mean_speed = mean(distance)) %>% 
  dplyr::ungroup() %>% 
  # Make date a factor
  dplyr::mutate(date = as.factor(date))
#> `summarise()` has grouped output by 'assay', 'date',
#> 'time', 'quadrant'. You can override using the `.groups`
#> argument.

# DGE data frame

df_dge = df %>% 
  # remove iCab when paired with a different test fish
  dplyr::filter(!(fish == "ref" & test_fish != "icab")) %>% 
  # Get individual
  tidyr::unite(date, time, quadrant, fish,
               col = "indiv", 
               remove = F) %>% 
  # add `line` %>% 
  dplyr::mutate(line = dplyr::case_when(fish == "ref" ~ ref_fish,
                                        fish == "test" ~ test_fish))  %>% 
  # Group by assay and individual to get mean speed
  dplyr::group_by(assay, date, time, quadrant, line, indiv) %>% 
  # Calculate mean speed
  dplyr::summarise(mean_speed = mean(distance)) %>% 
  dplyr::ungroup()
#> `summarise()` has grouped output by 'assay', 'date',
#> 'time', 'quadrant', 'line'. You can override using the
#> `.groups` argument.

# SGE data frame 

df_sge = df %>% 
  # keep only iCabs
  ## add `line`
  dplyr::mutate(line = dplyr::case_when(fish == "ref" ~ ref_fish,
                                        fish == "test" ~ test_fish)) %>% 
  # keep only iCab test fishes
  dplyr::filter(line == "icab") %>% 
  # Get individual
  tidyr::unite(date, time, quadrant, fish,
               col = "indiv", 
               remove = F) %>% 
  # Group by assay and individual to get mean speed
  dplyr::group_by(assay, date, time, quadrant, test_fish, indiv) %>% 
  # Calculate mean speed
  dplyr::summarise(mean_speed = mean(distance)) %>% 
  dplyr::ungroup()
#> `summarise()` has grouped output by 'assay', 'date',
#> 'time', 'quadrant', 'test_fish'. You can override using the
#> `.groups` argument.
```


## Controls


```r
model_fit = lm(mean_speed ~ assay + date + time + quadrant, data = df_control)

summary(model_fit)
#> 
#> Call:
#> lm(formula = mean_speed ~ assay + date + time + quadrant, data = df_control)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.30927 -0.36602  0.06531  0.35224  1.44668 
#> 
#> Coefficients:
#>                   Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)      0.9551959  0.2385450   4.004 8.11e-05 ***
#> assayopen_field  0.5160690  0.0663243   7.781 1.66e-13 ***
#> date20190612     0.1976414  0.1579453   1.251 0.211931    
#> date20190614     0.6531556  0.1675304   3.899 0.000123 ***
#> date20190615     0.3596363  0.1514893   2.374 0.018318 *  
#> date20190616     0.3963429  0.1478817   2.680 0.007826 ** 
#> time             0.0001001  0.0001381   0.725 0.469200    
#> quadrantq2       0.1635167  0.0937967   1.743 0.082453 .  
#> quadrantq3       0.1018693  0.0937967   1.086 0.278449    
#> quadrantq4      -0.2136390  0.0937967  -2.278 0.023552 *  
#> ---
#> Signif. codes:  
#> 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.5469 on 262 degrees of freedom
#> Multiple R-squared:  0.2773,	Adjusted R-squared:  0.2524 
#> F-statistic: 11.17 on 9 and 262 DF,  p-value: 9.486e-15
```


## DGE


```r
## All variables with interactions
model_dge = aov(mean_speed ~ assay * date * time * quadrant * line, data = df_dge)

dge_vars = summary(model_dge)[[1]] %>% 
  dplyr::mutate(TOT_SS = sum(`Sum Sq`),
                PROP_SS = `Sum Sq` / TOT_SS,
                VAR_EXP = PROP_SS * 100)

DT::datatable(dge_vars)
```

```{=html}
<div id="htmlwidget-a6e4dc872c0e594add90" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-a6e4dc872c0e594add90">{"x":{"filter":"none","vertical":false,"data":[["assay                   ","time                    ","quadrant                ","line                    ","assay:time              ","assay:quadrant          ","time:quadrant           ","assay:line              ","time:line               ","quadrant:line           ","assay:time:quadrant     ","assay:time:line         ","assay:quadrant:line     ","time:quadrant:line      ","assay:time:quadrant:line","Residuals               "],[1,1,3,4,1,3,3,4,4,12,3,4,12,12,12,670],[46.2641726229261,0.477162707390935,4.89124978680217,6.62925242717848,0.0891432273633267,0.271696770834394,0.816445406761689,2.19581781814295,3.74030287054715,5.36826726781128,0.652018061758932,0.33086008769393,2.30442928254236,8.3750637881556,1.55030026486377,229.943116954351],[46.2641726229261,0.477162707390935,1.63041659560072,1.65731310679462,0.0891432273633267,0.0905655902781312,0.272148468920563,0.548954454535737,0.935075717636786,0.44735560565094,0.217339353919644,0.0827150219234825,0.192035773545196,0.6979219823463,0.129191688738647,0.343198682021419],[134.802885461077,1.39033956826546,4.75064934981005,4.82901943863288,0.259742335950355,0.263886765953486,0.792976439529502,1.59952378401464,2.72459005998871,1.30348870518978,0.633275607702012,0.24101206169061,0.559546943520055,2.03358001911774,0.376434105101208,null],[1.59458952912767e-28,0.238766308499698,0.00275532899152628,0.00076102519572237,0.610465052956788,0.851431410613099,0.498026665407507,0.172692211764521,0.0285813902676841,0.211385318427479,0.593726980864489,0.915086857319263,0.874874182023666,0.0194238656102273,0.971724040662664,null],[313.899299345124,313.899299345124,313.899299345124,313.899299345124,313.899299345124,313.899299345124,313.899299345124,313.899299345124,313.899299345124,313.899299345124,313.899299345124,313.899299345124,313.899299345124,313.899299345124,313.899299345124,313.899299345124],[0.147385396270222,0.00152011396134499,0.0155822258826528,0.0211190418105706,0.000283986703854716,0.000865553925737407,0.0026009787484872,0.00699529378601354,0.0119156139511952,0.0171018771912231,0.00207715679238282,0.00105403257791333,0.00734130113494996,0.0266807342533996,0.00493884589133552,0.732537847118718],[14.7385396270222,0.152011396134499,1.55822258826528,2.11190418105706,0.0283986703854716,0.0865553925737407,0.26009787484872,0.699529378601354,1.19156139511952,1.71018771912231,0.207715679238282,0.105403257791333,0.734130113494996,2.66807342533996,0.493884589133552,73.2537847118718]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Df<\/th>\n      <th>Sum Sq<\/th>\n      <th>Mean Sq<\/th>\n      <th>F value<\/th>\n      <th>Pr(&gt;F)<\/th>\n      <th>TOT_SS<\/th>\n      <th>PROP_SS<\/th>\n      <th>VAR_EXP<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7,8]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
## No interactions
model_dge_noint = aov(mean_speed ~ assay + date + time + quadrant + line, data = df_dge)

dge_vars_noint = summary(model_dge_noint)[[1]] %>% 
  dplyr::mutate(TOT_SS = sum(`Sum Sq`),
                PROP_SS = `Sum Sq` / TOT_SS,
                VAR_EXP = PROP_SS * 100)

DT::datatable(dge_vars_noint)
```

```{=html}
<div id="htmlwidget-cdfe524b76ec02a532dd" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-cdfe524b76ec02a532dd">{"x":{"filter":"none","vertical":false,"data":[["assay      ","time       ","quadrant   ","line       ","Residuals  "],[1,1,3,4,740],[46.2641726229261,0.477162707390935,4.89124978680217,6.62925242717848,255.637461800826],[46.2641726229261,0.477162707390935,1.63041659560072,1.65731310679462,0.345456029460576],[133.922029657919,1.38125453516043,4.71960671274602,4.79746470016021,null],[1.38724660696488e-28,0.240265575045801,0.00285853241021431,0.000796397786962862,null],[313.899299345124,313.899299345124,313.899299345124,313.899299345124,313.899299345124],[0.147385396270222,0.00152011396134499,0.0155822258826528,0.0211190418105706,0.81439322207521],[14.7385396270222,0.152011396134499,1.55822258826528,2.11190418105706,81.439322207521]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Df<\/th>\n      <th>Sum Sq<\/th>\n      <th>Mean Sq<\/th>\n      <th>F value<\/th>\n      <th>Pr(&gt;F)<\/th>\n      <th>TOT_SS<\/th>\n      <th>PROP_SS<\/th>\n      <th>VAR_EXP<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7,8]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


## SGE


```r
model_sge = aov(mean_speed ~ assay * date * time * quadrant * test_fish, data = df_sge)

sge_vars = summary(model_sge)[[1]] %>% 
  dplyr::mutate(TOT_SS = sum(`Sum Sq`),
                PROP_SS = `Sum Sq` / TOT_SS,
                VAR_EXP = PROP_SS * 100)

DT::datatable(sge_vars)
```

```{=html}
<div id="htmlwidget-c31c2515472e149c0c66" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-c31c2515472e149c0c66">{"x":{"filter":"none","vertical":false,"data":[["assay                        ","time                         ","quadrant                     ","test_fish                    ","assay:time                   ","assay:quadrant               ","time:quadrant                ","assay:test_fish              ","time:test_fish               ","quadrant:test_fish           ","assay:time:quadrant          ","assay:time:test_fish         ","assay:quadrant:test_fish     ","time:quadrant:test_fish      ","assay:time:quadrant:test_fish","Residuals                    "],[1,1,3,4,1,3,3,4,4,12,3,4,12,12,12,670],[38.9131270547654,0.0990213326289756,4.06490447676763,6.72803013681188,0.00464892767486915,0.162822924911266,1.05877293498819,2.34412064412351,4.62349185479105,6.57458135510301,0.765926400983002,0.44823570739265,1.5756628872581,8.9963834273588,3.13499781922733,243.564222038123],[38.9131270547654,0.0990213326289756,1.35496815892254,1.68200753420297,0.00464892767486915,0.0542743083037553,0.35292431166273,0.586030161030878,1.15587296369776,0.547881779591917,0.255308800327668,0.112058926848163,0.131305240604842,0.749698618946567,0.261249818268944,0.363528689609138],[107.042795154914,0.272389320181145,3.72726609385188,4.62689075795212,0.012788337778423,0.149298555671384,0.970829323105667,1.61206027964662,3.17959213876776,1.50712115784037,0.7023071565608,0.308253323743573,0.361196363197688,2.06228185113157,0.718649795012978,null],[2.24508997832523e-23,0.601905973071706,0.0112087496035428,0.00108607009492853,0.909996618834323,0.930140639282038,0.405973894183611,0.169450257403212,0.0133017478807115,0.116350449657347,0.550864696828989,0.872516731595336,0.976201903841061,0.0174787030956122,0.733954755446139,null],[323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908],[0.120452094158206,0.000306511652602737,0.01258254717208,0.0208260137613194,1.43903385929364e-05,0.000504003758292783,0.00327733664472334,0.00725601517829142,0.0143116042935642,0.0203510268224171,0.0023708564680402,0.00138747342396678,0.00487732312518847,0.0278474978932037,0.00970410452945331,0.753931200780057],[12.0452094158206,0.0306511652602737,1.258254717208,2.08260137613194,0.00143903385929364,0.0504003758292783,0.327733664472334,0.725601517829142,1.43116042935642,2.03510268224171,0.23708564680402,0.138747342396678,0.487732312518847,2.78474978932037,0.970410452945331,75.3931200780057]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Df<\/th>\n      <th>Sum Sq<\/th>\n      <th>Mean Sq<\/th>\n      <th>F value<\/th>\n      <th>Pr(&gt;F)<\/th>\n      <th>TOT_SS<\/th>\n      <th>PROP_SS<\/th>\n      <th>VAR_EXP<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7,8]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

