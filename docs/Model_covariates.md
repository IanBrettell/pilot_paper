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
  dplyr::ungroup() %>% 
  # Convert `date` to factor
  dplyr::mutate(date = factor(date, levels = unique(df$date)))
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
# Separate by assay
df_control_of = df_control %>% 
  dplyr::filter(assay == "open_field")

df_control_no = df_control %>% 
  dplyr::filter(assay == "novel_object")

model_fit_of = lm(mean_speed ~ date + time + quadrant, data = df_control_of)

model_fit_no = lm(mean_speed ~ date + time + quadrant, data = df_control_no)

summary(model_fit_of)
#> 
#> Call:
#> lm(formula = mean_speed ~ date + time + quadrant, data = df_control_of)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.27802 -0.28921  0.07469  0.32925  1.29864 
#> 
#> Coefficients:
#>                Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)   1.3653659  0.3420269   3.992  0.00011 ***
#> date20190612  0.2601871  0.2286833   1.138  0.25736    
#> date20190614  0.8020814  0.2425613   3.307  0.00123 ** 
#> date20190615  0.4432846  0.2193359   2.021  0.04538 *  
#> date20190616  0.5258515  0.2141126   2.456  0.01540 *  
#> time          0.0001055  0.0002000   0.528  0.59876    
#> quadrantq2    0.1527986  0.1358049   1.125  0.26265    
#> quadrantq3    0.0925545  0.1358049   0.682  0.49678    
#> quadrantq4   -0.1949450  0.1358049  -1.435  0.15361    
#> ---
#> Signif. codes:  
#> 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.5599 on 127 degrees of freedom
#> Multiple R-squared:  0.1522,	Adjusted R-squared:  0.09878 
#> F-statistic:  2.85 on 8 and 127 DF,  p-value: 0.006037
summary(model_fit_no)
#> 
#> Call:
#> lm(formula = mean_speed ~ date + time + quadrant, data = df_control_no)
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.3553 -0.3706  0.0544  0.3735  1.4014 
#> 
#> Coefficients:
#>                Estimate Std. Error t value Pr(>|t|)   
#> (Intercept)   1.061e+00  3.346e-01   3.171   0.0019 **
#> date20190612  1.351e-01  2.237e-01   0.604   0.5470   
#> date20190614  5.042e-01  2.373e-01   2.125   0.0355 * 
#> date20190615  2.760e-01  2.146e-01   1.286   0.2007   
#> date20190616  2.668e-01  2.095e-01   1.274   0.2050   
#> time          9.476e-05  1.957e-04   0.484   0.6290   
#> quadrantq2    1.742e-01  1.329e-01   1.311   0.1921   
#> quadrantq3    1.112e-01  1.329e-01   0.837   0.4042   
#> quadrantq4   -2.323e-01  1.329e-01  -1.749   0.0828 . 
#> ---
#> Signif. codes:  
#> 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.5478 on 127 degrees of freedom
#> Multiple R-squared:  0.1211,	Adjusted R-squared:  0.06574 
#> F-statistic: 2.187 on 8 and 127 DF,  p-value: 0.03254
```

## DGE


```r
## All variables with interactions
model_dge_of = aov(mean_speed ~ date + time + quadrant + line,
                   data = df_dge %>% 
                     dplyr::filter(assay == "open_field"))

model_dge_no = aov(mean_speed ~ date + time + quadrant + line,
                   data = df_dge %>% 
                     dplyr::filter(assay == "novel_object"))

dge_vars_of = summary(model_dge_of)[[1]] %>% 
  dplyr::mutate(TOT_SS = sum(`Sum Sq`),
                PROP_SS = `Sum Sq` / TOT_SS,
                VAR_EXP = PROP_SS * 100)

dge_vars_no = summary(model_dge_no)[[1]] %>% 
  dplyr::mutate(TOT_SS = sum(`Sum Sq`),
                PROP_SS = `Sum Sq` / TOT_SS,
                VAR_EXP = PROP_SS * 100)

DT::datatable(dge_vars_of)
```

```{=html}
<div id="htmlwidget-3712586809b21a344f24" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-3712586809b21a344f24">{"x":{"filter":"none","vertical":false,"data":[["date       ","time       ","quadrant   ","line       ","Residuals  "],[5,1,3,4,361],[3.52206719652733,0.768244711192437,2.06070858371038,6.54617178682088,132.308667274163],[0.704413439305466,0.768244711192437,0.686902861236793,1.63654294670522,0.366506003529537],[1.92196971542567,2.09613131516009,1.87419265884258,4.46525549634913,null],[0.0899214087732978,0.148540253780587,0.133529547180944,0.0015580467801878,null],[145.205859552414,145.205859552414,145.205859552414,145.205859552414,145.205859552414],[0.0242556822939779,0.00529072802957465,0.014191635172729,0.0450820084464839,0.911179946057235],[2.42556822939779,0.529072802957465,1.4191635172729,4.50820084464839,91.1179946057235]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Df<\/th>\n      <th>Sum Sq<\/th>\n      <th>Mean Sq<\/th>\n      <th>F value<\/th>\n      <th>Pr(&gt;F)<\/th>\n      <th>TOT_SS<\/th>\n      <th>PROP_SS<\/th>\n      <th>VAR_EXP<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7,8]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
DT::datatable(dge_vars_no)
```

```{=html}
<div id="htmlwidget-86b60c3be7e7d6bf1b54" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-86b60c3be7e7d6bf1b54">{"x":{"filter":"none","vertical":false,"data":[["date       ","time       ","quadrant   ","line       ","Residuals  "],[5,1,3,4,361],[5.59676917785545,0.130591659750291,3.10898382355555,3.61834928856371,109.974573220059],[1.11935383557109,0.130591659750291,1.03632794118518,0.904587322140928,0.304638706980773],[3.67436510831087,0.428677173181847,3.40182622049598,2.96937750001027,null],[0.00295437320562422,0.513055689063896,0.0179164423999687,0.0195717513268915,null],[122.429267169784,122.429267169784,122.429267169784,122.429267169784,122.429267169784],[0.0457143075935748,0.0010666702723066,0.0253941226262838,0.0295546103657209,0.898270289142114],[4.57143075935748,0.10666702723066,2.53941226262838,2.95546103657209,89.8270289142114]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Df<\/th>\n      <th>Sum Sq<\/th>\n      <th>Mean Sq<\/th>\n      <th>F value<\/th>\n      <th>Pr(&gt;F)<\/th>\n      <th>TOT_SS<\/th>\n      <th>PROP_SS<\/th>\n      <th>VAR_EXP<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7,8]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-2121ff038ac523d747e8" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-2121ff038ac523d747e8">{"x":{"filter":"none","vertical":false,"data":[["assay                        ","time                         ","quadrant                     ","test_fish                    ","assay:time                   ","assay:quadrant               ","time:quadrant                ","assay:test_fish              ","time:test_fish               ","quadrant:test_fish           ","assay:time:quadrant          ","assay:time:test_fish         ","assay:quadrant:test_fish     ","time:quadrant:test_fish      ","assay:time:quadrant:test_fish","Residuals                    "],[1,1,3,4,1,3,3,4,4,12,3,4,12,12,12,670],[38.9131270547654,0.0990213326289756,4.06490447676763,6.72803013681188,0.00464892767486915,0.162822924911266,1.05877293498819,2.34412064412351,4.62349185479105,6.57458135510301,0.765926400983002,0.44823570739265,1.5756628872581,8.9963834273588,3.13499781922733,243.564222038123],[38.9131270547654,0.0990213326289756,1.35496815892254,1.68200753420297,0.00464892767486915,0.0542743083037553,0.35292431166273,0.586030161030878,1.15587296369776,0.547881779591917,0.255308800327668,0.112058926848163,0.131305240604842,0.749698618946567,0.261249818268944,0.363528689609138],[107.042795154914,0.272389320181145,3.72726609385188,4.62689075795212,0.012788337778423,0.149298555671384,0.970829323105667,1.61206027964662,3.17959213876776,1.50712115784037,0.7023071565608,0.308253323743573,0.361196363197688,2.06228185113157,0.718649795012978,null],[2.24508997832523e-23,0.601905973071706,0.0112087496035428,0.00108607009492853,0.909996618834323,0.930140639282038,0.405973894183611,0.169450257403212,0.0133017478807115,0.116350449657347,0.550864696828989,0.872516731595336,0.976201903841061,0.0174787030956122,0.733954755446139,null],[323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908,323.058949922908],[0.120452094158206,0.000306511652602737,0.01258254717208,0.0208260137613194,1.43903385929364e-05,0.000504003758292783,0.00327733664472334,0.00725601517829142,0.0143116042935642,0.0203510268224171,0.0023708564680402,0.00138747342396678,0.00487732312518847,0.0278474978932037,0.00970410452945331,0.753931200780057],[12.0452094158206,0.0306511652602737,1.258254717208,2.08260137613194,0.00143903385929364,0.0504003758292783,0.327733664472334,0.725601517829142,1.43116042935642,2.03510268224171,0.23708564680402,0.138747342396678,0.487732312518847,2.78474978932037,0.970410452945331,75.3931200780057]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Df<\/th>\n      <th>Sum Sq<\/th>\n      <th>Mean Sq<\/th>\n      <th>F value<\/th>\n      <th>Pr(&gt;F)<\/th>\n      <th>TOT_SS<\/th>\n      <th>PROP_SS<\/th>\n      <th>VAR_EXP<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7,8]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

