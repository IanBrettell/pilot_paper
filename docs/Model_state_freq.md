# Model proportion of time spent in HMM states


```r
library(tidyverse)
library(DT)
```

## Read in and clean data


```r
IN = "/hps/nobackup/birney/users/ian/pilot/hmm_out/0.08/dist_angle/14.csv"
N_STATES = 14
```


```r
# Read 

raw = readr::read_csv(IN)
#> Rows: 9152328 Columns: 15
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr (6): assay, ref_fish, test_fish, tank_side, quadrant...
#> dbl (9): date, time, frame, seconds, x, y, distance, ang...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

# Create line recode vector
line_vec = c("iCab", "HdrR", "HNI", "Kaga", "HO5")
names(line_vec) = c("icab", "hdr", "hni", "kaga", "ho5")

# Clean

df = raw %>% 
  # Get individual
  tidyr::unite(date, time, quadrant, fish,
               col = "indiv", 
               remove = F) %>% 
  # add `line` %>% 
  dplyr::mutate(line = dplyr::case_when(fish == "ref" ~ ref_fish,
                                        fish == "test" ~ test_fish)) %>% 
  # recode and order `assay` 
  dplyr::mutate(assay = stringr::str_replace(assay, pattern = "_", " "),
                assay = factor(assay, levels = c("open field", "novel object"))) %>% 
  # recode and order `line`
  dplyr::mutate(line = dplyr::recode(line, !!!line_vec),
                line = factor(line, levels = line_vec)) %>% 
  # convert `date` to factor
  dplyr::mutate(date = factor(date))

# Recode states by mean distance

rank_df = df %>% 
  dplyr::group_by(state) %>% 
  dplyr::summarise(mean_dist = mean(distance)) %>% 
  # rank
  dplyr::arrange(mean_dist) %>% 
  dplyr::mutate(rank = 1:nrow(.))

recode_vec = rank_df %>% 
  dplyr::pull(rank)
names(recode_vec) = rank_df %>% 
  dplyr::pull(state)

# Recode `state`

df = df %>% 
  dplyr::mutate(state_recode = dplyr::recode(state, !!!recode_vec),
                state_recode = factor(state_recode, levels = recode_vec))
```

## DGE


```r
# Get proportion of time each fish spent in each state
df_dge = df %>% 
  # remove iCab when paired with a different test fish
  dplyr::filter(!(fish == "ref" & test_fish != "icab")) %>% 
  ## count rows per fish per state
  dplyr::count(indiv, assay, line, date, time, quadrant, tank_side, state_recode) %>% 
  # add total row count per fish
  dplyr::add_count(indiv, assay, line, date, time, quadrant, tank_side, wt = n, name = "nn") %>% 
  # get proportion of time fish spent in each state
  dplyr::mutate(state_freq = n / nn)

# Split by assay

df_dge %>% 
  ggplot() + 
  geom_histogram(aes(state_freq, fill = state_recode),
                 bins = 40) +
  facet_grid(rows = vars(state_recode),
             cols = vars(assay)) +
  theme_bw() +
  scale_fill_viridis_d() +
  guides(fill = "none")
```

<img src="Model_state_freq_files/figure-html/unnamed-chunk-4-1.png" width="672" />


### Inverse-normalise


```r
# Add function
invnorm = function(x) {
  res = rank(x)
  res = qnorm(res/(length(res)+0.5))
  return(res)
}

df_dge = df_dge %>% 
  dplyr::group_by(assay, state_recode) %>% 
  dplyr::mutate(state_freq_invnorm = invnorm(state_freq)) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(indiv, assay, line, date, time, quadrant, tank_side, state_recode)

df_dge %>% 
  ggplot() + 
  geom_histogram(aes(state_freq_invnorm, fill = state_recode),
                 bins = 40) +
  facet_grid(rows = vars(state_recode),
             cols = vars(assay)) +
  theme_bw() +
  scale_fill_viridis_d() +
  guides(fill = "none")
```

<img src="Model_state_freq_files/figure-html/unnamed-chunk-5-1.png" width="672" />

### Calculate variance explained


```r
aov_dge = df_dge %>% 
  dplyr::group_by(assay, state_recode) %>% 
  tidyr::nest() %>%
  dplyr::mutate(model = purrr::map(data, ~aov(
    state_freq_invnorm ~ date + time + quadrant + tank_side + line,
    data = .))) %>%
  select(-data) %>% 
  dplyr::mutate(model_tidy = purrr::map(model, broom::tidy)) %>%
  tidyr::unnest(model_tidy) %>% 
  rstatix::adjust_pvalue(p.col = "p.value", method = "fdr") %>% 
  rstatix::add_significance(p.col = "p.value.adj")

DT::datatable(aov_dge %>% 
                dplyr::select(-model),
              options = list(pageLength = nrow(aov_dge)))
```

```{=html}
<div id="htmlwidget-e833a9afd9baa5a913ab" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-e833a9afd9baa5a913ab">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150","151","152","153","154","155","156","157","158","159","160","161","162","163","164","165","166","167","168"],["open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object"],["1","1","1","1","1","1","2","2","2","2","2","2","3","3","3","3","3","3","4","4","4","4","4","4","5","5","5","5","5","5","6","6","6","6","6","6","7","7","7","7","7","7","8","8","8","8","8","8","9","9","9","9","9","9","10","10","10","10","10","10","11","11","11","11","11","11","12","12","12","12","12","12","13","13","13","13","13","13","14","14","14","14","14","14","1","1","1","1","1","1","2","2","2","2","2","2","3","3","3","3","3","3","4","4","4","4","4","4","5","5","5","5","5","5","6","6","6","6","6","6","7","7","7","7","7","7","8","8","8","8","8","8","9","9","9","9","9","9","10","10","10","10","10","10","11","11","11","11","11","11","12","12","12","12","12","12","13","13","13","13","13","13","14","14","14","14","14","14"],["date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals","date","time","quadrant","tank_side","line","Residuals"],[5,1,3,1,4,348,5,1,3,1,4,350,5,1,3,1,4,350,5,1,3,1,4,360,5,1,3,1,4,360,5,1,3,1,4,360,5,1,3,1,4,359,5,1,3,1,4,360,5,1,3,1,4,358,5,1,3,1,4,358,5,1,3,1,4,354,5,1,3,1,4,353,5,1,3,1,4,346,5,1,3,1,4,325,5,1,3,1,4,359,5,1,3,1,4,360,5,1,3,1,4,360,5,1,3,1,4,360,5,1,3,1,4,359,5,1,3,1,4,359,5,1,3,1,4,356,5,1,3,1,4,359,5,1,3,1,4,352,5,1,3,1,4,353,5,1,3,1,4,350,5,1,3,1,4,347,5,1,3,1,4,337,5,1,3,1,4,291],[3.15442967690833,0.347478241091135,1.00647898909834,0.683697106829503,95.1007717682075,256.582483421782,6.6400463520967,0.348828078905289,4.04359358436237,0.000548872147699985,78.6931530701627,268.7735726784,4.30920014396655,0.845257349483785,1.54631339934977,0.340373353122194,100.767037299505,249.746025343984,28.1077761028714,1.44011492041859,1.83263110905378,2.06245096234054,65.3947472634475,270.436602013226,10.364079320235,0.776047260477553,8.62879497904958,0.571068804498988,27.5858894029797,321.348497963631,34.5948537007323,0.186227849095346,7.14701133731338,0.0864680425529636,20.6542798329875,306.599333269876,38.3490044049224,1.92172867922057,3.67628731453639,0.423356923305605,20.0798576382859,303.830260981124,47.5191343967888,1.02682938328681,4.10644797686908,0.254385651655601,21.4234642876771,294.73212275998,36.8054795022356,0.487847426289672,0.716305406420791,0.031048297980607,26.9560144337769,302.286916426198,21.2740896160196,0.184931756825796,11.0111736441928,0.0100786096025915,13.9144732582602,320.86406815149,2.57731434623361,0.000476284531867499,2.35988884195201,0.181892725278605,8.84152533541682,349.116318860553,37.0606026793912,0.366439176572604,1.77263126931216,0.00413436531444921,4.27598251192048,318.185234306218,19.4019777610247,0.238556323656118,3.99594597336709,1.1802351146164,16.4715889056354,313.4634696925,32.8854948134613,0.343406337181554,5.16726433151594,0.21496568434997,9.52017545849232,283.059424842503,7.67987852941367,0.113709521982583,6.54789652851052,0.206792138232191,63.2487006201052,290.481768517795,7.8731536929422,0.014475872934505,9.2106002540751,0.171710093453093,52.7052185386261,299.287170078477,7.39619436873009,0.245066119611779,6.89992727754472,0.100729358811276,71.2262664413859,283.407464764763,25.2524921621949,0.392147889628339,3.0908267682456,0.370063888581188,48.5013030063215,291.670181530093,7.71686385671121,5.1713097939688,4.78722547491769,1.75544859690215,7.06220770004084,341.788273621352,32.1889161719621,0.658117017399158,11.4700824744222,0.00237105985865739,15.2130736853189,308.744184898619,16.3408430075704,1.91459214532102,0.121962856262269,1.98257205948963,32.6474426266454,312.28138357608,44.0913821504822,0.539421011908874,6.23908834727133,0.565012255698005,7.4534746331961,309.354805905566,14.4164804250166,0.00609716869597615,0.272856051327491,0.0843342035794606,38.5647225887732,307.952807203204,31.766645255921,0.627853906900899,9.7433811298812,0.399089962974556,2.48161862669201,317.22345010558,10.2098878741739,0.961318962559193,3.85328081897689,1.40844124803531,5.16956384088122,337.635419089037,34.1567021915216,0.0788471548222035,6.69776976721109,1.34946688087623,2.43932554116829,311.01744672101,15.6584945523478,1.59962896287377,8.57633454400744,3.55237394547964,2.17558411781667,314.625818374483,14.0279227409415,1.20209483265935,2.08104310852102,5.0788514253136,6.87352447253534,268.398263814159],[0.630885935381665,0.347478241091135,0.335492996366113,0.683697106829503,23.7751929420519,0.737305986844201,1.32800927041934,0.348828078905289,1.34786452812079,0.000548872147699985,19.6732882675407,0.767924493366856,0.861840028793309,0.845257349483785,0.515437799783258,0.340373353122194,25.1917593248764,0.713560072411384,5.62155522057429,1.44011492041859,0.61087703635126,2.06245096234054,16.3486868158619,0.751212783370071,2.072815864047,0.776047260477553,2.87626499301653,0.571068804498988,6.89647235074492,0.892634716565642,6.91897074014645,0.186227849095346,2.38233711243779,0.0864680425529636,5.16356995824689,0.851664814638546,7.66980088098448,1.92172867922057,1.22542910484546,0.423356923305605,5.01996440957148,0.846323846744078,9.50382687935776,1.02682938328681,1.36881599228969,0.254385651655601,5.35586607191927,0.818700340999945,7.36109590044712,0.487847426289672,0.23876846880693,0.031048297980607,6.73900360844422,0.844376861525692,4.25481792320393,0.184931756825796,3.67039121473092,0.0100786096025915,3.47861831456506,0.89626834679187,0.515462869246723,0.000476284531867499,0.786629613984002,0.181892725278605,2.2103813338542,0.986204290566534,7.41212053587824,0.366439176572604,0.590877089770721,0.00413436531444921,1.06899562798012,0.901374601434046,3.88039555220495,0.238556323656118,1.33198199112236,1.1802351146164,4.11789722640885,0.905963785238438,6.57709896269226,0.343406337181554,1.72242144383865,0.21496568434997,2.38004386462308,0.870952076438472,1.53597570588273,0.113709521982583,2.18263217617017,0.206792138232191,15.8121751550263,0.80914141648411,1.57463073858844,0.014475872934505,3.0702000846917,0.171710093453093,13.1763046346565,0.831353250217992,1.47923887374602,0.245066119611779,2.29997575918157,0.100729358811276,17.8065666103465,0.787242957679897,5.05049843243897,0.392147889628339,1.0302755894152,0.370063888581188,12.1253257515804,0.810194948694704,1.54337277134224,5.1713097939688,1.59574182497256,1.75544859690215,1.76555192501021,0.952056472482874,6.43778323439243,0.658117017399158,3.8233608248074,0.00237105985865739,3.80326842132972,0.860011657099218,3.26816860151408,1.91459214532102,0.0406542854207563,1.98257205948963,8.16186065666134,0.877194897685619,8.81827643009645,0.539421011908874,2.07969611575711,0.565012255698005,1.86336865829903,0.861712551268985,2.88329608500331,0.00609716869597615,0.0909520171091637,0.0843342035794606,9.6411806471933,0.874865929554556,6.3533290511842,0.627853906900899,3.2477937099604,0.399089962974556,0.620404656673003,0.898650000299092,2.04197757483478,0.961318962559193,1.28442693965896,1.40844124803531,1.2923909602203,0.964672625968676,6.83134043830433,0.0788471548222035,2.2325899224037,1.34946688087623,0.609831385292073,0.896303881040376,3.13169891046957,1.59962896287377,2.85877818133581,3.55237394547964,0.543896029454168,0.933607769657219,2.80558454818829,1.20209483265935,0.693681036173674,5.0788514253136,1.71838111813384,0.922330803485083],[0.855663654762886,0.47128091632404,0.455025460734535,0.927290865704002,32.2460326733733,null,1.72934876005435,0.454247887544128,1.75520450221775,0.000714747546719773,25.6187794991223,null,1.20780304576296,1.18456368589591,0.72234675076676,0.477007285415993,35.30432867375,null,7.48330612180881,1.91705326679611,0.813187754354711,2.74549502883594,21.76305725592,null,2.32213225138949,0.869389511830044,3.22221950327317,0.639756435528455,7.72597370767595,null,8.12405376061347,0.218663312014813,2.79727079420192,0.101528255091366,6.06291333103664,null,9.06248938924649,2.27067769224951,1.44794349061515,0.500230408175684,5.93149351620418,null,11.6084315633116,1.25421882936149,1.67193773318556,0.310718878344303,6.54191259451256,null,8.71778495581519,0.577760296993676,0.282774765257664,0.0367706641374639,7.9810377516357,null,4.74725894140271,0.206335253819625,4.09519228015675,0.0112450803809675,3.88122410773127,null,0.522673521274796,0.00048294713014673,0.797633534459798,0.184437166841077,2.2413016805924,null,8.22313001063697,0.406533727475365,0.655528887579772,0.0045867337596063,1.18596155946639,null,4.28316850566348,0.263317725877236,1.47023756669457,1.30273983778036,4.54532211276533,null,7.55161981998776,0.394288441892031,1.97763056135308,0.246816891727286,2.7326921067295,null,1.89827844007401,0.140531085995666,2.69746688490396,0.255569834913094,19.5419179304077,null,1.89405735549305,0.017412421170795,3.69301507378079,0.206542878623579,15.8492249007284,null,1.8790118848513,0.311296680676598,2.9215577436982,0.127952060832831,22.6188960302988,null,6.23368294331602,0.484016705189442,1.27163911731993,0.456759066663393,14.9659360023354,null,1.62109372285162,5.43172589382489,1.676099970006,1.84384923335914,1.85446134345982,null,7.48569299177501,0.765241973136688,4.44570814040291,0.00275700897666303,4.4223451972244,null,3.72570407116683,2.18263028019481,0.0463457841900564,2.26012721314319,9.30450083350405,null,10.233431574265,0.625987182285445,2.41344531038161,0.655685303487747,2.16240166811423,null,3.29570050404336,0.00696926064897803,0.103961091678896,0.0963967171774536,11.0201807174068,null,7.06985928789814,0.698663447050503,3.61408079772932,0.444099441208178,0.69037406828745,null,2.1167570426126,0.996523521742814,1.33146406882771,1.46001991776332,1.3397197405934,null,7.62167896715444,0.0879692217004376,2.49088503311203,1.50559080399145,0.680384630918051,null,3.35440536406353,1.71338437281974,3.06207625327008,3.80499612464013,0.58257444628365,null,3.0418419699171,1.30332287300518,0.752095705307205,5.50653995954906,1.86308546959597,null],[0.511161782743186,0.492855633698571,0.713909669132647,0.336236524175467,7.1915000376114e-23,null,0.127121024258225,0.500769438628488,0.155453156334498,0.978686522845907,1.23514661809867e-18,null,0.304947397038901,0.277177742488119,0.53918537876444,0.490238795167534,8.89643816271305e-25,null,1.06864365312464e-06,0.16703942148241,0.487206896188696,0.0984000722900847,4.23519798622899e-16,null,0.0427382143164334,0.351749275063671,0.0227808488211416,0.424326468250441,5.55342083013867e-06,null,2.797734671354e-07,0.640343234447149,0.0400755714109849,0.750188099868159,9.9003995392974e-05,null,3.96486233929521e-08,0.1327216520085,0.228559436946264,0.479857753278543,0.000124406429256282,null,2.0631213389185e-10,0.263494833706004,0.172670997808013,0.577585806540874,4.31790313742151e-05,null,8.14832675482615e-08,0.447691827991959,0.837834080887352,0.848041439227003,3.58324624825054e-06,null,0.000326145833119661,0.64993035929542,0.0070605223355062,0.91560776919056,0.00422858023881566,null,0.759118912424825,0.982479435403555,0.495832743508238,0.667848683748236,0.0642095501994558,null,2.31428838322516e-07,0.524147452308506,0.579932645688118,0.946042549887349,0.316631977383493,null,0.000857658378511682,0.608177417031997,0.222393046258242,0.254502691594458,0.00136819664070412,null,1.00186107021785e-06,0.530494747673124,0.117167392155754,0.619661254196299,0.029102991604308,null,0.0938876809246076,0.707974967730609,0.0457303792575437,0.613490445281966,1.45512060062563e-14,null,0.0946001928212784,0.895092706037791,0.0121267633599868,0.649764562816888,5.83826202862662e-12,null,0.0972102355118656,0.577232371703986,0.0339902482767068,0.720774557785267,1.10562122184044e-16,null,1.46434270405006e-05,0.487057390382731,0.283873848101953,0.49957644756966,2.51693276662793e-11,null,0.153572321906318,0.0203274684905996,0.171771547070511,0.175353293370395,0.117921154075797,null,1.06545403602395e-06,0.382277873657228,0.00439744575358866,0.958153740782838,0.00167848986848748,null,0.00266750554790706,0.140459834178415,0.986749339633163,0.133629719105702,3.67228977504731e-07,null,3.49998411077103e-09,0.42935347165845,0.0664283107283512,0.418623719981101,0.0727706281689472,null,0.00636799195874394,0.933515666179548,0.957720288657745,0.756381083512856,1.97323524881579e-08,null,2.57264577013685e-06,0.403798572566568,0.0135001557551062,0.505585730858258,0.599012643951974,null,0.0629280135243458,0.318842827082929,0.263911302395327,0.227742700875489,0.254764438586312,null,8.22595757226226e-07,0.766952365111646,0.0600783911281052,0.220644971079411,0.605943275672154,null,0.00569285601990138,0.191438709990963,0.0282888782636979,0.0519284813933455,0.675466134801132,null,0.0107937363484691,0.254545700499051,0.521887728326401,0.0196161248176425,0.116948069656098,null],[0.638952228428983,0.638952228428983,0.713909669132647,0.638952228428983,3.5957500188057e-22,null,0.25908859389083,0.62596179828561,0.25908859389083,0.978686522845907,6.17573309049337e-18,null,0.508245661731501,0.508245661731501,0.53918537876444,0.53918537876444,4.44821908135653e-24,null,2.67160913281159e-06,0.208799276853013,0.487206896188696,0.164000120483475,2.11759899311449e-15,null,0.0712303571940557,0.424326468250441,0.0569521220528541,0.424326468250441,2.77671041506934e-05,null,1.398867335677e-06,0.750188099868159,0.0667926190183081,0.750188099868159,0.000247509988482435,null,1.9824311696476e-07,0.2212027533475,0.28569929618283,0.479857753278543,0.000311016073140705,null,1.03156066945925e-09,0.329368542132505,0.287784996346688,0.577585806540874,0.000107947578435538,null,4.07416337741308e-07,0.746153046653265,0.848041439227003,0.848041439227003,8.95811562062634e-06,null,0.00163072916559831,0.812412949119276,0.0117675372258437,0.91560776919056,0.0105714505970391,null,0.948898640531032,0.982479435403555,0.948898640531032,0.948898640531032,0.321047750997279,null,1.15714419161258e-06,0.724915807110147,0.724915807110147,0.946042549887349,0.724915807110147,null,0.00342049160176031,0.608177417031997,0.318128364493072,0.318128364493072,0.00342049160176031,null,5.00930535108927e-06,0.619661254196299,0.195278986926257,0.619661254196299,0.0727574790107701,null,0.156479468207679,0.707974967730609,0.114325948143859,0.707974967730609,7.27560300312813e-14,null,0.157666988035464,0.895092706037791,0.030316908399967,0.81220570352111,2.91913101431331e-11,null,0.162017059186443,0.720774557785267,0.0849756206917669,0.720774557785267,5.52810610920219e-16,null,3.66085676012516e-05,0.49957644756966,0.473123080169922,0.49957644756966,1.25846638331396e-10,null,0.175353293370395,0.101637342452998,0.175353293370395,0.175353293370395,0.175353293370395,null,5.32727018011975e-06,0.477847342071535,0.00732907625598111,0.958153740782838,0.0041962246712187,null,0.00666876386976765,0.175574792723019,0.986749339633163,0.175574792723019,1.83614488752365e-06,null,1.74999205538551e-08,0.42935347165845,0.121284380281579,0.42935347165845,0.121284380281579,null,0.0159199798968599,0.957720288657745,0.957720288657745,0.957720288657745,9.86617624407894e-08,null,1.28632288506842e-05,0.599012643951974,0.0337503893877654,0.599012643951974,0.599012643951974,null,0.314640067621729,0.318842827082929,0.318842827082929,0.318842827082929,0.318842827082929,null,4.11297878613113e-06,0.766952365111646,0.150195977820263,0.367741618465684,0.757429094590193,null,0.0284642800995069,0.239298387488704,0.0707221956592448,0.0865474689889092,0.675466134801132,null,0.0490403120441062,0.318182125623813,0.521887728326401,0.0490403120441062,0.19491344942683,null],["ns","ns","ns","ns","****","","ns","ns","ns","ns","****","","ns","ns","ns","ns","****","","****","ns","ns","ns","****","","ns","ns","ns","ns","****","","****","ns","ns","ns","***","","****","ns","ns","ns","***","","****","ns","ns","ns","***","","****","ns","ns","ns","****","","**","ns","*","ns","*","","ns","ns","ns","ns","ns","","****","ns","ns","ns","ns","","**","ns","ns","ns","**","","****","ns","ns","ns","ns","","ns","ns","ns","ns","****","","ns","ns","*","ns","****","","ns","ns","ns","ns","****","","****","ns","ns","ns","****","","ns","ns","ns","ns","ns","","****","ns","**","ns","**","","**","ns","ns","ns","****","","****","ns","ns","ns","ns","","*","ns","ns","ns","****","","****","ns","*","ns","ns","","ns","ns","ns","ns","ns","","****","ns","ns","ns","ns","","*","ns","ns","ns","ns","","*","ns","ns","*","ns",""]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>assay<\/th>\n      <th>state_recode<\/th>\n      <th>term<\/th>\n      <th>df<\/th>\n      <th>sumsq<\/th>\n      <th>meansq<\/th>\n      <th>statistic<\/th>\n      <th>p.value<\/th>\n      <th>p.value.adj<\/th>\n      <th>p.value.adj.signif<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"pageLength":168,"columnDefs":[{"className":"dt-right","targets":[4,5,6,7,8,9]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false,"lengthMenu":[10,25,50,100,168]}},"evals":[],"jsHooks":[]}</script>
```

### Save to Google Drive

To `aov_state_freq` here: https://docs.google.com/spreadsheets/d/1_l72BZkmWyNAOfCUI8WGP4UfQuIPQtPZZmlRjQffvEs


```r
# Open field
dge_tbl_of = aov_dge %>% 
  dplyr::filter(assay == "open field") %>% 
  # add variance explained
  dplyr::group_by(assay, state_recode) %>% 
  dplyr::mutate(tot_ss = sum(sumsq)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(var_expl_perc = (sumsq / tot_ss) * 100 ) %>%
  # select and rename key columns
  dplyr::select(State = state_recode,
                Variable = term,
                Statistic = statistic,
                `p-value` = p.value,
                `p-value FDR-adj` = p.value.adj,
                `Significance (p-value FDR-adj)` = p.value.adj.signif,
                `Variance explained (%)` = var_expl_perc) %>% 
  # show only 2 decimals
  dplyr::mutate(dplyr::across(c("Statistic", 
                                #"p-value", 
                                "Variance explained (%)"),
                              ~ format(round(.x, 2), nsmall = 2)))

## Write to google sheet
googlesheets4::write_sheet(
  data = dge_tbl_of,
  ss = "https://docs.google.com/spreadsheets/d/1_l72BZkmWyNAOfCUI8WGP4UfQuIPQtPZZmlRjQffvEs",
  sheet = "DGE_OF")
#> ! Using an auto-discovered, cached token.
#>   To suppress this message, modify your code or options to
#>   clearly consent to the use of a cached token.
#>   See gargle's "Non-interactive auth" vignette for more
#>   details:
#>   <https://gargle.r-lib.org/articles/non-interactive-auth.html>
#> ℹ The googlesheets4 package is using a cached token for
#>   'brettell@ebi.ac.uk'.
#> Auto-refreshing stale OAuth token.
#> ✔ Writing to "aov_state_freq".
#> ✔ Writing to sheet 'DGE_OF'.

# Novel object
dge_tbl_no = aov_dge %>% 
  dplyr::filter(assay == "novel object") %>% 
  # add variance explained
  dplyr::group_by(assay, state_recode) %>% 
  dplyr::mutate(tot_ss = sum(sumsq)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(var_expl_perc = (sumsq / tot_ss) * 100 ) %>%
  # select and rename key columns
  dplyr::select(State = state_recode,
                Variable = term,
                Statistic = statistic,
                `p-value` = p.value,
                `p-value FDR-adj` = p.value.adj,
                `Significance (p-value FDR-adj)` = p.value.adj.signif,
                `Variance explained (%)` = var_expl_perc) %>% 
  # show only 2 decimals
  dplyr::mutate(dplyr::across(c("Statistic", 
                                #"p-value", 
                                "Variance explained (%)"),
                              ~ format(round(.x, 2), nsmall = 2)))

## Write to google sheet
googlesheets4::write_sheet(
  data = dge_tbl_no,
  ss = "https://docs.google.com/spreadsheets/d/1_l72BZkmWyNAOfCUI8WGP4UfQuIPQtPZZmlRjQffvEs",
  sheet = "DGE_NO")
#> ✔ Writing to "aov_state_freq".
#> ✔ Writing to sheet 'DGE_NO'.
```

## SGE


```r
# Get proportion of time each fish spent in each state
df_sge = df %>% 
  # take all iCab fishes
  dplyr::filter(line == "iCab") %>% 
  ## count rows per fish per state
  dplyr::count(indiv, assay, test_fish, date, time, quadrant, tank_side, state_recode) %>% 
  # add total row count per fish
  dplyr::add_count(indiv, assay, test_fish, date, time, quadrant, tank_side, wt = n, name = "nn") %>% 
  # get proportion of time fish spent in each state
  dplyr::mutate(state_freq = n / nn)

# Split by assay

df_sge %>% 
  ggplot() + 
  geom_histogram(aes(state_freq, fill = state_recode),
                 bins = 40) +
  facet_grid(rows = vars(state_recode),
             cols = vars(assay)) +
  theme_bw() +
  scale_fill_viridis_d(option = "inferno") +
  guides(fill = "none")
```

<img src="Model_state_freq_files/figure-html/unnamed-chunk-8-1.png" width="672" />


### Inverse-normalise


```r
df_sge = df_sge %>% 
  dplyr::group_by(assay, state_recode) %>% 
  dplyr::mutate(state_freq_invnorm = invnorm(state_freq)) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(indiv, assay, test_fish, date, time, quadrant, tank_side, state_recode)

df_sge %>% 
  ggplot() + 
  geom_histogram(aes(state_freq_invnorm, fill = state_recode),
                 bins = 40) +
  facet_grid(rows = vars(state_recode),
             cols = vars(assay)) +
  theme_bw() +
  scale_fill_viridis_d(option = "inferno") +
  guides(fill = "none")
```

<img src="Model_state_freq_files/figure-html/unnamed-chunk-9-1.png" width="672" />

### Calculate variance explained


```r
aov_sge = df_sge %>% 
  dplyr::group_by(assay, state_recode) %>% 
  tidyr::nest() %>%
  dplyr::mutate(model = purrr::map(data, ~aov(
    state_freq_invnorm ~ date + time + quadrant + tank_side + test_fish,
    data = .))) %>%
  select(-data) %>% 
  dplyr::mutate(model_tidy = purrr::map(model, broom::tidy)) %>%
  tidyr::unnest(model_tidy) %>% 
  rstatix::adjust_pvalue(p.col = "p.value", method = "fdr") %>% 
  rstatix::add_significance(p.col = "p.value.adj")

DT::datatable(aov_sge %>% 
                dplyr::select(-model),
              options = list(pageLength = nrow(aov_sge)))
```

```{=html}
<div id="htmlwidget-8a3e16949627605f642f" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-8a3e16949627605f642f">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150","151","152","153","154","155","156","157","158","159","160","161","162","163","164","165","166","167","168"],["open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","open field","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object","novel object"],["1","1","1","1","1","1","2","2","2","2","2","2","3","3","3","3","3","3","4","4","4","4","4","4","5","5","5","5","5","5","6","6","6","6","6","6","7","7","7","7","7","7","8","8","8","8","8","8","9","9","9","9","9","9","10","10","10","10","10","10","11","11","11","11","11","11","12","12","12","12","12","12","13","13","13","13","13","13","14","14","14","14","14","14","1","1","1","1","1","1","2","2","2","2","2","2","3","3","3","3","3","3","4","4","4","4","4","4","5","5","5","5","5","5","6","6","6","6","6","6","7","7","7","7","7","7","8","8","8","8","8","8","9","9","9","9","9","9","10","10","10","10","10","10","11","11","11","11","11","11","12","12","12","12","12","12","13","13","13","13","13","13","14","14","14","14","14","14"],["date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals","date","time","quadrant","tank_side","test_fish","Residuals"],[5,1,3,1,4,358,5,1,3,1,4,359,5,1,3,1,4,359,5,1,3,1,4,360,5,1,3,1,4,360,5,1,3,1,4,360,5,1,3,1,4,360,5,1,3,1,4,360,5,1,3,1,4,359,5,1,3,1,4,360,5,1,3,1,4,358,5,1,3,1,4,357,5,1,3,1,4,351,5,1,3,1,4,334,5,1,3,1,4,360,5,1,3,1,4,360,5,1,3,1,4,360,5,1,3,1,4,360,5,1,3,1,4,360,5,1,3,1,4,360,5,1,3,1,4,359,5,1,3,1,4,359,5,1,3,1,4,358,5,1,3,1,4,358,5,1,3,1,4,354,5,1,3,1,4,352,5,1,3,1,4,348,5,1,3,1,4,311],[3.25290968899676,0.0652265457284949,4.90851195912554,2.87726687289756,31.3084569942097,324.832284218777,2.42362437740144,0.0950013391639979,4.48586447554095,3.20735649285021,33.3936902087185,324.3724902352,2.58850216319511,0.0316382338210911,4.10620081227298,2.50458474784877,29.2896835813389,329.146724776042,21.8345615695885,3.80184424175251,0.917112083735977,3.14564315985691,1.60211785719192,337.97556607135,23.816256809334,12.3148362440912,4.13214816969663,0.844136563423729,14.5386995614366,313.62988216836,20.701088351733,2.84660177498955,7.76183379190151,2.11636997031362,2.02349482834192,333.825462962468,46.3141528152969,12.6198476759238,0.793695671254986,0.912976988729216,5.09519842365965,303.542448316773,43.8255566142029,6.31627948417343,4.66278800031271,0.619956992907118,5.87158526738869,307.965911275883,34.6228450387267,2.58336053990702,0.673333616472155,0.253963158681997,3.79816373931059,326.348779653772,24.6707509786805,0.984891177928369,4.49912806939216,0.203602436502357,7.46389355822884,331.43212206475,14.4896250467408,3.04886355557776,4.71764053551041,1.20582853623104,9.91734078285153,333.829585442238,42.9153747191465,9.22615666702598,1.5299400240728,0.0246525485180906,10.917867118947,301.513480624104,14.3705060329228,7.90887120780442,2.47740880901984,1.03528062425861,7.87139712769898,326.358042667001,22.4718314190128,3.01530382022284,2.22210946179224,0.00978610530341082,6.7229615629175,306.621728954857,16.199331113407,0.296253753519537,6.69966572837077,2.7059836455035,29.295387469833,314.053265681763,15.1960337664672,1.12051363831315,5.7617579497087,4.18764608609857,34.5630365448587,308.443154714321,18.9440652194448,0.162720986090302,5.71041133392056,3.07335604198609,32.2487000624419,309.138793812209,26.255890417228,0.235428115488087,2.76500277920823,1.047304490053,1.06564461370872,337.907706592159,12.7854442814239,13.6744908214727,0.277833003150079,0.783120687324545,5.23979547152759,336.516023214044,37.080247086131,1.21410832342639,5.06093341020619,0.858148584933949,17.5531462409295,307.50759619116,10.9958756913357,8.06181582966964,0.575280548525423,0.688918050987581,1.58842590497573,346.370833097009,51.9165478726221,2.32628356956149,2.01414837210212,0.854482105593943,10.3295857251616,300.810339313009,5.95341462572355,0.517202139192336,1.8686667767198,0.0432084151656896,5.28731676886602,353.613419138529,35.7784059153713,2.73904640937633,3.20043237794087,0.583008358915127,13.9503223939526,311.019931839673,8.40646364119297,3.34323738360616,6.05686347921536,0.472700028076471,9.59498058514185,335.41611649869,38.1875726794605,3.78846948415161,2.47449804605547,1.8988771980739,8.29132613267755,306.110496881954,13.4345083558294,5.91153375248292,8.87512923661305,1.08417824273781,6.45421664673391,321.322523240678,15.3670266106918,4.14672572243424,1.24550071759405,6.65509928569545,3.83745592687805,285.995277630058],[0.650581937799352,0.0652265457284949,1.63617065304185,2.87726687289756,7.82711424855242,0.90735274921446,0.484724875480289,0.0950013391639979,1.49528815851365,3.20735649285021,8.34842255217962,0.903544541045125,0.517700432639023,0.0316382338210911,1.36873360409099,2.50458474784877,7.32242089533473,0.916843244501508,4.3669123139177,3.80184424175251,0.305704027911992,3.14564315985691,0.40052946429798,0.938821016864861,4.7632513618668,12.3148362440912,1.37738272323221,0.844136563423729,3.63467489035914,0.871194117134332,4.14021767034659,2.84660177498955,2.58727793063384,2.11636997031362,0.50587370708548,0.927292952673523,9.26283056305938,12.6198476759238,0.264565223751662,0.912976988729216,1.27379960591491,0.843173467546591,8.76511132284059,6.31627948417343,1.5542626667709,0.619956992907118,1.46789631684717,0.855460864655229,6.92456900774534,2.58336053990702,0.224444538824052,0.253963158681997,0.949540934827647,0.909049525497972,4.9341501957361,0.984891177928369,1.49970935646405,0.203602436502357,1.86597338955721,0.920644783513194,2.89792500934816,3.04886355557776,1.57254684517014,1.20582853623104,2.47933519571288,0.932484875536977,8.58307494382929,9.22615666702598,0.509980008024267,0.0246525485180906,2.72946677973676,0.844575575977883,2.87410120658456,7.90887120780442,0.825802936339946,1.03528062425861,1.96784928192475,0.929794993353279,4.49436628380256,3.01530382022284,0.740703153930748,0.00978610530341082,1.68074039072938,0.918029128607358,3.23986622268141,0.296253753519537,2.23322190945692,2.7059836455035,7.32384686745825,0.87237018244934,3.03920675329343,1.12051363831315,1.92058598323623,4.18764608609857,8.64075913621468,0.856786540873113,3.78881304388895,0.162720986090302,1.90347044464019,3.07335604198609,8.06217501561046,0.858718871700581,5.2511780834456,0.235428115488087,0.92166759306941,1.047304490053,0.26641115342718,0.938632518311553,2.55708885628478,13.6744908214727,0.0926110010500265,0.783120687324545,1.3099488678819,0.934766731150122,7.41604941722621,1.21410832342639,1.68697780340206,0.858148584933949,4.38828656023238,0.854187767197666,2.19917513826713,8.06181582966964,0.191760182841808,0.688918050987581,0.397106476243932,0.96482126210866,10.3833095745244,2.32628356956149,0.671382790700708,0.854482105593943,2.58239643129041,0.837911808671335,1.19068292514471,0.517202139192336,0.622888925573268,0.0432084151656896,1.3218291922165,0.987746980833881,7.15568118307425,2.73904640937633,1.06681079264696,0.583008358915127,3.48758059848815,0.868770759328696,1.68129272823859,3.34323738360616,2.01895449307179,0.472700028076471,2.39874514628546,0.947503153951103,7.63751453589209,3.78846948415161,0.824832682018488,1.8988771980739,2.07283153316939,0.869632093414641,2.68690167116588,5.91153375248292,2.95837641220435,1.08417824273781,1.61355416168348,0.923340584024937,3.07340532213837,4.14672572243424,0.415166905864683,6.65509928569545,0.959363981719513,0.919598963440702],[0.717011039380873,0.0718866458331281,1.80323546102479,3.17105654376267,8.62631898710699,null,0.536470371366098,0.1051429507328,1.65491361032856,3.54974917909445,9.23963587066006,null,0.564655338569353,0.0345077896476108,1.49287635841739,2.73174805275522,7.98655706877784,null,4.65148546471696,4.04959430334074,0.325625462596559,3.35063138058158,0.426630270416746,null,5.46749716071871,14.135582417153,1.58102849427279,0.968941992170924,4.17206087469332,null,4.46484323903221,3.06979770177523,2.79014083216564,2.28230999083063,0.545538177149919,null,10.9856760436399,14.967083478853,0.313773184207843,1.08278666712051,1.5107206938346,null,10.2460693235489,7.38348151872388,1.81687173661335,0.724705265339023,1.71591288099283,null,7.61737266619452,2.84182595936329,0.246900232087028,0.279372192117781,1.04454257792775,null,5.35945055475936,1.06978412908615,1.62897719437582,0.22115200145425,2.02681145103177,null,3.10774478532895,3.26961180343225,1.68640466609668,1.29313468546785,2.65884762397368,null,10.1625895751147,10.9240154812002,0.603829926568493,0.029189274730739,3.23176144014877,null,3.09111280134903,8.50603763662063,0.88815592925674,1.11345041827436,2.11643351060405,null,4.89566849651108,3.2845404641975,0.806840579290101,0.0106599071842701,1.8308137926724,null,3.7138663010981,0.339596377179869,2.55994754793973,3.10187544226461,8.39534295738444,null,3.54721579799364,1.30780957083113,2.24161549185758,4.88761889493634,10.0850780492061,null,4.41216930097938,0.18949273324814,2.2166398193516,3.5790013976282,9.38860817119854,null,5.5944983590507,0.250820327332771,0.981925913591125,1.11577690908997,0.283829025981766,null,2.7355368682608,14.6287735386644,0.0990739164797611,0.837771244127405,1.4013644519314,null,8.68198972405803,1.42135999841056,1.97494961668268,1.00463694036415,5.13737930786463,null,2.27936015160024,8.35576095415868,0.198752028352596,0.71403697041452,0.411585535932364,null,12.3918883432244,2.77628689020417,0.801257105763087,1.0197757052128,3.08194299753966,null,1.20545336837123,0.523618041085482,0.6306158739634,0.0437444163374834,1.33822650725855,null,8.23655850089095,3.15278383850397,1.22795430347044,0.671072722757867,4.01438533817946,null,1.77444552160864,3.52847098151051,2.13081559111726,0.498890189552726,2.53164871935535,null,8.78246627939307,4.35640486688578,0.948484638808296,2.18354084727702,2.38357294868263,null,2.90997895863452,6.40233284961215,3.20399261484693,1.17419104228232,1.74751786025675,null,3.34211481778878,4.50927620331277,0.451465173809382,7.23695822882971,1.04324169541256,null],[0.610997870240196,0.788763793803917,0.146200801497491,0.0758018236459194,1.1758933162741e-06,null,0.748637569947325,0.745931906115061,0.176425623124173,0.0603622990754019,4.07977106367334e-07,null,0.727104203495294,0.852735750538648,0.21614815525178,0.0992457809581984,3.54370041919834e-06,null,0.000397158660072426,0.0449263238393462,0.806839669996746,0.0680058070094789,0.789407541613645,null,7.27101080613777e-05,0.000198425266954733,0.193597585792333,0.325605270352408,0.00257446142656329,null,0.000584346506291467,0.0806107999972839,0.0404554169180545,0.131734076050591,0.702399751487864,null,7.39888376716735e-10,0.000129845069771199,0.815428491559777,0.298773324215287,0.198452859771479,null,3.39703623107868e-09,0.00690061369766515,0.143672013091979,0.395170718579796,0.145831640094812,null,8.08850555823182e-07,0.0927086792451851,0.863528411686472,0.597439600026599,0.38405958066864,null,9.11015086985351e-05,0.301688502249001,0.182279586855347,0.638448315417117,0.0901417975682084,null,0.00925218353554147,0.0714133371234412,0.169555712725325,0.2562316141058,0.0326645008420339,null,4.08187572932548e-09,0.00104550589147936,0.612883711328382,0.864439531806072,0.0126518339772367,null,0.00958356238132449,0.0037673780970686,0.447381021993354,0.292058900309825,0.0783239397723372,null,0.000244840496045301,0.0708330175040772,0.490774942096932,0.917829039823973,0.122493845754579,null,0.0027275882536164,0.560427121241477,0.0548083272435089,0.0790510458248486,1.74568408613279e-06,null,0.00382505588636675,0.253551498908996,0.0831192146483545,0.0276780975483724,9.52931885175952e-08,null,0.000651502743447812,0.663600008889325,0.0858625730181005,0.0593166427334099,3.15077873765563e-07,null,5.57698924249621e-05,0.616804107775373,0.401305804779578,0.291537946170223,0.888410803197771,null,0.0192733647769105,0.00015425633047318,0.960495864792288,0.360647648384339,0.232985726654627,null,8.73160643455464e-08,0.233964955148826,0.117353335785165,0.316864326212498,0.000490490230618667,null,0.0463479793345519,0.00407850327333165,0.897212534526697,0.398669677291914,0.800294708805012,null,4.19610193907146e-11,0.0965419623048581,0.493803941966109,0.313252353189194,0.0162389132547677,null,0.306018516512677,0.469774939819327,0.595655522916514,0.834448807317285,0.255259369613487,null,2.2225970925064e-07,0.0766469852154807,0.299355824972164,0.413222855161023,0.00337140471421031,null,0.117318199310516,0.0611450512632606,0.0960076829328991,0.480452772524604,0.0402167055862224,null,7.24361062046331e-08,0.0375889109745523,0.417277251174342,0.140387375655113,0.0511194755982198,null,0.0137271451315802,0.011837430442711,0.0233824364435104,0.279291949872022,0.139059617875067,null,0.00590007211284791,0.0345007456448566,0.716440611669289,0.00752752093948407,0.384978551654001,null],[0.763747337800245,0.788763793803917,0.243668002495819,0.189504559114798,5.87946658137052e-06,null,0.748637569947325,0.748637569947325,0.294042705206955,0.150905747688505,2.03988553183667e-06,null,0.852735750538648,0.852735750538648,0.360246925419634,0.248114452395496,1.77185020959917e-05,null,0.00198579330036213,0.112315809598366,0.806839669996746,0.113343011682465,0.806839669996746,null,0.000363550540306889,0.000496063167386832,0.241996982240416,0.325605270352408,0.00429076904427215,null,0.00292173253145733,0.134351333328807,0.101138542295136,0.164667595063239,0.702399751487864,null,3.69944188358368e-09,0.000324612674427998,0.815428491559777,0.373466655269108,0.330754766285798,null,1.69851811553934e-08,0.0172515342441629,0.182289550118515,0.395170718579796,0.182289550118515,null,4.04425277911591e-06,0.231771698112963,0.863528411686472,0.746799500033249,0.640099301114399,null,0.000455507543492676,0.377110627811252,0.303799311425578,0.638448315417117,0.225354493920521,null,0.0462609176777074,0.119022228539069,0.211944640906656,0.2562316141058,0.0816612521050846,null,2.04093786466274e-08,0.0026137647286984,0.766104639160477,0.864439531806072,0.0210863899620611,null,0.0239589059533112,0.018836890485343,0.447381021993354,0.365073625387281,0.130539899620562,null,0.00122420248022651,0.177082543760193,0.613468677621165,0.917829039823973,0.204156409590965,null,0.006818970634041,0.560427121241477,0.0913472120725149,0.0988138072810608,8.72842043066393e-06,null,0.00956263971591688,0.253551498908996,0.103899018310443,0.0461301625806207,4.76465942587976e-07,null,0.00162875685861953,0.663600008889325,0.107328216272626,0.0988610712223498,1.57538936882782e-06,null,0.00027884946212481,0.771005134719217,0.668843007965963,0.668843007965963,0.888410803197771,null,0.0481834119422763,0.000771281652365899,0.960495864792288,0.450809560480424,0.388309544424378,null,4.36580321727732e-07,0.292456193936033,0.195588892975275,0.316864326212498,0.00122622557654667,null,0.11586994833638,0.0203925163666583,0.897212534526697,0.66444946215319,0.897212534526697,null,2.09805096953573e-10,0.160903270508097,0.493803941966109,0.391565441486492,0.0405972831369191,null,0.744569403645643,0.744569403645643,0.744569403645643,0.834448807317285,0.744569403645643,null,1.1112985462532e-06,0.127744975359135,0.374194781215205,0.413222855161023,0.00842851178552577,null,0.146647749138145,0.146647749138145,0.146647749138145,0.480452772524604,0.146647749138145,null,3.62180531023165e-07,0.085199125997033,0.417277251174342,0.175484219568891,0.085199125997033,null,0.0343178628289505,0.0343178628289505,0.0389707274058507,0.279291949872022,0.173824522343833,null,0.0188188023487102,0.0575012427414276,0.716440611669289,0.0188188023487102,0.481223189567501,null],["ns","ns","ns","ns","****","","ns","ns","ns","ns","****","","ns","ns","ns","ns","****","","**","ns","ns","ns","ns","","***","***","ns","ns","**","","**","ns","ns","ns","ns","","****","***","ns","ns","ns","","****","*","ns","ns","ns","","****","ns","ns","ns","ns","","***","ns","ns","ns","ns","","*","ns","ns","ns","ns","","****","**","ns","ns","*","","*","*","ns","ns","ns","","**","ns","ns","ns","ns","","**","ns","ns","ns","****","","**","ns","ns","*","****","","**","ns","ns","ns","****","","***","ns","ns","ns","ns","","*","***","ns","ns","ns","","****","ns","ns","ns","**","","ns","*","ns","ns","ns","","****","ns","ns","ns","*","","ns","ns","ns","ns","ns","","****","ns","ns","ns","**","","ns","ns","ns","ns","ns","","****","ns","ns","ns","ns","","*","*","*","ns","ns","","*","ns","ns","*","ns",""]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>assay<\/th>\n      <th>state_recode<\/th>\n      <th>term<\/th>\n      <th>df<\/th>\n      <th>sumsq<\/th>\n      <th>meansq<\/th>\n      <th>statistic<\/th>\n      <th>p.value<\/th>\n      <th>p.value.adj<\/th>\n      <th>p.value.adj.signif<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"pageLength":168,"columnDefs":[{"className":"dt-right","targets":[4,5,6,7,8,9]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false,"lengthMenu":[10,25,50,100,168]}},"evals":[],"jsHooks":[]}</script>
```

### Save to Google Drive

To `aov_state_freq` here: https://docs.google.com/spreadsheets/d/1_l72BZkmWyNAOfCUI8WGP4UfQuIPQtPZZmlRjQffvEs


```r
# Open field
sge_tbl_of = aov_sge %>% 
  dplyr::filter(assay == "open field") %>% 
  dplyr::select(-model) %>% 
  # add variance explained
  dplyr::group_by(assay, state_recode) %>% 
  dplyr::mutate(tot_ss = sum(sumsq)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(var_expl_perc = (sumsq / tot_ss) * 100 ) %>%
  # select and rename key columns
  dplyr::select(State = state_recode,
                Variable = term,
                Statistic = statistic,
                `p-value` = p.value,
                `p-value FDR-adj` = p.value.adj,
                `Significance (p-value FDR-adj)` = p.value.adj.signif,
                `Variance explained (%)` = var_expl_perc) %>% 
  # show only 2 decimals
  dplyr::mutate(dplyr::across(c("Statistic", 
                                #"p-value", 
                                "Variance explained (%)"),
                              ~ format(round(.x, 2), nsmall = 2)))

## Write to Google sheet
googlesheets4::write_sheet(
  data = sge_tbl_of,
  ss = "https://docs.google.com/spreadsheets/d/1_l72BZkmWyNAOfCUI8WGP4UfQuIPQtPZZmlRjQffvEs",
  sheet = "SGE_OF")
#> ✔ Writing to "aov_state_freq".
#> ✔ Writing to sheet 'SGE_OF'.

# Novel object
sge_tbl_no = aov_sge %>% 
  dplyr::filter(assay == "novel object") %>% 
  dplyr::select(-model) %>% 
  # add variance explained
  dplyr::group_by(assay, state_recode) %>% 
  dplyr::mutate(tot_ss = sum(sumsq)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(var_expl_perc = (sumsq / tot_ss) * 100 ) %>%
  # select and rename key columns
  dplyr::select(State = state_recode,
                Variable = term,
                Statistic = statistic,
                `p-value` = p.value,
                `p-value FDR-adj` = p.value.adj,
                `Significance (p-value FDR-adj)` = p.value.adj.signif,
                `Variance explained (%)` = var_expl_perc) %>% 
  # show only 2 decimals
  dplyr::mutate(dplyr::across(c("Statistic", 
                                #"p-value", 
                                "Variance explained (%)"),
                              ~ format(round(.x, 2), nsmall = 2)))

## Write to Google sheet
googlesheets4::write_sheet(
  data = sge_tbl_no,
  ss = "https://docs.google.com/spreadsheets/d/1_l72BZkmWyNAOfCUI8WGP4UfQuIPQtPZZmlRjQffvEs",
  sheet = "SGE_NO")
#> ✔ Writing to "aov_state_freq".
#> ✔ Writing to sheet 'SGE_NO'.
```

## Write final table with only significant variables

### DGE


```r
final_dge = dplyr::bind_rows(
  list(
    "open field" = dge_tbl_of,
    "novel_object" = dge_tbl_no
    ),
  .id = "Assay") %>% 
  # filter for significant rows
  dplyr::filter(`p-value FDR-adj` < 0.05) %>% 
  # remove p-value
  dplyr::select(-`p-value`) %>% 
  # convert adjusted p-value to character in scientific notation
  dplyr::mutate(`p-value FDR-adj` = as.character(scales::scientific(`p-value FDR-adj`, digits = 3))) %>% 
  # remove underscores from values
  dplyr::mutate(dplyr::across(c("Assay", "Variable"),
                              ~stringr::str_replace(., pattern = "_", " "))) %>% 
  # rename columns
  dplyr::rename("p-value (FDR-adjusted)" = "p-value FDR-adj",
                "Significance" = "Significance (p-value FDR-adj)")

## Write to Google sheet
googlesheets4::write_sheet(
  data = final_dge,
  ss = "https://docs.google.com/spreadsheets/d/1_l72BZkmWyNAOfCUI8WGP4UfQuIPQtPZZmlRjQffvEs",
  sheet = "DGE_FINAL")
#> ✔ Writing to "aov_state_freq".
#> ✔ Writing to sheet 'DGE_FINAL'.
```

### SGE


```r
final_sge = dplyr::bind_rows(
  list(
    "open field" = sge_tbl_of,
    "novel_object" = sge_tbl_no
    ),
  .id = "Assay") %>% 
  # filter for significant rows
  dplyr::filter(`p-value FDR-adj` < 0.05) %>% 
  # remove p-value
  dplyr::select(-`p-value`) %>% 
  # convert adjusted p-value to character in scientific notation
  dplyr::mutate(`p-value FDR-adj` = as.character(scales::scientific(`p-value FDR-adj`, digits = 3))) %>% 
  # remove underscores from values
  dplyr::mutate(dplyr::across(c("Assay", "Variable"),
                              ~stringr::str_replace(., pattern = "_", " "))) %>% 
  # rename columns
  dplyr::rename("p-value (FDR-adjusted)" = "p-value FDR-adj",
                "Significance" = "Significance (p-value FDR-adj)")

## Write to Google sheet
googlesheets4::write_sheet(
  data = final_sge,
  ss = "https://docs.google.com/spreadsheets/d/1_l72BZkmWyNAOfCUI8WGP4UfQuIPQtPZZmlRjQffvEs",
  sheet = "SGE_FINAL")
#> ✔ Writing to "aov_state_freq".
#> ✔ Writing to sheet 'SGE_FINAL'.
```

