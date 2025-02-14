Survival Analysis
================
Polly Wu (rw3031)
2025-01-10

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(survival)
```

``` r
seer = read_csv("./SEER_RECODE_ses.csv")
```

    ## Rows: 69311 Columns: 26
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (11): age, Sex, race, stage, month_diagnosis_treatment, radiation_recode...
    ## dbl (15): year_of_diagnosis, surgery, survival_month, death_s, death_o, sex_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
seer=
  seer|>
  janitor::clean_names()|>
  mutate(race_time = survival_month*race_c)|>
  mutate(yost_c = case_when(
    yost_ses == "Group 1 (lowest SES)" ~ 1,
    yost_ses == "Group 2" ~ 2,
    yost_ses == "Group 3" ~3,
    yost_ses == "Group 4" ~ 4,
    yost_ses == "Group 5 (highest SES)"~ 5)
  )
```

``` r
ph_test = function(covariate){
cox_model= coxph(Surv(survival_month, death_s) ~ covariate, data = seer)

result = cox.zph(cox_model)

print(result)
}
```

# ph test for persistent poverty

``` r
cox_model= coxph(Surv(survival_month, death_s) ~ persistent_poverty_c, data =  seer)

result = cox.zph(cox_model)

print(result)
```

    ##                         chisq df p
    ## persistent_poverty_c 3.29e-05  1 1
    ## GLOBAL               3.29e-05  1 1

``` r
cox_model= coxph(Surv(survival_month, death_o) ~ persistent_poverty_c, data =  seer)

result = cox.zph(cox_model)

print(result)
```

    ##                      chisq df    p
    ## persistent_poverty_c 0.197  1 0.66
    ## GLOBAL               0.197  1 0.66

# ph test for overall death

## ph test for the stratified analysis of different race groups

``` r
seer_w = 
seer|>
  filter(race_c == 1)

cox_model= coxph(Surv(survival_month, death_o) ~ persistent_poverty_c, data =  seer_w)

result = cox.zph(cox_model)

print(result)
```

    ##                      chisq df    p
    ## persistent_poverty_c 0.135  1 0.71
    ## GLOBAL               0.135  1 0.71

``` r
seer_w = 
seer|>
  filter(race_c == 2)

cox_model= coxph(Surv(survival_month, death_o) ~ persistent_poverty_c, data =  seer_w)

result = cox.zph(cox_model)

print(result)
```

    ##                      chisq df    p
    ## persistent_poverty_c 0.132  1 0.72
    ## GLOBAL               0.132  1 0.72

``` r
seer_w = 
seer|>
  filter(race_c == 3)

cox_model= coxph(Surv(survival_month, death_o) ~ persistent_poverty_c, data =  seer_w)

result = cox.zph(cox_model)

print(result)
```

    ##                       chisq df   p
    ## persistent_poverty_c 0.0143  1 0.9
    ## GLOBAL               0.0143  1 0.9

``` r
seer_w = 
seer|>
  filter(race_c == 4)

cox_model= coxph(Surv(survival_month, death_o) ~ persistent_poverty_c, data =  seer_w)

result = cox.zph(cox_model)

print(result)
```

    ##                      chisq df    p
    ## persistent_poverty_c  0.13  1 0.72
    ## GLOBAL                0.13  1 0.72

``` r
seer_w = 
seer|>
  filter(race_c == 5)

cox_model= coxph(Surv(survival_month, death_o) ~ persistent_poverty_c, data =  seer_w)

result = cox.zph(cox_model)

print(result)
```

    ##                      chisq df    p
    ## persistent_poverty_c 0.298  1 0.58
    ## GLOBAL               0.298  1 0.58

## ph analysis for stages

``` r
seer_w = 
seer|>
  filter(stage_c == 0)

cox_model= coxph(Surv(survival_month, death_o) ~ persistent_poverty_c, data =  seer_w)

result = cox.zph(cox_model)

print(result)
```

    ##                       chisq df    p
    ## persistent_poverty_c 0.0847  1 0.77
    ## GLOBAL               0.0847  1 0.77

``` r
seer_w = 
seer|>
  filter(stage_c == 1)

cox_model= coxph(Surv(survival_month, death_o) ~ persistent_poverty_c, data =  seer_w)

result = cox.zph(cox_model)

print(result)
```

    ##                      chisq df    p
    ## persistent_poverty_c 0.406  1 0.52
    ## GLOBAL               0.406  1 0.52

``` r
seer_w = 
seer|>
  filter(stage_c == 2)

cox_model= coxph(Surv(survival_month, death_o) ~ persistent_poverty_c, data =  seer_w)

result = cox.zph(cox_model)

print(result)
```

    ##                      chisq df    p
    ## persistent_poverty_c  2.08  1 0.15
    ## GLOBAL                2.08  1 0.15

# ph test for cancer specific death

## ph test for the stratified analysis of different race groups

``` r
seer_w = 
seer|>
  filter(race_c == 1)

cox_model= coxph(Surv(survival_month, death_s) ~ persistent_poverty_c, data =  seer_w)

result = cox.zph(cox_model)

print(result)
```

    ##                       chisq df    p
    ## persistent_poverty_c 0.0798  1 0.78
    ## GLOBAL               0.0798  1 0.78

``` r
seer_w = 
seer|>
  filter(race_c == 2)

cox_model= coxph(Surv(survival_month, death_s) ~ persistent_poverty_c, data =  seer_w)

result = cox.zph(cox_model)

print(result)
```

    ##                        chisq df    p
    ## persistent_poverty_c 0.00723  1 0.93
    ## GLOBAL               0.00723  1 0.93

``` r
seer_w = 
seer|>
  filter(race_c == 3)

cox_model= coxph(Surv(survival_month, death_s) ~ persistent_poverty_c, data =  seer_w)

result = cox.zph(cox_model)

print(result)
```

    ##                      chisq df    p
    ## persistent_poverty_c 0.592  1 0.44
    ## GLOBAL               0.592  1 0.44

``` r
seer_w = 
seer|>
  filter(race_c == 4)

cox_model= coxph(Surv(survival_month, death_s) ~ persistent_poverty_c, data =  seer_w)

result = cox.zph(cox_model)

print(result)
```

    ##                       chisq df   p
    ## persistent_poverty_c 0.0636  1 0.8
    ## GLOBAL               0.0636  1 0.8

``` r
seer_w = 
seer|>
  filter(race_c == 5)

cox_model= coxph(Surv(survival_month, death_s) ~ persistent_poverty_c, data =  seer_w)

result = cox.zph(cox_model)

print(result)
```

    ##                      chisq df    p
    ## persistent_poverty_c  2.09  1 0.15
    ## GLOBAL                2.09  1 0.15

## ph analysis for stages

``` r
seer_w = 
seer|>
  filter(stage_c == 0)

cox_model= coxph(Surv(survival_month, death_s) ~ persistent_poverty_c, data =  seer_w)

result = cox.zph(cox_model)

print(result)
```

    ##                      chisq df    p
    ## persistent_poverty_c 0.544  1 0.46
    ## GLOBAL               0.544  1 0.46

``` r
seer_w = 
seer|>
  filter(stage_c == 1)

cox_model= coxph(Surv(survival_month, death_s) ~ persistent_poverty_c, data =  seer_w)

result = cox.zph(cox_model)

print(result)
```

    ##                      chisq df    p
    ## persistent_poverty_c 0.422  1 0.52
    ## GLOBAL               0.422  1 0.52

``` r
seer_w = 
seer|>
  filter(stage_c == 2)

cox_model= coxph(Surv(survival_month, death_s) ~ persistent_poverty_c, data =  seer_w)

result = cox.zph(cox_model)

print(result)
```

    ##                      chisq df    p
    ## persistent_poverty_c  1.88  1 0.17
    ## GLOBAL                1.88  1 0.17

``` r
cox_model= coxph(Surv(survival_month, death_o) ~ persistent_poverty_c+sex+age_c+marital+race_c+rural+race_time, data =  seer)

result = cox.zph(cox_model)

print(result)
```

    ##                         chisq df       p
    ## persistent_poverty_c   282.50  1 < 2e-16
    ## sex                      9.69  1  0.0019
    ## age_c                   16.27  1 5.5e-05
    ## marital                 53.55  5 2.6e-10
    ## race_c               48648.66  1 < 2e-16
    ## rural                 1055.60  3 < 2e-16
    ## race_time            73591.70  1 < 2e-16
    ## GLOBAL               93871.88 13 < 2e-16
