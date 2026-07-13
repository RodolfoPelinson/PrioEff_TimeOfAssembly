Effect of first survey on the effect of time
================
Rodolfo Pelinson
2026-07-13

## Loading packages and functions

``` r
source(paste(sep = "/",dir,"functions/pairwise_gllvm.R"))
source(paste(sep = "/",dir,"functions/pairwise_mvabund.R"))
source(paste(sep = "/",dir,"functions/my.anova.gllvm.R"))
source(paste(sep = "/",dir,"functions/My_coefplot.R"))
source(paste(sep = "/",dir,"functions/remove_sp.R"))
source(paste(sep = "/",dir,"functions/extract_mles.R"))
source(paste(sep = "/",dir,"functions/my_ordiplot.R"))
source(paste(sep = "/",dir,"functions/get_scaled_lvs.R"))
source(paste(sep = "/",dir,"functions/text_contour.R"))
source(paste(sep = "/",dir,"functions/letters.R"))
```

``` r
library(vegan)
library(gllvm)
library(mvabund)
library(DHARMa)
library(glmmTMB)
library(vioplot)
library(yarrr)
library(colorspace)
```

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31 ucrt)
    ## Platform: x86_64-w64-mingw32/x64
    ## Running under: Windows 11 x64 (build 26200)
    ## 
    ## Matrix products: default
    ##   LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] LC_COLLATE=Portuguese_Brazil.utf8  LC_CTYPE=Portuguese_Brazil.utf8   
    ## [3] LC_MONETARY=Portuguese_Brazil.utf8 LC_NUMERIC=C                      
    ## [5] LC_TIME=Portuguese_Brazil.utf8    
    ## 
    ## time zone: Europe/Berlin
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] colorspace_2.1-2       yarrr_0.1.14           circlize_0.4.17       
    ##  [4] BayesFactor_0.9.12-4.8 Matrix_1.7-4           coda_0.19-4.1         
    ##  [7] jpeg_0.1-11            vioplot_0.5.1          zoo_1.8-15            
    ## [10] sm_2.2-6.0             glmmTMB_1.1.14         DHARMa_0.4.7          
    ## [13] mvabund_4.2.8          gllvm_2.0.5            TMB_1.9.20            
    ## [16] vegan_2.7-3            permute_0.9-10        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] sandwich_3.1-1      shape_1.4.6.1       stringi_1.8.7      
    ##  [4] lattice_0.22-7      magrittr_2.0.4      lme4_2.0-1         
    ##  [7] digest_0.6.39       evaluate_1.0.5      grid_4.5.2         
    ## [10] estimability_1.5.1  mvtnorm_1.3-6       fastmap_1.2.0      
    ## [13] GlobalOptions_0.1.3 mgcv_1.9-3          pbapply_1.7-4      
    ## [16] numDeriv_2016.8-1.1 reformulas_0.4.4    Rdpack_2.6.6       
    ## [19] cli_3.6.5           rlang_1.1.7         rbibutils_2.4.1    
    ## [22] splines_4.5.2       yaml_2.3.12         otel_0.2.0         
    ## [25] tools_4.5.2         parallel_4.5.2      MatrixModels_0.5-4 
    ## [28] nloptr_2.2.1        minqa_1.2.8         boot_1.3-32        
    ## [31] lifecycle_1.0.5     emmeans_2.0.2       stringr_1.6.0      
    ## [34] tweedie_3.0.17      MASS_7.3-65         cluster_2.1.8.1    
    ## [37] glue_1.8.0          Rcpp_1.1.1          statmod_1.5.1      
    ## [40] xfun_0.57           rstudioapi_0.18.0   knitr_1.51         
    ## [43] xtable_1.8-8        htmltools_0.5.9     nlme_3.1-168       
    ## [46] rmarkdown_2.31      compiler_4.5.2      alabama_2025.1.0

## Loading and preparing data

``` r
source(paste(sep = "/",dir,"ajeitando_planilhas.R"))
```

``` r
comm_all_drop_atrasado <- comm_all[Exp_design_all$treatments != "atrasado",]
Exp_design_all_drop_atrasado <- Exp_design_all[Exp_design_all$treatments != "atrasado",]

Exp_design_drop_atrasado <- Exp_design[Exp_design$treatments != "atrasado",]

ncol(comm_all_drop_atrasado)
```

    ## [1] 26

``` r
comm_all_drop_atrasado_rm <- remove_sp(comm_all_drop_atrasado, 10)
ncol(comm_all_drop_atrasado_rm)
```

    ## [1] 11

``` r
Exp_design_all_drop_atrasado$treatments_AM <- paste(Exp_design_all_drop_atrasado$treatments, Exp_design_all_drop_atrasado$AM, sep="_")

Exp_design_all_drop_atrasado$AM_numeric <- as.numeric(Exp_design_all_drop_atrasado$AM)

Exp_design_all_drop_atrasado$treatments <- as.factor(Exp_design_all_drop_atrasado$treatments)

AM_days <- as.numeric(Exp_design_all_drop_atrasado$AM)
AM_days[AM_days == 1] <- 32
AM_days[AM_days == 2] <- 80
AM_days[AM_days == 3] <- 116
AM_days[AM_days == 4] <- 158

AM_days_sc <- scale(AM_days)

Exp_design_all_drop_atrasado$AM_days_sc <- c(AM_days_sc)
Exp_design_all_drop_atrasado$AM_days_sc_squared <- c(AM_days_sc)^2

Exp_design_all_drop_atrasado$AM_days <- c(AM_days)
Exp_design_all_drop_atrasado$AM_days_squared <- c(AM_days)^2


comm_AM1 <- comm_all_drop_atrasado_rm[Exp_design_all_drop_atrasado$AM == 1,]
comm_AM1_rm <- remove_sp(comm_AM1, 3)



comm_AM1_pred <- rbind(comm_AM1_rm, comm_AM1_rm, comm_AM1_rm, comm_AM1_rm)

nrow(comm_AM1_pred)
```

    ## [1] 96

``` r
nrow(Exp_design_all_drop_atrasado)
```

    ## [1] 96

``` r
nrow(comm_all_drop_atrasado_rm)
```

    ## [1] 96

``` r
comm_AM1_pred_st <- decostand(comm_AM1_pred, "stand")
colnames(comm_AM1_pred_st) <- gsub(" ", "_", colnames(comm_AM1_pred_st))

predictors <- data.frame(Exp_design_all_drop_atrasado, comm_AM1_pred_st)
```

## Using NMDS to describe initial communities.

``` r
comm_AM1_rm_total <- data.frame(decostand(comm_AM1_rm, method = "total", MARGIN = 2))
nmds_comm_AM1_rm_total <- metaMDS(comm_AM1_rm_total, distance  = "bray", k = 2)
```

    ## Run 0 stress 0.1917598 
    ## Run 1 stress 0.1917598 
    ## ... Procrustes: rmse 1.066959e-05  max resid 3.324409e-05 
    ## ... Similar to previous best
    ## Run 2 stress 0.2346589 
    ## Run 3 stress 0.2284806 
    ## Run 4 stress 0.1917598 
    ## ... New best solution
    ## ... Procrustes: rmse 5.771373e-06  max resid 1.865517e-05 
    ## ... Similar to previous best
    ## Run 5 stress 0.2527162 
    ## Run 6 stress 0.2288855 
    ## Run 7 stress 0.1917598 
    ## ... Procrustes: rmse 8.776373e-06  max resid 3.1457e-05 
    ## ... Similar to previous best
    ## Run 8 stress 0.2379068 
    ## Run 9 stress 0.1917598 
    ## ... Procrustes: rmse 4.42891e-06  max resid 1.582488e-05 
    ## ... Similar to previous best
    ## Run 10 stress 0.1917598 
    ## ... Procrustes: rmse 4.824192e-06  max resid 1.678249e-05 
    ## ... Similar to previous best
    ## Run 11 stress 0.1917598 
    ## ... Procrustes: rmse 2.559498e-06  max resid 6.907674e-06 
    ## ... Similar to previous best
    ## Run 12 stress 0.1917598 
    ## ... Procrustes: rmse 3.908659e-06  max resid 1.45046e-05 
    ## ... Similar to previous best
    ## Run 13 stress 0.2355307 
    ## Run 14 stress 0.2297681 
    ## Run 15 stress 0.1917598 
    ## ... Procrustes: rmse 4.590371e-06  max resid 1.619287e-05 
    ## ... Similar to previous best
    ## Run 16 stress 0.2309966 
    ## Run 17 stress 0.2316998 
    ## Run 18 stress 0.2283595 
    ## Run 19 stress 0.1917598 
    ## ... Procrustes: rmse 2.185291e-06  max resid 7.707593e-06 
    ## ... Similar to previous best
    ## Run 20 stress 0.1917598 
    ## ... Procrustes: rmse 6.325846e-06  max resid 2.272347e-05 
    ## ... Similar to previous best
    ## *** Best solution repeated 9 times

``` r
comm_AM1_nmds <- data.frame(nmds_comm_AM1_rm_total$points)
comm_AM1_nmds_st <- decostand(comm_AM1_nmds, method = "stand")

comm_AM1_nmds_st_pred <- rbind(comm_AM1_nmds_st, comm_AM1_nmds_st, comm_AM1_nmds_st, comm_AM1_nmds_st)

predictors_nmds <- data.frame(Exp_design_all_drop_atrasado, comm_AM1_nmds_st_pred)
```

### Testing for effects

``` r
nrow(predictors_nmds)
```

    ## [1] 96

``` r
nrow(comm_all_drop_atrasado_rm)
```

    ## [1] 96

``` r
set.seed(1)
control <- permute::how(within = permute::Within(type = 'free'),
                        plots = Plots(strata = predictors_nmds$sites, type = 'free'),
                        nperm = 999)
permutations <- shuffleSet(nrow(comm_all_drop_atrasado_rm), control = control)


comm_all_drop_atrasado_mvabund <- mvabund(comm_all_drop_atrasado_rm)


mod_full_lin <- manyglm(comm_all_drop_atrasado_mvabund ~ block2 + ( MDS1 + MDS2 ) + AM_days_sc + AM_days_sc:( MDS1 + MDS2 ), family="negative.binomial",  data = predictors_nmds, composition = FALSE)

mod_full_quad <- manyglm(comm_all_drop_atrasado_mvabund ~ block2 + ( MDS1 + MDS2 ) + AM_days_sc + AM_days_sc_squared + AM_days_sc:( MDS1 + MDS2 ) + AM_days_sc_squared:( MDS1 + MDS2 ), family="negative.binomial",  data = predictors_nmds, composition = FALSE)

anova_quadratic_test <- anova(mod_full_lin, mod_full_quad, bootID = permutations, show.time = "all", cor.type = "I", test = "LR", resamp = "pit.trap")
```

    ## Using <int> bootID matrix from input. 
    ## Resampling begins for test 1.
    ##  Resampling run 0 finished. Time elapsed: 0.00 minutes...
    ##  Resampling run 100 finished. Time elapsed: 0.05 minutes...
    ##  Resampling run 200 finished. Time elapsed: 0.10 minutes...
    ##  Resampling run 300 finished. Time elapsed: 0.16 minutes...
    ##  Resampling run 400 finished. Time elapsed: 0.21 minutes...
    ##  Resampling run 500 finished. Time elapsed: 0.26 minutes...
    ##  Resampling run 600 finished. Time elapsed: 0.31 minutes...
    ##  Resampling run 700 finished. Time elapsed: 0.36 minutes...
    ##  Resampling run 800 finished. Time elapsed: 0.41 minutes...
    ##  Resampling run 900 finished. Time elapsed: 0.46 minutes...
    ## Time elapsed: 0 hr 0 min 30 sec

``` r
anova_quadratic_test
```

    ## Analysis of Deviance Table
    ## 
    ## mod_full_lin: comm_all_drop_atrasado_mvabund ~ block2 + (MDS1 + MDS2) + AM_days_sc + AM_days_sc:(MDS1 + MDS2)
    ## mod_full_quad: comm_all_drop_atrasado_mvabund ~ block2 + (MDS1 + MDS2) + AM_days_sc + AM_days_sc_squared + AM_days_sc:(MDS1 + MDS2) + AM_days_sc_squared:(MDS1 + MDS2)
    ## 
    ## Multivariate test:
    ##               Res.Df Df.diff   Dev Pr(>Dev)   
    ## mod_full_lin      88                          
    ## mod_full_quad     85       3 126.8    0.007 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Arguments:
    ##  Test statistics calculated assuming uncorrelated response (for faster computation) 
    ##  P-value calculated using 999 iterations via PIT-trap resampling.

``` r
########################################################
mod_block <- manyglm(comm_all_drop_atrasado_mvabund ~ block2 + ( MDS1 + MDS2 ) + AM_days_sc + AM_days_sc_squared, family="negative.binomial",  data = predictors_nmds, composition = FALSE)

mod_no_block <- manyglm(comm_all_drop_atrasado_mvabund ~ ( MDS1 + MDS2 ) + AM_days_sc + AM_days_sc_squared, family="negative.binomial",  data = predictors_nmds, composition = FALSE)

anova_block <- anova(mod_no_block, mod_block, bootID = permutations, show.block = "all", cor.type = "I", test = "LR", resamp = "pit.trap")
```

    ## Warning in anova.manyglm(mod_no_block, mod_block, bootID = permutations, :
    ## show.block is not a manyglm object nor a valid argument- removed from input,
    ## default value is used instead.

    ## Using <int> bootID matrix from input. 
    ## Time elapsed: 0 hr 0 min 20 sec

``` r
anova_block
```

    ## Analysis of Deviance Table
    ## 
    ## mod_no_block: comm_all_drop_atrasado_mvabund ~ (MDS1 + MDS2) + AM_days_sc + AM_days_sc_squared
    ## mod_block: comm_all_drop_atrasado_mvabund ~ block2 + (MDS1 + MDS2) + AM_days_sc + AM_days_sc_squared
    ## 
    ## Multivariate test:
    ##              Res.Df Df.diff   Dev Pr(>Dev)   
    ## mod_no_block     91                          
    ## mod_block        89       2 70.02    0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Arguments:
    ##  Test statistics calculated assuming uncorrelated response (for faster computation) 
    ##  P-value calculated using 999 iterations via PIT-trap resampling.

``` r
########################################################


mod_time <- manyglm(comm_all_drop_atrasado_mvabund ~ block2 + ( MDS1 + MDS2 ) + AM_days_sc + AM_days_sc_squared, family="negative.binomial",  data = predictors_nmds, composition = FALSE)

mod_no_time <- manyglm(comm_all_drop_atrasado_mvabund ~ block2 + ( MDS1 + MDS2 ), family="negative.binomial",  data = predictors_nmds, composition = FALSE)

anova_time <- anova(mod_no_time, mod_time, bootID = permutations, show.time = "all", cor.type = "I", test = "LR", resamp = "pit.trap")
```

    ## Using <int> bootID matrix from input. 
    ## Resampling begins for test 1.
    ##  Resampling run 0 finished. Time elapsed: 0.00 minutes...
    ##  Resampling run 100 finished. Time elapsed: 0.03 minutes...
    ##  Resampling run 200 finished. Time elapsed: 0.06 minutes...
    ##  Resampling run 300 finished. Time elapsed: 0.10 minutes...
    ##  Resampling run 400 finished. Time elapsed: 0.13 minutes...
    ##  Resampling run 500 finished. Time elapsed: 0.17 minutes...
    ##  Resampling run 600 finished. Time elapsed: 0.23 minutes...
    ##  Resampling run 700 finished. Time elapsed: 0.27 minutes...
    ##  Resampling run 800 finished. Time elapsed: 0.30 minutes...
    ##  Resampling run 900 finished. Time elapsed: 0.33 minutes...
    ## Time elapsed: 0 hr 0 min 21 sec

``` r
anova_time
```

    ## Analysis of Deviance Table
    ## 
    ## mod_no_time: comm_all_drop_atrasado_mvabund ~ block2 + (MDS1 + MDS2)
    ## mod_time: comm_all_drop_atrasado_mvabund ~ block2 + (MDS1 + MDS2) + AM_days_sc + AM_days_sc_squared
    ## 
    ## Multivariate test:
    ##             Res.Df Df.diff   Dev Pr(>Dev)    
    ## mod_no_time     91                           
    ## mod_time        89       2 258.7    0.001 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Arguments:
    ##  Test statistics calculated assuming uncorrelated response (for faster computation) 
    ##  P-value calculated using 999 iterations via PIT-trap resampling.

``` r
########################################################



mod_prior <- manyglm(comm_all_drop_atrasado_mvabund ~ block2 + ( MDS1 + MDS2 ) + AM_days_sc + AM_days_sc_squared, family="negative.binomial",  data = predictors_nmds, composition = FALSE)

mod_no_prior <- manyglm(comm_all_drop_atrasado_mvabund ~ block2 + AM_days_sc + AM_days_sc_squared, family="negative.binomial",  data = predictors_nmds, composition = FALSE)

anova_prior <- anova(mod_no_prior, mod_prior, bootID = permutations, show.time = "all", cor.type = "I", test = "LR", resamp = "pit.trap")
```

    ## Using <int> bootID matrix from input. 
    ## Resampling begins for test 1.
    ##  Resampling run 0 finished. Time elapsed: 0.00 minutes...
    ##  Resampling run 100 finished. Time elapsed: 0.03 minutes...
    ##  Resampling run 200 finished. Time elapsed: 0.06 minutes...
    ##  Resampling run 300 finished. Time elapsed: 0.10 minutes...
    ##  Resampling run 400 finished. Time elapsed: 0.13 minutes...
    ##  Resampling run 500 finished. Time elapsed: 0.16 minutes...
    ##  Resampling run 600 finished. Time elapsed: 0.19 minutes...
    ##  Resampling run 700 finished. Time elapsed: 0.22 minutes...
    ##  Resampling run 800 finished. Time elapsed: 0.25 minutes...
    ##  Resampling run 900 finished. Time elapsed: 0.29 minutes...
    ## Time elapsed: 0 hr 0 min 18 sec

``` r
anova_prior
```

    ## Analysis of Deviance Table
    ## 
    ## mod_no_prior: comm_all_drop_atrasado_mvabund ~ block2 + AM_days_sc + AM_days_sc_squared
    ## mod_prior: comm_all_drop_atrasado_mvabund ~ block2 + (MDS1 + MDS2) + AM_days_sc + AM_days_sc_squared
    ## 
    ## Multivariate test:
    ##              Res.Df Df.diff   Dev Pr(>Dev)   
    ## mod_no_prior     91                          
    ## mod_prior        89       2 59.44    0.007 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Arguments:
    ##  Test statistics calculated assuming uncorrelated response (for faster computation) 
    ##  P-value calculated using 999 iterations via PIT-trap resampling.

``` r
########################################################
mod_no_int <- manyglm(comm_all_drop_atrasado_mvabund ~ block2 + ( MDS1 + MDS2 ) + AM_days_sc + AM_days_sc_squared, family="negative.binomial",  data = predictors_nmds, composition = FALSE)

mod_int <- manyglm(comm_all_drop_atrasado_mvabund ~ block2 + ( MDS1 + MDS2 ) + AM_days_sc + AM_days_sc_squared + AM_days_sc:( MDS1 + MDS2 ) + AM_days_sc_squared:( MDS1 + MDS2 ), family="negative.binomial",  data = predictors_nmds, composition = FALSE)

anova_interaction <- anova(mod_no_int, mod_int, bootID = permutations, show.time = "all", cor.type = "I", test = "LR", resamp = "pit.trap")
```

    ## Using <int> bootID matrix from input. 
    ## Resampling begins for test 1.
    ##  Resampling run 0 finished. Time elapsed: 0.00 minutes...
    ##  Resampling run 100 finished. Time elapsed: 0.05 minutes...
    ##  Resampling run 200 finished. Time elapsed: 0.09 minutes...
    ##  Resampling run 300 finished. Time elapsed: 0.14 minutes...
    ##  Resampling run 400 finished. Time elapsed: 0.19 minutes...
    ##  Resampling run 500 finished. Time elapsed: 0.23 minutes...
    ##  Resampling run 600 finished. Time elapsed: 0.28 minutes...
    ##  Resampling run 700 finished. Time elapsed: 0.33 minutes...
    ##  Resampling run 800 finished. Time elapsed: 0.38 minutes...
    ##  Resampling run 900 finished. Time elapsed: 0.42 minutes...
    ## Time elapsed: 0 hr 0 min 28 sec

``` r
anova_interaction
```

    ## Analysis of Deviance Table
    ## 
    ## mod_no_int: comm_all_drop_atrasado_mvabund ~ block2 + (MDS1 + MDS2) + AM_days_sc + AM_days_sc_squared
    ## mod_int: comm_all_drop_atrasado_mvabund ~ block2 + (MDS1 + MDS2) + AM_days_sc + AM_days_sc_squared + AM_days_sc:(MDS1 + MDS2) + AM_days_sc_squared:(MDS1 + MDS2)
    ## 
    ## Multivariate test:
    ##            Res.Df Df.diff   Dev Pr(>Dev)    
    ## mod_no_int     89                           
    ## mod_int        85       4 81.74    0.001 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Arguments:
    ##  Test statistics calculated assuming uncorrelated response (for faster computation) 
    ##  P-value calculated using 999 iterations via PIT-trap resampling.

Which are the species who have at least one of the coefficients for the
effect of time in interaction with the effect of the first survey whose
confidence interval do not include zero?

``` r
coefs <- mod_int$coefficients
LCL <- coefs + mod_int$stderr.coefficients * qnorm(0.025) 
UCL <- coefs + mod_int$stderr.coefficients * qnorm(0.975) 

coefs[LCL < 0 &  UCL > 0] <- NA

coefs <- coefs[-c(1:7),]

coefs <- coefs[, colSums(coefs, na.rm = TRUE) != 0]
```

### Making predictions for plotting

``` r
comm_all_drop_atrasado_mvabund2 <- comm_all_drop_atrasado_mvabund + 0.001

mod_int <- manyglm(comm_all_drop_atrasado_mvabund2 ~ block2 + ( MDS1 + MDS2 ) + AM_days_sc + AM_days_sc_squared + AM_days_sc:( MDS1 + MDS2 ) + AM_days_sc_squared:( MDS1 + MDS2 ), family="negative.binomial",  data = predictors_nmds, composition = FALSE)
```

    ## Warning in manyglm(comm_all_drop_atrasado_mvabund2 ~ block2 + (MDS1 + MDS2) + :
    ## Non-integer data are fitted to the negative.binomial model.

``` r
#High abundances of 
#Aedes
#Culex
#Dendropsophus_minutus
#Tanypodinae
#Callibaetis
#Culex

new_data_AB_high_MDS1_high_MDS2 <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("AB", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
             MDS1 = rep(quantile(predictors_nmds$MDS1, 0.9), 100),
             MDS2 = rep(quantile(predictors_nmds$MDS2, 0.9), 100))

new_data_CD_high_MDS1_high_MDS2 <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("CD", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
             MDS1 = rep(quantile(predictors_nmds$MDS1, 0.9), 100),
             MDS2 = rep(quantile(predictors_nmds$MDS2, 0.9), 100))


new_data_EF_high_MDS1_high_MDS2 <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("EF", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
             MDS1 = rep(quantile(predictors_nmds$MDS1, 0.9), 100),
             MDS2 = rep(quantile(predictors_nmds$MDS2, 0.9), 100))


###########################################################################################



new_data_AB_low_MDS1_high_MDS2 <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("AB", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
             MDS1 = rep(quantile(predictors_nmds$MDS1, 0.1), 100),
             MDS2 = rep(quantile(predictors_nmds$MDS2, 0.9), 100))

new_data_CD_low_MDS1_high_MDS2 <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("CD", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
             MDS1 = rep(quantile(predictors_nmds$MDS1, 0.1), 100),
             MDS2 = rep(quantile(predictors_nmds$MDS2, 0.9), 100))


new_data_EF_low_MDS1_high_MDS2 <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("EF", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
             MDS1 = rep(quantile(predictors_nmds$MDS1, 0.1), 100),
             MDS2 = rep(quantile(predictors_nmds$MDS2, 0.9), 100))


###########################################################################################




new_data_AB_low_MDS1_low_MDS2 <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("AB", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
             MDS1 = rep(quantile(predictors_nmds$MDS1, 0.1), 100),
             MDS2 = rep(quantile(predictors_nmds$MDS2, 0.1), 100))

new_data_CD_low_MDS1_low_MDS2 <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("CD", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
             MDS1 = rep(quantile(predictors_nmds$MDS1, 0.1), 100),
             MDS2 = rep(quantile(predictors_nmds$MDS2, 0.1), 100))


new_data_EF_low_MDS1_low_MDS2 <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("EF", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
             MDS1 = rep(quantile(predictors_nmds$MDS1, 0.1), 100),
             MDS2 = rep(quantile(predictors_nmds$MDS2, 0.1), 100))


###########################################################################################



new_data_AB_high_MDS1_low_MDS2 <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("AB", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
             MDS1 = rep(quantile(predictors_nmds$MDS1, 0.9), 100),
             MDS2 = rep(quantile(predictors_nmds$MDS2, 0.1), 100))

new_data_CD_high_MDS1_low_MDS2 <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("CD", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
             MDS1 = rep(quantile(predictors_nmds$MDS1, 0.9), 100),
             MDS2 = rep(quantile(predictors_nmds$MDS2, 0.1), 100))


new_data_EF_high_MDS1_low_MDS2 <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("EF", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
             MDS1 = rep(quantile(predictors_nmds$MDS1, 0.9), 100),
             MDS2 = rep(quantile(predictors_nmds$MDS2, 0.1), 100))


###########################################################################################





###########################################################################################


AB_high_MDS1_high_MDS2 <- predict.manyglm(mod_int, newdata = new_data_AB_high_MDS1_high_MDS2, type = "response")
```

    ## Warning: glm.fit: algoritmo não convergiu

    ## Warning: glm.fit: algoritmo não convergiu

``` r
CD_high_MDS1_high_MDS2 <- predict.manyglm(mod_int, newdata = new_data_CD_high_MDS1_high_MDS2, type = "response")
```

    ## Warning: glm.fit: algoritmo não convergiu
    ## Warning: glm.fit: algoritmo não convergiu

``` r
EF_high_MDS1_high_MDS2 <- predict.manyglm(mod_int, newdata = new_data_EF_high_MDS1_high_MDS2, type = "response")
```

    ## Warning: glm.fit: algoritmo não convergiu
    ## Warning: glm.fit: algoritmo não convergiu

``` r
high_MDS1_high_MDS2 <- array(NA, dim = c(nrow(AB_high_MDS1_high_MDS2), ncol(AB_high_MDS1_high_MDS2), 3))
high_MDS1_high_MDS2[,,1] <- AB_high_MDS1_high_MDS2
high_MDS1_high_MDS2[,,2] <- CD_high_MDS1_high_MDS2
high_MDS1_high_MDS2[,,3] <- EF_high_MDS1_high_MDS2

high_MDS1_high_MDS2_mean <- as.data.frame(apply(high_MDS1_high_MDS2, MARGIN = c(1,2), FUN = mean))
colnames(high_MDS1_high_MDS2_mean) <- colnames(AB_high_MDS1_high_MDS2)
#high_MDS1_high_MDS2_mean






AB_low_MDS1_high_MDS2 <- predict.manyglm(mod_int, newdata = new_data_AB_low_MDS1_high_MDS2, type = "response")
```

    ## Warning: glm.fit: algoritmo não convergiu
    ## Warning: glm.fit: algoritmo não convergiu

``` r
CD_low_MDS1_high_MDS2 <- predict.manyglm(mod_int, newdata = new_data_CD_low_MDS1_high_MDS2, type = "response")
```

    ## Warning: glm.fit: algoritmo não convergiu
    ## Warning: glm.fit: algoritmo não convergiu

``` r
EF_low_MDS1_high_MDS2 <- predict.manyglm(mod_int, newdata = new_data_EF_low_MDS1_high_MDS2, type = "response")
```

    ## Warning: glm.fit: algoritmo não convergiu
    ## Warning: glm.fit: algoritmo não convergiu

``` r
low_MDS1_high_MDS2 <- array(NA, dim = c(nrow(AB_low_MDS1_high_MDS2), ncol(AB_low_MDS1_high_MDS2), 3))
low_MDS1_high_MDS2[,,1] <- AB_low_MDS1_high_MDS2
low_MDS1_high_MDS2[,,2] <- CD_low_MDS1_high_MDS2
low_MDS1_high_MDS2[,,3] <- EF_low_MDS1_high_MDS2

low_MDS1_high_MDS2_mean <- as.data.frame(apply(low_MDS1_high_MDS2, MARGIN = c(1,2), FUN = mean))
colnames(low_MDS1_high_MDS2_mean) <- colnames(AB_low_MDS1_high_MDS2)
#low_MDS1_high_MDS2_mean





AB_high_MDS1_low_MDS2 <- predict.manyglm(mod_int, newdata = new_data_AB_high_MDS1_low_MDS2, type = "response")
```

    ## Warning: glm.fit: algoritmo não convergiu
    ## Warning: glm.fit: algoritmo não convergiu

``` r
CD_high_MDS1_low_MDS2 <- predict.manyglm(mod_int, newdata = new_data_CD_high_MDS1_low_MDS2, type = "response")
```

    ## Warning: glm.fit: algoritmo não convergiu
    ## Warning: glm.fit: algoritmo não convergiu

``` r
EF_high_MDS1_low_MDS2 <- predict.manyglm(mod_int, newdata = new_data_EF_high_MDS1_low_MDS2, type = "response")
```

    ## Warning: glm.fit: algoritmo não convergiu
    ## Warning: glm.fit: algoritmo não convergiu

``` r
high_MDS1_low_MDS2 <- array(NA, dim = c(nrow(AB_high_MDS1_low_MDS2), ncol(AB_high_MDS1_low_MDS2), 3))
high_MDS1_low_MDS2[,,1] <- AB_high_MDS1_low_MDS2
high_MDS1_low_MDS2[,,2] <- CD_high_MDS1_low_MDS2
high_MDS1_low_MDS2[,,3] <- EF_high_MDS1_low_MDS2

high_MDS1_low_MDS2_mean <- as.data.frame(apply(high_MDS1_low_MDS2, MARGIN = c(1,2), FUN = mean))
colnames(high_MDS1_low_MDS2_mean) <- colnames(AB_high_MDS1_low_MDS2)
#high_MDS1_low_MDS2_mean




AB_low_MDS1_low_MDS2 <- predict.manyglm(mod_int, newdata = new_data_AB_low_MDS1_low_MDS2, type = "response")
```

    ## Warning: glm.fit: algoritmo não convergiu
    ## Warning: glm.fit: algoritmo não convergiu

``` r
CD_low_MDS1_low_MDS2 <- predict.manyglm(mod_int, newdata = new_data_CD_low_MDS1_low_MDS2, type = "response")
```

    ## Warning: glm.fit: algoritmo não convergiu
    ## Warning: glm.fit: algoritmo não convergiu

``` r
EF_low_MDS1_low_MDS2 <- predict.manyglm(mod_int, newdata = new_data_EF_low_MDS1_low_MDS2, type = "response")
```

    ## Warning: glm.fit: algoritmo não convergiu
    ## Warning: glm.fit: algoritmo não convergiu

``` r
low_MDS1_low_MDS2 <- array(NA, dim = c(nrow(AB_low_MDS1_low_MDS2), ncol(AB_low_MDS1_low_MDS2), 3))
low_MDS1_low_MDS2[,,1] <- AB_low_MDS1_low_MDS2
low_MDS1_low_MDS2[,,2] <- CD_low_MDS1_low_MDS2
low_MDS1_low_MDS2[,,3] <- EF_low_MDS1_low_MDS2

low_MDS1_low_MDS2_mean <- as.data.frame(apply(low_MDS1_low_MDS2, MARGIN = c(1,2), FUN = mean))
colnames(low_MDS1_low_MDS2_mean) <- colnames(AB_low_MDS1_low_MDS2)
#low_MDS1_low_MDS2_mean
```

### Ploting predicted effects

``` r
col_control <- "#56B4E9"
col_closed <- "#009E73"
col_exc_drag <- "#CC79A7"
col_exc_amph <- "#E69F00"

col_high_1_high_2 <- "#290AD8" #blue dark 
col_low_1_high_2 <- "#D82632" #red dark
col_high_1_low_2 <- "#72D9FF" #blue light 
col_low_1_low_2 <- "#FFAD72" #red light 

sp_names <- colnames(high_MDS1_high_MDS2_mean)
sp_names <- gsub("\\.", " ", sp_names)

sp_fonts <- rep(3, ncol(high_MDS1_high_MDS2_mean))
sp_fonts[sp_names == "Chironominae" | sp_names == "Tanypodinae"] <- 1

sp_letters <- c("a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)", "i)", "j)", "k)")


#svg(filename = "Plots/Community structure analysis/predicted_species_trajectories.svg", width = 10, height = 12, pointsize = 13)


close.screen(all.screens = TRUE)
```

    ## [1] FALSE

``` r
split.screen(matrix(c(0,0.3333333,       0.75,1,
                      0.3333333,0.666666,0.75,1,
                      0.666666,1,        0.75,1,
                      0,0.3333333,       0.5,0.75,
                      0.3333333,0.666666,0.5,0.75,
                      0.666666,1,        0.5,0.75,
                      0,0.3333333,       0.25,0.5,
                      0.3333333,0.666666,0.25,0.5,
                      0.666666,1,        0.25,0.5,
                      0,0.3333333,       0,0.25,
                      0.3333333,0.666666,0,0.25,
                      0.666666,1,        0,0.25), ncol = 4, nrow = 12, byrow = TRUE))
```

    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12

``` r
for(i in 1:11){
  
screen(i)
  
par(mar = c(4,4,1,1), bty = "l")
  
ymax <- max(c(high_MDS1_high_MDS2_mean[,i], high_MDS1_low_MDS2_mean[,i], low_MDS1_high_MDS2_mean[,i], low_MDS1_low_MDS2_mean[,i]))*1.3

plot(NA, xlim = c(32 - 10, 158 + 10), ylim = c(0, ymax), xaxt = "n", yaxt = "n", xlab = "", ylab = "")

lines(y = high_MDS1_high_MDS2_mean[,i],                                                                                                                         x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= "white", lwd = 3)
lines(y = high_MDS1_high_MDS2_mean[,i],
          x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= col_high_1_high_2, lwd = 2)


lines(y = low_MDS1_high_MDS2_mean[,i],                                                                                                                    x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= "white", lwd = 3)
lines(y = low_MDS1_high_MDS2_mean[,i],
          x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= col_low_1_high_2, lwd = 2)


lines(y = high_MDS1_low_MDS2_mean[,i],                                                                                                                         x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= "white", lwd = 3)
lines(y = high_MDS1_low_MDS2_mean[,i],
          x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= col_high_1_low_2, lwd = 2)


lines(y = low_MDS1_low_MDS2_mean[,i],                                                                                                                         x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= "white", lwd = 3)
lines(y = low_MDS1_low_MDS2_mean[,i],
          x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= col_low_1_low_2, lwd = 2)


axis(1, at = c(32, 80, 116, 158))
title(xlab = "Days", line = 2.5, cex.lab = 1.25)

axis(2, gap.axis = -10, las = 2)
#title(ylab = "Abundance of", line = 3.75, cex.lab = 1.25)
title(ylab = sp_names[i], line = 2.5, cex.lab = 1, font.lab = sp_fonts[i])

letters(x = 10, y = 97, sp_letters[i], cex = 1.5)

}

screen(12)

par(mar = c(0.1,0.1,0.1,0.1), bty = "n")
plot(NA, xlim = c(0, 100), ylim = c(0, 100), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend(legend = c("High nMDS1, high nMDS2",
                 "High nMDS1, low nMDS2",
                  "Low nMDS1, high nMDS2",
                  "Low nMDS1, low nMDS2"), x = 50, y = 50, bty = "n", lty = 1, lwd = 3, col = c(col_high_1_high_2, col_high_1_low_2, col_low_1_high_2,col_low_1_low_2), xjust = 0.5, yjust = 0.5)
```

![](Effect_first_survey_on_time_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
#dev.off()
```

### Plotting coefficients

``` r
sp_names <- colnames(mod_int$coefficients)
sp_names <- gsub("\\.", " ", sp_names)

sp_names[sp_names == "Dendropsophus minutus"] <- "D. minutus"
sp_names[sp_names == "Scinax fuscovarius"] <- "S. fuscovarius"

sp_font <- rep(3, length(sp_names))
sp_font[sp_names == "Chironominae"] <- 1
sp_font[sp_names == "Tanypodinae"] <- 1

interaction_coefs <- data.frame(t(mod_int$coefficients[(nrow(mod_int$coefficients)-3):nrow(mod_int$coefficients),]))
interaction_coefs_se <- data.frame(t(mod_int$stderr.coefficients[(nrow(mod_int$stderr.coefficients)-3):nrow(mod_int$stderr.coefficients),]))

interaction_upper <- interaction_coefs + (interaction_coefs_se * qnorm(0.975))
interaction_lower <- interaction_coefs + (interaction_coefs_se * qnorm(0.025))


#svg(filename = "Plots/Community structure analysis/community_trajectories_coefs.svg", width = 10, height = 4.5, pointsize = 12)


close.screen(all.screens = TRUE)
split.screen(matrix(c(0 , 0.1, 0, 1  ,
                      0.1, 0.325, 0, 1,
                      0.325, 0.55, 0, 1,
                      0.55, 0.775, 0, 1,
                      0.775, 1, 0, 1), ncol = 4, nrow = 5, byrow = TRUE))
```

    ## [1] 1 2 3 4 5

``` r
screen(2)
par(mar = c(4,2,2,0.5))

My_coefplot(mles = interaction_coefs$MDS1.AM_days_sc, upper = interaction_upper$MDS1.AM_days_sc,
            lower = interaction_lower$MDS1.AM_days_sc, col_sig = "red",
            cex_sig = 1.5, species_labels = sp_names, yaxis_font = sp_font, cex.axis = 0.9)
title(xlab = "nMDS1 : Time", cex.lab = 1, line = 2.5)
#screen(2, new = FALSE)
letters(x = 12, y = 94, "a)", cex = 1.5)

screen(3)
par(mar = c(4,2,2,0.5))

My_coefplot(mles = interaction_coefs$MDS2.AM_days_sc, upper = interaction_upper$MDS2.AM_days_sc,
            lower = interaction_lower$MDS2.AM_days_sc, col_sig = "red",
            cex_sig = 1.5, species_labels = FALSE, yaxis_font = sp_font, cex.axis = 0.9)
title(xlab = "nMDS2 : Time", cex.lab = 1, line = 2.5)
#screen(2, new = FALSE)
letters(x = 12, y = 94, "b)", cex = 1.5)

screen(4)
par(mar = c(4,2,2,0.5))

My_coefplot(mles = interaction_coefs$MDS1.AM_days_sc_squared, upper = interaction_upper$MDS1.AM_days_sc_squared,
            lower = interaction_lower$MDS1.AM_days_sc_squared, col_sig = "red",
            cex_sig = 1.5, species_labels = FALSE, yaxis_font = sp_font, cex.axis = 0.9, xlim = c(-2,5))
title(xlab = "nMDS1 : Time²", cex.lab = 1, line = 2.5)
#screen(2, new = FALSE)
letters(x = 12, y = 94, "c)", cex = 1.5)

screen(5)
par(mar = c(4,2,2,0.5))

My_coefplot(mles = interaction_coefs$MDS2.AM_days_sc_squared, upper = interaction_upper$MDS2.AM_days_sc_squared,
            lower = interaction_lower$MDS2.AM_days_sc_squared, col_sig = "red",
            cex_sig = 1.5, species_labels = FALSE, yaxis_font = sp_font, cex.axis = 0.9, xlim = c(-6,2))
title(xlab = "nMDS2 : Time²", cex.lab = 1, line = 2.5)
#screen(2, new = FALSE)
letters(x = 12, y = 94, "d)", cex = 1.5)
```

![](Effect_first_survey_on_time_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
#dev.off()
```

``` r
comm_trajectories <- rbind(high_MDS1_high_MDS2_mean, low_MDS1_high_MDS2_mean, high_MDS1_low_MDS2_mean, low_MDS1_low_MDS2_mean)

comm_trajectories_relative <- decostand(comm_trajectories, method = "total", MARGIN = 2)

set.seed(5); nmds_comm_trajectories <- metaMDS(comm_trajectories_relative, distance = "bray", k = 2, trymax  = 50, try = 100, engine = "monoMDS", autotransform  = FALSE)
```

    ## Run 0 stress 0.1592961 
    ## Run 1 stress 0.1592961 
    ## ... New best solution
    ## ... Procrustes: rmse 1.696257e-05  max resid 0.0001040694 
    ## ... Similar to previous best
    ## Run 2 stress 0.1592967 
    ## ... Procrustes: rmse 9.751421e-05  max resid 0.001726844 
    ## ... Similar to previous best
    ## Run 3 stress 0.1592962 
    ## ... Procrustes: rmse 0.0001644513  max resid 0.002403504 
    ## ... Similar to previous best
    ## Run 4 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001700941  max resid 0.002425061 
    ## ... Similar to previous best
    ## Run 5 stress 0.1592963 
    ## ... Procrustes: rmse 0.000165495  max resid 0.002406413 
    ## ... Similar to previous best
    ## Run 6 stress 0.1733057 
    ## Run 7 stress 0.1592961 
    ## ... Procrustes: rmse 1.860249e-05  max resid 0.0002066058 
    ## ... Similar to previous best
    ## Run 8 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001743537  max resid 0.002441999 
    ## ... Similar to previous best
    ## Run 9 stress 0.1592961 
    ## ... Procrustes: rmse 1.545236e-05  max resid 8.839851e-05 
    ## ... Similar to previous best
    ## Run 10 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001698159  max resid 0.002423374 
    ## ... Similar to previous best
    ## Run 11 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001691222  max resid 0.002414183 
    ## ... Similar to previous best
    ## Run 12 stress 0.1592962 
    ## ... Procrustes: rmse 0.0001663904  max resid 0.002411799 
    ## ... Similar to previous best
    ## Run 13 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001755855  max resid 0.002458383 
    ## ... Similar to previous best
    ## Run 14 stress 0.1592962 
    ## ... Procrustes: rmse 0.00016031  max resid 0.00238173 
    ## ... Similar to previous best
    ## Run 15 stress 0.1732934 
    ## Run 16 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001713499  max resid 0.002432471 
    ## ... Similar to previous best
    ## Run 17 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001725873  max resid 0.002430775 
    ## ... Similar to previous best
    ## Run 18 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001740137  max resid 0.002437892 
    ## ... Similar to previous best
    ## Run 19 stress 0.1592961 
    ## ... Procrustes: rmse 3.151641e-05  max resid 0.0002162954 
    ## ... Similar to previous best
    ## Run 20 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001708145  max resid 0.002425823 
    ## ... Similar to previous best
    ## Run 21 stress 0.1732823 
    ## Run 22 stress 0.1732833 
    ## Run 23 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001706946  max resid 0.002428621 
    ## ... Similar to previous best
    ## Run 24 stress 0.4189868 
    ## Run 25 stress 0.1592961 
    ## ... Procrustes: rmse 2.195086e-05  max resid 0.0002022978 
    ## ... Similar to previous best
    ## Run 26 stress 0.1733567 
    ## Run 27 stress 0.1732952 
    ## Run 28 stress 0.1592962 
    ## ... Procrustes: rmse 0.0001557557  max resid 0.002347389 
    ## ... Similar to previous best
    ## Run 29 stress 0.1592962 
    ## ... Procrustes: rmse 0.0001617494  max resid 0.002388669 
    ## ... Similar to previous best
    ## Run 30 stress 0.1592961 
    ## ... Procrustes: rmse 2.530078e-05  max resid 0.0001445069 
    ## ... Similar to previous best
    ## Run 31 stress 0.1732927 
    ## Run 32 stress 0.1592962 
    ## ... Procrustes: rmse 0.0001646565  max resid 0.00240281 
    ## ... Similar to previous best
    ## Run 33 stress 0.1592967 
    ## ... Procrustes: rmse 9.480252e-05  max resid 0.001713659 
    ## ... Similar to previous best
    ## Run 34 stress 0.1592962 
    ## ... Procrustes: rmse 0.000159705  max resid 0.002380452 
    ## ... Similar to previous best
    ## Run 35 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001728705  max resid 0.002433273 
    ## ... Similar to previous best
    ## Run 36 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001641068  max resid 0.002400116 
    ## ... Similar to previous best
    ## Run 37 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001677178  max resid 0.002412935 
    ## ... Similar to previous best
    ## Run 38 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001708569  max resid 0.002448093 
    ## ... Similar to previous best
    ## Run 39 stress 0.1592961 
    ## ... New best solution
    ## ... Procrustes: rmse 9.105846e-06  max resid 6.860157e-05 
    ## ... Similar to previous best
    ## Run 40 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001742708  max resid 0.002446581 
    ## ... Similar to previous best
    ## Run 41 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001638928  max resid 0.002459091 
    ## ... Similar to previous best
    ## Run 42 stress 0.17328 
    ## Run 43 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001727211  max resid 0.002439813 
    ## ... Similar to previous best
    ## Run 44 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001751888  max resid 0.00244938 
    ## ... Similar to previous best
    ## Run 45 stress 0.1592961 
    ## ... Procrustes: rmse 8.908776e-06  max resid 9.451516e-05 
    ## ... Similar to previous best
    ## Run 46 stress 0.1592962 
    ## ... Procrustes: rmse 0.0001695503  max resid 0.002426309 
    ## ... Similar to previous best
    ## Run 47 stress 0.1592962 
    ## ... Procrustes: rmse 0.0001686835  max resid 0.002428108 
    ## ... Similar to previous best
    ## Run 48 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001728935  max resid 0.002440225 
    ## ... Similar to previous best
    ## Run 49 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001715283  max resid 0.002434706 
    ## ... Similar to previous best
    ## Run 50 stress 0.1592963 
    ## ... Procrustes: rmse 0.0001777584  max resid 0.002458186 
    ## ... Similar to previous best
    ## *** Best solution repeated 11 times

``` r
nmds_comm_trajectories
```

    ## 
    ## Call:
    ## metaMDS(comm = comm_trajectories_relative, distance = "bray",      k = 2, try = 100, trymax = 50, engine = "monoMDS", autotransform = FALSE) 
    ## 
    ## global Multidimensional Scaling using monoMDS
    ## 
    ## Data:     comm_trajectories_relative 
    ## Distance: bray 
    ## 
    ## Dimensions: 2 
    ## Stress:     0.1592961 
    ## Stress type 1, weak ties
    ## Best solution was repeated 11 times in 50 tries
    ## The best solution was from try 39 (random start)
    ## Scaling: centring, PC rotation, halfchange scaling 
    ## Species: expanded scores based on 'comm_trajectories_relative'

``` r
points_high_MDS1_high_MDS2 <- nmds_comm_trajectories$points[1:100,]
points_low_MDS1_high_MDS2 <- nmds_comm_trajectories$points[101:200,]
points_high_MDS1_low_MDS2 <- nmds_comm_trajectories$points[201:300,]
points_low_MDS1_low_MDS2 <- nmds_comm_trajectories$points[301:400,]

species_pred <- nmds_comm_trajectories$species
```

``` r
find_day <- function(day, min, max, length = 100){
  
  new_day <- day - min
  range <- max - min
  
  result <- round((new_day*length)/range)
  
  if(result == 0){result <- 1}
  
  return(result)
  
}
```

``` r
species <- nmds_comm_AM1_rm_total$species

sp_names <- rownames(species)
sp_names <- gsub("\\.", " ", sp_names)

sp_font <- rep(3, length(sp_names))
sp_font[sp_names == "Chironominae"] <- 1
sp_font[sp_names == "Tanypodinae"] <- 1

col_control <- "#56B4E9"
col_closed <- "#009E73"
col_exc_drag <- "#CC79A7"
col_exc_amph <- "#E69F00"

col_high_1_high_2 <- "#290AD8" #blue dark 
col_low_1_high_2 <- "#D82632" #red dark
col_high_1_low_2 <- "#72D9FF" #blue light 
col_low_1_low_2 <- "#FFAD72" #red light 

#svg(filename = "Plots/Community structure analysis/community_trajectories.svg", width = 10, height = 4.5, pointsize = 12)




close.screen(all.screens = TRUE)
split.screen(matrix(c(0  , 0.5, 0, 1  ,
                      0.5, 1  , 0, 1  ), ncol = 4, nrow = 2, byrow = TRUE))
```

    ## [1] 1 2

``` r
screen(1)
par(mar = c(4,4,1,1), bty = "o")


xmin <- min(c(comm_AM1_nmds[,1], species[,1]))*1.2
xmax <- max(c(comm_AM1_nmds[,1], species[,1]))*1.2
ymin <- min(c(comm_AM1_nmds[,2], species[,2]))*1.2
ymax <- max(c(comm_AM1_nmds[,2], species[,2]))*1.2




plot(comm_AM1_nmds[,1], comm_AM1_nmds[,2], xlim = c(xmin,xmax), ylim = c(ymin, ymax),
     type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")



rect(xleft = 0, ybottom = 0, xright = xmax*1.2, ytop = ymax*1.2 ,col = transparent(col_high_1_high_2, trans.val = 0.6), border = "white", lwd = 4)
rect(xleft = xmin*1.2, ybottom = ymin*1.2, xright = 0, ytop = 0,col = transparent(col_low_1_low_2, trans.val = 0.6), border = "white", lwd = 4)
rect(xleft = 0, ybottom = ymin*1.2, xright = xmax*1.2, ytop = 0,col = transparent(col_high_1_low_2, trans.val = 0.6), border = "white", lwd = 4)
rect(xleft = xmin*1.2, ybottom = 0, xright = 0, ytop = ymax*1.2,col = transparent(col_low_1_high_2, trans.val = 0.6), border = "white", lwd = 4)


abline(h = 0, v = 0, lty = 2)

library(scales)
library(shape)


#pal <- col_numeric(palette = c("white", "black"), domain = urb, na.color = "grey50", alpha = FALSE, reverse = FALSE)
#col <-pal(urb)

#points(comm_AM1_nmds[,1],comm_AM1_nmds[,2], col = "grey70", bg = "white", pch = 16, cex = 1.5)


points(comm_AM1_nmds[Exp_design_drop_atrasado$treatments == "controle",1],comm_AM1_nmds[Exp_design_drop_atrasado$treatments == "controle",2], col = "white", bg = col_control, pch = 21, cex = 1.5)

points(comm_AM1_nmds[Exp_design_drop_atrasado$treatments == "fechado",1],comm_AM1_nmds[Exp_design_drop_atrasado$treatments == "fechado",2], col = "white", bg = col_closed, pch = 21, cex = 1.5)

points(comm_AM1_nmds[Exp_design_drop_atrasado$treatments == "fechado_dia",1],comm_AM1_nmds[Exp_design_drop_atrasado$treatments == "fechado_dia",2], col = "white", bg = col_exc_drag, pch = 21, cex = 1.5)

points(comm_AM1_nmds[Exp_design_drop_atrasado$treatments == "fechado_noite",1],comm_AM1_nmds[Exp_design_drop_atrasado$treatments == "fechado_noite",2], col = "white", bg = col_exc_amph, pch = 21, cex = 1.5)

#Arrows(x0 <- rep(0, nrow(species)),
#       y0 <- rep(0, nrow(species)),
#       x1 <- species[,1],
#       y1 <- species[,2], arr.type = "triangle", arr.length = 0.4, col = "#C49C94", lwd = 1.5)


for(i in 1:length(sp_names)){
  text_contour(x = species[i,1],
       y = species[i,2], labels = sp_names[i], cex = 0.8, col = "white", thick = 0.004, font = sp_font[i])
  text(x = species[i,1], y = species[i,2], labels = sp_names[i], cex = 0.8, col = "black", font = sp_font[i])
}


axis(1, cex.axis = 1)
axis(2, cex.axis = 1, las = 2)
title(xlab = "nMDS 1", cex.lab = 1.2, line = 3)
title(ylab = "nMDS 2", cex.lab = 1.2, line = 3)
#title(main = "Morphological traits", line = 0.5, adj = 0, cex.main = 1.5)

par(new= TRUE, bty = "n")
plot(NA, ylab = "", xlab = "", xlim = c(0,100), ylim = c(0,100), xaxt = "n", yaxt = "n")

legend(x = 105, y = -5, xjust = 1, yjust = 0, legend = c("Control", "Both excluded", "Dragonfly exc.", "Amphibian exc."), bty = "n", bg = "transparent", pch = 21, pt.bg = c(col_control, col_closed, col_exc_drag, col_exc_amph), col = "white", cex = 0.8, pt.cex = 1.2)

letters(x = 0, y = 97, "a)", cex = 1.5)




screen(2)

#set.seed(1);species_pred <- jitter(species_pred, amount = 0.25)

xmin <- min(c(points_high_MDS1_high_MDS2[,1], points_low_MDS1_high_MDS2[,1], points_high_MDS1_low_MDS2[,1], points_low_MDS1_low_MDS2[,1],species_pred[,1]))*1.1
xmax <- max(c(points_high_MDS1_high_MDS2[,1], points_low_MDS1_high_MDS2[,1], points_high_MDS1_low_MDS2[,1], points_low_MDS1_low_MDS2[,1],species_pred[,1]))*1.1

ymin <- min(c(points_high_MDS1_high_MDS2[,2], points_low_MDS1_high_MDS2[,2], points_high_MDS1_low_MDS2[,2], points_low_MDS1_low_MDS2[,2],species_pred[,2]))*1.1
ymax <- max(c(points_high_MDS1_high_MDS2[,2], points_low_MDS1_high_MDS2[,2], points_high_MDS1_low_MDS2[,2], points_low_MDS1_low_MDS2[,2],species_pred[,2]))*1.1

par(mar = c(4, 4, 1, 1), bty = "o")
  
plot(NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = "", ylab = "", xaxt = "n", yaxt = "n")

lines(x = points_high_MDS1_high_MDS2[,1], y = points_high_MDS1_high_MDS2[,2], lwd = 3, lty = 1, col = col_high_1_high_2) #controle

lines(x = points_low_MDS1_high_MDS2[,1], y = points_low_MDS1_high_MDS2[,2], lwd = 3, lty = 1, col = col_low_1_high_2) #atrasado

lines(x = points_high_MDS1_low_MDS2[,1], y = points_high_MDS1_low_MDS2[,2], lwd = 3, lty = 1, col = col_high_1_low_2) #fechado

lines(x = points_low_MDS1_low_MDS2[,1], y = points_low_MDS1_low_MDS2[,2], lwd = 3, lty = 1, col = col_low_1_low_2) #fechado




#####################################################################



##################################################

points(x = points_high_MDS1_high_MDS2[find_day(day = 32, min = 32, max = 158),1],
       y = points_high_MDS1_high_MDS2[find_day(day = 32, min = 32, max = 158),2],
       pch = 21, bg = "white", col = col_high_1_high_2, lwd = 3, cex = 3.5)

text(x = points_high_MDS1_high_MDS2[find_day(day = 32, min = 32, max = 158),1],
       y = points_high_MDS1_high_MDS2[find_day(day = 32, min = 32, max = 158),2], labels = "32", cex = 0.7, col = col_high_1_high_2, font = 2)

points(x = points_high_MDS1_high_MDS2[find_day(day = 80, min = 32, max = 158),1],
       y = points_high_MDS1_high_MDS2[find_day(day = 80, min = 32, max = 158),2],
       pch = 21, bg = "white", col = col_high_1_high_2, lwd = 3, cex = 3.5)

text(x = points_high_MDS1_high_MDS2[find_day(day = 80, min = 32, max = 158),1],
       y = points_high_MDS1_high_MDS2[find_day(day = 80, min = 32, max = 158),2], labels = "80", cex = 0.7, col = col_high_1_high_2, font = 2)

points(x = points_high_MDS1_high_MDS2[find_day(day = 116, min = 32, max = 158),1],
       y = points_high_MDS1_high_MDS2[find_day(day = 116, min = 32, max = 158),2],
       pch = 21, bg = "white", col = col_high_1_high_2, lwd = 3, cex = 3.5)

text(x = points_high_MDS1_high_MDS2[find_day(day = 116, min = 32, max = 158),1],
       y = points_high_MDS1_high_MDS2[find_day(day = 116, min = 32, max = 158),2], labels = "116", cex = 0.7, col = col_high_1_high_2, font = 2)

points(x = points_high_MDS1_high_MDS2[find_day(day = 158, min = 32, max = 158),1],
       y = points_high_MDS1_high_MDS2[find_day(day = 158, min = 32, max = 158),2],
       pch = 21, bg = "white", col = col_high_1_high_2, lwd = 3, cex = 3.5)

text(x = points_high_MDS1_high_MDS2[find_day(day = 158, min = 32, max = 158),1],
       y = points_high_MDS1_high_MDS2[find_day(day = 158, min = 32, max = 158),2], labels = "158", cex = 0.7, col = col_high_1_high_2, font = 2)
#################################################


##################################################

points(x = points_low_MDS1_high_MDS2[find_day(day = 32, min = 32, max = 158),1],
       y = points_low_MDS1_high_MDS2[find_day(day = 32, min = 32, max = 158),2],
       pch = 21, bg = "white", col = col_low_1_high_2, lwd = 3, cex = 3.5)

text(x = points_low_MDS1_high_MDS2[find_day(day = 32, min = 32, max = 158),1],
       y = points_low_MDS1_high_MDS2[find_day(day = 32, min = 32, max = 158),2], labels = "32", cex = 0.7, font = 2, col = col_low_1_high_2)

points(x = points_low_MDS1_high_MDS2[find_day(day = 80, min = 32, max = 158),1],
       y = points_low_MDS1_high_MDS2[find_day(day = 80, min = 32, max = 158),2],
       pch = 21, bg = "white", col = col_low_1_high_2, lwd = 3, cex = 3.5)

text(x = points_low_MDS1_high_MDS2[find_day(day = 80, min = 32, max = 158),1],
       y = points_low_MDS1_high_MDS2[find_day(day = 80, min = 32, max = 158),2], labels = "80", cex = 0.7, font = 2, col = col_low_1_high_2)

points(x = points_low_MDS1_high_MDS2[find_day(day = 116, min = 32, max = 158),1],
       y = points_low_MDS1_high_MDS2[find_day(day = 116, min = 32, max = 158),2],
       pch = 21, bg = "white", col = col_low_1_high_2, lwd = 3, cex = 3.5)

text(x = points_low_MDS1_high_MDS2[find_day(day = 116, min = 32, max = 158),1],
       y = points_low_MDS1_high_MDS2[find_day(day = 116, min = 32, max = 158),2], labels = "116", cex = 0.7, font = 2, col = col_low_1_high_2)

points(x = points_low_MDS1_high_MDS2[find_day(day = 158, min = 32, max = 158),1],
       y = points_low_MDS1_high_MDS2[find_day(day = 158, min = 32, max = 158),2],
       pch = 21, bg = "white", col = col_low_1_high_2, lwd = 3, cex = 3.5)

text(x = points_low_MDS1_high_MDS2[find_day(day = 158, min = 32, max = 158),1],
       y = points_low_MDS1_high_MDS2[find_day(day = 158, min = 32, max = 158),2], labels = "158", cex = 0.7, font = 2, col = col_low_1_high_2)
#################################################



##################################################

points(x = points_high_MDS1_low_MDS2[find_day(day = 32, min = 32, max = 158),1],
       y = points_high_MDS1_low_MDS2[find_day(day = 32, min = 32, max = 158),2],
       pch = 21, bg = "white", col = col_high_1_low_2, lwd = 3, cex = 3.5)

text(x = points_high_MDS1_low_MDS2[find_day(day = 32, min = 32, max = 158),1],
       y = points_high_MDS1_low_MDS2[find_day(day = 32, min = 32, max = 158),2], labels = "32", cex = 0.7, font = 2, col = col_high_1_low_2)

points(x = points_high_MDS1_low_MDS2[find_day(day = 80, min = 32, max = 158),1],
       y = points_high_MDS1_low_MDS2[find_day(day = 80, min = 32, max = 158),2],
       pch = 21, bg = "white", col = col_high_1_low_2, lwd = 3, cex = 3.5)

text(x = points_high_MDS1_low_MDS2[find_day(day = 80, min = 32, max = 158),1],
       y = points_high_MDS1_low_MDS2[find_day(day = 80, min = 32, max = 158),2], labels = "80", cex = 0.7, font = 2, col = col_high_1_low_2)

points(x = points_high_MDS1_low_MDS2[find_day(day = 116, min = 32, max = 158),1],
       y = points_high_MDS1_low_MDS2[find_day(day = 116, min = 32, max = 158),2],
       pch = 21, bg = "white", col = col_high_1_low_2, lwd = 3, cex = 3.5)

text(x = points_high_MDS1_low_MDS2[find_day(day = 116, min = 32, max = 158),1],
       y = points_high_MDS1_low_MDS2[find_day(day = 116, min = 32, max = 158),2], labels = "116", cex = 0.7, font = 2, col = col_high_1_low_2)

points(x = points_high_MDS1_low_MDS2[find_day(day = 158, min = 32, max = 158),1],
       y = points_high_MDS1_low_MDS2[find_day(day = 158, min = 32, max = 158),2],
       pch = 21, bg = "white", col = col_high_1_low_2, lwd = 3, cex = 3.5)

text(x = points_high_MDS1_low_MDS2[find_day(day = 158, min = 32, max = 158),1],
       y = points_high_MDS1_low_MDS2[find_day(day = 158, min = 32, max = 158),2], labels = "158", cex = 0.7, font = 2, col = col_high_1_low_2)
#################################################





##################################################

points(x = points_low_MDS1_low_MDS2[find_day(day = 32, min = 32, max = 158),1],
       y = points_low_MDS1_low_MDS2[find_day(day = 32, min = 32, max = 158),2],
       pch = 21, bg = "white", col = col_low_1_low_2, lwd = 3, cex = 3.5)

text(x = points_low_MDS1_low_MDS2[find_day(day = 32, min = 32, max = 158),1],
       y = points_low_MDS1_low_MDS2[find_day(day = 32, min = 32, max = 158),2], labels = "32", cex = 0.7, font = 2, col = col_low_1_low_2)

points(x = points_low_MDS1_low_MDS2[find_day(day = 80, min = 32, max = 158),1],
       y = points_low_MDS1_low_MDS2[find_day(day = 80, min = 32, max = 158),2],
       pch = 21, bg = "white", col = col_low_1_low_2, lwd = 3, cex = 3.5)

text(x = points_low_MDS1_low_MDS2[find_day(day = 80, min = 32, max = 158),1],
       y = points_low_MDS1_low_MDS2[find_day(day = 80, min = 32, max = 158),2], labels = "80", cex = 0.7, font = 2, col = col_low_1_low_2)

points(x = points_low_MDS1_low_MDS2[find_day(day = 116, min = 32, max = 158),1],
       y = points_low_MDS1_low_MDS2[find_day(day = 116, min = 32, max = 158),2],
       pch = 21, bg = "white", col = col_low_1_low_2, lwd = 3, cex = 3.5)

text(x = points_low_MDS1_low_MDS2[find_day(day = 116, min = 32, max = 158),1],
       y = points_low_MDS1_low_MDS2[find_day(day = 116, min = 32, max = 158),2], labels = "116", cex = 0.7, font = 2, col = col_low_1_low_2)

points(x = points_low_MDS1_low_MDS2[find_day(day = 158, min = 32, max = 158),1],
       y = points_low_MDS1_low_MDS2[find_day(day = 158, min = 32, max = 158),2],
       pch = 21, bg = "white", col = col_low_1_low_2, lwd = 3, cex = 3.5)

text(x = points_low_MDS1_low_MDS2[find_day(day = 158, min = 32, max = 158),1],
       y = points_low_MDS1_low_MDS2[find_day(day = 158, min = 32, max = 158),2], labels = "158", cex = 0.7, font = 2, col = col_low_1_low_2)
#################################################



########################## species_pred #################################


sp_names <- rownames(species_pred)
sp_names <- substr(sp_names, start = 1, stop = 3)


for(i in 1:nrow(species_pred)){
  #points(x = species_pred[i,1],
   #    y = species_pred[i,2],
    #   pch = 21, bg = "white", col = "grey50", lwd = 2, cex = 4)

text_contour(x = species_pred[i,1],
       y = species_pred[i,2], labels = sp_names[i], cex = 0.8, col = "white", thick = 0.004)
  
text(x = species_pred[i,1],
       y = species_pred[i,2], labels = sp_names[i], cex = 0.8, col = "grey30")
}




axis(1)
axis(2, cex.axis = 1, las = 2)

title(xlab = "nMDS 1", line = 3, cex.lab = 1.2)
title(ylab = "nMDS 2", line = 3, cex.lab = 1.2)

letters(x = 0, y = 97, "b)", cex = 1.5)
```

![](Effect_first_survey_on_time_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
#par(new = TRUE)
#plot(NA, xlim = c(0, 100), ylim = c(0, 100), xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")

#legend(x = 100, y = 100, xjust = 1, yjust = 1, lty = c(1,1,1), lwd = 3, legend = c("Control", "Closed", "Late"), col = c("#009E73", "#56B4E9","#D55E00"), bty = "n")

#dev.off()
```

## Using treatments to describe initial communities.

### Testing for effects

``` r
nrow(predictors_nmds)
```

    ## [1] 96

``` r
nrow(comm_all_drop_atrasado_rm)
```

    ## [1] 96

``` r
#comm_all_drop_atrasado_rm2 <- comm_all_drop_atrasado_rm[,colnames(comm_all_drop_atrasado_rm) != "Microvelia"]


set.seed(3)
control <- permute::how(within = permute::Within(type = 'free'),
                        plots = Plots(strata = predictors_nmds$sites, type = 'free'),
                        nperm = 999)
permutations <- shuffleSet(nrow(comm_all_drop_atrasado_rm), control = control)


comm_all_drop_atrasado_mvabund <- mvabund(comm_all_drop_atrasado_rm)


mod_full_lin <- manyglm(comm_all_drop_atrasado_mvabund ~ block2 + treatments  + AM_days_sc + AM_days_sc:treatments, family="negative.binomial",  data = predictors_nmds, composition = FALSE)


mod_full_quad <- manyglm(comm_all_drop_atrasado_mvabund ~ block2 + treatments + AM_days_sc + AM_days_sc_squared + AM_days_sc:treatments + AM_days_sc_squared:treatments, family="negative.binomial",  data = predictors_nmds, composition = FALSE)

anova_quadratic_test <- anova(mod_full_lin, mod_full_quad, bootID = permutations, show.time = "all", cor.type = "I", test = "LR", resamp = "pit.trap")
```

    ## Using <int> bootID matrix from input. 
    ## Resampling begins for test 1.
    ##  Resampling run 0 finished. Time elapsed: 0.00 minutes...
    ##  Resampling run 100 finished. Time elapsed: 0.04 minutes...
    ##  Resampling run 200 finished. Time elapsed: 0.08 minutes...
    ##  Resampling run 300 finished. Time elapsed: 0.11 minutes...
    ##  Resampling run 400 finished. Time elapsed: 0.15 minutes...
    ##  Resampling run 500 finished. Time elapsed: 0.19 minutes...
    ##  Resampling run 600 finished. Time elapsed: 0.23 minutes...
    ##  Resampling run 700 finished. Time elapsed: 0.27 minutes...
    ##  Resampling run 800 finished. Time elapsed: 0.31 minutes...
    ##  Resampling run 900 finished. Time elapsed: 0.34 minutes...
    ## Time elapsed: 0 hr 0 min 22 sec

``` r
anova_quadratic_test
```

    ## Analysis of Deviance Table
    ## 
    ## mod_full_lin: comm_all_drop_atrasado_mvabund ~ block2 + treatments + AM_days_sc + AM_days_sc:treatments
    ## mod_full_quad: comm_all_drop_atrasado_mvabund ~ block2 + treatments + AM_days_sc + AM_days_sc_squared + AM_days_sc:treatments + AM_days_sc_squared:treatments
    ## 
    ## Multivariate test:
    ##               Res.Df Df.diff   Dev Pr(>Dev)    
    ## mod_full_lin      86                           
    ## mod_full_quad     82       4 118.3    0.001 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Arguments:
    ##  Test statistics calculated assuming uncorrelated response (for faster computation) 
    ##  P-value calculated using 999 iterations via PIT-trap resampling.

``` r
########################################################
mod_block <- manyglm(comm_all_drop_atrasado_mvabund ~ block2 + treatments + AM_days_sc + AM_days_sc_squared, family="negative.binomial",  data = predictors_nmds, composition = FALSE)

mod_no_block <- manyglm(comm_all_drop_atrasado_mvabund ~ treatments + AM_days_sc + AM_days_sc_squared, family="negative.binomial",  data = predictors_nmds, composition = FALSE)

anova_block <- anova(mod_no_block, mod_block, bootID = permutations, show.time = "all", cor.type = "I", test = "LR", resamp = "pit.trap")
```

    ## Using <int> bootID matrix from input. 
    ## Resampling begins for test 1.
    ##  Resampling run 0 finished. Time elapsed: 0.00 minutes...
    ##  Resampling run 100 finished. Time elapsed: 0.03 minutes...
    ##  Resampling run 200 finished. Time elapsed: 0.06 minutes...
    ##  Resampling run 300 finished. Time elapsed: 0.10 minutes...
    ##  Resampling run 400 finished. Time elapsed: 0.13 minutes...
    ##  Resampling run 500 finished. Time elapsed: 0.16 minutes...
    ##  Resampling run 600 finished. Time elapsed: 0.19 minutes...
    ##  Resampling run 700 finished. Time elapsed: 0.23 minutes...
    ##  Resampling run 800 finished. Time elapsed: 0.26 minutes...
    ##  Resampling run 900 finished. Time elapsed: 0.29 minutes...
    ## Time elapsed: 0 hr 0 min 19 sec

``` r
anova_block
```

    ## Analysis of Deviance Table
    ## 
    ## mod_no_block: comm_all_drop_atrasado_mvabund ~ treatments + AM_days_sc + AM_days_sc_squared
    ## mod_block: comm_all_drop_atrasado_mvabund ~ block2 + treatments + AM_days_sc + AM_days_sc_squared
    ## 
    ## Multivariate test:
    ##              Res.Df Df.diff   Dev Pr(>Dev)   
    ## mod_no_block     90                          
    ## mod_block        88       2 79.58    0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Arguments:
    ##  Test statistics calculated assuming uncorrelated response (for faster computation) 
    ##  P-value calculated using 999 iterations via PIT-trap resampling.

``` r
########################################################






mod_time <- manyglm(comm_all_drop_atrasado_mvabund ~ block2 + treatments + AM_days_sc + AM_days_sc_squared, family="negative.binomial",  data = predictors_nmds, composition = FALSE)

mod_no_time <- manyglm(comm_all_drop_atrasado_mvabund ~ block2 + treatments, family="negative.binomial",  data = predictors_nmds, composition = FALSE)

anova_time <- anova(mod_no_time, mod_time, bootID = permutations, show.time = "all", cor.type = "I", test = "LR", resamp = "pit.trap")
```

    ## Using <int> bootID matrix from input. 
    ## Resampling begins for test 1.
    ##  Resampling run 0 finished. Time elapsed: 0.00 minutes...
    ##  Resampling run 100 finished. Time elapsed: 0.03 minutes...
    ##  Resampling run 200 finished. Time elapsed: 0.10 minutes...
    ##  Resampling run 300 finished. Time elapsed: 0.16 minutes...
    ##  Resampling run 400 finished. Time elapsed: 0.21 minutes...
    ##  Resampling run 500 finished. Time elapsed: 0.26 minutes...
    ##  Resampling run 600 finished. Time elapsed: 0.31 minutes...
    ##  Resampling run 700 finished. Time elapsed: 0.37 minutes...
    ##  Resampling run 800 finished. Time elapsed: 0.42 minutes...
    ##  Resampling run 900 finished. Time elapsed: 0.47 minutes...
    ## Time elapsed: 0 hr 0 min 31 sec

``` r
anova_time
```

    ## Analysis of Deviance Table
    ## 
    ## mod_no_time: comm_all_drop_atrasado_mvabund ~ block2 + treatments
    ## mod_time: comm_all_drop_atrasado_mvabund ~ block2 + treatments + AM_days_sc + AM_days_sc_squared
    ## 
    ## Multivariate test:
    ##             Res.Df Df.diff   Dev Pr(>Dev)    
    ## mod_no_time     90                           
    ## mod_time        88       2 254.5    0.001 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Arguments:
    ##  Test statistics calculated assuming uncorrelated response (for faster computation) 
    ##  P-value calculated using 999 iterations via PIT-trap resampling.

``` r
########################################################



mod_prior <- manyglm(comm_all_drop_atrasado_mvabund ~ block2 + treatments + AM_days_sc + AM_days_sc_squared, family="negative.binomial",  data = predictors_nmds, composition = FALSE)

mod_no_prior <- manyglm(comm_all_drop_atrasado_mvabund ~ block2 + AM_days_sc + AM_days_sc_squared, family="negative.binomial",  data = predictors_nmds, composition = FALSE)

anova_prior <- anova(mod_no_prior, mod_prior, bootID = permutations, show.time = "all", cor.type = "I", test = "LR", resamp = "pit.trap")
```

    ## Using <int> bootID matrix from input. 
    ## Resampling begins for test 1.
    ##  Resampling run 0 finished. Time elapsed: 0.00 minutes...
    ##  Resampling run 100 finished. Time elapsed: 0.05 minutes...
    ##  Resampling run 200 finished. Time elapsed: 0.10 minutes...
    ##  Resampling run 300 finished. Time elapsed: 0.15 minutes...
    ##  Resampling run 400 finished. Time elapsed: 0.20 minutes...
    ##  Resampling run 500 finished. Time elapsed: 0.25 minutes...
    ##  Resampling run 600 finished. Time elapsed: 0.30 minutes...
    ##  Resampling run 700 finished. Time elapsed: 0.34 minutes...
    ##  Resampling run 800 finished. Time elapsed: 0.39 minutes...
    ##  Resampling run 900 finished. Time elapsed: 0.43 minutes...
    ## Time elapsed: 0 hr 0 min 28 sec

``` r
anova_prior
```

    ## Analysis of Deviance Table
    ## 
    ## mod_no_prior: comm_all_drop_atrasado_mvabund ~ block2 + AM_days_sc + AM_days_sc_squared
    ## mod_prior: comm_all_drop_atrasado_mvabund ~ block2 + treatments + AM_days_sc + AM_days_sc_squared
    ## 
    ## Multivariate test:
    ##              Res.Df Df.diff   Dev Pr(>Dev)  
    ## mod_no_prior     91                         
    ## mod_prior        88       3 68.05    0.046 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Arguments:
    ##  Test statistics calculated assuming uncorrelated response (for faster computation) 
    ##  P-value calculated using 999 iterations via PIT-trap resampling.

``` r
########################################################
mod_no_int <- manyglm(comm_all_drop_atrasado_mvabund ~ block2 + treatments + AM_days_sc + AM_days_sc_squared, family="negative.binomial",  data = predictors_nmds, composition = FALSE)

mod_int <- manyglm(comm_all_drop_atrasado_mvabund ~ block2 + treatments + AM_days_sc + AM_days_sc_squared + AM_days_sc:treatments + AM_days_sc_squared:treatments, family="negative.binomial",  data = predictors_nmds, composition = FALSE)

anova_interaction <- anova(mod_no_int, mod_int, bootID = permutations, show.time = "all", cor.type = "I", test = "LR", resamp = "pit.trap")
```

    ## Using <int> bootID matrix from input. 
    ## Resampling begins for test 1.
    ##  Resampling run 0 finished. Time elapsed: 0.00 minutes...
    ##  Resampling run 100 finished. Time elapsed: 0.05 minutes...
    ##  Resampling run 200 finished. Time elapsed: 0.11 minutes...
    ##  Resampling run 300 finished. Time elapsed: 0.17 minutes...
    ##  Resampling run 400 finished. Time elapsed: 0.23 minutes...
    ##  Resampling run 500 finished. Time elapsed: 0.28 minutes...
    ##  Resampling run 600 finished. Time elapsed: 0.33 minutes...
    ##  Resampling run 700 finished. Time elapsed: 0.38 minutes...
    ##  Resampling run 800 finished. Time elapsed: 0.42 minutes...
    ##  Resampling run 900 finished. Time elapsed: 0.46 minutes...
    ## Time elapsed: 0 hr 0 min 30 sec

``` r
anova_interaction
```

    ## Analysis of Deviance Table
    ## 
    ## mod_no_int: comm_all_drop_atrasado_mvabund ~ block2 + treatments + AM_days_sc + AM_days_sc_squared
    ## mod_int: comm_all_drop_atrasado_mvabund ~ block2 + treatments + AM_days_sc + AM_days_sc_squared + AM_days_sc:treatments + AM_days_sc_squared:treatments
    ## 
    ## Multivariate test:
    ##            Res.Df Df.diff   Dev Pr(>Dev)
    ## mod_no_int     88                       
    ## mod_int        82       6 95.53    0.108
    ## Arguments:
    ##  Test statistics calculated assuming uncorrelated response (for faster computation) 
    ##  P-value calculated using 999 iterations via PIT-trap resampling.

``` r
########################################################


## Just out of curiosity, lets see if there is any interaction effec of time effects with blocks
########################################################
mod_block <- manyglm(comm_all_drop_atrasado_mvabund ~ block2 + treatments + AM_days_sc + AM_days_sc_squared, family="negative.binomial",  data = predictors_nmds, composition = FALSE)

mod_block_int <- manyglm(comm_all_drop_atrasado_mvabund ~ block2 + treatments + AM_days_sc + AM_days_sc_squared + AM_days_sc:block2 + AM_days_sc_squared:block2, family="negative.binomial",  data = predictors_nmds, composition = FALSE)

anova_block <- anova(mod_block, mod_block_int, bootID = permutations, show.time = "all", cor.type = "I", test = "LR", resamp = "pit.trap")
```

    ## Using <int> bootID matrix from input. 
    ## Resampling begins for test 1.
    ##  Resampling run 0 finished. Time elapsed: 0.00 minutes...
    ##  Resampling run 100 finished. Time elapsed: 0.05 minutes...
    ##  Resampling run 200 finished. Time elapsed: 0.09 minutes...
    ##  Resampling run 300 finished. Time elapsed: 0.13 minutes...
    ##  Resampling run 400 finished. Time elapsed: 0.18 minutes...
    ##  Resampling run 500 finished. Time elapsed: 0.23 minutes...
    ##  Resampling run 600 finished. Time elapsed: 0.29 minutes...
    ##  Resampling run 700 finished. Time elapsed: 0.36 minutes...
    ##  Resampling run 800 finished. Time elapsed: 0.43 minutes...
    ##  Resampling run 900 finished. Time elapsed: 0.48 minutes...
    ## Time elapsed: 0 hr 0 min 30 sec

``` r
anova_block
```

    ## Analysis of Deviance Table
    ## 
    ## mod_block: comm_all_drop_atrasado_mvabund ~ block2 + treatments + AM_days_sc + AM_days_sc_squared
    ## mod_block_int: comm_all_drop_atrasado_mvabund ~ block2 + treatments + AM_days_sc + AM_days_sc_squared + AM_days_sc:block2 + AM_days_sc_squared:block2
    ## 
    ## Multivariate test:
    ##               Res.Df Df.diff   Dev Pr(>Dev)
    ## mod_block         88                       
    ## mod_block_int     84       4 53.03    0.435
    ## Arguments:
    ##  Test statistics calculated assuming uncorrelated response (for faster computation) 
    ##  P-value calculated using 999 iterations via PIT-trap resampling.

``` r
########################################################
```

### Making predictions for plotting

``` r
comm_all_drop_atrasado_mvabund2 <- comm_all_drop_atrasado_mvabund + 0.001

mod_int <- manyglm(comm_all_drop_atrasado_mvabund2 ~ block2 + treatments + AM_days_sc + AM_days_sc_squared + AM_days_sc:treatments + AM_days_sc_squared:treatments, family="negative.binomial",  data = predictors_nmds, composition = FALSE)
```

    ## Warning in manyglm(comm_all_drop_atrasado_mvabund2 ~ block2 + treatments + :
    ## Non-integer data are fitted to the negative.binomial model.

``` r
#High abundances of 
#Aedes
#Culex
#Dendropsophus_minutus
#Tanypodinae
#Callibaetis
#Culex

new_data_AB_control <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("AB", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
                                treatments = factor(rep("controle", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$treatments))))

new_data_CD_control <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("CD", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
             treatments = factor(rep("controle", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$treatments))))


new_data_EF_control <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("EF", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
             treatments = factor(rep("controle", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$treatments))))


###########################################################################################



new_data_AB_fechado <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("AB", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
                                treatments = factor(rep("fechado", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$treatments))))

new_data_CD_fechado <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("CD", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
             treatments = factor(rep("fechado", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$treatments))))


new_data_EF_fechado <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("EF", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
             treatments = factor(rep("fechado", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$treatments))))



###########################################################################################


new_data_AB_fechado_dia <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("AB", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
                                treatments = factor(rep("fechado_dia", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$treatments))))

new_data_CD_fechado_dia <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("CD", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
             treatments = factor(rep("fechado_dia", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$treatments))))


new_data_EF_fechado_dia <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("EF", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
             treatments = factor(rep("fechado_dia", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$treatments))))


###########################################################################################




new_data_AB_fechado_noite <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("AB", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
                                treatments = factor(rep("fechado_noite", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$treatments))))

new_data_CD_fechado_noite <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("CD", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
             treatments = factor(rep("fechado_noite", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$treatments))))


new_data_EF_fechado_noite <- 
  data.frame(AM_days_sc = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100),
                                AM_days_sc_squared = seq(from = min(Exp_design_all_drop_atrasado$AM_days_sc), to = max(Exp_design_all_drop_atrasado$AM_days_sc), length.out = 100)^2,
                                block2 = factor(rep("EF", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$block2))),
             treatments = factor(rep("fechado_noite", 100), levels = levels(as.factor(Exp_design_all_drop_atrasado$treatments))))


###########################################################################################





###########################################################################################


AB_control <- predict.manyglm(mod_int, newdata = new_data_AB_control, type = "response")
CD_control <- predict.manyglm(mod_int, newdata = new_data_CD_control, type = "response")
EF_control <- predict.manyglm(mod_int, newdata = new_data_EF_control, type = "response")

control <- array(NA, dim = c(nrow(AB_control), ncol(AB_control), 3))
control[,,1] <- AB_control
control[,,2] <- CD_control
control[,,3] <- EF_control

control_mean <- as.data.frame(apply(control, MARGIN = c(1,2), FUN = mean))
colnames(control_mean) <- colnames(AB_control)
control_mean
```

    ##           Aedes  Callibaetis Chironominae     Culex Dendropsophus.minutus
    ## 1   29.42870720 0.7988902358     13.06247 60.233506             0.2363311
    ## 2   30.16782366 0.8423019713     12.90390 55.680777             0.2556978
    ## 3   30.86306143 0.8858778542     12.75327 51.542587             0.2763616
    ## 4   31.51056942 0.9294054118     12.61033 47.777226             0.2983821
    ## 5   32.10670406 0.9726618292     12.47488 44.347530             0.3218195
    ## 6   32.64806334 1.0154156972     12.34669 41.220355             0.3467340
    ## 7   33.13151912 1.0574289512     12.22559 38.366113             0.3731857
    ## 8   33.55424720 1.0984589825     12.11138 35.758366             0.4012344
    ## 9   33.91375481 1.1382608966     12.00388 33.373466             0.4309389
    ## 10  34.20790503 1.1765898901     11.90296 31.190242             0.4623574
    ## 11  34.43493798 1.2132037170     11.80844 29.189723             0.4955464
    ## 12  34.59348839 1.2478652076     11.72019 27.354891             0.5305611
    ## 13  34.68259918 1.2803448058     11.63809 25.670468             0.5674544
    ## 14  34.70173118 1.3104230867     11.56201 24.122725             0.6062770
    ## 15  34.65076850 1.3378932160     11.49184 22.699314             0.6470765
    ## 16  34.53001974 1.3625633131     11.42748 21.389119             0.6898977
    ## 17  34.34021483 1.3842586793     11.36884 20.182123             0.7347816
    ## 18  34.08249773 1.4028238568     11.31583 19.069292             0.7817652
    ## 19  33.75841493 1.4181244826     11.26837 18.042474             0.8308811
    ## 20  33.36989998 1.4300489076     11.22641 17.094303             0.8821571
    ## 21  32.91925428 1.4385095529     11.18987 16.218119             0.9356156
    ## 22  32.40912440 1.4434439816     11.15870 15.407897             0.9912735
    ## 23  31.84247618 1.4448156654     11.13287 14.658180             1.0491415
    ## 24  31.22256606 1.4426144352     11.11233 13.964021             1.1092236
    ## 25  30.55290990 1.4368566069     11.09706 13.320936             1.1715171
    ## 26  29.83724985 1.4275847795     11.08704 12.724852             1.2360119
    ## 27  29.07951965 1.4148673107     11.08224 12.172074             1.3026903
    ## 28  28.28380871 1.3987974773     11.08267 11.659238             1.3715265
    ## 29  27.45432560 1.3794923362     11.08832 11.183289             1.4424864
    ## 30  26.59536127 1.3570913044     11.09920 10.741446             1.5155273
    ## 31  25.71125254 1.3317544831     11.11533 10.331174             1.5905974
    ## 32  24.80634614 1.3036607527     11.13673  9.950169             1.6676361
    ## 33  23.88496397 1.2730056721     11.16343  9.596326             1.7465732
    ## 34  22.95136963 1.2399992154     11.19547  9.267728             1.8273293
    ## 35  22.00973682 1.2048633831     11.23289  8.962629             1.9098151
    ## 36  21.06411987 1.1678297260     11.27575  8.679432             1.9939319
    ## 37  20.11842652 1.1291368191     11.32411  8.416683             2.0795713
    ## 38  19.17639339 1.0890277260     11.37804  8.173056             2.1666154
    ## 39  18.24156417 1.0477474882     11.43761  7.947338             2.2549365
    ## 40  17.31727069 1.0055406775     11.50291  7.738428             2.3443978
    ## 41  16.40661707 0.9626490432     11.57404  7.545318             2.4348533
    ## 42  15.51246674 0.9193092862     11.65110  7.367094             2.5261479
    ## 43  14.63743268 0.8757509853     11.73419  7.202920             2.6181183
    ## 44  13.78387048 0.8321947011     11.82346  7.052040             2.7105925
    ## 45  12.95387441 0.7888502757     11.91901  6.913767             2.8033913
    ## 46  12.14927623 0.7459153444     12.02100  6.787479             2.8963276
    ## 47  11.37164664 0.7035740690     12.12958  6.672615             2.9892082
    ## 48  10.62229923 0.6619961020     12.24491  6.568669             3.0818332
    ## 49   9.90229665 0.6213357803     12.36716  6.475190             3.1739976
    ## 50   9.21245881 0.5817315512     12.49651  6.391775             3.2654915
    ## 51   8.55337300 0.5433056215     12.63318  6.318066             3.3561009
    ## 52   7.92540545 0.5061638235     12.77735  6.253752             3.4456088
    ## 53   7.32871432 0.4703956853     12.92927  6.198561             3.5337956
    ## 54   6.76326372 0.4360746902     13.08915  6.152264             3.6204401
    ## 55   6.22883859 0.4032587099     13.25726  6.114667             3.7053208
    ## 56   5.72506019 0.3719905916     13.43386  6.085614             3.7882162
    ## 57   5.25140198 0.3422988808     13.61923  6.064987             3.8689062
    ## 58   4.80720558 0.3141986586     13.81366  6.052699             3.9471728
    ## 59   4.39169683 0.2876924720     14.01747  6.048700             4.0228012
    ## 60   4.00400153 0.2627713385     14.23099  6.052974             4.0955809
    ## 61   3.64316083 0.2394158027     14.45458  6.065539             4.1653062
    ## 62   3.30814607 0.2175970264     14.68859  6.086446             4.2317779
    ## 63   2.99787305 0.1972778942     14.93343  6.115780             4.2948034
    ## 64   2.71121540 0.1784141169     15.18951  6.153664             4.3541983
    ## 65   2.44701727 0.1609553179     15.45726  6.200255             4.4097871
    ## 66   2.20410495 0.1448460877     15.73715  6.255745             4.4614039
    ## 67   1.98129769 0.1300269970     16.02966  6.320367             4.5088932
    ## 68   1.77741737 0.1164355556     16.33529  6.394394             4.5521112
    ## 69   1.59129733 0.1040071105     16.65461  6.478139             4.5909257
    ## 70   1.42179006 0.0926756768     16.98816  6.571960             4.6252176
    ## 71   1.26777399 0.0823746962     17.33657  6.676261             4.6548808
    ## 72   1.12815930 0.0730377203     17.70046  6.791498             4.6798234
    ## 73   1.00189277 0.0645990180     18.08050  6.918175             4.6999675
    ## 74   0.88796181 0.0569941062     18.47741  7.056858             4.7152503
    ## 75   0.78539757 0.0501602050     18.89194  7.208169             4.7256237
    ## 76   0.69327740 0.0440366204     19.32486  7.372797             4.7310553
    ## 77   0.61072650 0.0385650577     19.77703  7.551504             4.7315279
    ## 78   0.53691892 0.0336898682     20.24932  7.745124             4.7270401
    ## 79   0.47107806 0.0293582364     20.74265  7.954578             4.7176059
    ## 80   0.41247656 0.0255203096     21.25802  8.180873             4.7032551
    ## 81   0.36043578 0.0221292777     21.79646  8.425117             4.6840325
    ## 82   0.31432487 0.0191414071     22.35908  8.688525             4.6599984
    ## 83   0.27355951 0.0165160355     22.94703  8.972427             4.6312277
    ## 84   0.23760037 0.0142155308     23.56154  9.278282             4.5978097
    ## 85   0.20595135 0.0122052218     24.20390  9.607691             4.5598479
    ## 86   0.17815761 0.0104533043     24.87550  9.962407             4.5174589
    ## 87   0.15380355 0.0089307275     25.57779 10.344352             4.4707724
    ## 88   0.13251058 0.0076110645     26.31230 10.755636             4.4199302
    ## 89   0.11393495 0.0064703727     27.08066 11.198573             4.3650855
    ## 90   0.09776550 0.0054870448     27.88459 11.675704             4.3064023
    ## 91   0.08372141 0.0046416570     28.72592 12.189819             4.2440544
    ## 92   0.07155000 0.0039168137     29.60659 12.743984             4.1782246
    ## 93   0.06102460 0.0032969936     30.52864 13.341570             4.1091038
    ## 94   0.05194246 0.0027683985     31.49424 13.986289             4.0368903
    ## 95   0.04412272 0.0023188062     32.50569 14.682223             3.9617884
    ## 96   0.03740454 0.0019374281     33.56544 15.433872             3.8840079
    ## 97   0.03164525 0.0016147753     34.67608 16.246200             3.8037628
    ## 98   0.02671868 0.0013425297     35.84036 17.124680             3.7212705
    ## 99   0.02251353 0.0011134251     37.06118 18.075359             3.6367509
    ## 100  0.01893192 0.0009211354     38.34165 19.104919             3.5504252
    ##       Microvelia Scinax.fuscovarius Tanypodinae    Dasyhelea Erythrodiplax
    ## 1   0.2027767629         73.1628564    1.787848 0.0001321673    0.05899589
    ## 2   0.1578266813         72.1563950    1.846936 0.0001690972    0.06146980
    ## 3   0.1235033767         71.1203910    1.906979 0.0002156636    0.06403111
    ## 4   0.0971657949         70.0565227    1.967947 0.0002741863    0.06668212
    ## 5   0.0768571281         68.9664944    2.029803 0.0003474905    0.06942515
    ## 6   0.0611210865         67.8520321    2.092511 0.0004390038    0.07226256
    ## 7   0.0488690783         66.7148785    2.156029 0.0005528686    0.07519674
    ## 8   0.0392837916         65.5567890    2.220315 0.0006940709    0.07823009
    ## 9   0.0317489099         64.3795268    2.285324 0.0008685884    0.08136501
    ## 10  0.0257976654         63.1848590    2.351008 0.0010835587    0.08460396
    ## 11  0.0210750281         61.9745517    2.417318 0.0013474701    0.08794939
    ## 12  0.0173098015         60.7503662    2.484199 0.0016703754    0.09140376
    ## 13  0.0142939475         59.5140547    2.551598 0.0020641313    0.09496955
    ## 14  0.0118672050         58.2673561    2.619457 0.0025426630    0.09864926
    ## 15  0.0099056020         57.0119923    2.687716 0.0031222564    0.10244539
    ## 16  0.0083128410         55.7496643    2.756314 0.0038218759    0.10636043
    ## 17  0.0070138140         54.4820485    2.825187 0.0046635094    0.11039690
    ## 18  0.0059497016         53.2107933    2.894269 0.0056725378    0.11455730
    ## 19  0.0050742550         51.9375157    2.963491 0.0068781272    0.11884415
    ## 20  0.0043509647         50.6637980    3.032785 0.0083136412    0.12325993
    ## 21  0.0037508959         49.3911852    3.102078 0.0100170677    0.12780715
    ## 22  0.0032510273         48.1211819    3.171298 0.0120314565    0.13248830
    ## 23  0.0028329728         46.8552497    3.240369 0.0144053583    0.13730583
    ## 24  0.0024819918         45.5948048    3.309216 0.0171932584    0.14226221
    ## 25  0.0021862230         44.3412159    3.377761 0.0204559935    0.14735987
    ## 26  0.0019360865         43.0958020    3.445925 0.0242611413    0.15260123
    ## 27  0.0017238172         41.8598310    3.513629 0.0286833678    0.15798868
    ## 28  0.0015430992         40.6345175    3.580793 0.0338047182    0.16352456
    ## 29  0.0013887774         39.4210218    3.647334 0.0397148348    0.16921121
    ## 30  0.0012566305         38.2204486    3.713173 0.0465110844    0.17505091
    ## 31  0.0011431908         37.0338459    3.778226 0.0542985752    0.18104592
    ## 32  0.0010456011         35.8622047    3.842411 0.0631900451    0.18719843
    ## 33  0.0009615004         34.7064575    3.905645 0.0733056013    0.19351062
    ## 34  0.0008889332         33.5674787    3.967848 0.0847722908    0.19998457
    ## 35  0.0008262756         32.4460842    4.028936 0.0977234839    0.20662235
    ## 36  0.0007721771         31.3430310    4.088829 0.1122980532    0.21342595
    ## 37  0.0007255127         30.2590178    4.147445 0.1286393321    0.22039730
    ## 38  0.0006853452         29.1946850    4.204704 0.1468938418    0.22753825
    ## 39  0.0006508933         28.1506155    4.260529 0.1672097762    0.23485060
    ## 40  0.0006215077         27.1273349    4.314840 0.1897352419    0.24233606
    ## 41  0.0005966495         26.1253126    4.367561 0.2146162532    0.24999627
    ## 42  0.0005758751         25.1449628    4.418618 0.2419944888    0.25783278
    ## 43  0.0005588219         24.1866450    4.467938 0.2720048237    0.26584705
    ## 44  0.0005451986         23.2506660    4.515449 0.3047726558    0.27404045
    ## 45  0.0005347764         22.3372804    4.561082 0.3404110560    0.28241425
    ## 46  0.0005273827         21.4466927    4.604771 0.3790177749    0.29096962
    ## 47  0.0005228965         20.5790582    4.646451 0.4206721508    0.29970764
    ## 48  0.0005212448         19.7344849    4.686060 0.4654319684    0.30862925
    ## 49  0.0005224009         18.9130350    4.723538 0.5133303251    0.31773530
    ## 50  0.0005263835         18.1147269    4.758830 0.5643725675    0.32702651
    ## 51  0.0005332572         17.3395368    4.791882 0.6185333690    0.33650349
    ## 52  0.0005431345         16.5874002    4.822644 0.6757540169    0.34616671
    ## 53  0.0005561786         15.8582145    4.851069 0.7359399878    0.35601650
    ## 54  0.0005726078         15.1518403    4.877114 0.7989588837    0.36605309
    ## 55  0.0005927021         14.4681038    4.900738 0.8646388043    0.37627653
    ## 56  0.0006168105         13.8067983    4.921905 0.9327672261    0.38668676
    ## 57  0.0006453619         13.1676865    4.940582 1.0030904543    0.39728356
    ## 58  0.0006788768         12.5505023    4.956740 1.0753137036    0.40806654
    ## 59  0.0007179841         11.9549530    4.970355 1.1491018577    0.41903519
    ## 60  0.0007634399         11.3807208    4.981404 1.2240809426    0.43018881
    ## 61  0.0008161520         10.8274652    4.989871 1.2998403368    0.44152657
    ## 62  0.0008772097         10.2948248    4.995742 1.3759357247    0.45304745
    ## 63  0.0009479206          9.7824189    4.999008 1.4518927852    0.46475027
    ## 64  0.0010298563          9.2898497    4.999664 1.5272115890    0.47663367
    ## 65  0.0011249093          8.8167038    4.997710 1.6013716598    0.48869613
    ## 66  0.0012353629          8.3625542    4.993147 1.6738376380    0.50093594
    ## 67  0.0013639792          7.9269619    4.985984 1.7440654692    0.51335122
    ## 68  0.0015141089          7.5094777    4.976231 1.8115090225    0.52593990
    ## 69  0.0016898286          7.1096436    4.963904 1.8756270300    0.53869972
    ## 70  0.0018961136          6.7269943    4.949021 1.9358902284    0.55162824
    ## 71  0.0021390564          6.3610589    4.931607 1.9917885711    0.56472284
    ## 72  0.0024261424          6.0113624    4.911688 2.0428383769    0.57798068
    ## 73  0.0027666008          5.6774268    4.889295 2.0885892739    0.59139875
    ## 74  0.0031718517          5.3587723    4.864462 2.1286308018    0.60497384
    ## 75  0.0036560778          5.0549189    4.837229 2.1625985358    0.61870254
    ## 76  0.0042369579          4.7653875    4.807636 2.1901796051    0.63258125
    ## 77  0.0049366125          4.4897006    4.775729 2.2111174901    0.64660617
    ## 78  0.0057828256          4.2273838    4.741556 2.2252159939    0.66077330
    ## 79  0.0068106306          3.9779663    4.705170 2.2323423028    0.67507845
    ## 80  0.0080643751          3.7409824    4.666625 2.2324290699    0.68951721
    ## 81  0.0096004211          3.5159716    4.625979 2.2254754743    0.70408499
    ## 82  0.0114906876          3.3024799    4.583292 2.2115472341    0.71877701
    ## 83  0.0138273169          3.1000603    4.538628 2.1907755726    0.73358828
    ## 84  0.0167288461          2.9082734    4.492053 2.1633551598    0.74851361
    ## 85  0.0203483981          2.7266880    4.443634 2.1295410748    0.76354762
    ## 86  0.0248845982          2.5548818    4.393441 2.0896448557    0.77868475
    ## 87  0.0305961801          2.3924414    4.341547 2.0440297221    0.79391923
    ## 88  0.0378216048          2.2389631    4.288025 1.9931050728    0.80924512
    ## 89  0.0470055205          2.0940531    4.232952 1.9373203747    0.82465627
    ## 90  0.0587345909          1.9573279    4.176404 1.8771585716    0.84014637
    ## 91  0.0737862091          1.8284143    4.118460 1.8131291464    0.85570892
    ## 92  0.0931950010          1.7069498    4.059199 1.7457609762    0.87133724
    ## 93  0.1183439894          1.5925829    3.998702 1.6755951195    0.88702447
    ## 94  0.1510900837          1.4849727    3.937050 1.6031776721    0.90276360
    ## 95  0.1939375420          1.3837895    3.874323 1.5290528203    0.91854744
    ## 96  0.2502787521          1.2887145    3.810606 1.4537562140    0.93436865
    ## 97  0.3247298735          1.1994400    3.745979 1.3778087661    0.95021971
    ## 98  0.4236007054          1.1156693    3.680525 1.3017109751    0.96609297
    ## 99  0.5555552736          1.0371166    3.614327 1.2259378478    0.98198064
    ## 100 0.7325445524          0.9635068    3.547466 1.1509344853    0.99787476
    ##        Pantala
    ## 1   0.03903966
    ## 2   0.04342746
    ## 3   0.04824344
    ## 4   0.05352139
    ## 5   0.05929688
    ## 6   0.06560722
    ## 7   0.07249146
    ## 8   0.07999031
    ## 9   0.08814614
    ## 10  0.09700286
    ## 11  0.10660587
    ## 12  0.11700194
    ## 13  0.12823908
    ## 14  0.14036636
    ## 15  0.15343381
    ## 16  0.16749214
    ## 17  0.18259259
    ## 18  0.19878666
    ## 19  0.21612583
    ## 20  0.23466129
    ## 21  0.25444364
    ## 22  0.27552252
    ## 23  0.29794626
    ## 24  0.32176155
    ## 25  0.34701297
    ## 26  0.37374261
    ## 27  0.40198966
    ## 28  0.43178993
    ## 29  0.46317540
    ## 30  0.49617378
    ## 31  0.53080804
    ## 32  0.56709593
    ## 33  0.60504952
    ## 34  0.64467476
    ## 35  0.68597102
    ## 36  0.72893067
    ## 37  0.77353867
    ## 38  0.81977221
    ## 39  0.86760033
    ## 40  0.91698363
    ## 41  0.96787396
    ## 42  1.02021425
    ## 43  1.07393827
    ## 44  1.12897053
    ## 45  1.18522621
    ## 46  1.24261115
    ## 47  1.30102186
    ## 48  1.36034573
    ## 49  1.42046114
    ## 50  1.48123774
    ## 51  1.54253681
    ## 52  1.60421162
    ## 53  1.66610794
    ## 54  1.72806456
    ## 55  1.78991392
    ## 56  1.85148281
    ## 57  1.91259306
    ## 58  1.97306240
    ## 59  2.03270531
    ## 60  2.09133389
    ## 61  2.14875888
    ## 62  2.20479059
    ## 63  2.25923998
    ## 64  2.31191964
    ## 65  2.36264494
    ## 66  2.41123501
    ## 67  2.45751385
    ## 68  2.50131139
    ## 69  2.54246453
    ## 70  2.58081812
    ## 71  2.61622595
    ## 72  2.64855170
    ## 73  2.67766975
    ## 74  2.70346609
    ## 75  2.72583895
    ## 76  2.74469957
    ## 77  2.75997272
    ## 78  2.77159723
    ## 79  2.77952640
    ## 80  2.78372827
    ## 81  2.78418592
    ## 82  2.78089748
    ## 83  2.77387623
    ## 84  2.76315046
    ## 85  2.74876329
    ## 86  2.73077240
    ## 87  2.70924963
    ## 88  2.68428048
    ## 89  2.65596360
    ## 90  2.62441009
    ## 91  2.58974279
    ## 92  2.55209550
    ## 93  2.51161208
    ## 94  2.46844559
    ## 95  2.42275730
    ## 96  2.37471565
    ## 97  2.32449531
    ## 98  2.27227603
    ## 99  2.21824164
    ## 100 2.16257897

``` r
AB_fechado <- predict.manyglm(mod_int, newdata = new_data_AB_fechado, type = "response")
CD_fechado <- predict.manyglm(mod_int, newdata = new_data_CD_fechado, type = "response")
EF_fechado <- predict.manyglm(mod_int, newdata = new_data_EF_fechado, type = "response")

fechado <- array(NA, dim = c(nrow(AB_fechado), ncol(AB_fechado), 3))
fechado[,,1] <- AB_fechado
fechado[,,2] <- CD_fechado
fechado[,,3] <- EF_fechado

fechado_mean <- as.data.frame(apply(fechado, MARGIN = c(1,2), FUN = mean))
colnames(fechado_mean) <- colnames(AB_fechado)
fechado_mean
```

    ##            Aedes Callibaetis Chironominae     Culex Dendropsophus.minutus
    ## 1   34.039024561   6.2292594     38.84541 59.001946             0.8127767
    ## 2   31.001825403   5.8218683     36.98143 54.454454             0.8313861
    ## 3   28.241690447   5.4461193     35.24730 50.314528             0.8500890
    ## 4   25.732819488   5.0993020     33.63305 46.542138             0.8688725
    ## 5   23.451861523   4.7789571     32.12956 43.101479             0.8877237
    ## 6   21.377678005   4.4828514     30.72851 39.960503             0.9066292
    ## 7   19.491129375   4.2089558     29.42228 37.090496             0.9255750
    ## 8   17.774882544   3.9554254     28.20392 34.465712             0.9445471
    ## 9   16.213237243   3.7205817     27.06704 32.063048             0.9635311
    ## 10  14.791969350   3.5028965     26.00580 29.861751             0.9825121
    ## 11  13.498189520   3.3009776     25.01484 27.843169             1.0014751
    ## 12  12.320215584   3.1135558     24.08927 25.990521             1.0204049
    ## 13  11.247457370   2.9394735     23.22457 24.288697             1.0392857
    ## 14  10.270312702   2.7776739     22.41660 22.724084             1.0581018
    ## 15   9.380073479   2.6271918     21.66158 21.284403             1.0768372
    ## 16   8.568840850   2.4871450     20.95602 19.958573             1.0954757
    ## 17   7.829448568   2.3567268     20.29670 18.736585             1.1140008
    ## 18   7.155393745   2.2351990     19.68070 17.609390             1.1323961
    ## 19   6.540774265   2.1218856     19.10529 16.568801             1.1506448
    ## 20   5.980232209   2.0161672     18.56799 15.607409             1.1687303
    ## 21   5.468902703   1.9174760     18.06651 14.718496             1.1866356
    ## 22   5.002367669   1.8252911     17.59876 13.895974             1.2043440
    ## 23   4.576613984   1.7391345     17.16279 13.134317             1.2218385
    ## 24   4.187995635   1.6585669     16.75683 12.428505             1.2391023
    ## 25   3.833199475   1.5831849     16.37926 11.773979             1.2561184
    ## 26   3.509214234   1.5126174     16.02856 11.166588             1.2728701
    ## 27   3.213302464   1.4465230     15.70338 10.602559             1.2893407
    ## 28   2.942975147   1.3845875     15.40246 10.078451             1.3055135
    ## 29   2.695968693   1.3265215     15.12464  9.591131             1.3213721
    ## 30   2.470224116   1.2720583     14.86888  9.137740             1.3369001
    ## 31   2.263868166   1.2209518     14.63422  8.715668             1.3520816
    ## 32   2.075196233   1.1729752     14.41979  8.322532             1.3669005
    ## 33   1.902656858   1.1279192     14.22482  7.956155             1.3813413
    ## 34   1.744837691   1.0855902     14.04859  7.614544             1.3953886
    ## 35   1.600452760   1.0458098     13.89046  7.295876             1.4090273
    ## 36   1.468330933   1.0084126     13.74988  6.998484             1.4222428
    ## 37   1.347405447   0.9732460     13.62635  6.720838             1.4350207
    ## 38   1.236704418   0.9401688     13.51942  6.461536             1.4473470
    ## 39   1.135342220   0.9090502     13.42873  6.219293             1.4592082
    ## 40   1.042511675   0.8797690     13.35395  5.992931             1.4705911
    ## 41   0.957476952   0.8522133     13.29484  5.781365             1.4814831
    ## 42   0.879567132   0.8262791     13.25118  5.583602             1.4918719
    ## 43   0.808170361   0.8018701     13.22282  5.398727             1.5017459
    ## 44   0.742728532   0.7788971     13.20966  5.225903             1.5110939
    ## 45   0.682732473   0.7572773     13.21167  5.064355             1.5199053
    ## 46   0.627717554   0.7369341     13.22884  4.913375             1.5281701
    ## 47   0.577259715   0.7177962     13.26124  4.772309             1.5358788
    ## 48   0.530971842   0.6997977     13.30898  4.640558             1.5430226
    ## 49   0.488500477   0.6828772     13.37221  4.517568             1.5495932
    ## 50   0.449522834   0.6669781     13.45117  4.402832             1.5555831
    ## 51   0.413744072   0.6520476     13.54613  4.295884             1.5609852
    ## 52   0.380894827   0.6380370     13.65741  4.196293             1.5657934
    ## 53   0.350728959   0.6249011     13.78542  4.103666             1.5700020
    ## 54   0.323021509   0.6125978     13.93059  4.017642             1.5736061
    ## 55   0.297566831   0.6010885     14.09345  3.937887             1.5766015
    ## 56   0.274176899   0.5903374     14.27458  3.864099             1.5789848
    ## 57   0.252679764   0.5803111     14.47463  3.796000             1.5807530
    ## 58   0.232918148   0.5709793     14.69433  3.733336             1.5819042
    ## 59   0.214748163   0.5623136     14.93449  3.675876             1.5824370
    ## 60   0.198038144   0.5542882     15.19599  3.623410             1.5823508
    ## 61   0.182667591   0.5468793     15.47982  3.575750             1.5816456
    ## 62   0.168526195   0.5400652     15.78705  3.532724             1.5803223
    ## 63   0.155512963   0.5338260     16.11886  3.494179             1.5783825
    ## 64   0.143535404   0.5281436     16.47653  3.459980             1.5758284
    ## 65   0.132508806   0.5230018     16.86147  3.430006             1.5726630
    ## 66   0.122355559   0.5183858     17.27521  3.404153             1.5688900
    ## 67   0.113004551   0.5142826     17.71941  3.382332             1.5645138
    ## 68   0.104390608   0.5106806     18.19590  3.364468             1.5595396
    ## 69   0.096453986   0.5075698     18.70664  3.350498             1.5539730
    ## 70   0.089139911   0.5049414     19.25380  3.340376             1.5478206
    ## 71   0.082398153   0.5027880     19.83970  3.334066             1.5410895
    ## 72   0.076182640   0.5011039     20.46690  3.331547             1.5337875
    ## 73   0.070451107   0.4998842     21.13816  3.332811             1.5259230
    ## 74   0.065164774   0.4991256     21.85650  3.337862             1.5175049
    ## 75   0.060288048   0.4988260     22.62518  3.346716             1.5085429
    ## 76   0.055788260   0.4989846     23.44778  3.359405             1.4990473
    ## 77   0.051635414   0.4996018     24.32818  3.375972             1.4890287
    ## 78   0.047801969   0.5006793     25.27061  3.396473             1.4784985
    ## 79   0.044262625   0.5022202     26.27967  3.420979             1.4674685
    ## 80   0.040994143   0.5042286     27.36039  3.449575             1.4559511
    ## 81   0.037975169   0.5067101     28.51826  3.482360             1.4439590
    ## 82   0.035186080   0.5096717     29.75924  3.519449             1.4315055
    ## 83   0.032608837   0.5131215     31.08986  3.560973             1.4186042
    ## 84   0.030226858   0.5170694     32.51726  3.607078             1.4052694
    ## 85   0.028024893   0.5215262     34.04923  3.657930             1.3915153
    ## 86   0.025988917   0.5265048     35.69429  3.713711             1.3773568
    ## 87   0.024106028   0.5320192     37.46179  3.774625             1.3628091
    ## 88   0.022364356   0.5380853     39.36193  3.840894             1.3478876
    ## 89   0.020752977   0.5447205     41.40593  3.912766             1.3326080
    ## 90   0.019261836   0.5519442     43.60606  3.990509             1.3169863
    ## 91   0.017881676   0.5597775     45.97580  4.074419             1.3010385
    ## 92   0.016603972   0.5682435     48.52996  4.164817             1.2847811
    ## 93   0.015420877   0.5773675     51.28481  4.262056             1.2682305
    ## 94   0.014325157   0.5871770     54.25825  4.366518             1.2514035
    ## 95   0.013310150   0.5977017     57.46997  4.478621             1.2343166
    ## 96   0.012369718   0.6089741     60.94167  4.598819             1.2169869
    ## 97   0.011498201   0.6210290     64.69727  4.727605             1.1994311
    ## 98   0.010690383   0.6339045     68.76314  4.865518             1.1816661
    ## 99   0.009941454   0.6476413     73.16842  5.013140             1.1637089
    ## 100  0.009246978   0.6622838     77.94527  5.171106             1.1455762
    ##       Microvelia Scinax.fuscovarius Tanypodinae    Dasyhelea Erythrodiplax
    ## 1   1.241708e-03         0.82591293    5.118157 2.840211e-08  1.818516e-03
    ## 2   1.054813e-03         0.86380856    5.272859 4.196292e-08  2.450084e-03
    ## 3   9.007502e-04         0.90179259    5.428258 6.175895e-08  3.272164e-03
    ## 4   7.732258e-04         0.93972708    5.584143 9.054263e-08  4.331907e-03
    ## 5   6.672388e-04         0.97746841    5.740297 1.322286e-07  5.684774e-03
    ## 6   5.788009e-04         1.01486819    5.896496 1.923608e-07  7.394985e-03
    ## 7   5.047196e-04         1.05177409    6.052509 2.787576e-07  9.535674e-03
    ## 8   4.424296e-04         1.08803086    6.208100 4.023980e-07  1.218865e-02
    ## 9   3.898622e-04         1.12348137    6.363026 5.786341e-07  1.544364e-02
    ## 10  3.453433e-04         1.15796772    6.517041 8.288409e-07  1.939697e-02
    ## 11  3.075134e-04         1.19133238    6.669896 1.182653e-06  2.414950e-02
    ## 12  2.752644e-04         1.22341939    6.821335 1.680980e-06  2.980385e-02
    ## 13  2.476903e-04         1.25407551    6.971103 2.380053e-06  3.646084e-02
    ## 14  2.240480e-04         1.28315149    7.118942 3.356833e-06  4.421514e-02
    ## 15  2.037258e-04         1.31050323    7.264590 4.716197e-06  5.315025e-02
    ## 16  1.862190e-04         1.33599298    7.407789 6.600443e-06  6.333295e-02
    ## 17  1.711099e-04         1.35949049    7.548278 9.201809e-06  7.480731e-02
    ## 18  1.580516e-04         1.38087411    7.685797 1.277887e-05  8.758878e-02
    ## 19  1.467561e-04         1.40003186    7.820090 1.767789e-05  1.016583e-01
    ## 20  1.369828e-04         1.41686237    7.950902 2.436058e-05  1.169573e-01
    ## 21  1.285313e-04         1.43127581    8.077980 3.343981e-05  1.333834e-01
    ## 22  1.212342e-04         1.44319466    8.201078 4.572554e-05  1.507878e-01
    ## 23  1.149514e-04         1.45255442    8.319954 6.228350e-05  1.689744e-01
    ## 24  1.095661e-04         1.45930419    8.434371 8.450961e-05  1.877005e-01
    ## 25  1.049811e-04         1.46340711    8.544098 1.142242e-04  2.066807e-01
    ## 26  1.011158e-04         1.46484074    8.648913 1.537904e-04  2.255925e-01
    ## 27  9.790396e-05         1.46359723    8.748602 2.062621e-04  2.440840e-01
    ## 28  9.529153e-05         1.45968338    8.842957 2.755678e-04  2.617845e-01
    ## 29  9.323551e-05         1.45312061    8.931784 3.667386e-04  2.783163e-01
    ## 30  9.170255e-05         1.44394477    9.014894 4.861874e-04  2.933076e-01
    ## 31  9.066809e-05         1.43220576    9.092114 6.420514e-04  3.064066e-01
    ## 32  9.011572e-05         1.41796715    9.163279 8.446074e-04  3.172947e-01
    ## 33  9.003671e-05         1.40130554    9.228236 1.106774e-03  3.256999e-01
    ## 34  9.042983e-05         1.38230992    9.286847 1.444714e-03  3.314076e-01
    ## 35  9.130126e-05         1.36108087    9.338984 1.878556e-03  3.342699e-01
    ## 36  9.266481e-05         1.33772964    9.384536 2.433241e-03  3.342121e-01
    ## 37  9.454225e-05         1.31237723    9.423402 3.139534e-03  3.312356e-01
    ## 38  9.696389e-05         1.28515332    9.455499 4.035193e-03  3.254183e-01
    ## 39  9.996941e-05         1.25619515    9.480756 5.166333e-03  3.169107e-01
    ## 40  1.036089e-04         1.22564642    9.499117 6.588999e-03  3.059298e-01
    ## 41  1.079445e-04         1.19365606    9.510543 8.370964e-03  2.927500e-01
    ## 42  1.130515e-04         1.16037705    9.515008 1.059377e-02  2.776910e-01
    ## 43  1.190215e-04         1.12596520    9.512503 1.335502e-02  2.611060e-01
    ## 44  1.259644e-04         1.09057798    9.503033 1.677095e-02  2.433672e-01
    ## 45  1.340117e-04         1.05437330    9.486618 2.097923e-02  2.248522e-01
    ## 46  1.433214e-04         1.01750837    9.463296 2.614211e-02  2.059312e-01
    ## 47  1.540820e-04         0.98013861    9.433116 3.244970e-02  1.869551e-01
    ## 48  1.665199e-04         0.94241659    9.396146 4.012358e-02  1.682451e-01
    ## 49  1.809061e-04         0.90449103    9.352465 4.942057e-02  1.500851e-01
    ## 50  1.975666e-04         0.86650590    9.302168 6.063659e-02  1.327159e-01
    ## 51  2.168935e-04         0.82859956    9.245366 7.411067e-02  1.163317e-01
    ## 52  2.393606e-04         0.79090403    9.182180 9.022892e-02  1.010796e-01
    ## 53  2.655412e-04         0.75354432    9.112746 1.094283e-01  8.706002e-02
    ## 54  2.961311e-04         0.71663783    9.037213 1.322004e-01  7.432999e-02
    ## 55  3.319779e-04         0.68029389    8.955741 1.590944e-01  6.290706e-02
    ## 56  3.741169e-04         0.64461340    8.868504 1.907198e-01  5.277458e-02
    ## 57  4.238172e-04         0.60968851    8.775684 2.277486e-01  4.388743e-02
    ## 58  4.826393e-04         0.57560241    8.677475 2.709161e-01  3.617809e-02
    ## 59  5.525097e-04         0.54242927    8.574080 3.210205e-01  2.956249e-02
    ## 60  6.358140e-04         0.51023418    8.465712 3.789219e-01  2.394565e-02
    ## 61  7.355179e-04         0.47907322    8.352591 4.455390e-01  1.922658e-02
    ## 62  8.553216e-04         0.44899360    8.234946 5.218440e-01  1.530268e-02
    ## 63  9.998588e-04         0.42003389    8.113011 6.088561e-01  1.207322e-02
    ## 64  1.174954e-03         0.39222424    7.987027 7.076322e-01  9.442103e-03
    ## 65  1.387957e-03         0.36558675    7.857240 8.192559e-01  7.319886e-03
    ## 66  1.648179e-03         0.34013582    7.723900 9.448231e-01  5.625096e-03
    ## 67  1.967458e-03         0.31587861    7.587262 1.085427e+00  4.284948e-03
    ## 68  2.360911e-03         0.29281545    7.447583 1.242137e+00  3.235574e-03
    ## 69  2.847914e-03         0.27094034    7.305120 1.415981e+00  2.421849e-03
    ## 70  3.453402e-03         0.25024146    7.160134 1.607920e+00  1.796936e-03
    ## 71  4.209596e-03         0.23070170    7.012885 1.818822e+00  1.321626e-03
    ## 72  5.158300e-03         0.21229913    6.863633 2.049440e+00  9.635500e-04
    ## 73  6.353980e-03         0.19500761    6.712637 2.300378e+00  6.963540e-04
    ## 74  7.867887e-03         0.17879725    6.560155 2.572066e+00  4.988568e-04
    ## 75  9.793622e-03         0.16363493    6.406441 2.864732e+00  3.542516e-04
    ## 76  1.225467e-02         0.14948483    6.251746 3.178373e+00  2.493662e-04
    ## 77  1.541462e-02         0.13630888    6.096318 3.512730e+00  1.740017e-04
    ## 78  1.949113e-02         0.12406723    5.940399 3.867263e+00  1.203537e-04
    ## 79  2.477504e-02         0.11271869    5.784229 4.241130e+00  8.251929e-05
    ## 80  3.165663e-02         0.10222114    5.628039 4.633173e+00  5.608432e-05
    ## 81  4.066192e-02         0.09253188    5.472056 5.041902e+00  3.778483e-05
    ## 82  5.250301e-02         0.08360803    5.316499 5.465493e+00  2.523385e-05
    ## 83  6.814805e-02         0.07540680    5.161580 5.901783e+00  1.670473e-05
    ## 84  8.891924e-02         0.06788580    5.007506 6.348281e+00  1.096190e-05
    ## 85  1.166302e-01         0.06100330    4.854472 6.802178e+00  7.130529e-06
    ## 86  1.537798e-01         0.05471842    4.702668 7.260372e+00  4.597777e-06
    ## 87  2.038265e-01         0.04899139    4.552274 7.719493e+00  2.938760e-06
    ## 88  2.715782e-01         0.04378364    4.403463 8.175939e+00  1.861960e-06
    ## 89  3.637494e-01         0.03905799    4.256396 8.625922e+00  1.169409e-06
    ## 90  4.897591e-01         0.03477874    4.111227 9.065513e+00  7.280358e-07
    ## 91  6.628814e-01         0.03091176    3.968101 9.490700e+00  4.492923e-07
    ## 92  9.019079e-01         0.02742455    3.827152 9.897446e+00  2.748497e-07
    ## 93  1.233563e+00         0.02428629    3.688506 1.028175e+01  1.666678e-07
    ## 94  1.696031e+00         0.02146786    3.552279 1.063971e+01  1.001839e-07
    ## 95  2.344116e+00         0.01894185    3.418578 1.096760e+01  5.969450e-08
    ## 96  3.256848e+00         0.01668252    3.287499 1.126192e+01  3.525825e-08
    ## 97  4.548716e+00         0.01466584    3.159131 1.151947e+01  2.064321e-08
    ## 98  6.386357e+00         0.01286940    3.033551 1.173738e+01  1.198074e-08
    ## 99  9.013439e+00         0.01127238    2.910830 1.191322e+01  6.892552e-09
    ## 100 1.278795e+01         0.00985550    2.791028 1.204497e+01  3.930670e-09
    ##        Pantala
    ## 1   0.08966493
    ## 2   0.09774335
    ## 3   0.10640092
    ## 4   0.11566371
    ## 5   0.12555744
    ## 6   0.13610728
    ## 7   0.14733769
    ## 8   0.15927219
    ## 9   0.17193314
    ## 10  0.18534158
    ## 11  0.19951690
    ## 12  0.21447668
    ## 13  0.23023645
    ## 14  0.24680938
    ## 15  0.26420609
    ## 16  0.28243438
    ## 17  0.30149901
    ## 18  0.32140142
    ## 19  0.34213956
    ## 20  0.36370759
    ## 21  0.38609575
    ## 22  0.40929012
    ## 23  0.43327246
    ## 24  0.45802006
    ## 25  0.48350557
    ## 26  0.50969698
    ## 27  0.53655743
    ## 28  0.56404526
    ## 29  0.59211393
    ## 30  0.62071206
    ## 31  0.64978349
    ## 32  0.67926735
    ## 33  0.70909822
    ## 34  0.73920626
    ## 35  0.76951742
    ## 36  0.79995371
    ## 37  0.83043347
    ## 38  0.86087167
    ## 39  0.89118028
    ## 40  0.92126867
    ## 41  0.95104403
    ## 42  0.98041179
    ## 43  1.00927615
    ## 44  1.03754055
    ## 45  1.06510820
    ## 46  1.09188264
    ## 47  1.11776828
    ## 48  1.14267093
    ## 49  1.16649844
    ## 50  1.18916120
    ## 51  1.21057271
    ## 52  1.23065016
    ## 53  1.24931492
    ## 54  1.26649309
    ## 55  1.28211596
    ## 56  1.29612047
    ## 57  1.30844966
    ## 58  1.31905301
    ## 59  1.32788683
    ## 60  1.33491453
    ## 61  1.34010689
    ## 62  1.34344225
    ## 63  1.34490667
    ## 64  1.34449402
    ## 65  1.34220604
    ## 66  1.33805228
    ## 67  1.33205010
    ## 68  1.32422451
    ## 69  1.31460798
    ## 70  1.30324027
    ## 71  1.29016811
    ## 72  1.27544488
    ## 73  1.25913029
    ## 74  1.24128993
    ## 75  1.22199484
    ## 76  1.20132108
    ## 77  1.17934917
    ## 78  1.15616362
    ## 79  1.13185234
    ## 80  1.10650614
    ## 81  1.08021815
    ## 82  1.05308322
    ## 83  1.02519741
    ## 84  0.99665739
    ## 85  0.96755990
    ## 86  0.93800125
    ## 87  0.90807675
    ## 88  0.87788025
    ## 89  0.84750366
    ## 90  0.81703652
    ## 91  0.78656559
    ## 92  0.75617444
    ## 93  0.72594318
    ## 94  0.69594810
    ## 95  0.66626140
    ## 96  0.63695103
    ## 97  0.60808041
    ## 98  0.57970836
    ## 99  0.55188895
    ## 100 0.52467143

``` r
AB_fechado_dia <- predict.manyglm(mod_int, newdata = new_data_AB_fechado_dia, type = "response")
CD_fechado_dia <- predict.manyglm(mod_int, newdata = new_data_CD_fechado_dia, type = "response")
EF_fechado_dia <- predict.manyglm(mod_int, newdata = new_data_EF_fechado_dia, type = "response")

fechado_dia <- array(NA, dim = c(nrow(AB_fechado_dia), ncol(AB_fechado_dia), 3))
fechado_dia[,,1] <- AB_fechado_dia
fechado_dia[,,2] <- CD_fechado_dia
fechado_dia[,,3] <- EF_fechado_dia

fechado_dia_mean <- as.data.frame(apply(fechado_dia, MARGIN = c(1,2), FUN = mean))
colnames(fechado_dia_mean) <- colnames(AB_fechado_dia)
fechado_dia_mean
```

    ##          Aedes Callibaetis Chironominae     Culex Dendropsophus.minutus
    ## 1   96.9410183   8.5809797     28.17368 97.831014              2.972319
    ## 2   92.3357069   7.7959967     27.54490 93.899762              3.184734
    ## 3   87.9326308   7.0975407     26.94314 90.140043              3.406852
    ## 4   83.7237635   6.4750872     26.36723 86.543883              3.638612
    ## 5   79.7013550   5.9194973     25.81607 83.103694              3.879900
    ## 6   75.8579243   5.4228241     25.28861 79.812262              4.130547
    ## 7   72.1862518   4.9781467     24.78386 76.662724              4.390328
    ## 8   68.6793722   4.5794291     24.30090 73.648553              4.658958
    ## 9   65.3305671   4.2213996     23.83883 70.763536              4.936087
    ## 10  62.1333582   3.8994475     23.39682 68.001763              5.221307
    ## 11  59.0815000   3.6095342     22.97407 65.357610              5.514142
    ## 12  56.1689734   3.3481178     22.56983 62.825722              5.814053
    ## 13  53.3899786   3.1120873     22.18339 60.401004              6.120436
    ## 14  50.7389288   2.8987069     21.81408 58.078603              6.432623
    ## 15  48.2104437   2.7055671     21.46125 55.853900              6.749881
    ## 16  45.7993433   2.5305435     21.12431 53.722497              7.071417
    ## 17  43.5006413   2.3717602     20.80268 51.680203              7.396379
    ## 18  41.3095397   2.2275590     20.49582 49.723029              7.723856
    ## 19  39.2214224   2.0964724     20.20321 47.847172              8.052885
    ## 20  37.2318495   1.9771998     19.92439 46.049012              8.382454
    ## 21  35.3365518   1.8685875     19.65888 44.325097              8.711504
    ## 22  33.5314253   1.7696110     19.40626 42.672139              9.038939
    ## 23  31.8125254   1.6793595     19.16612 41.087004              9.363627
    ## 24  30.1760622   1.5970224     18.93807 39.566704              9.684408
    ## 25  28.6183949   1.5218779     18.72175 38.108391             10.000101
    ## 26  27.1360272   1.4532827     18.51683 36.709349             10.309510
    ## 27  25.7256021   1.3906629     18.32297 35.366990             10.611431
    ## 28  24.3838971   1.3335064     18.13989 34.078844             10.904663
    ## 29  23.1078198   1.2813561     17.96729 32.842556             11.188010
    ## 30  21.8944033   1.2338036     17.80490 31.655880             11.460294
    ## 31  20.7408019   1.1904845     17.65249 30.516671             11.720362
    ## 32  19.6442866   1.1510731     17.50982 29.422886             11.967091
    ## 33  18.6022410   1.1152791     17.37668 28.372573             12.199401
    ## 34  17.6121574   1.0828435     17.25285 27.363870             12.416259
    ## 35  16.6716328   1.0535358     17.13817 26.394999             12.616686
    ## 36  15.7783652   1.0271512     17.03245 25.464263             12.799771
    ## 37  14.9301495   1.0035083     16.93554 24.570043             12.964668
    ## 38  14.1248744   0.9824467     16.84730 23.710793             13.110612
    ## 39  13.3605188   0.9638258     16.76760 22.885034             13.236916
    ## 40  12.6351482   0.9475225     16.69632 22.091357             13.342986
    ## 41  11.9469113   0.9334306     16.63335 21.328414             13.428316
    ## 42  11.2940376   0.9214590     16.57861 20.594918             13.492499
    ## 43  10.6748334   0.9115310     16.53201 19.889639             13.535227
    ## 44  10.0876794   0.9035837     16.49349 19.211403             13.556296
    ## 45   9.5310275   0.8975668     16.46299 18.559087             13.555604
    ## 46   9.0033982   0.8934427     16.44046 17.931618             13.533154
    ## 47   8.5033780   0.8911854     16.42588 17.327969             13.489055
    ## 48   8.0296164   0.8907809     16.41923 16.747162             13.423517
    ## 49   7.5808237   0.8922267     16.42048 16.188257             13.336856
    ## 50   7.1557685   0.8955319     16.42965 15.650359             13.229484
    ## 51   6.7532753   0.9007169     16.44675 15.132611             13.101913
    ## 52   6.3722223   0.9078145     16.47180 14.634193             12.954743
    ## 53   6.0115390   0.9168691     16.50484 14.154320             12.788666
    ## 54   5.6702042   0.9279383     16.54591 13.692243             12.604453
    ## 55   5.3472443   0.9410925     16.59508 13.247243             12.402953
    ## 56   5.0417306   0.9564163     16.65243 12.818634             12.185084
    ## 57   4.7527780   0.9740094     16.71802 12.405759             11.951826
    ## 58   4.4795431   0.9939872     16.79196 12.007989             11.704217
    ## 59   4.2212220   1.0164826     16.87436 11.624722             11.443339
    ## 60   3.9770490   1.0416469     16.96533 11.255380             11.170316
    ## 61   3.7462952   1.0696523     17.06501 10.899413             10.886306
    ## 62   3.5282661   1.1006930     17.17356 10.556292             10.592486
    ## 63   3.3223008   1.1349880     17.29112 10.225512             10.290053
    ## 64   3.1277703   1.1727833     17.41788  9.906586              9.980209
    ## 65   2.9440762   1.2143554     17.55402  9.599051              9.664157
    ## 66   2.7706491   1.2600138     17.69976  9.302463              9.343093
    ## 67   2.6069476   1.3101055     17.85531  9.016396              9.018196
    ## 68   2.4524567   1.3650190     18.02090  8.740440              8.690625
    ## 69   2.3066872   1.4251896     18.19680  8.474205              8.361509
    ## 70   2.1691737   1.4911044     18.38327  8.217315              8.031944
    ## 71   2.0394744   1.5633094     18.58060  7.969412              7.702984
    ## 72   1.9171693   1.6424165     18.78910  7.730151              7.375639
    ## 73   1.8018596   1.7291121     19.00910  7.499201              7.050868
    ## 74   1.6931667   1.8241666     19.24094  7.276246              6.729579
    ## 75   1.5907312   1.9284452     19.48500  7.060981              6.412620
    ## 76   1.4942118   2.0429211     19.74166  6.853116              6.100781
    ## 77   1.4032847   2.1686894     20.01135  6.652371              5.794791
    ## 78   1.3176429   2.3069842     20.29449  6.458478              5.495312
    ## 79   1.2369950   2.4591971     20.59156  6.271180              5.202945
    ## 80   1.1610647   2.6269001     20.90305  6.090230              4.918227
    ## 81   1.0895903   2.8118701     21.22948  5.915390              4.641626
    ## 82   1.0223234   3.0161186     21.57139  5.746435              4.373549
    ## 83   0.9590288   3.2419258     21.92938  5.583145              4.114341
    ## 84   0.8994837   3.4918791     22.30405  5.425311              3.864282
    ## 85   0.8434770   3.7689191     22.69606  5.272733              3.623595
    ## 86   0.7908088   4.0763916     23.10609  5.125216              3.392446
    ## 87   0.7412898   4.4181095     23.53487  4.982576              3.170943
    ## 88   0.6947408   4.7984229     23.98315  4.844635              2.959146
    ## 89   0.6509924   5.2223029     24.45175  4.711221              2.757062
    ## 90   0.6098840   5.6954372     24.94153  4.582171              2.564656
    ## 91   0.5712641   6.2243436     25.45337  4.457326              2.381847
    ## 92   0.5349890   6.8165014     25.98824  4.336535              2.208519
    ## 93   0.5009232   7.4805061     26.54714  4.219653              2.044516
    ## 94   0.4689382   8.2262502     27.13112  4.106538              1.889654
    ## 95   0.4389130   9.0651360     27.74132  3.997057              1.743719
    ## 96   0.4107329  10.0103259     28.37891  3.891080              1.606471
    ## 97   0.3842899  11.0770365     29.04515  3.788483              1.477651
    ## 98   0.3594815  12.2828865     29.74136  3.689146              1.356978
    ## 99   0.3362115  13.6483066     30.46893  3.592954              1.244161
    ## 100  0.3143886  15.1970252     31.22934  3.499797              1.138891
    ##       Microvelia Scinax.fuscovarius Tanypodinae    Dasyhelea Erythrodiplax
    ## 1   0.6542149831         86.1691260   0.5400349 5.953076e-07    0.04732806
    ## 2   0.4914014160         84.7546542   0.6059726 7.899954e-07    0.05164622
    ## 3   0.3712833100         83.2809161   0.6786456 1.045820e-06    0.05626843
    ## 4   0.2821808628         81.7518335   0.7585636 1.381142e-06    0.06120648
    ## 5   0.2157261502         80.1714206   0.8462521 1.819570e-06    0.06647165
    ## 6   0.1658941745         78.5437670   0.9422506 2.391378e-06    0.07207454
    ## 7   0.1283253875         76.8730201   1.0471092 3.135285e-06    0.07802498
    ## 8   0.0998497862         75.1633677   1.1613854 4.100671e-06    0.08433190
    ## 9   0.0781510470         73.4190206   1.2856407 5.350346e-06    0.09100317
    ## 10  0.0615283924         71.6441960   1.4204362 6.963987e-06    0.09804547
    ## 11  0.0487269740         69.8431004   1.5663280 9.042388e-06    0.10546416
    ## 12  0.0388165044         68.0199136   1.7238622 1.171271e-05    0.11326317
    ## 13  0.0311040190         66.1787731   1.8935696 1.513495e-05    0.12144479
    ## 14  0.0250708871         64.3237592   2.0759593 1.950983e-05    0.13000962
    ## 15  0.0203271270         62.4588801   2.2715133 2.508853e-05    0.13895637
    ## 16  0.0165781248         60.5880588   2.4806791 3.218444e-05    0.14828179
    ## 17  0.0136002821         58.7151201   2.7038635 4.118754e-05    0.15798054
    ## 18  0.0112231182         56.8437785   2.9414251 5.258173e-05    0.16804507
    ## 19  0.0093160596         54.9776273   3.1936675 6.696580e-05    0.17846552
    ## 20  0.0077786475         53.1201286   3.4608315 8.507860e-05    0.18922970
    ## 21  0.0065332462         51.2746037   3.7430884 1.078293e-04    0.20032292
    ## 22  0.0055195932         49.4442252   4.0405323 1.363334e-04    0.21172805
    ## 23  0.0046907061         47.6320102   4.3531731 1.719558e-04    0.22342540
    ## 24  0.0040097978         45.8408133   4.6809301 2.163618e-04    0.23539275
    ## 25  0.0034479411         44.0733221   5.0236254 2.715773e-04    0.24760535
    ## 26  0.0029822930         42.3320530   5.3809778 3.400599e-04    0.26003592
    ## 27  0.0025947402         40.6193476   5.7525979 4.247824e-04    0.27265476
    ## 28  0.0022708610         38.9373712   6.1379832 5.293303e-04    0.28542973
    ## 29  0.0019991268         37.2881108   6.5365147 6.580155e-04    0.29832642
    ## 30  0.0017702852         35.6733753   6.9474535 8.160083e-04    0.31130825
    ## 31  0.0015768822         34.0947959   7.3699396 1.009490e-03    0.32433659
    ## 32  0.0014128900         32.5538275   7.8029903 1.245830e-03    0.33737092
    ## 33  0.0012734167         31.0517506   8.2455014 1.533786e-03    0.35036906
    ## 34  0.0011544786         29.5896747   8.6962484 1.883734e-03    0.36328732
    ## 35  0.0010528204         28.1685416   9.1538896 2.307935e-03    0.37608078
    ## 36  0.0009657747         26.7891297   9.6169702 2.820829e-03    0.38870348
    ## 37  0.0008911493         25.4520587  10.0839279 3.439372e-03    0.40110873
    ## 38  0.0008271385         24.1577953  10.5531002 4.183411e-03    0.41324936
    ## 39  0.0007722520         22.9066588  11.0227321 5.076110e-03    0.42507803
    ## 40  0.0007252588         21.6988276  11.4909863 6.144417e-03    0.43654753
    ## 41  0.0006851412         20.5343452  11.9559538 7.419582e-03    0.44761105
    ## 42  0.0006510589         19.4131281  12.4156658 8.937732e-03    0.45822255
    ## 43  0.0006223197         18.3349719  12.8681071 1.074050e-02    0.46833705
    ## 44  0.0005983563         17.2995597  13.3112298 1.287569e-02    0.47791094
    ## 45  0.0005787079         16.3064686  13.7429688 1.539805e-02    0.48690229
    ## 46  0.0005630047         15.3551781  14.1612570 1.837004e-02    0.49527119
    ## 47  0.0005509570         14.4450772  14.5640417 2.186269e-02    0.50297998
    ## 48  0.0005423461         13.5754723  14.9493010 2.595650e-02    0.50999361
    ## 49  0.0005370175         12.7455945  15.3150608 3.074241e-02    0.51627984
    ## 50  0.0005348764         11.9546073  15.6594111 3.632275e-02    0.52180952
    ## 51  0.0005358849         11.2016140  15.9805232 4.281231e-02    0.52655679
    ## 52  0.0005400610         10.4856646  16.2766650 5.033936e-02    0.53049932
    ## 53  0.0005474786          9.8057629  16.5462173 5.904674e-02    0.53361846
    ## 54  0.0005582703          9.1608735  16.7876878 6.909287e-02    0.53589938
    ## 55  0.0005726313          8.5499280  16.9997255 8.065284e-02    0.53733120
    ## 56  0.0005908249          7.9718312  17.1811330 9.391938e-02    0.53790708
    ## 57  0.0006131907          7.4254675  17.3308776 1.091038e-01    0.53762426
    ## 58  0.0006401554          6.9097061  17.4481013 1.264368e-01    0.53648409
    ## 59  0.0006722462          6.4234068  17.5321290 1.461694e-01    0.53449203
    ## 60  0.0007101081          5.9654243  17.5824750 1.685732e-01    0.53165759
    ## 61  0.0007545250          5.5346139  17.5988472 1.939410e-01    0.52799425
    ## 62  0.0008064472          5.1298349  17.5811508 2.225871e-01    0.52351939
    ## 63  0.0008670244          4.7499552  17.5294884 2.548469e-01    0.51825409
    ## 64  0.0009376480          4.3938549  17.4441595 2.910770e-01    0.51222303
    ## 65  0.0010200030          4.0604295  17.3256573 3.316542e-01    0.50545426
    ## 66  0.0011161336          3.7485932  17.1746643 3.769747e-01    0.49797899
    ## 67  0.0012285250          3.4572813  16.9920454 4.274527e-01    0.48983134
    ## 68  0.0013602068          3.1854529  16.7788397 4.835184e-01    0.48104811
    ## 69  0.0015148826          2.9320930  16.5362510 5.456160e-01    0.47166849
    ## 70  0.0016970949          2.6962139  16.2656360 6.142007e-01    0.46173373
    ## 71  0.0019124336          2.4768576  15.9684922 6.897356e-01    0.45128691
    ## 72  0.0021678026          2.2730961  15.6464434 7.726879e-01    0.44037257
    ## 73  0.0024717594          2.0840331  15.3012258 8.635245e-01    0.42903645
    ## 74  0.0028349523          1.9088048  14.9346716 9.627075e-01    0.41732510
    ## 75  0.0032706829          1.7465800  14.5486935 1.070689e+00    0.40528565
    ## 76  0.0037956330          1.5965610  14.1452679 1.187903e+00    0.39296542
    ## 77  0.0044308098          1.4579835  13.7264182 1.314765e+00    0.38041168
    ## 78  0.0052027757          1.3301168  13.2941978 1.451658e+00    0.36767131
    ## 79  0.0061452591          1.2122635  12.8506742 1.598930e+00    0.35479055
    ## 80  0.0073012700          1.1037592  12.3979122 1.756887e+00    0.34181470
    ## 81  0.0087258898          1.0039723  11.9379584 1.925783e+00    0.32878790
    ## 82  0.0104899675          0.9123032  11.4728266 2.105814e+00    0.31575288
    ## 83  0.0126850350          0.8281837  11.0044836 2.297110e+00    0.30275073
    ## 84  0.0154298712          0.7510767  10.5348357 2.499728e+00    0.28982075
    ## 85  0.0188793066          0.6804747  10.0657173 2.713643e+00    0.27700024
    ## 86  0.0232360811          0.6158993   9.5988799 2.938744e+00    0.26432438
    ## 87  0.0287668849          0.5569003   9.1359824 3.174827e+00    0.25182606
    ## 88  0.0358241499          0.5030549   8.6785828 3.421585e+00    0.23953586
    ## 89  0.0448757848          0.4539660   8.2281317 3.678611e+00    0.22748187
    ## 90  0.0565459268          0.4092619   7.7859663 3.945385e+00    0.21568972
    ## 91  0.0716710395          0.3685950   7.3533066 4.221280e+00    0.20418249
    ## 92  0.0913774819          0.3316405   6.9312520 4.505551e+00    0.19298072
    ## 93  0.1171892501          0.2980957   6.5207804 4.797344e+00    0.18210244
    ## 94  0.1511783045          0.2676788   6.1227471 5.095689e+00    0.17156315
    ## 95  0.1961752758          0.2401278   5.7378862 5.399507e+00    0.16137588
    ## 96  0.2560661546          0.2151993   5.3668123 5.707612e+00    0.15155130
    ## 97  0.3362119873          0.1926679   5.0100233 6.018716e+00    0.14209771
    ## 98  0.4440453372          0.1723249   4.6679043 6.331438e+00    0.13302121
    ## 99  0.5899219294          0.1539772   4.3407324 6.644312e+00    0.12432575
    ## 100 0.7883423790          0.1374470   4.0286817 6.955795e+00    0.11601328
    ##        Pantala
    ## 1   0.03476592
    ## 2   0.03662018
    ## 3   0.03856088
    ## 4   0.04059133
    ## 5   0.04271491
    ## 6   0.04493507
    ## 7   0.04725537
    ## 8   0.04967945
    ## 9   0.05221102
    ## 10  0.05485389
    ## 11  0.05761194
    ## 12  0.06048913
    ## 13  0.06348951
    ## 14  0.06661722
    ## 15  0.06987644
    ## 16  0.07327147
    ## 17  0.07680665
    ## 18  0.08048642
    ## 19  0.08431526
    ## 20  0.08829773
    ## 21  0.09243847
    ## 22  0.09674216
    ## 23  0.10121354
    ## 24  0.10585740
    ## 25  0.11067861
    ## 26  0.11568205
    ## 27  0.12087266
    ## 28  0.12625541
    ## 29  0.13183530
    ## 30  0.13761737
    ## 31  0.14360667
    ## 32  0.14980827
    ## 33  0.15622725
    ## 34  0.16286869
    ## 35  0.16973767
    ## 36  0.17683925
    ## 37  0.18417850
    ## 38  0.19176043
    ## 39  0.19959005
    ## 40  0.20767231
    ## 41  0.21601211
    ## 42  0.22461432
    ## 43  0.23348371
    ## 44  0.24262501
    ## 45  0.25204283
    ## 46  0.26174172
    ## 47  0.27172610
    ## 48  0.28200032
    ## 49  0.29256855
    ## 50  0.30343488
    ## 51  0.31460324
    ## 52  0.32607739
    ## 53  0.33786095
    ## 54  0.34995735
    ## 55  0.36236986
    ## 56  0.37510152
    ## 57  0.38815519
    ## 58  0.40153351
    ## 59  0.41523887
    ## 60  0.42927345
    ## 61  0.44363915
    ## 62  0.45833764
    ## 63  0.47337028
    ## 64  0.48873819
    ## 65  0.50444217
    ## 66  0.52048271
    ## 67  0.53685999
    ## 68  0.55357388
    ## 69  0.57062390
    ## 70  0.58800923
    ## 71  0.60572869
    ## 72  0.62378074
    ## 73  0.64216346
    ## 74  0.66087456
    ## 75  0.67991137
    ## 76  0.69927078
    ## 77  0.71894932
    ## 78  0.73894308
    ## 79  0.75924775
    ## 80  0.77985859
    ## 81  0.80077041
    ## 82  0.82197762
    ## 83  0.84347417
    ## 84  0.86525356
    ## 85  0.88730886
    ## 86  0.90963269
    ## 87  0.93221720
    ## 88  0.95505413
    ## 89  0.97813471
    ## 90  1.00144978
    ## 91  1.02498969
    ## 92  1.04874434
    ## 93  1.07270322
    ## 94  1.09685533
    ## 95  1.12118927
    ## 96  1.14569320
    ## 97  1.17035482
    ## 98  1.19516146
    ## 99  1.22010000
    ## 100 1.24515694

``` r
AB_fechado_noite <- predict.manyglm(mod_int, newdata = new_data_AB_fechado_noite, type = "response")
CD_fechado_noite <- predict.manyglm(mod_int, newdata = new_data_CD_fechado_noite, type = "response")
EF_fechado_noite <- predict.manyglm(mod_int, newdata = new_data_EF_fechado_noite, type = "response")

fechado_noite <- array(NA, dim = c(nrow(AB_fechado_noite), ncol(AB_fechado_noite), 3))
fechado_noite[,,1] <- AB_fechado_noite
fechado_noite[,,2] <- CD_fechado_noite
fechado_noite[,,3] <- EF_fechado_noite

fechado_noite_mean <- as.data.frame(apply(fechado_noite, MARGIN = c(1,2), FUN = mean))
colnames(fechado_noite_mean) <- colnames(AB_fechado_noite)
```

### Ploting predicted effects

``` r
sp_names <- colnames(control_mean)
sp_names <- gsub("\\.", " ", sp_names)

sp_fonts <- rep(3, ncol(control_mean))
sp_fonts[sp_names == "Chironominae" | sp_names == "Tanypodinae"] <- 1

sp_letters <- c("a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)", "i)", "j)", "k)")


#svg(filename = "Plots/Community structure analysis/predicted_species_trajectories.svg", width = 10, height = 12, pointsize = 13)


close.screen(all.screens = TRUE)
split.screen(matrix(c(0,0.3333333,       0.75,1,
                      0.3333333,0.666666,0.75,1,
                      0.666666,1,        0.75,1,
                      0,0.3333333,       0.5,0.75,
                      0.3333333,0.666666,0.5,0.75,
                      0.666666,1,        0.5,0.75,
                      0,0.3333333,       0.25,0.5,
                      0.3333333,0.666666,0.25,0.5,
                      0.666666,1,        0.25,0.5,
                      0,0.3333333,       0,0.25,
                      0.3333333,0.666666,0,0.25,
                      0.666666,1,        0,0.25), ncol = 4, nrow = 12, byrow = TRUE))
```

    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12

``` r
for(i in 1:11){
  
screen(i)
  
par(mar = c(4,4,1,1), bty = "l")
  
ymax <- max(c(control_mean[,i], fechado_mean[,i], fechado_dia_mean[,i], fechado_noite_mean[,i]))*1.3

plot(NA, xlim = c(32 - 10, 158 + 10), ylim = c(0, ymax), xaxt = "n", yaxt = "n", xlab = "", ylab = "")

lines(y = control_mean[,i],                                                                                                                         x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= "white", lwd = 3)
lines(y = control_mean[,i],
          x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= "#56B4E9", lwd = 2)


lines(y = fechado_mean[,i],                                                                                                                    x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= "white", lwd = 3)
lines(y = fechado_mean[,i],
          x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= "#009E73", lwd = 2)


lines(y = fechado_dia_mean[,i],                                                                                                                         x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= "white", lwd = 3)
lines(y = fechado_dia_mean[,i],
          x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= "#CC79A7", lwd = 2)


lines(y = fechado_noite_mean[,i],                                                                                                                         x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= "white", lwd = 3)
lines(y = fechado_noite_mean[,i],
          x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= "#E69F00", lwd = 2)


axis(1, at = c(32, 80, 116, 158))
title(xlab = "Days", line = 2.5, cex.lab = 1.25)

axis(2, gap.axis = -10, las = 2)
#title(ylab = "Abundance of", line = 3.75, cex.lab = 1.25)
title(ylab = sp_names[i], line = 2.5, cex.lab = 1, font.lab = sp_fonts[i])

letters(x = 10, y = 97, sp_letters[i], cex = 1.5)

}

screen(12)

par(mar = c(0.1,0.1,0.1,0.1), bty = "n")
plot(NA, xlim = c(0, 100), ylim = c(0, 100), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend(legend = c("Control",
                 "Both excl.",
                  "Dragonfly excl.",
                  "Amphibian excl."), x = 50, y = 50, bty = "n", lty = 1, lwd = 3, col = c("#56B4E9", "#009E73","#CC79A7","#E69F00"), xjust = 0.5, yjust = 0.5)
```

![](Effect_first_survey_on_time_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
#dev.off()
```

``` r
comm_trajectories <- rbind(control_mean, fechado_mean, fechado_dia_mean, fechado_noite_mean)

comm_trajectories_relative <- decostand(comm_trajectories, method = "total", MARGIN = 2)

#comm_trajectories_relative <- comm_trajectories_relative[,colnames(comm_trajectories_relative) != "Microvelia"]

set.seed(5); nmds_comm_trajectories <- metaMDS(comm_trajectories_relative, distance = "bray", k = 2, trymax  = 100, try = 50, engine = "monoMDS")
```

    ## Run 0 stress 0.1567881 
    ## Run 1 stress 0.1569253 
    ## ... Procrustes: rmse 0.00176793  max resid 0.03371151 
    ## Run 2 stress 0.1567881 
    ## ... Procrustes: rmse 9.80751e-05  max resid 0.001589101 
    ## ... Similar to previous best
    ## Run 3 stress 0.1602076 
    ## Run 4 stress 0.1600284 
    ## Run 5 stress 0.1570587 
    ## ... Procrustes: rmse 0.003642055  max resid 0.04541914 
    ## Run 6 stress 0.1567879 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002093835  max resid 0.002811186 
    ## ... Similar to previous best
    ## Run 7 stress 0.1602067 
    ## Run 8 stress 0.1567882 
    ## ... Procrustes: rmse 0.0001690787  max resid 0.001721471 
    ## ... Similar to previous best
    ## Run 9 stress 0.1600357 
    ## Run 10 stress 0.1602067 
    ## Run 11 stress 0.1567895 
    ## ... Procrustes: rmse 0.0001704459  max resid 0.002881909 
    ## ... Similar to previous best
    ## Run 12 stress 0.161833 
    ## Run 13 stress 0.1567975 
    ## ... Procrustes: rmse 0.000575074  max resid 0.003300781 
    ## ... Similar to previous best
    ## Run 14 stress 0.156788 
    ## ... Procrustes: rmse 0.0001773843  max resid 0.001788087 
    ## ... Similar to previous best
    ## Run 15 stress 0.1567881 
    ## ... Procrustes: rmse 0.0001934232  max resid 0.001829891 
    ## ... Similar to previous best
    ## Run 16 stress 0.165371 
    ## Run 17 stress 0.1567881 
    ## ... Procrustes: rmse 0.0002019545  max resid 0.001849717 
    ## ... Similar to previous best
    ## Run 18 stress 0.157309 
    ## Run 19 stress 0.1570588 
    ## ... Procrustes: rmse 0.003601611  max resid 0.04547346 
    ## Run 20 stress 0.1567941 
    ## ... Procrustes: rmse 0.0004746372  max resid 0.002943548 
    ## ... Similar to previous best
    ## Run 21 stress 0.1567882 
    ## ... Procrustes: rmse 0.0002120144  max resid 0.00189447 
    ## ... Similar to previous best
    ## Run 22 stress 0.1567881 
    ## ... Procrustes: rmse 0.0002017798  max resid 0.001872537 
    ## ... Similar to previous best
    ## Run 23 stress 0.1570488 
    ## ... Procrustes: rmse 0.003085657  max resid 0.04457712 
    ## Run 24 stress 0.1567883 
    ## ... Procrustes: rmse 0.0002131162  max resid 0.001967099 
    ## ... Similar to previous best
    ## Run 25 stress 0.1600353 
    ## Run 26 stress 0.1567895 
    ## ... Procrustes: rmse 0.0002829833  max resid 0.002218092 
    ## ... Similar to previous best
    ## Run 27 stress 0.1570586 
    ## ... Procrustes: rmse 0.003592124  max resid 0.04543725 
    ## Run 28 stress 0.1569842 
    ## ... Procrustes: rmse 0.002191826  max resid 0.04268784 
    ## Run 29 stress 0.1569117 
    ## ... Procrustes: rmse 0.001917473  max resid 0.00766943 
    ## ... Similar to previous best
    ## Run 30 stress 0.1567892 
    ## ... Procrustes: rmse 0.000357264  max resid 0.006422666 
    ## ... Similar to previous best
    ## Run 31 stress 0.1663415 
    ## Run 32 stress 0.1567878 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001557473  max resid 0.001719447 
    ## ... Similar to previous best
    ## Run 33 stress 0.1567881 
    ## ... Procrustes: rmse 0.0001397673  max resid 0.002411939 
    ## ... Similar to previous best
    ## Run 34 stress 0.1600837 
    ## Run 35 stress 0.1717095 
    ## Run 36 stress 0.1570587 
    ## ... Procrustes: rmse 0.003636336  max resid 0.04544818 
    ## Run 37 stress 0.1567886 
    ## ... Procrustes: rmse 9.71305e-05  max resid 0.001013222 
    ## ... Similar to previous best
    ## Run 38 stress 0.156788 
    ## ... Procrustes: rmse 3.434837e-05  max resid 0.0004766276 
    ## ... Similar to previous best
    ## Run 39 stress 0.156788 
    ## ... Procrustes: rmse 2.224975e-05  max resid 0.0003357224 
    ## ... Similar to previous best
    ## Run 40 stress 0.1567879 
    ## ... Procrustes: rmse 2.472701e-05  max resid 0.0003738466 
    ## ... Similar to previous best
    ## Run 41 stress 0.1569859 
    ## ... Procrustes: rmse 0.002230339  max resid 0.04270595 
    ## Run 42 stress 0.156788 
    ## ... Procrustes: rmse 3.218856e-05  max resid 0.0004508234 
    ## ... Similar to previous best
    ## Run 43 stress 0.1567904 
    ## ... Procrustes: rmse 0.0003307045  max resid 0.003121873 
    ## ... Similar to previous best
    ## Run 44 stress 0.1743174 
    ## Run 45 stress 0.1567881 
    ## ... Procrustes: rmse 4.455981e-05  max resid 0.0006738972 
    ## ... Similar to previous best
    ## Run 46 stress 0.1567881 
    ## ... Procrustes: rmse 4.952399e-05  max resid 0.0006876095 
    ## ... Similar to previous best
    ## Run 47 stress 0.1570896 
    ## ... Procrustes: rmse 0.003850066  max resid 0.04604649 
    ## Run 48 stress 0.1567895 
    ## ... Procrustes: rmse 0.0001489089  max resid 0.00151336 
    ## ... Similar to previous best
    ## Run 49 stress 0.1569267 
    ## ... Procrustes: rmse 0.001719062  max resid 0.0334077 
    ## Run 50 stress 0.1567901 
    ## ... Procrustes: rmse 0.0003038527  max resid 0.003052737 
    ## ... Similar to previous best
    ## *** Best solution repeated 12 times

``` r
points_control <- nmds_comm_trajectories$points[1:100,]
points_fechado <- nmds_comm_trajectories$points[101:200,]
points_fechado_dia <- nmds_comm_trajectories$points[201:300,]
points_fechado_noite <- nmds_comm_trajectories$points[301:400,]

species_pred <- nmds_comm_trajectories$species
```

``` r
find_day <- function(day, min, max, length = 100){
  
  new_day <- day - min
  range <- max - min
  
  result <- round((new_day*length)/range)
  
  if(result == 0){result <- 1}
  
  return(result)
  
}
```

``` r
#screen(2)

#set.seed(1);species_pred <- jitter(species_pred, amount = 0.25)

xmin <- min(c(points_control[,1], points_fechado[,1], points_fechado_dia[,1], points_fechado_noite[,1],species_pred[,1]))*1.1
xmax <- max(c(points_control[,1], points_fechado[,1], points_fechado_dia[,1], points_fechado_noite[,1],species_pred[,1]))*1.1

ymin <- min(c(points_control[,2], points_fechado[,2], points_fechado_dia[,2], points_fechado_noite[,2],species_pred[,2]))*1.1
ymax <- max(c(points_control[,2], points_fechado[,2], points_fechado_dia[,2], points_fechado_noite[,2],species_pred[,2]))*1.1

par(mar = c(4, 4, 1, 1))
  
plot(NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = "", ylab = "", xaxt = "n", yaxt = "n")

lines(x = points_control[,1], y = points_control[,2], lwd = 3, lty = 1, col = "#56B4E9") #controle

lines(x = points_fechado[,1], y = points_fechado[,2], lwd = 3, lty = 1, col = "#009E73") #atrasado

lines(x = points_fechado_dia[,1], y = points_fechado_dia[,2], lwd = 3, lty = 1, col = "#CC79A7") #fechado

lines(x = points_fechado_noite[,1], y = points_fechado_noite[,2], lwd = 3, lty = 1, col = "#E69F00") #fechado



#####################################################################



##################################################

points(x = points_control[find_day(day = 32, min = 32, max = 158),1],
       y = points_control[find_day(day = 32, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#56B4E9", lwd = 3, cex = 3.5)

text(x = points_control[find_day(day = 32, min = 32, max = 158),1],
       y = points_control[find_day(day = 32, min = 32, max = 158),2], labels = "32", cex = 0.7, col = "#56B4E9", font = 2)

points(x = points_control[find_day(day = 80, min = 32, max = 158),1],
       y = points_control[find_day(day = 80, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#56B4E9", lwd = 3, cex = 3.5)

text(x = points_control[find_day(day = 80, min = 32, max = 158),1],
       y = points_control[find_day(day = 80, min = 32, max = 158),2], labels = "80", cex = 0.7, col = "#56B4E9", font = 2)

points(x = points_control[find_day(day = 116, min = 32, max = 158),1],
       y = points_control[find_day(day = 116, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#56B4E9", lwd = 3, cex = 3.5)

text(x = points_control[find_day(day = 116, min = 32, max = 158),1],
       y = points_control[find_day(day = 116, min = 32, max = 158),2], labels = "116", cex = 0.7, col = "#56B4E9", font = 2)

points(x = points_control[find_day(day = 158, min = 32, max = 158),1],
       y = points_control[find_day(day = 158, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#56B4E9", lwd = 3, cex = 3.5)

text(x = points_control[find_day(day = 158, min = 32, max = 158),1],
       y = points_control[find_day(day = 158, min = 32, max = 158),2], labels = "158", cex = 0.7, col = "#56B4E9", font = 2)
#################################################


##################################################

points(x = points_fechado[find_day(day = 32, min = 32, max = 158),1],
       y = points_fechado[find_day(day = 32, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#009E73", lwd = 3, cex = 3.5)

text(x = points_fechado[find_day(day = 32, min = 32, max = 158),1],
       y = points_fechado[find_day(day = 32, min = 32, max = 158),2], labels = "32", cex = 0.7, font = 2, col = "#009E73")

points(x = points_fechado[find_day(day = 80, min = 32, max = 158),1],
       y = points_fechado[find_day(day = 80, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#009E73", lwd = 3, cex = 3.5)

text(x = points_fechado[find_day(day = 80, min = 32, max = 158),1],
       y = points_fechado[find_day(day = 80, min = 32, max = 158),2], labels = "80", cex = 0.7, font = 2, col = "#009E73")

points(x = points_fechado[find_day(day = 116, min = 32, max = 158),1],
       y = points_fechado[find_day(day = 116, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#009E73", lwd = 3, cex = 3.5)

text(x = points_fechado[find_day(day = 116, min = 32, max = 158),1],
       y = points_fechado[find_day(day = 116, min = 32, max = 158),2], labels = "116", cex = 0.7, font = 2, col = "#009E73")

points(x = points_fechado[find_day(day = 158, min = 32, max = 158),1],
       y = points_fechado[find_day(day = 158, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#009E73", lwd = 3, cex = 3.5)

text(x = points_fechado[find_day(day = 158, min = 32, max = 158),1],
       y = points_fechado[find_day(day = 158, min = 32, max = 158),2], labels = "158", cex = 0.7, font = 2, col = "#009E73")
#################################################



##################################################

points(x = points_fechado_dia[find_day(day = 32, min = 32, max = 158),1],
       y = points_fechado_dia[find_day(day = 32, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#CC79A7", lwd = 3, cex = 3.5)

text(x = points_fechado_dia[find_day(day = 32, min = 32, max = 158),1],
       y = points_fechado_dia[find_day(day = 32, min = 32, max = 158),2], labels = "32", cex = 0.7, font = 2, col = "#CC79A7")

points(x = points_fechado_dia[find_day(day = 80, min = 32, max = 158),1],
       y = points_fechado_dia[find_day(day = 80, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#CC79A7", lwd = 3, cex = 3.5)

text(x = points_fechado_dia[find_day(day = 80, min = 32, max = 158),1],
       y = points_fechado_dia[find_day(day = 80, min = 32, max = 158),2], labels = "80", cex = 0.7, font = 2, col = "#CC79A7")

points(x = points_fechado_dia[find_day(day = 116, min = 32, max = 158),1],
       y = points_fechado_dia[find_day(day = 116, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#CC79A7", lwd = 3, cex = 3.5)

text(x = points_fechado_dia[find_day(day = 116, min = 32, max = 158),1],
       y = points_fechado_dia[find_day(day = 116, min = 32, max = 158),2], labels = "116", cex = 0.7, font = 2, col = "#CC79A7")

points(x = points_fechado_dia[find_day(day = 158, min = 32, max = 158),1],
       y = points_fechado_dia[find_day(day = 158, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#CC79A7", lwd = 3, cex = 3.5)

text(x = points_fechado_dia[find_day(day = 158, min = 32, max = 158),1],
       y = points_fechado_dia[find_day(day = 158, min = 32, max = 158),2], labels = "158", cex = 0.7, font = 2, col = "#CC79A7")
#################################################





##################################################

points(x = points_fechado_noite[find_day(day = 32, min = 32, max = 158),1],
       y = points_fechado_noite[find_day(day = 32, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#E69F00", lwd = 3, cex = 3.5)

text(x = points_fechado_noite[find_day(day = 32, min = 32, max = 158),1],
       y = points_fechado_noite[find_day(day = 32, min = 32, max = 158),2], labels = "32", cex = 0.7, font = 2, col = "#E69F00")

points(x = points_fechado_noite[find_day(day = 80, min = 32, max = 158),1],
       y = points_fechado_noite[find_day(day = 80, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#E69F00", lwd = 3, cex = 3.5)

text(x = points_fechado_noite[find_day(day = 80, min = 32, max = 158),1],
       y = points_fechado_noite[find_day(day = 80, min = 32, max = 158),2], labels = "80", cex = 0.7, font = 2, col = "#E69F00")

points(x = points_fechado_noite[find_day(day = 116, min = 32, max = 158),1],
       y = points_fechado_noite[find_day(day = 116, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#E69F00", lwd = 3, cex = 3.5)

text(x = points_fechado_noite[find_day(day = 116, min = 32, max = 158),1],
       y = points_fechado_noite[find_day(day = 116, min = 32, max = 158),2], labels = "116", cex = 0.7, font = 2, col = "#E69F00")

points(x = points_fechado_noite[find_day(day = 158, min = 32, max = 158),1],
       y = points_fechado_noite[find_day(day = 158, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#E69F00", lwd = 3, cex = 3.5)

text(x = points_fechado_noite[find_day(day = 158, min = 32, max = 158),1],
       y = points_fechado_noite[find_day(day = 158, min = 32, max = 158),2], labels = "158", cex = 0.7, font = 2, col = "#E69F00")
#################################################



########################## species_pred #################################


sp_names <- rownames(species_pred)
sp_names <- substr(sp_names, start = 1, stop = 3)


for(i in 1:nrow(species_pred)){
  #points(x = species_pred[i,1],
   #    y = species_pred[i,2],
    #   pch = 21, bg = "white", col = "grey50", lwd = 2, cex = 4)

text_contour(x = species_pred[i,1],
       y = species_pred[i,2], labels = sp_names[i], cex = 0.8, col = "white", thick = 0.004)
  
text(x = species_pred[i,1],
       y = species_pred[i,2], labels = sp_names[i], cex = 0.8, col = "grey30")
}




axis(1)
axis(2, cex.axis = 1, las = 2)

title(xlab = "nMDS 1", line = 3, cex.lab = 1.2)
title(ylab = "nMDS 2", line = 3, cex.lab = 1.2)

#letters(x = 10, y = 97, "b)", cex = 1.5)


par(new = TRUE)
plot(NA, xlim = c(0, 100), ylim = c(0, 100), xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")

legend(legend = c("Control",
                 "Both excl.",
                  "Dragonfly excl.",
                  "Amphibian excl."), x = 0, y = 0, bty = "n", lty = 1, lwd = 3, col = c("#56B4E9", "#009E73","#CC79A7","#E69F00"), xjust = 0, yjust = 0)
```

![](Effect_first_survey_on_time_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
#dev.off()
```
