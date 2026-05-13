Effect of first survey on the effect of time
================
Rodolfo Pelinson
2026-05-13

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
    ## time zone: Europe/London
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
    ## Run 1 stress 0.231817 
    ## Run 2 stress 0.230996 
    ## Run 3 stress 0.2309964 
    ## Run 4 stress 0.234659 
    ## Run 5 stress 0.1917598 
    ## ... New best solution
    ## ... Procrustes: rmse 5.5457e-06  max resid 1.838722e-05 
    ## ... Similar to previous best
    ## Run 6 stress 0.1917598 
    ## ... Procrustes: rmse 1.487994e-06  max resid 4.033407e-06 
    ## ... Similar to previous best
    ## Run 7 stress 0.1917598 
    ## ... Procrustes: rmse 1.75744e-06  max resid 5.838785e-06 
    ## ... Similar to previous best
    ## Run 8 stress 0.2329255 
    ## Run 9 stress 0.2355307 
    ## Run 10 stress 0.2382912 
    ## Run 11 stress 0.1917598 
    ## ... Procrustes: rmse 6.779448e-06  max resid 2.309183e-05 
    ## ... Similar to previous best
    ## Run 12 stress 0.2365435 
    ## Run 13 stress 0.1917598 
    ## ... Procrustes: rmse 2.178956e-06  max resid 7.282448e-06 
    ## ... Similar to previous best
    ## Run 14 stress 0.234434 
    ## Run 15 stress 0.1917598 
    ## ... Procrustes: rmse 1.567165e-06  max resid 4.483632e-06 
    ## ... Similar to previous best
    ## Run 16 stress 0.1917598 
    ## ... Procrustes: rmse 1.163158e-06  max resid 3.712417e-06 
    ## ... Similar to previous best
    ## Run 17 stress 0.2284807 
    ## Run 18 stress 0.2344339 
    ## Run 19 stress 0.1917598 
    ## ... New best solution
    ## ... Procrustes: rmse 1.30027e-06  max resid 3.084484e-06 
    ## ... Similar to previous best
    ## Run 20 stress 0.1917598 
    ## ... Procrustes: rmse 1.292555e-05  max resid 4.683998e-05 
    ## ... Similar to previous best
    ## *** Best solution repeated 2 times

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
    ##  Resampling run 200 finished. Time elapsed: 0.11 minutes...
    ##  Resampling run 300 finished. Time elapsed: 0.16 minutes...
    ##  Resampling run 400 finished. Time elapsed: 0.21 minutes...
    ##  Resampling run 500 finished. Time elapsed: 0.27 minutes...
    ##  Resampling run 600 finished. Time elapsed: 0.32 minutes...
    ##  Resampling run 700 finished. Time elapsed: 0.37 minutes...
    ##  Resampling run 800 finished. Time elapsed: 0.42 minutes...
    ##  Resampling run 900 finished. Time elapsed: 0.47 minutes...
    ## Time elapsed: 0 hr 0 min 31 sec

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
    ## Time elapsed: 0 hr 0 min 21 sec

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
    ##  Resampling run 200 finished. Time elapsed: 0.07 minutes...
    ##  Resampling run 300 finished. Time elapsed: 0.10 minutes...
    ##  Resampling run 400 finished. Time elapsed: 0.13 minutes...
    ##  Resampling run 500 finished. Time elapsed: 0.17 minutes...
    ##  Resampling run 600 finished. Time elapsed: 0.20 minutes...
    ##  Resampling run 700 finished. Time elapsed: 0.23 minutes...
    ##  Resampling run 800 finished. Time elapsed: 0.27 minutes...
    ##  Resampling run 900 finished. Time elapsed: 0.30 minutes...
    ## Time elapsed: 0 hr 0 min 20 sec

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
    ##  Resampling run 200 finished. Time elapsed: 0.07 minutes...
    ##  Resampling run 300 finished. Time elapsed: 0.10 minutes...
    ##  Resampling run 400 finished. Time elapsed: 0.13 minutes...
    ##  Resampling run 500 finished. Time elapsed: 0.16 minutes...
    ##  Resampling run 600 finished. Time elapsed: 0.20 minutes...
    ##  Resampling run 700 finished. Time elapsed: 0.23 minutes...
    ##  Resampling run 800 finished. Time elapsed: 0.26 minutes...
    ##  Resampling run 900 finished. Time elapsed: 0.29 minutes...
    ## Time elapsed: 0 hr 0 min 19 sec

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
    ##  Resampling run 200 finished. Time elapsed: 0.10 minutes...
    ##  Resampling run 300 finished. Time elapsed: 0.14 minutes...
    ##  Resampling run 400 finished. Time elapsed: 0.19 minutes...
    ##  Resampling run 500 finished. Time elapsed: 0.24 minutes...
    ##  Resampling run 600 finished. Time elapsed: 0.29 minutes...
    ##  Resampling run 700 finished. Time elapsed: 0.34 minutes...
    ##  Resampling run 800 finished. Time elapsed: 0.39 minutes...
    ##  Resampling run 900 finished. Time elapsed: 0.43 minutes...
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
high_MDS1_high_MDS2_mean
```

    ##           Aedes Callibaetis Chironominae      Culex Dendropsophus.minutus
    ## 1   11.99895854  0.04850235     8.173712 143.308578              8.797558
    ## 2   11.95167625  0.04609302     8.343942 125.917785              8.812179
    ## 3   11.88964602  0.04390795     8.515753 110.857952              8.826591
    ## 4   11.81309962  0.04192632     8.689095  97.793844              8.840794
    ## 5   11.72232194  0.04012969     8.863920  86.441256              8.854786
    ## 6   11.61764923  0.03850175     9.040177  76.558869              8.868566
    ## 7   11.49946701  0.03702803     9.217811  67.941454              8.882134
    ## 8   11.36820771  0.03569573     9.396766  60.414203              8.895487
    ## 9   11.22434807  0.03449352     9.576986  53.827989              8.908625
    ## 10  11.06840618  0.03341137     9.758410  48.055395              8.921547
    ## 11  10.90093845  0.03244044     9.940976  42.987385              8.934252
    ## 12  10.72253626  0.03157291    10.124621  38.530513              8.946739
    ## 13  10.53382251  0.03080194    10.309280  34.604570              8.959007
    ## 14  10.33544800  0.03012153    10.494885  31.140601              8.971055
    ## 15  10.12808769  0.02952647    10.681366  28.079245              8.982882
    ## 16   9.91243691  0.02901226    10.868652  25.369315              8.994488
    ## 17   9.68920749  0.02857506    11.056670  22.966613              9.005871
    ## 18   9.45912389  0.02821164    11.245345  20.832915              9.017030
    ## 19   9.22291927  0.02791933    11.434602  18.935119              9.027965
    ## 20   8.98133174  0.02769601    11.624360  17.244512              9.038674
    ## 21   8.73510049  0.02754007    11.814541  15.736156              9.049158
    ## 22   8.48496216  0.02745038    12.005063  14.388361              9.059414
    ## 23   8.23164725  0.02742630    12.195842  13.182229              9.069443
    ## 24   7.97587666  0.02746766    12.386795  12.101279              9.079243
    ## 25   7.71835846  0.02757475    12.577834  11.131113              9.088814
    ## 26   7.45978475  0.02774835    12.768873  10.259136              9.098155
    ## 27   7.20082878  0.02798970    12.959823   9.474316              9.107266
    ## 28   6.94214227  0.02830054    13.150593   8.766976              9.116144
    ## 29   6.68435293  0.02868316    13.341092   8.128618              9.124791
    ## 30   6.42806221  0.02914034    13.531228   7.551765              9.133205
    ## 31   6.17384337  0.02967549    13.720907   7.029834              9.141385
    ## 32   5.92223968  0.03029260    13.910035   6.557022              9.149332
    ## 33   5.67376298  0.03099638    14.098515   6.128201              9.157043
    ## 34   5.42889242  0.03179221    14.286253   5.738843              9.164519
    ## 35   5.18807349  0.03268633    14.473151   5.384936              9.171759
    ## 36   4.95171727  0.03368582    14.659110   5.062926              9.178763
    ## 37   4.72019995  0.03479875    14.844033   4.769662              9.185530
    ## 38   4.49386260  0.03603426    15.027820   4.502341              9.192059
    ## 39   4.27301110  0.03740272    15.210372   4.258476              9.198350
    ## 40   4.05791633  0.03891584    15.391590   4.035848              9.204402
    ## 41   3.84881461  0.04058682    15.571373   3.832484              9.210215
    ## 42   3.64590823  0.04243061    15.749622   3.646622              9.215789
    ## 43   3.44936626  0.04446405    15.926235   3.476690              9.221122
    ## 44   3.25932544  0.04670618    16.101113   3.321285              9.226216
    ## 45   3.07589126  0.04917849    16.274155   3.179152              9.231068
    ## 46   2.89913918  0.05190528    16.445261   3.049167              9.235680
    ## 47   2.72911597  0.05491405    16.614332   2.930327              9.240049
    ## 48   2.56584108  0.05823592    16.781268   2.821732              9.244177
    ## 49   2.40930816  0.06190618    16.945971   2.722579              9.248063
    ## 50   2.25948669  0.06596484    17.108342   2.632146              9.251706
    ## 51   2.11632353  0.07045741    17.268283   2.549790              9.255107
    ## 52   1.97974461  0.07543559    17.425698   2.474934              9.258264
    ## 53   1.84965664  0.08095832    17.580490   2.407065              9.261178
    ## 54   1.72594876  0.08709280    17.732565   2.345724              9.263848
    ## 55   1.60849427  0.09391577    17.881829   2.290503              9.266275
    ## 56   1.49715228  0.10151504    18.028188   2.241040              9.268457
    ## 57   1.39176935  0.10999115    18.171552   2.197017              9.270396
    ## 58   1.29218115  0.11945949    18.311829   2.158152              9.272090
    ## 59   1.19821394  0.13005263    18.448931   2.124200              9.273540
    ## 60   1.10968616  0.14192312    18.582771   2.094951              9.274745
    ## 61   1.02640985  0.15524681    18.713264   2.070223              9.275706
    ## 62   0.94819201  0.17022674    18.840324   2.049865              9.276422
    ## 63   0.87483592  0.18709769    18.963871   2.033754              9.276893
    ## 64   0.80614239  0.20613162    19.083823   2.021791              9.277119
    ## 65   0.74191089  0.22764407    19.200103   2.013906              9.277101
    ## 66   0.68194063  0.25200178    19.312635   2.010050              9.276838
    ## 67   0.62603157  0.27963172    19.421344   2.010201              9.276330
    ## 68   0.57398525  0.31103180    19.526158   2.014359              9.275577
    ## 69   0.52560570  0.34678372    19.627008   2.022550              9.274579
    ## 70   0.48070013  0.38756821    19.723827   2.034823              9.273337
    ## 71   0.43907959  0.43418333    19.816550   2.051251              9.271850
    ## 72   0.40055955  0.48756630    19.905116   2.071933              9.270119
    ## 73   0.36496043  0.54881979    19.989463   2.096997              9.268144
    ## 74   0.33210797  0.61924340    20.069537   2.126594              9.265924
    ## 75   0.30183365  0.70037165    20.145282   2.160908              9.263461
    ## 76   0.27397494  0.79401969    20.216647   2.200153              9.260754
    ## 77   0.24837557  0.90233859    20.283583   2.244576              9.257803
    ## 78   0.22488565  1.02788218    20.346046   2.294461              9.254609
    ## 79   0.20336183  1.17368805    20.403991   2.350130              9.251172
    ## 80   0.18366736  1.34337597    20.457381   2.411949              9.247493
    ## 81   0.16567210  1.54126746    20.506178   2.480328              9.243570
    ## 82   0.14925249  1.77253162    20.550348   2.555731              9.239406
    ## 83   0.13429153  2.04336295    20.589861   2.638675              9.235000
    ## 84   0.12067867  2.36119896    20.624691   2.729742              9.230352
    ## 85   0.10830967  2.73498655    20.654812   2.829581              9.225463
    ## 86   0.09708649  3.17550903    20.680205   2.938919              9.220334
    ## 87   0.08691710  3.69578801    20.700852   3.058567              9.214964
    ## 88   0.07771528  4.31157850    20.716738   3.189431              9.209354
    ## 89   0.06940049  5.04197963    20.727852   3.332525              9.203504
    ## 90   0.06189754  5.91018969    20.734187   3.488979              9.197416
    ## 91   0.05513650  6.94444107    20.735739   3.660061              9.191089
    ## 92   0.04905234  8.17916028    20.732505   3.847186              9.184524
    ## 93   0.04358481  9.65640972    20.724490   4.051938              9.177722
    ## 94   0.03867813 11.42768326    20.711697   4.276096              9.170683
    ## 95   0.03428077 13.55614635    20.694136   4.521649              9.163407
    ## 96   0.03034523 16.11943618    20.671820   4.790835              9.155895
    ## 97   0.02682781 19.21316841    20.644763   5.086165              9.148148
    ## 98   0.02368835 22.95533718    20.612984   5.410465              9.140167
    ## 99   0.02089004 27.49184630    20.576505   5.766916              9.131951
    ## 100  0.01839918 33.00347593    20.535351   6.159104              9.123502
    ##       Microvelia Scinax.fuscovarius Tanypodinae    Dasyhelea Erythrodiplax
    ## 1   0.9686609015        207.6948917    1.276219 0.0004320865    0.04184581
    ## 2   0.6816058315        190.6704577    1.359761 0.0005073910    0.04499056
    ## 3   0.4830239943        175.1076403    1.446932 0.0005949583    0.04832285
    ## 4   0.3447291378        160.8758576    1.537734 0.0006966298    0.05184961
    ## 5   0.2477770979        147.8566108    1.632159 0.0008144966    0.05557765
    ## 6   0.1793569506        135.9423307    1.730180 0.0009509295    0.05951365
    ## 7   0.1307522411        125.0353370    1.831759 0.0011086107    0.06366411
    ## 8   0.0959961802        115.0468994    1.936837 0.0012905699    0.06803535
    ## 9   0.0709794522        105.8963915    2.045341 0.0015002229    0.07263337
    ## 10  0.0528548951         97.5105254    2.157180 0.0017414129    0.07746394
    ## 11  0.0396379940         89.8226625    2.272243 0.0020184569    0.08253244
    ## 12  0.0299372567         82.7721898    2.390404 0.0023361941    0.08784389
    ## 13  0.0227712156         76.3039566    2.511513 0.0027000396    0.09340285
    ## 14  0.0174435273         70.3677660    2.635407 0.0031160404    0.09921343
    ## 15  0.0134572489         64.9179149    2.761898 0.0035909368    0.10527919
    ## 16  0.0104556766         59.9127774    2.890784 0.0041322272    0.11160312
    ## 17  0.0081812917         55.3144294    3.021840 0.0047482367    0.11818760
    ## 18  0.0064471161         51.0883073    3.154824 0.0054481904    0.12503431
    ## 19  0.0051166180         47.2029001    3.289476 0.0062422899    0.13214424
    ## 20  0.0040895391         43.6294707    3.425518 0.0071417941    0.13951761
    ## 21  0.0032918469         40.3418024    3.562654 0.0081591040    0.14715382
    ## 22  0.0026685710         37.3159708    3.700572 0.0093078494    0.15505144
    ## 23  0.0021786717         34.5301357    3.838947 0.0106029808    0.16320813
    ## 24  0.0017913427         31.9643529    3.977435 0.0120608617    0.17162065
    ## 25  0.0014833358         29.6004042    4.115685 0.0136993653    0.18028476
    ## 26  0.0012370126         27.4216422    4.253328 0.0155379707    0.18919526
    ## 27  0.0010389214         25.4128498    4.389991 0.0175978613    0.19834590
    ## 28  0.0008787495         23.5601132    4.525288 0.0199020230    0.20772938
    ## 29  0.0007485510         21.8507057    4.658829 0.0224753416    0.21733735
    ## 30  0.0006421723         20.2729830    4.790217 0.0253446980    0.22716037
    ## 31  0.0005548245         18.8162873    4.919053 0.0285390611    0.23718788
    ## 32  0.0004827626         17.4708611    5.044936 0.0320895765    0.24740824
    ## 33  0.0004230439         16.2277675    5.167467 0.0360296490    0.25780870
    ## 34  0.0003733457         15.0788192    5.286250 0.0403950192    0.26837541
    ## 35  0.0003318262         14.0165131    5.400892 0.0452238309    0.27909343
    ## 36  0.0002970190         13.0339704    5.511011 0.0505566887    0.28994676
    ## 37  0.0002677513         12.1248833    5.616230 0.0564367031    0.30091832
    ## 38  0.0002430820         11.2834653    5.716187 0.0629095229    0.31199005
    ## 39  0.0002222531         10.5044065    5.810531 0.0700233510    0.32314289
    ## 40  0.0002046524          9.7828328    5.898929 0.0778289428    0.33435683
    ## 41  0.0001897841          9.1142689    5.981062 0.0863795852    0.34561099
    ## 42  0.0001772460          8.4946040    6.056635 0.0957310538    0.35688363
    ## 43  0.0001667121          7.9200611    6.125371 0.1059415462    0.36815224
    ## 44  0.0001579180          7.3871687    6.187016 0.1170715899    0.37939361
    ## 45  0.0001506503          6.8927352    6.241343 0.1291839228    0.39058388
    ## 46  0.0001447379          6.4338253    6.288147 0.1423433433    0.40169864
    ## 47  0.0001400453          6.0077386    6.327253 0.1566165303    0.41271300
    ## 48  0.0001364673          5.6119899    6.358514 0.1720718296    0.42360168
    ## 49  0.0001339252          5.2442916    6.381811 0.1887790069    0.43433912
    ## 50  0.0001323641          4.9025369    6.397056 0.2068089643    0.44489954
    ## 51  0.0001317504          4.5847853    6.404191 0.2262334219    0.45525708
    ## 52  0.0001320710          4.2892488    6.403188 0.2471245614    0.46538587
    ## 53  0.0001333328          4.0142790    6.394052 0.2695546330    0.47526018
    ## 54  0.0001355627          3.7583564    6.376817 0.2935955253    0.48485447
    ## 55  0.0001388089          3.5200795    6.351548 0.3193182979    0.49414354
    ## 56  0.0001431425          3.2981550    6.318342 0.3467926796    0.50310262
    ## 57  0.0001486598          3.0913897    6.277325 0.3760865303    0.51170749
    ## 58  0.0001554864          2.8986818    6.228650 0.4072652725    0.51993459
    ## 59  0.0001637816          2.7190138    6.172501 0.4403912907    0.52776112
    ## 60  0.0001737447          2.5514460    6.109086 0.4755233049    0.53516513
    ## 61  0.0001856231          2.3951099    6.038642 0.5127157194    0.54212567
    ## 62  0.0001997223          2.2492027    5.961426 0.5520179520    0.54862283
    ## 63  0.0002164187          2.1129822    5.877721 0.5934737474    0.55463787
    ## 64  0.0002361766          1.9857619    5.787828 0.6371204802    0.56015330
    ## 65  0.0002595690          1.8669067    5.692070 0.6829884525    0.56515299
    ## 66  0.0002873047          1.7558286    5.590783 0.7311001930    0.56962218
    ## 67  0.0003202628          1.6519837    5.484323 0.7814697621    0.57354764
    ## 68  0.0003595375          1.5548677    5.373055 0.8341020722    0.57691768
    ## 69  0.0004064954          1.4640141    5.257356 0.8889922280    0.57972221
    ## 70  0.0004628509          1.3789901    5.137613 0.9461248958    0.58195281
    ## 71  0.0005307627          1.2993948    5.014219 1.0054737081    0.58360276
    ## 72  0.0006129620          1.2248565    4.887571 1.0670007125    0.58466710
    ## 73  0.0007129197          1.1550303    4.758070 1.1306558710    0.58514259
    ## 74  0.0008350674          1.0895964    4.626115 1.1963766195    0.58502780
    ## 75  0.0009850910          1.0282578    4.492105 1.2640874927    0.58432308
    ## 76  0.0011703211          0.9707390    4.356435 1.3336998232    0.58303055
    ## 77  0.0014002565          0.9167840    4.219495 1.4051115221    0.58115414
    ## 78  0.0016872679          0.8661550    4.081668 1.4782069463    0.57869948
    ## 79  0.0020475494          0.8186313    3.943326 1.5528568587    0.57567397
    ## 80  0.0025024108          0.7740075    3.804833 1.6289184884    0.57208668
    ## 81  0.0030800424          0.7320927    3.666540 1.7062356932    0.56794832
    ## 82  0.0038179363          0.6927094    3.528784 1.7846392291    0.56327119
    ## 83  0.0047662249          0.6556925    3.391890 1.8639471307    0.55806913
    ## 84  0.0059923100          0.6208882    3.256164 1.9439652019    0.55235744
    ## 85  0.0075873111          0.5881535    3.121897 2.0244876205    0.54615279
    ## 86  0.0096750984          0.5573552    2.989364 2.1052976542    0.53947316
    ## 87  0.0124250107          0.5283692    2.858821 2.1861684873    0.53233776
    ## 88  0.0160698577          0.5010800    2.730506 2.2668641555    0.52476691
    ## 89  0.0209315395          0.4753798    2.604636 2.3471405847    0.51678197
    ## 90  0.0274577024          0.4511682    2.481412 2.4267467282    0.50840522
    ## 91  0.0362744717          0.4283515    2.361014 2.5054257968    0.49965978
    ## 92  0.0482627295          0.4068424    2.243604 2.5829165733    0.49056948
    ## 93  0.0646690593          0.3865594    2.129324 2.6589548033    0.48115877
    ## 94  0.0872680130          0.3674264    2.018297 2.7332746523    0.47145258
    ## 95  0.1186007711          0.3493724    1.910629 2.8056102178    0.46147627
    ## 96  0.1623281422          0.3323310    1.806407 2.8756970867    0.45125547
    ## 97  0.2237556471          0.3162404    1.705700 2.9432739220    0.44081596
    ## 98  0.3106190389          0.3010425    1.608561 3.0080840701    0.43018364
    ## 99  0.4342661906          0.2866833    1.515027 3.0698771707    0.41938435
    ## 100 0.6114456325          0.2731122    1.425119 3.1284107590    0.40844377
    ##        Pantala
    ## 1   0.01920183
    ## 2   0.02119671
    ## 3   0.02337237
    ## 4   0.02574221
    ## 5   0.02832026
    ## 6   0.03112127
    ## 7   0.03416064
    ## 8   0.03745443
    ## 9   0.04101937
    ## 10  0.04487281
    ## 11  0.04903274
    ## 12  0.05351773
    ## 13  0.05834688
    ## 14  0.06353986
    ## 15  0.06911676
    ## 16  0.07509813
    ## 17  0.08150484
    ## 18  0.08835807
    ## 19  0.09567921
    ## 20  0.10348980
    ## 21  0.11181140
    ## 22  0.12066551
    ## 23  0.13007350
    ## 24  0.14005642
    ## 25  0.15063497
    ## 26  0.16182929
    ## 27  0.17365889
    ## 28  0.18614248
    ## 29  0.19929780
    ## 30  0.21314154
    ## 31  0.22768910
    ## 32  0.24295451
    ## 33  0.25895019
    ## 34  0.27568687
    ## 35  0.29317334
    ## 36  0.31141637
    ## 37  0.33042049
    ## 38  0.35018784
    ## 39  0.37071803
    ## 40  0.39200800
    ## 41  0.41405183
    ## 42  0.43684065
    ## 43  0.46036251
    ## 44  0.48460224
    ## 45  0.50954137
    ## 46  0.53515803
    ## 47  0.56142688
    ## 48  0.58831906
    ## 49  0.61580215
    ## 50  0.64384013
    ## 51  0.67239342
    ## 52  0.70141884
    ## 53  0.73086970
    ## 54  0.76069586
    ## 55  0.79084380
    ## 56  0.82125672
    ## 57  0.85187470
    ## 58  0.88263485
    ## 59  0.91347145
    ## 60  0.94431623
    ## 61  0.97509850
    ## 62  1.00574548
    ## 63  1.03618250
    ## 64  1.06633331
    ## 65  1.09612041
    ## 66  1.12546532
    ## 67  1.15428893
    ## 68  1.18251187
    ## 69  1.21005483
    ## 70  1.23683894
    ## 71  1.26278617
    ## 72  1.28781964
    ## 73  1.31186407
    ## 74  1.33484608
    ## 75  1.35669464
    ## 76  1.37734136
    ## 77  1.39672090
    ## 78  1.41477129
    ## 79  1.43143426
    ## 80  1.44665555
    ## 81  1.46038523
    ## 82  1.47257794
    ## 83  1.48319315
    ## 84  1.49219539
    ## 85  1.49955445
    ## 86  1.50524554
    ## 87  1.50924943
    ## 88  1.51155257
    ## 89  1.51214714
    ## 90  1.51103114
    ## 91  1.50820834
    ## 92  1.50368831
    ## 93  1.49748635
    ## 94  1.48962340
    ## 95  1.48012591
    ## 96  1.46902572
    ## 97  1.45635986
    ## 98  1.44217035
    ## 99  1.42650399
    ## 100 1.40941205

``` r
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
low_MDS1_high_MDS2_mean
```

    ##            Aedes  Callibaetis Chironominae     Culex Dendropsophus.minutus
    ## 1   17.087666368  0.001319290     29.72936 26.608207              1.990859
    ## 2   16.798348755  0.001713952     29.23765 24.617315              1.931771
    ## 3   16.494231017  0.002218600     28.76352 22.809805              1.875960
    ## 4   16.176300062  0.002861421     28.30638 21.166948              1.823242
    ## 5   15.845573324  0.003677109     27.86567 19.672099              1.773447
    ## 6   15.503093363  0.004708184     27.44084 18.310448              1.726413
    ## 7   15.149922436  0.006006514     27.03137 17.068802              1.681993
    ## 8   14.787137053  0.007635084     26.63677 15.935397              1.640048
    ## 9   14.415822582  0.009670018     26.25655 14.899735              1.600449
    ## 10  14.037067917  0.012202895     25.89028 13.952435              1.563075
    ## 11  13.651960264  0.015343367     25.53750 13.085106              1.527815
    ## 12  13.261580069  0.019222092     25.19781 12.290238              1.494564
    ## 13  12.866996115  0.023994007     24.87082 11.561100              1.463225
    ## 14  12.469260836  0.029841942     24.55614 10.891653              1.433708
    ## 15  12.069405858  0.036980566     24.25342 10.276476              1.405928
    ## 16  11.668437793  0.045660661     23.96230  9.710698              1.379806
    ## 17  11.267334320  0.056173696     23.68246  9.189936              1.355271
    ## 18  10.867040563  0.068856664     23.41358  8.710244              1.332255
    ## 19  10.468465771  0.084097121     23.15536  8.268066              1.310693
    ## 20  10.072480342  0.102338359     22.90753  7.860196              1.290529
    ## 21   9.679913166  0.124084615     22.66979  7.483738              1.271708
    ## 22   9.291549316  0.149906200     22.44190  7.136078              1.254180
    ## 23   8.908128087  0.180444405     22.22361  6.814852              1.237898
    ## 24   8.530341376  0.216416026     22.01467  6.517920              1.222822
    ## 25   8.158832402  0.258617307     21.81487  6.243346              1.208910
    ## 26   7.794194771  0.307927106     21.62399  5.989377              1.196129
    ## 27   7.436971864  0.365309045     21.44182  5.754421              1.184445
    ## 28   7.087656548  0.431812391     21.26818  5.537037              1.173828
    ## 29   6.746691195  0.508571413     21.10289  5.335916              1.164251
    ## 30   6.414467988  0.596802946     20.94576  5.149872              1.155692
    ## 31   6.091329516  0.697801891     20.79663  4.977825              1.148128
    ## 32   5.777569616  0.812934387     20.65536  4.818796              1.141541
    ## 33   5.473434465  0.943628430     20.52179  4.671898              1.135913
    ## 34   5.179123887  1.091361701     20.39579  4.536323              1.131233
    ## 35   4.894792870  1.257646453     20.27723  4.411338              1.127487
    ## 36   4.620553250  1.444011332     20.16598  4.296279              1.124668
    ## 37   4.356475565  1.651980063     20.06194  4.190544              1.122767
    ## 38   4.102591039  1.883047051     19.96500  4.093589              1.121780
    ## 39   3.858893685  2.138649976     19.87506  4.004919              1.121706
    ## 40   3.625342498  2.420139613     19.79203  3.924091              1.122543
    ## 41   3.401863725  2.728747173     19.71582  3.850705              1.124293
    ## 42   3.188353190  3.065549597     19.64637  3.784402              1.126962
    ## 43   2.984678645  3.431433323     19.58359  3.724860              1.130556
    ## 44   2.790682148  3.827057184     19.52744  3.671796              1.135082
    ## 45   2.606182426  4.252815167     19.47785  3.624957              1.140553
    ## 46   2.430977232  4.708799886     19.43477  3.584124              1.146982
    ## 47   2.264845664  5.194767680     19.39816  3.549106              1.154385
    ## 48   2.107550434  5.710106323     19.36799  3.519741              1.162780
    ## 49   1.958840079  6.253806357     19.34422  3.495893              1.172188
    ## 50   1.818451106  6.824437076     19.32684  3.477455              1.182633
    ## 51   1.686110051  7.420128177     19.31582  3.464341              1.194140
    ## 52   1.561535449  8.038558018     19.31115  3.456492              1.206740
    ## 53   1.444439714  8.676949370     19.31283  3.453872              1.220464
    ## 54   1.334530903  9.332073391     19.32086  3.456470              1.235347
    ## 55   1.231514391 10.000262436     19.33525  3.464297              1.251428
    ## 56   1.135094422 10.677432106     19.35601  3.477388              1.268749
    ## 57   1.044975551 11.359112724     19.38316  3.495804              1.287355
    ## 58   0.960863976 12.040490230     19.41673  3.519628              1.307296
    ## 59   0.882468750 12.716456191     19.45676  3.548970              1.328624
    ## 60   0.809502885 13.381666397     19.50328  3.583964              1.351398
    ## 61   0.741684338 14.030607241     19.55634  3.624772              1.375680
    ## 62   0.678736896 14.657668856     19.61599  3.671585              1.401537
    ## 63   0.620390943 15.257223698     19.68229  3.724622              1.429040
    ## 64   0.566384138 15.823709131     19.75531  3.784136              1.458267
    ## 65   0.516461981 16.351712329     19.83511  3.850410              1.489301
    ## 66   0.470378295 16.836055724     19.92179  3.923766              1.522232
    ## 67   0.427895612 17.271881132     20.01543  4.004561              1.557156
    ## 68   0.388785476 17.654730640     20.11611  4.093197              1.594176
    ## 69   0.352828676 17.980622403     20.22396  4.190116              1.633403
    ## 70   0.319815391 18.246119518     20.33906  4.295813              1.674955
    ## 71   0.289545277 18.448390350     20.46154  4.410831              1.718961
    ## 72   0.261827497 18.585258791     20.59153  4.535772              1.765556
    ## 73   0.236480679 18.655243231     20.72916  4.671301              1.814889
    ## 74   0.213332840 18.657583261     20.87457  4.818150              1.867116
    ## 75   0.192221255 18.592253423     21.02791  4.977125              1.922408
    ## 76   0.172992289 18.459963670     21.18934  5.149115              1.980946
    ## 77   0.155501192 18.262146529     21.35903  5.335099              2.042925
    ## 78   0.139611869 18.000931268     21.53716  5.536153              2.108557
    ## 79   0.125196619 17.679105757     21.72391  5.753465              2.178065
    ## 80   0.112135857 17.300066934     21.91949  5.988344              2.251694
    ## 81   0.100317812 16.867761144     22.12410  6.242230              2.329704
    ## 82   0.089638223 16.386615794     22.33796  6.516713              2.412376
    ## 83   0.080000014 15.861463983     22.56130  6.813546              2.500012
    ## 84   0.071312972 15.297463903     22.79437  7.134665              2.592938
    ## 85   0.063493409 14.700014878     23.03742  7.482209              2.691503
    ## 86   0.056463839 14.074671945     23.29071  7.858539              2.796087
    ## 87   0.050152640 13.427060855     23.55453  8.266271              2.907095
    ## 88   0.044493732 12.762795279     23.82917  8.708297              3.024967
    ## 89   0.039426254 12.087397888     24.11494  9.187823              3.150176
    ## 90   0.034894247 11.406226807     24.41216  9.708404              3.283235
    ## 91   0.030846351 10.724408726     24.72117 10.273983              3.424696
    ## 92   0.027235503 10.046779741     25.04232 10.888940              3.575155
    ## 93   0.024018653  9.377834726     25.37598 11.558147              3.735258
    ## 94   0.021156486  8.721685805     25.72254 12.287020              3.905703
    ## 95   0.018613158  8.082030215     26.08241 13.081597              4.087245
    ## 96   0.016356043  7.462127609     26.45601 13.948604              4.280702
    ## 97   0.014355491  6.864786628     26.84378 14.895549              4.486960
    ## 98   0.012584603  6.292360323     27.24620 15.930818              4.706979
    ## 99   0.011019011  5.746749865     27.66373 17.063788              4.941800
    ## 100  0.009636678  5.229415795     28.09690 18.304952              5.192553
    ##     Microvelia Scinax.fuscovarius Tanypodinae   Dasyhelea Erythrodiplax
    ## 1   0.07663103          2.2398563   13.116005 0.002019766    0.08156934
    ## 2   0.07725306          2.2437210   13.290883 0.002381627    0.08644849
    ## 3   0.07794350          2.2462265   13.458838 0.002802735    0.09153452
    ## 4   0.07870408          2.2473681   13.619551 0.003291744    0.09682991
    ## 5   0.07953672          2.2471439   13.772713 0.003858386    0.10233664
    ## 6   0.08044357          2.2455542   13.918028 0.004513580    0.10805625
    ## 7   0.08142694          2.2426019   14.055211 0.005269534    0.11398971
    ## 8   0.08248938          2.2382924   14.183993 0.006139867    0.12013748
    ## 9   0.08363367          2.2326335   14.304120 0.007139725    0.12649939
    ## 10  0.08486280          2.2256356   14.415353 0.008285899    0.13307468
    ## 11  0.08618005          2.2173113   14.517468 0.009596956    0.13986191
    ## 12  0.08758894          2.2076756   14.610260 0.011093358    0.14685900
    ## 13  0.08909328          2.1967461   14.693543 0.012797592    0.15406315
    ## 14  0.09069718          2.1845422   14.767146 0.014734290    0.16147080
    ## 15  0.09240506          2.1710859   14.830920 0.016930346    0.16907770
    ## 16  0.09422169          2.1564012   14.884735 0.019415033    0.17687877
    ## 17  0.09615218          2.1405142   14.928480 0.022220107    0.18486817
    ## 18  0.09820206          2.1234530   14.962066 0.025379896    0.19303927
    ## 19  0.10037722          2.1052476   14.985424 0.028931386    0.20138459
    ## 20  0.10268402          2.0859299   14.998505 0.032914277    0.20989587
    ## 21  0.10512929          2.0655334   15.001283 0.037371032    0.21856398
    ## 22  0.10772035          2.0440933   14.993751 0.042346894    0.22737901
    ## 23  0.11046506          2.0216465   14.975926 0.047889876    0.23633020
    ## 24  0.11337186          1.9982310   14.947844 0.054050730    0.24540597
    ## 25  0.11644980          1.9738864   14.909563 0.060882868    0.25459397
    ## 26  0.11970860          1.9486535   14.861161 0.068442259    0.26388102
    ## 27  0.12315871          1.9225741   14.802739 0.076787276    0.27325320
    ## 28  0.12681133          1.8956909   14.734414 0.085978503    0.28269585
    ## 29  0.13067849          1.8680477   14.656327 0.096078494    0.29219358
    ## 30  0.13477314          1.8396889   14.568637 0.107151483    0.30173034
    ## 31  0.13910915          1.8106595   14.471520 0.119263040    0.31128942
    ## 32  0.14370147          1.7810053   14.365174 0.132479677    0.32085351
    ## 33  0.14856615          1.7507720   14.249811 0.146868390    0.33040477
    ## 34  0.15372047          1.7200061   14.125661 0.162496159    0.33992482
    ## 35  0.15918299          1.6887539   13.992972 0.179429381    0.34939486
    ## 36  0.16497372          1.6570619   13.852005 0.197733252    0.35879569
    ## 37  0.17111420          1.6249765   13.703035 0.217471100    0.36810778
    ## 38  0.17762760          1.5925440   13.546353 0.238703664    0.37731132
    ## 39  0.18453894          1.5598103   13.382261 0.261488329    0.38638631
    ## 40  0.19187515          1.5268209   13.211073 0.285878324    0.39531264
    ## 41  0.19966530          1.4936210   13.033113 0.311921877    0.40407012
    ## 42  0.20794075          1.4602551   12.848715 0.339661354    0.41263859
    ## 43  0.21673535          1.4267669   12.658223 0.369132368    0.42099798
    ## 44  0.22608569          1.3931995   12.461986 0.400362888    0.42912839
    ## 45  0.23603126          1.3595951   12.260361 0.433372340    0.43701017
    ## 46  0.24661480          1.3259949   12.053710 0.468170729    0.44462401
    ## 47  0.25788251          1.2924392   11.842399 0.504757774    0.45195100
    ## 48  0.26988441          1.2589671   11.626799 0.543122091    0.45897269
    ## 49  0.28267465          1.2256166   11.407279 0.583240415    0.46567124
    ## 50  0.29631188          1.1924244   11.184214 0.625076894    0.47202940
    ## 51  0.31085970          1.1594261   10.957977 0.668582449    0.47803065
    ## 52  0.32638705          1.1266559   10.728938 0.713694237    0.48365926
    ## 53  0.34296877          1.0941465   10.497469 0.760335206    0.48890034
    ## 54  0.36068608          1.0619295   10.263936 0.808413769    0.49373992
    ## 55  0.37962721          1.0300347   10.028703 0.857823610    0.49816499
    ## 56  0.39988806          0.9984907    9.792128 0.908443623    0.50216360
    ## 57  0.42157290          0.9673245    9.554564 0.960137998    0.50572489
    ## 58  0.44479519          0.9365616    9.316357 1.012756470    0.50883912
    ## 59  0.46967845          0.9062260    9.077848 1.066134720    0.51149774
    ## 60  0.49635721          0.8763400    8.839366 1.120094950    0.51369343
    ## 61  0.52497810          0.8469246    8.601235 1.174446622    0.51542013
    ## 62  0.55570101          0.8179992    8.363769 1.228987368    0.51667304
    ## 63  0.58870040          0.7895816    8.127271 1.283504059    0.51744870
    ## 64  0.62416675          0.7616880    7.892033 1.337774046    0.51774493
    ## 65  0.66230811          0.7343332    7.658338 1.391566538    0.51756093
    ## 66  0.70335191          0.7075306    7.426457 1.444644138    0.51689719
    ## 67  0.74754685          0.6812921    7.196648 1.496764501    0.51575558
    ## 68  0.79516510          0.6556278    6.969159 1.547682103    0.51413925
    ## 69  0.84650466          0.6305470    6.744223 1.597150121    0.51205269
    ## 70  0.90189203          0.6060570    6.522063 1.644922377    0.50950166
    ## 71  0.96168511          0.5821642    6.302887 1.690755350    0.50649321
    ## 72  1.02627650          0.5588735    6.086891 1.734410218    0.50303560
    ## 73  1.09609709          0.5361885    5.874258 1.775654919    0.49913828
    ## 74  1.17162011          0.5141118    5.665158 1.814266183    0.49481186
    ## 75  1.25336556          0.4926444    5.459747 1.850031547    0.49006804
    ## 76  1.34190523          0.4717865    5.258168 1.882751293    0.48491960
    ## 77  1.43786822          0.4515371    5.060552 1.912240301    0.47938028
    ## 78  1.54194710          0.4318943    4.867016 1.938329791    0.47346475
    ## 79  1.65490479          0.4128549    4.677666 1.960868935    0.46718856
    ## 80  1.77758223          0.3944149    4.492593 1.979726305    0.46056806
    ## 81  1.91090694          0.3765696    4.311878 1.994791151    0.45362033
    ## 82  2.05590254          0.3593132    4.135588 2.005974481    0.44636307
    ## 83  2.21369946          0.3426392    3.963781 2.013209942    0.43881461
    ## 84  2.38554681          0.3265405    3.796500 2.016454472    0.43099374
    ## 85  2.57282577          0.3110090    3.633781 2.015688728    0.42291970
    ## 86  2.77706447          0.2960362    3.475646 2.010917276    0.41461206
    ## 87  2.99995468          0.2816130    3.322108 2.002168548    0.40609067
    ## 88  3.24337061          0.2677297    3.173171 1.989494562    0.39737557
    ## 89  3.50938980          0.2543762    3.028829 1.972970399    0.38848690
    ## 90  3.80031670          0.2415418    2.889066 1.952693473    0.37944485
    ## 91  4.11870912          0.2292155    2.753858 1.928782563    0.37026956
    ## 92  4.46740783          0.2173862    2.623175 1.901376665    0.36098106
    ## 93  4.84956998          0.2060420    2.496976 1.870633640    0.35159922
    ## 94  5.26870654          0.1951711    2.375216 1.836728699    0.34214362
    ## 95  5.72872450          0.1847614    2.257840 1.799852749    0.33263356
    ## 96  6.23397441          0.1748007    2.144790 1.760210597    0.32308794
    ## 97  6.78930379          0.1652764    2.036000 1.718019064    0.31352522
    ## 98  7.40011755          0.1561761    1.931401 1.673505013    0.30396340
    ## 99  8.07244592          0.1474872    1.830916 1.626903331    0.29441991
    ## 100 8.81302129          0.1391971    1.734467 1.578454876    0.28491159
    ##        Pantala
    ## 1   0.01846631
    ## 2   0.02105388
    ## 3   0.02395111
    ## 4   0.02718694
    ## 5   0.03079189
    ## 6   0.03479796
    ## 7   0.03923852
    ## 8   0.04414817
    ## 9   0.04956261
    ## 10  0.05551841
    ## 11  0.06205277
    ## 12  0.06920328
    ## 13  0.07700758
    ## 14  0.08550307
    ## 15  0.09472644
    ## 16  0.10471336
    ## 17  0.11549797
    ## 18  0.12711240
    ## 19  0.13958631
    ## 20  0.15294634
    ## 21  0.16721557
    ## 22  0.18241295
    ## 23  0.19855278
    ## 24  0.21564412
    ## 25  0.23369026
    ## 26  0.25268819
    ## 27  0.27262810
    ## 28  0.29349294
    ## 29  0.31525794
    ## 30  0.33789032
    ## 31  0.36134896
    ## 32  0.38558419
    ## 33  0.41053764
    ## 34  0.43614217
    ## 35  0.46232196
    ## 36  0.48899263
    ## 37  0.51606149
    ## 38  0.54342791
    ## 39  0.57098377
    ## 40  0.59861410
    ## 41  0.62619769
    ## 42  0.65360795
    ## 43  0.68071377
    ## 44  0.70738051
    ## 45  0.73347108
    ## 46  0.75884703
    ## 47  0.78336981
    ## 48  0.80690196
    ## 49  0.82930838
    ## 50  0.85045763
    ## 51  0.87022320
    ## 52  0.88848475
    ## 53  0.90512934
    ## 54  0.92005259
    ## 55  0.93315977
    ## 56  0.94436680
    ## 57  0.95360113
    ## 58  0.96080255
    ## 59  0.96592384
    ## 60  0.96893127
    ## 61  0.96980495
    ## 62  0.96853913
    ## 63  0.96514217
    ## 64  0.95963650
    ## 65  0.95205836
    ## 66  0.94245739
    ## 67  0.93089613
    ## 68  0.91744930
    ## 69  0.90220300
    ## 70  0.88525381
    ## 71  0.86670777
    ## 72  0.84667925
    ## 73  0.82528983
    ## 74  0.80266700
    ## 75  0.77894300
    ## 76  0.75425342
    ## 77  0.72873603
    ## 78  0.70252947
    ## 79  0.67577200
    ## 80  0.64860036
    ## 81  0.62114862
    ## 82  0.59354713
    ## 83  0.56592156
    ## 84  0.53839202
    ## 85  0.51107229
    ## 86  0.48406914
    ## 87  0.45748178
    ## 88  0.43140141
    ## 89  0.40591085
    ## 90  0.38108435
    ## 91  0.35698741
    ## 92  0.33367682
    ## 93  0.31120067
    ## 94  0.28959853
    ## 95  0.26890169
    ## 96  0.24913346
    ## 97  0.23030953
    ## 98  0.21243845
    ## 99  0.19552202
    ## 100 0.17955586

``` r
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
high_MDS1_low_MDS2_mean
```

    ##           Aedes Callibaetis Chironominae      Culex Dendropsophus.minutus
    ## 1   127.8799048 53.49990992    23.602488 141.140987            0.01189648
    ## 2   122.7854079 42.76933913    22.360115 131.029776            0.01499642
    ## 3   117.8338585 34.32108503    21.210610 121.738821            0.01882407
    ## 4   113.0244309 27.64639187    20.146293 113.195829            0.02352863
    ## 5   108.3561205 22.35449641    19.160199 105.335316            0.02928444
    ## 6   103.8277527 18.14430354    18.246004  98.097926            0.03629394
    ## 7    99.4379928 14.78307211    17.397962  91.429826            0.04479075
    ## 8    95.1853548 12.09032760    16.610850  85.282162            0.05504269
    ## 9    91.0682109  9.92568351    15.879917  79.610574            0.06735472
    ## 10   87.0848013  8.17959420    15.200836  74.374756            0.08207170
    ## 11   83.2332426  6.76631285    14.569665  69.538064            0.09958086
    ## 12   79.5115375  5.61851275    13.982813  65.067165            0.12031377
    ## 13   75.9175837  4.68316703    13.437003  60.931716            0.14474777
    ## 14   72.4491827  3.91838320    12.929244  57.104086            0.17340653
    ## 15   69.1040488  3.29096410    12.456807  53.559091            0.20685978
    ## 16   65.8798170  2.77452322    12.017198  50.273771            0.24572182
    ## 17   62.7740520  2.34802418    11.608138  47.227175            0.29064873
    ## 18   59.7842562  1.99464561    11.227545  44.400178            0.34233411
    ## 19   56.9078773  1.70089658    10.873514  41.775312            0.40150316
    ## 20   54.1423162  1.45592516    10.544303  39.336611            0.46890494
    ## 21   51.4849344  1.25097657    10.238321  37.069474            0.54530275
    ## 22   48.9330611  1.07896729     9.954110  34.960541            0.63146256
    ## 23   46.4840002  0.93414945     9.690340  32.997581            0.72813947
    ## 24   44.1350368  0.81184558     9.445795  31.169391            0.83606220
    ## 25   41.8834433  0.70823841     9.219361  29.465701            0.95591580
    ## 26   39.7264862  0.62020392     9.010026  27.877092            1.08832283
    ## 27   37.6614309  0.54517824     8.816864  26.394924            1.23382306
    ## 28   35.6855478  0.48105142     8.639032  25.011262            1.39285231
    ## 29   33.7961173  0.42608228     8.475766  23.718817            1.56572072
    ## 30   31.9904344  0.37883005     8.326369  22.510892            1.75259099
    ## 31   30.2658136  0.33809936     8.190214  21.381325            1.95345716
    ## 32   28.6195934  0.30289581     8.066733  20.324449            2.16812469
    ## 33   27.0491394  0.27239000     7.955418  19.335045            2.39619226
    ## 34   25.5518491  0.24588839     7.855814  18.408306            2.63703631
    ## 35   24.1251546  0.22280959     7.767518  17.539804            2.88979878
    ## 36   22.7665259  0.20266498     7.690175  16.725452            3.15337878
    ## 37   21.4734740  0.18504295     7.623477  15.961483            3.42642888
    ## 38   20.2435533  0.16959589     7.567157  15.244419            3.70735648
    ## 39   19.0743642  0.15602964     7.520996  14.571046            3.99433068
    ## 40   17.9635548  0.14409465     7.484810  13.938397            4.28529505
    ## 41   16.9088232  0.13357882     7.458459  13.343729            4.57798627
    ## 42   15.9079191  0.12430149     7.441839  12.784501            4.86995879
    ## 43   14.9586448  0.11610850     7.434887  12.258368            5.15861517
    ## 44   14.0588571  0.10886812     7.437574  11.763152            5.44124172
    ## 45   13.2064676  0.10246755     7.449911  11.296842            5.71504891
    ## 46   12.3994441  0.09681018     7.471947  10.857570            5.97721568
    ## 47   11.6358106  0.09181310     7.503767  10.443605            6.22493681
    ## 48   10.9136486  0.08740519     7.545496  10.053343            6.45547213
    ## 49   10.2310963  0.08352545     7.597297   9.685294            6.66619661
    ## 50    9.5863497  0.08012156     7.659374   9.338075            6.85464977
    ## 51    8.9776621  0.07714876     7.731973   9.010402            7.01858340
    ## 52    8.4033437  0.07456885     7.815383   8.701081            7.15600613
    ## 53    7.8617619  0.07234941     7.909938   8.409002            7.26522363
    ## 54    7.3513405  0.07046305     8.016019   8.133135            7.34487339
    ## 55    6.8705591  0.06888695     8.134058   7.872520            7.39395298
    ## 56    6.4179526  0.06760229     8.264540   7.626263            7.41184108
    ## 57    5.9921106  0.06659396     8.408006   7.393533            7.39831056
    ## 58    5.5916763  0.06585022     8.565055   7.173557            7.35353324
    ## 59    5.2153459  0.06536249     8.736354   6.965612            7.27807633
    ## 60    4.8618672  0.06512519     8.922635   6.769027            7.17289042
    ## 61    4.5300392  0.06513559     9.124707   6.583176            7.03928955
    ## 62    4.2187105  0.06539382     9.343457   6.407476            6.87892394
    ## 63    3.9267782  0.06590282     9.579859   6.241381            6.69374601
    ## 64    3.6531870  0.06666844     9.834981   6.084384            6.48597088
    ## 65    3.3969279  0.06769952    10.109992   5.936013            6.25803235
    ## 66    3.1570367  0.06900806    10.406171   5.795825            6.01253554
    ## 67    2.9325933  0.07060949    10.724919   5.663410            5.75220763
    ## 68    2.7227196  0.07252292    11.067765   5.538382            5.47984778
    ## 69    2.5265791  0.07477157    11.436383   5.420385            5.19827773
    ## 70    2.3433748  0.07738320    11.832604   5.309083            4.91029414
    ## 71    2.1723486  0.08039070    12.258430   5.204167            4.61862371
    ## 72    2.0127793  0.08383280    12.716051   5.105345            4.32588228
    ## 73    1.8639819  0.08775484    13.207863   5.012349            4.03453843
    ## 74    1.7253059  0.09220983    13.736487   4.924926            3.74688238
    ## 75    1.5961342  0.09725956    14.304797   4.842842            3.46500056
    ## 76    1.4758819  0.10297608    14.915939   4.765881            3.19075604
    ## 77    1.3639947  0.10944336    15.573362   4.693841            2.92577508
    ## 78    1.2599481  0.11675929    16.280848   4.626533            2.67143933
    ## 79    1.1632459  0.12503812    17.042548   4.563787            2.42888380
    ## 80    1.0734190  0.13441336    17.863021   4.505440            2.19899988
    ## 81    0.9900244  0.14504121    18.747276   4.451345            1.98244303
    ## 82    0.9126441  0.15710476    19.700821   4.401368            1.77964456
    ## 83    0.8408835  0.17081904    20.729715   4.355382            1.59082670
    ## 84    0.7743711  0.18643704    21.840632   4.313274            1.41602041
    ## 85    0.7127567  0.20425707    23.040928   4.274941            1.25508517
    ## 86    0.6557109  0.22463166    24.338712   4.240289            1.10773003
    ## 87    0.6029237  0.24797838    25.742937   4.209234            0.97353530
    ## 88    0.5541039  0.27479299    27.263492   4.181700            0.85197435
    ## 89    0.5089780  0.30566553    28.911307   4.157621            0.74243491
    ## 90    0.4672891  0.34129998    30.698479   4.136940            0.64423941
    ## 91    0.4287964  0.38253841    32.638399   4.119607            0.55666413
    ## 92    0.3932743  0.43039062    34.745913   4.105580            0.47895668
    ## 93    0.3605113  0.48607081    37.037483   4.094827            0.41035173
    ## 94    0.3303095  0.55104271    39.531390   4.087322            0.35008488
    ## 95    0.3024838  0.62707568    42.247944   4.083047            0.29740446
    ## 96    0.2768612  0.71631434    45.209733   4.081992            0.25158147
    ## 97    0.2532800  0.82136521    48.441902   4.084154            0.21191752
    ## 98    0.2315894  0.94540512    51.972463   4.089539            0.17775102
    ## 99    0.2116485  1.09231670    55.832656   4.098160            0.14846169
    ## 100   0.1933262  1.26685869    60.057348   4.110036            0.12347350
    ##       Microvelia Scinax.fuscovarius Tanypodinae    Dasyhelea Erythrodiplax
    ## 1   2.337724e-02          61.327303   0.1766115  0.002984816    0.03668307
    ## 2   1.158347e-02          62.046680   0.2002408  0.003238410    0.03989417
    ## 3   5.826041e-03          62.707816   0.2266121  0.003513922    0.04333495
    ## 4   2.974390e-03          63.308678   0.2559830  0.003813276    0.04701671
    ## 5   1.541387e-03          63.847406   0.2886266  0.004138569    0.05095083
    ## 6   8.108025e-04          64.322323   0.3248319  0.004492086    0.05514872
    ## 7   4.329202e-04          64.731941   0.3649038  0.004876316    0.05962175
    ## 8   2.346336e-04          65.074971   0.4091620  0.005293971    0.06438120
    ## 9   1.290809e-04          65.350329   0.4579409  0.005748004    0.06943822
    ## 10  7.208139e-05          65.557143   0.5115885  0.006241637    0.07480371
    ## 11  4.085770e-05          65.694757   0.5704655  0.006778378    0.08048830
    ## 12  2.350792e-05          65.762732   0.6349438  0.007362054    0.08650228
    ## 13  1.372916e-05          65.760851   0.7054047  0.007996834    0.09285545
    ## 14  8.138849e-06          65.689121   0.7822376  0.008687264    0.09955713
    ## 15  4.897468e-06          65.547770   0.8658372  0.009438302    0.10661601
    ## 16  2.991368e-06          65.337248   0.9566015  0.010255353    0.11404010
    ## 17  1.854631e-06          65.058223   1.0549286  0.011144311    0.12183662
    ## 18  1.167171e-06          64.711580   1.1612143  0.012111605    0.13001192
    ## 19  7.455918e-07          64.298414   1.2758478  0.013164249    0.13857141
    ## 20  4.834562e-07          63.820023   1.3992092  0.014309893    0.14751942
    ## 21  3.182019e-07          63.277906   1.5316645  0.015556881    0.15685915
    ## 22  2.125875e-07          62.673751   1.6735623  0.016914321    0.16659256
    ## 23  1.441658e-07          62.009428   1.8252289  0.018392150    0.17672030
    ## 24  9.923760e-08          61.286978   1.9869643  0.020001211    0.18724162
    ## 25  6.933935e-08          60.508604   2.1590366  0.021753342    0.19815428
    ## 26  4.917821e-08          59.676659   2.3416781  0.023661462    0.20945445
    ## 27  3.540423e-08          58.793636   2.5350797  0.025739674    0.22113671
    ## 28  2.587183e-08          57.862152   2.7393864  0.028003376    0.23319391
    ## 29  1.919060e-08          56.884938   2.9546919  0.030469380    0.24561712
    ## 30  1.444905e-08          55.864825   3.1810342  0.033156046    0.25839564
    ## 31  1.104281e-08          54.804730   3.4183909  0.036083424    0.27151688
    ## 32  8.566622e-09          53.707642   3.6666743  0.039273412    0.28496635
    ## 33  6.745730e-09          52.576610   3.9257279  0.042749929    0.29872765
    ## 34  5.391849e-09          51.414724   4.1953219  0.046539108    0.31278244
    ## 35  4.374577e-09          50.225110   4.4751503  0.050669498    0.32711044
    ## 36  3.602665e-09          49.010905   4.7648277  0.055172292    0.34168943
    ## 37  3.011626e-09          47.775253   5.0638870  0.060081579    0.35649529
    ## 38  2.555453e-09          46.521286   5.3717780  0.065434613    0.37150198
    ## 39  2.201020e-09          45.252114   5.6878659  0.071272110    0.38668166
    ## 40  1.924287e-09          43.970812   6.0114311  0.077638580    0.40200468
    ## 41  1.707674e-09          42.680405   6.3416700  0.084582679    0.41743970
    ## 42  1.538260e-09          41.383864   6.6776960  0.092157603    0.43295373
    ## 43  1.406513e-09          40.084086   7.0185417  0.100421520    0.44851226
    ## 44  1.305411e-09          38.783890   7.3631622  0.109438037    0.46407936
    ## 45  1.229817e-09          37.486009   7.7104386  0.119276719    0.47961780
    ## 46  1.176043e-09          36.193075   8.0591835  0.130013653    0.49508917
    ## 47  1.141551e-09          34.907618   8.4081459  0.141732066    0.51045405
    ## 48  1.124752e-09          33.632054   8.7560186  0.154523011    0.52567217
    ## 49  1.124884e-09          32.368681   9.1014449  0.168486105    0.54070255
    ## 50  1.141953e-09          31.119676   9.4430276  0.183730351    0.55550369
    ## 51  1.176734e-09          29.887086   9.7793369  0.200375032    0.57003377
    ## 52  1.230829e-09          28.672828  10.1089208  0.218550693    0.58425079
    ## 53  1.306793e-09          27.478684  10.4303149  0.238400222    0.59811286
    ## 54  1.408333e-09          26.306301  10.7420529  0.260080031    0.61157830
    ## 55  1.540612e-09          25.157187  11.0426775  0.283761353    0.62460592
    ## 56  1.710687e-09          24.032713  11.3307519  0.309631665    0.63715518
    ## 57  1.928136e-09          22.934115  11.6048709  0.337896252    0.64918645
    ## 58  2.205941e-09          21.862490  11.8636720  0.368779921    0.66066115
    ## 59  2.561768e-09          20.818800  12.1058472  0.402528879    0.67154202
    ## 60  3.019779e-09          19.803876  12.3301535  0.439412802    0.68179328
    ## 61  3.613267e-09          18.818421  12.5354237  0.479727099    0.69138084
    ## 62  4.388482e-09          17.863008  12.7205764  0.523795403    0.70027248
    ## 63  5.410260e-09          16.938090  12.8846255  0.571972302    0.70843805
    ## 64  6.770354e-09          16.044003  13.0266891  0.624646339    0.71584962
    ## 65  8.599915e-09          15.180969  13.1459973  0.682243307    0.72248165
    ## 66  1.108834e-08          14.349100  13.2418988  0.745229871    0.72831112
    ## 67  1.451203e-08          13.548409  13.3138674  0.814117535    0.73331768
    ## 68  1.927877e-08          12.778809  13.3615066  0.889467018    0.73748378
    ## 69  2.599681e-08          12.040123  13.3845530  0.971893042    0.74079474
    ## 70  3.558362e-08          11.332087  13.3828791  1.062069611    0.74323885
    ## 71  4.943902e-08          10.654359  13.3564939  1.160735792    0.74480746
    ## 72  6.972346e-08          10.006523  13.3055438  1.268702089    0.74549500
    ## 73  9.981080e-08           9.388095  13.2303101  1.386857433    0.74529901
    ## 74  1.450326e-07           8.798533  13.1312073  1.516176870    0.74422021
    ## 75  2.139160e-07           8.237235  13.0087788  1.657730012    0.74226241
    ## 76  3.202656e-07           7.703553  12.8636922  1.812690330    0.73943258
    ## 77  4.867062e-07           7.196796  12.6967331  1.982345369    0.73574071
    ## 78  7.507803e-07           6.716233  12.5087979  2.168107980    0.73119986
    ## 79  1.175570e-06           6.261101  12.3008863  2.371528674    0.72582598
    ## 80  1.868415e-06           5.830612  12.0740916  2.594309208    0.71963788
    ## 81  3.014309e-06           5.423954  11.8295921  2.838317525    0.71265711
    ## 82  4.936189e-06           5.040299  11.5686401  3.105604199    0.70490782
    ## 83  8.205125e-06           4.678807  11.2925518  3.398420513    0.69641663
    ## 84  1.384421e-05           4.338627  11.0026959  3.719238361    0.68721248
    ## 85  2.371049e-05           4.018907  10.7004826  4.070772139    0.67732646
    ## 86  4.121947e-05           3.718794  10.3873523  4.456002841    0.66679164
    ## 87  7.273673e-05           3.437436  10.0647640  4.878204578    0.65564288
    ## 88  1.302851e-04           3.173990   9.7341843  5.340973763    0.64391664
    ## 89  2.368781e-04           2.927622   9.3970768  5.848261244    0.63165080
    ## 90  4.371644e-04           2.697509   9.0548910  6.404407674    0.61888442
    ## 91  8.189437e-04           2.482843   8.7090526  7.014182455    0.60565757
    ## 92  1.557230e-03           2.282832   8.3609541  7.682826628    0.59201109
    ## 93  3.005668e-03           2.096705   8.0119458  8.416100094    0.57798643
    ## 94  5.888691e-03           1.923707   7.6633278  9.220333627    0.56362538
    ## 95  1.171079e-02           1.763108   7.3163428 10.102486165    0.54896990
    ## 96  2.363974e-02           1.614201   6.9721693 11.070207917    0.53406194
    ## 97  4.843830e-02           1.476300   6.6319164 12.131909889    0.51894320
    ## 98  1.007452e-01           1.348746   6.2966188 13.296840485    0.50365496
    ## 99  2.126913e-01           1.230903   5.9672330 14.575169917    0.48823793
    ## 100 4.557895e-01           1.122164   5.6446345 15.978083224    0.47273200
    ##        Pantala
    ## 1   0.05402980
    ## 2   0.05875851
    ## 3   0.06384320
    ## 4   0.06930506
    ## 5   0.07516604
    ## 6   0.08144883
    ## 7   0.08817683
    ## 8   0.09537411
    ## 9   0.10306542
    ## 10  0.11127609
    ## 11  0.12003205
    ## 12  0.12935970
    ## 13  0.13928591
    ## 14  0.14983795
    ## 15  0.16104339
    ## 16  0.17293003
    ## 17  0.18552581
    ## 18  0.19885876
    ## 19  0.21295681
    ## 20  0.22784776
    ## 21  0.24355914
    ## 22  0.26011808
    ## 23  0.27755118
    ## 24  0.29588439
    ## 25  0.31514286
    ## 26  0.33535077
    ## 27  0.35653123
    ## 28  0.37870609
    ## 29  0.40189577
    ## 30  0.42611910
    ## 31  0.45139320
    ## 32  0.47773323
    ## 33  0.50515230
    ## 34  0.53366121
    ## 35  0.56326839
    ## 36  0.59397964
    ## 37  0.62579799
    ## 38  0.65872358
    ## 39  0.69275343
    ## 40  0.72788136
    ## 41  0.76409778
    ## 42  0.80138963
    ## 43  0.83974017
    ## 44  0.87912893
    ## 45  0.91953158
    ## 46  0.96091984
    ## 47  1.00326140
    ## 48  1.04651988
    ## 49  1.09065475
    ## 50  1.13562132
    ## 51  1.18137076
    ## 52  1.22785004
    ## 53  1.27500203
    ## 54  1.32276549
    ## 55  1.37107517
    ## 56  1.41986192
    ## 57  1.46905276
    ## 58  1.51857101
    ## 59  1.56833651
    ## 60  1.61826571
    ## 61  1.66827194
    ## 62  1.71826559
    ## 63  1.76815435
    ## 64  1.81784348
    ## 65  1.86723610
    ## 66  1.91623345
    ## 67  1.96473522
    ## 68  2.01263991
    ## 69  2.05984509
    ## 70  2.10624786
    ## 71  2.15174511
    ## 72  2.19623398
    ## 73  2.23961219
    ## 74  2.28177844
    ## 75  2.32263280
    ## 76  2.36207711
    ## 77  2.40001536
    ## 78  2.43635407
    ## 79  2.47100271
    ## 80  2.50387400
    ## 81  2.53488438
    ## 82  2.56395425
    ## 83  2.59100841
    ## 84  2.61597631
    ## 85  2.63879240
    ## 86  2.65939640
    ## 87  2.67773356
    ## 88  2.69375491
    ## 89  2.70741750
    ## 90  2.71868454
    ## 91  2.72752560
    ## 92  2.73391676
    ## 93  2.73784069
    ## 94  2.73928673
    ## 95  2.73825094
    ## 96  2.73473616
    ## 97  2.72875191
    ## 98  2.72031443
    ## 99  2.70944658
    ## 100 2.69617770

``` r
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
low_MDS1_low_MDS2_mean
```

    ##           Aedes Callibaetis Chironominae    Culex Dendropsophus.minutus
    ## 1   182.1132344   1.4552268    85.846779 26.20575           0.002692134
    ## 2   172.5776420   1.5903620    78.351109 25.61673           0.003287454
    ## 3   163.4681874   1.7341909    71.642734 25.04862           0.004000775
    ## 4   154.7703116   1.8868329    65.630393 24.50062           0.004852324
    ## 5   146.4696893   2.0483566    60.234279 23.97197           0.005865121
    ## 6   138.5522418   2.2187752    55.384505 23.46191           0.007065216
    ## 7   131.0041481   2.3980412    51.019794 22.96974           0.008481941
    ## 8   123.8118551   2.5860421    47.086345 22.49480           0.010148141
    ## 9   116.9620867   2.7825959    43.536861 22.03643           0.012100382
    ## 10  110.4418514   2.9874476    40.329710 21.59402           0.014379143
    ## 11  104.2384493   3.2002659    37.428200 21.16698           0.017028968
    ## 12   98.3394782   3.4206405    34.799952 20.75474           0.020098564
    ## 13   92.7328378   3.6480805    32.416356 20.35678           0.023640850
    ## 14   87.4067343   3.8820131    30.252105 19.97257           0.027712934
    ## 15   82.3496830   4.1217836    28.284784 19.60162           0.032376009
    ## 16   77.5505109   4.3666560    26.494516 19.24346           0.037695150
    ## 17   72.9983573   4.6158149    24.863655 18.89764           0.043739015
    ## 18   68.6826756   4.8683679    23.376517 18.56372           0.050579427
    ## 19   64.5932322   5.1233498    22.019147 18.24129           0.058290821
    ## 20   60.7201061   5.3797269    20.779114 17.92996           0.066949573
    ## 21   57.0536876   5.6364032    19.645335 17.62935           0.076633175
    ## 22   53.5846763   5.8922272    18.607914 17.33910           0.087419274
    ## 23   50.3040783   6.1459996    17.658012 17.05885           0.099384564
    ## 24   47.2032036   6.3964818    16.787722 16.78828           0.112603542
    ## 25   44.2736622   6.6424061    15.989968 16.52706           0.127147121
    ## 26   41.5073600   6.8824851    15.258409 16.27490           0.143081137
    ## 27   38.8964952   7.1154233    14.587363 16.03150           0.160464734
    ## 28   36.4335528   7.3399282    13.971730 15.79658           0.179348681
    ## 29   34.1112998   7.5547220    13.406932 15.56988           0.199773622
    ## 30   31.9227802   7.7585537    12.888860 15.35114           0.221768300
    ## 31   29.8613089   7.9502110    12.413821 15.14011           0.245347796
    ## 32   27.9204662   8.1285326    11.978496 14.93657           0.270511805
    ## 33   26.0940918   8.2924192    11.579904 14.74027           0.297243011
    ## 34   24.3762782   8.4408456    11.215366 14.55102           0.325505590
    ## 35   22.7613651   8.5728706    10.882479 14.36860           0.355243891
    ## 36   21.2439321   8.6876474    10.579084 14.19282           0.386381349
    ## 37   19.8187928   8.7844325    10.303248 14.02349           0.418819659
    ## 38   18.4809879   8.8625937    10.053241 13.86043           0.452438269
    ## 39   17.2257787   8.9216170     9.827519 13.70346           0.487094217
    ## 40   16.0486401   8.9611123     9.624708 13.55243           0.522622356
    ## 41   14.9452541   8.9808175     9.443589 13.40717           0.558835984
    ## 42   13.9115033   8.9806015     9.283087 13.26754           0.595527897
    ## 43   12.9434639   8.9604658     9.142262 13.13338           0.632471887
    ## 44   12.0373992   8.9205442     9.020295 13.00457           0.669424667
    ## 45   11.1897531   8.8611018     8.916485 12.88097           0.706128214
    ## 46   10.3971435   8.7825313     8.830239 12.76246           0.742312522
    ## 47    9.6563560   8.6853493     8.761069 12.64892           0.777698707
    ## 48    8.9643373   8.5701906     8.708585 12.54023           0.812002444
    ## 49    8.3181893   8.4378010     8.672493 12.43628           0.844937653
    ## 50    7.7151631   8.2890295     8.652591 12.33698           0.876220400
    ## 51    7.1526522   8.1248190     8.648767 12.24223           0.905572916
    ## 52    6.6281878   7.9461967     8.661002 12.15193           0.932727666
    ## 53    6.1394320   7.7542631     8.689362 12.06599           0.957431388
    ## 54    5.6841728   7.5501809     8.734007 11.98433           0.979449011
    ## 55    5.2603186   7.3351634     8.795186 11.90688           0.998567372
    ## 56    4.8658926   7.1104622     8.873243 11.83356           1.014598648
    ## 57    4.4990278   6.8773556     8.968619 11.76429           1.027383427
    ## 58    4.1579622   6.6371363     9.081856 11.69902           1.036793345
    ## 59    3.8410334   6.3910997     9.213604 11.63767           1.042733232
    ## 60    3.5466744   6.1405328     9.364623 11.58020           1.045142714
    ## 61    3.2734089   5.8867030     9.535795 11.52654           1.043997234
    ## 62    3.0198466   5.6308480     9.728131 11.47665           1.039308473
    ## 63    2.7846795   5.3741661     9.942778 11.43048           1.031124144
    ## 64    2.5666770   5.1178079    10.181035 11.38799           1.019527191
    ## 65    2.3646830   4.8628679    10.444363 11.34913           1.004634385
    ## 66    2.1776112   4.6103784    10.734402 11.31388           0.986594369
    ## 67    2.0044417   4.3613032    11.052986 11.28219           0.965585188
    ## 68    1.8442179   4.1165330    11.402162 11.25404           0.941811363
    ## 69    1.6960424   3.8768814    11.784216 11.22941           0.915500585
    ## 70    1.5590746   3.6430828    12.201690 11.20826           0.886900081
    ## 71    1.4325268   3.4157901    12.657420 11.19059           0.856272766
    ## 72    1.3156620   3.1955743    13.154556 11.17636           0.823893235
    ## 73    1.2077904   2.9829245    13.696610 11.16558           0.790043695
    ## 74    1.1082673   2.7782493    14.287487 11.15823           0.755009914
    ## 75    1.0164901   2.5818783    14.931534 11.15431           0.719077278
    ## 76    0.9318961   2.3940649    15.633596 11.15380           0.682527015
    ## 77    0.8539600   2.2149896    16.399069 11.15672           0.645632664
    ## 78    0.7821918   2.0447634    17.233971 11.16306           0.608656847
    ## 79    0.7161346   1.8834325    18.145017 11.17282           0.571848383
    ## 80    0.6553628   1.7309824    19.139706 11.18603           0.535439803
    ## 81    0.5994798   1.5873432    20.226418 11.20268           0.499645264
    ## 82    0.5481168   1.4523946    21.414531 11.22280           0.464658907
    ## 83    0.5009303   1.3259710    22.714543 11.24640           0.430653645
    ## 84    0.4576012   1.2078668    24.138224 11.27351           0.397780378
    ## 85    0.4178330   1.0978416    25.698777 11.30415           0.366167634
    ## 86    0.3813502   0.9956252    27.411037 11.33835           0.335921592
    ## 87    0.3478972   0.9009231    29.291683 11.37613           0.307126474
    ## 88    0.3172369   0.8134206    31.359494 11.41755           0.279845261
    ## 89    0.2891492   0.7327877    33.635635 11.46263           0.254120692
    ## 90    0.2634305   0.6586836    36.143987 11.51141           0.229976505
    ## 91    0.2398920   0.5907600    38.911533 11.56395           0.207418877
    ## 92    0.2183591   0.5286655    41.968794 11.62029           0.186438009
    ## 93    0.1986700   0.4720483    45.350333 11.68049           0.167009818
    ## 94    0.1806754   0.4205596    49.095347 11.74459           0.149097691
    ## 95    0.1642372   0.3738559    53.248336 11.81268           0.132654260
    ## 96    0.1492278   0.3316015    57.859890 11.88479           0.117623163
    ## 97    0.1355295   0.2934704    62.987595 11.96102           0.103940757
    ## 98    0.1230335   0.2591480    68.697084 12.04142           0.091537756
    ## 99    0.1116397   0.2283321    75.063261 12.12609           0.080340776
    ## 100   0.1012557   0.2007343    82.171727 12.21509           0.070273753
    ##       Microvelia Scinax.fuscovarius Tanypodinae  Dasyhelea Erythrodiplax
    ## 1   1.849380e-03          0.6613757    1.815079 0.01395237    0.07150570
    ## 2   1.312868e-03          0.7301364    1.957238 0.01520067    0.07665587
    ## 3   9.401231e-04          0.8043964    2.107865 0.01655341    0.08208629
    ## 4   6.790741e-04          0.8843956    2.267215 0.01801865    0.08780439
    ## 5   4.947870e-04          0.9703611    2.435530 0.01960499    0.09381717
    ## 6   3.636539e-04          1.0625039    2.613034 0.02132166    0.10013105
    ## 7   2.696043e-04          1.1610156    2.799932 0.02317848    0.10675192
    ## 8   2.016203e-04          1.2660646    2.996406 0.02518599    0.11368496
    ## 9   1.520934e-04          1.3777933    3.202615 0.02735538    0.12093465
    ## 10  1.157325e-04          1.4963135    3.418690 0.02969863    0.12850468
    ## 11  8.883191e-05          1.6217035    3.644730 0.03222848    0.13639786
    ## 12  6.877829e-05          1.7540044    3.880806 0.03495852    0.14461608
    ## 13  5.371586e-05          1.8932162    4.126951 0.03790323    0.15316023
    ## 14  4.231774e-05          2.0392953    4.383163 0.04107799    0.16203017
    ## 15  3.362878e-05          2.1921505    4.649397 0.04449917    0.17122461
    ## 16  2.695681e-05          2.3516406    4.925570 0.04818419    0.18074110
    ## 17  2.179690e-05          2.5175719    5.211554 0.05215152    0.19057594
    ## 18  1.777827e-05          2.6896957    5.507174 0.05642081    0.20072416
    ## 19  1.462693e-05          2.8677069    5.812208 0.06101286    0.21117944
    ## 20  1.213908e-05          3.0512424    6.126386 0.06594978    0.22193411
    ## 21  1.016218e-05          3.2398807    6.449386 0.07125497    0.23297907
    ## 22  8.581373e-06          3.4331412    6.780836 0.07695322    0.24430377
    ## 23  7.309631e-06          3.6304850    7.120311 0.08307077    0.25589622
    ## 24  6.280625e-06          3.8313161    7.467333 0.08963539    0.26774292
    ## 25  5.443510e-06          4.0349825    7.821370 0.09667644    0.27982888
    ## 26  4.759090e-06          4.2407793    8.181841 0.10422493    0.29213763
    ## 27  4.196987e-06          4.4479514    8.548109 0.11231362    0.30465120
    ## 28  3.733533e-06          4.6556973    8.919487 0.12097706    0.31735015
    ## 29  3.350204e-06          4.8631737    9.295239 0.13025173    0.33021359
    ## 30  3.032432e-06          5.0695005    9.674579 0.14017605    0.34321923
    ## 31  2.768725e-06          5.2737666   10.056675 0.15079048    0.35634339
    ## 32  2.549983e-06          5.4750360   10.440651 0.16213766    0.36956107
    ## 33  2.368991e-06          5.6723550   10.825590 0.17426241    0.38284604
    ## 34  2.220027e-06          5.8647589   11.210537 0.18721185    0.39617085
    ## 35  2.098563e-06          6.0512803   11.594501 0.20103553    0.40950698
    ## 36  2.001034e-06          6.2309565   11.976463 0.21578543    0.42282486
    ## 37  1.924667e-06          6.4028380   12.355373 0.23151613    0.43609404
    ## 38  1.867349e-06          6.5659965   12.730165 0.24828486    0.44928324
    ## 39  1.827529e-06          6.7195335   13.099750 0.26615157    0.46236047
    ## 40  1.804146e-06          6.8625883   13.463030 0.28517909    0.47529321
    ## 41  1.796585e-06          6.9943460   13.818899 0.30543314    0.48804845
    ## 42  1.804649e-06          7.1140452   14.166251 0.32698247    0.50059292
    ## 43  1.828548e-06          7.2209855   14.503981 0.34989893    0.51289313
    ## 44  1.868912e-06          7.3145341   14.830997 0.37425757    0.52491561
    ## 45  1.926815e-06          7.3941322   15.146222 0.40013672    0.53662700
    ## 46  2.003826e-06          7.4593003   15.448599 0.42761808    0.54799422
    ## 47  2.102077e-06          7.5096434   15.737101 0.45678679    0.55898462
    ## 48  2.224365e-06          7.5448547   16.010732 0.48773155    0.56956614
    ## 49  2.374282e-06          7.5647189   16.268535 0.52054467    0.57970746
    ## 50  2.556391e-06          7.5691140   16.509601 0.55532214    0.58937817
    ## 51  2.776457e-06          7.5580132   16.733065 0.59216374    0.59854888
    ## 52  3.041749e-06          7.5314845   16.938123 0.63117308    0.60719142
    ## 53  3.361434e-06          7.4896905   17.124026 0.67245767    0.61527894
    ## 54  3.747093e-06          7.4328865   17.290092 0.71612903    0.62278609
    ## 55  4.213405e-06          7.3614178   17.435708 0.76230266    0.62968910
    ## 56  4.779039e-06          7.2757168   17.560330 0.81109818    0.63596597
    ## 57  5.467853e-06          7.1762976   17.663493 0.86263933    0.64159652
    ## 58  6.310471e-06          7.0637520   17.744810 0.91705401    0.64656256
    ## 59  7.346415e-06          6.9387427   17.803973 0.97447434    0.65084792
    ## 60  8.626962e-06          6.8019976   17.840760 1.03503667    0.65443862
    ## 61  1.021902e-05          6.6543019   17.855030 1.09888160    0.65732287
    ## 62  1.221038e-05          6.4964915   17.846730 1.16615398    0.65949118
    ## 63  1.471695e-05          6.3294444   17.815892 1.23700294    0.66093639
    ## 64  1.789267e-05          6.1540731   17.762631 1.31158185    0.66165372
    ## 65  2.194327e-05          5.9713160   17.687150 1.39004833    0.66164080
    ## 66  2.714541e-05          5.7821293   17.589733 1.47256419    0.66089767
    ## 67  3.387350e-05          5.5874788   17.470747 1.55929543    0.65942679
    ## 68  4.263757e-05          5.3883317   17.330637 1.65041214    0.65723303
    ## 69  5.413694e-05          5.1856488   17.169926 1.74608848    0.65432362
    ## 70  6.933677e-05          4.9803770   16.989208 1.84650259    0.65070814
    ## 71  8.957820e-05          4.7734423   16.789149 1.95183647    0.64639845
    ## 72  1.167373e-04          4.5657433   16.570479 2.06227591    0.64140863
    ## 73  1.534567e-04          4.3581446   16.333989 2.17801038    0.63575490
    ## 74  2.034843e-04          4.1514721   16.080528 2.29923285    0.62945553
    ## 75  2.721728e-04          3.9465078   15.810994 2.42613967    0.62253076
    ## 76  3.672207e-04          3.7439854   15.526332 2.55893043    0.61500267
    ## 77  4.997794e-04          3.5445872   15.227527 2.69780771    0.60689508
    ## 78  6.861172e-04          3.3489412   14.915601 2.84297696    0.59823340
    ## 79  9.501387e-04          3.1576191   14.591601 2.99464621    0.58904451
    ## 80  1.327225e-03          2.9711343   14.256600 3.15302590    0.57935666
    ## 81  1.870125e-03          2.7899420   13.911687 3.31832859    0.56919924
    ## 82  2.658065e-03          2.6144384   13.557964 3.49076871    0.55860272
    ## 83  3.810916e-03          2.4449614   13.196537 3.67056225    0.54759844
    ## 84  5.511398e-03          2.2817913   12.828513 3.85792648    0.53621849
    ## 85  8.040129e-03          2.1251530   12.454993 4.05307962    0.52449554
    ## 86  1.183131e-02          1.9752170   12.077067 4.25624048    0.51246266
    ## 87  1.756191e-02          1.8321025   11.695810 4.46762811    0.50015324
    ## 88  2.629536e-02          1.6958798   11.312276 4.68746141    0.48760075
    ## 89  3.971507e-02          1.5665732   10.927491 4.91595876    0.47483866
    ## 90  6.050626e-02          1.4441646   10.542456 5.15333756    0.46190026
    ## 91  9.298525e-02          1.3285962   10.158134 5.39981381    0.44881851
    ## 92  1.441440e-01          1.2197749    9.775454 5.65560167    0.43562594
    ## 93  2.253968e-01          1.1175750    9.395302 5.92091296    0.42235451
    ## 94  3.555230e-01          1.0218427    9.018522 6.19595670    0.40903547
    ## 95  5.656613e-01          0.9323989    8.645913 6.48093858    0.39569924
    ## 96  9.078496e-01          0.8490432    8.278223 6.77606045    0.38237536
    ## 97  1.469739e+00          0.7715573    7.916155 7.08151977    0.36909231
    ## 98  2.400132e+00          0.6997081    7.560356 7.39750908    0.35587749
    ## 99  3.953656e+00          0.6332511    7.211424 7.72421539    0.34275711
    ## 100 6.569485e+00          0.5719333    6.869904 8.06181966    0.32975611
    ##        Pantala
    ## 1   0.05196021
    ## 2   0.05836259
    ## 3   0.06542405
    ## 4   0.07319468
    ## 5   0.08172610
    ## 6   0.09107125
    ## 7   0.10128405
    ## 8   0.11241907
    ## 9   0.12453120
    ## 10  0.13767515
    ## 11  0.15190503
    ## 12  0.16727382
    ## 13  0.18383281
    ## 14  0.20163098
    ## 15  0.22071443
    ## 16  0.24112565
    ## 17  0.26290286
    ## 18  0.28607929
    ## 19  0.31068248
    ## 20  0.33673348
    ## 21  0.36424623
    ## 22  0.39322674
    ## 23  0.42367246
    ## 24  0.45557161
    ## 25  0.48890252
    ## 26  0.52363313
    ## 27  0.55972046
    ## 28  0.59711015
    ## 29  0.63573621
    ## 30  0.67552069
    ## 31  0.71637361
    ## 32  0.75819290
    ## 33  0.80086455
    ## 34  0.84426278
    ## 35  0.88825044
    ## 36  0.93267951
    ## 37  0.97739172
    ## 38  1.02221933
    ## 39  1.06698605
    ## 40  1.11150804
    ## 41  1.15559511
    ## 42  1.19905194
    ## 43  1.24167950
    ## 44  1.28327651
    ## 45  1.32364094
    ## 46  1.36257167
    ## 47  1.39987010
    ## 48  1.43534180
    ## 49  1.46879825
    ## 50  1.50005843
    ## 51  1.52895048
    ## 52  1.55531328
    ## 53  1.57899792
    ## 54  1.59986911
    ## 55  1.61780644
    ## 56  1.63270562
    ## 57  1.64447936
    ## 58  1.65305836
    ## 59  1.65839187
    ## 60  1.66044826
    ## 61  1.65921534
    ## 62  1.65470042
    ## 63  1.64693028
    ## 64  1.63595091
    ## 65  1.62182704
    ## 66  1.60464152
    ## 67  1.58449447
    ## 68  1.56150236
    ## 69  1.53579687
    ## 70  1.50752363
    ## 71  1.47684085
    ## 72  1.44391783
    ## 73  1.40893344
    ## 74  1.37207450
    ## 75  1.33353409
    ## 76  1.29350994
    ## 77  1.25220269
    ## 78  1.20981430
    ## 79  1.16654637
    ## 80  1.12259866
    ## 81  1.07816752
    ## 82  1.03344457
    ## 83  0.98861535
    ## 84  0.94385814
    ## 85  0.89934291
    ## 86  0.85523039
    ## 87  0.81167122
    ## 88  0.76880533
    ## 89  0.72676138
    ## 90  0.68565637
    ## 91  0.64559536
    ## 92  0.60667138
    ## 93  0.56896536
    ## 94  0.53254629
    ## 95  0.49747140
    ## 96  0.46378648
    ## 97  0.43152630
    ## 98  0.40071506
    ## 99  0.37136698
    ## 100 0.34348685

### Ploting predicted effects

``` r
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
          col= "#56B4E9", lwd = 2)


lines(y = low_MDS1_high_MDS2_mean[,i],                                                                                                                    x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= "white", lwd = 3)
lines(y = low_MDS1_high_MDS2_mean[,i],
          x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= "#009E73", lwd = 2)


lines(y = high_MDS1_low_MDS2_mean[,i],                                                                                                                         x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= "white", lwd = 3)
lines(y = high_MDS1_low_MDS2_mean[,i],
          x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= "#CC79A7", lwd = 2)


lines(y = low_MDS1_low_MDS2_mean[,i],                                                                                                                         x = seq(from = min(Exp_design_all_drop_atrasado$AM_days), to = max(Exp_design_all_drop_atrasado$AM_days), length.out = 100),
          col= "white", lwd = 3)
lines(y = low_MDS1_low_MDS2_mean[,i],
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
legend(legend = c("High nMDS1, high nMDS2",
                 "High nMDS1, low nMDS2",
                  "Low nMDS1, high nMDS2",
                  "Low nMDS1, low nMDS2"), x = 50, y = 50, bty = "n", lty = 1, lwd = 3, col = c("#56B4E9", "#CC79A7", "#009E73","#E69F00"), xjust = 0.5, yjust = 0.5)
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

set.seed(5); nmds_comm_trajectories <- metaMDS(comm_trajectories_relative, distance = "bray", k = 2, trymax  = 100, try = 100, engine = "monoMDS", autotransform  = FALSE)
```

    ## Run 0 stress 0.1592978 
    ## Run 1 stress 0.1592978 
    ## ... New best solution
    ## ... Procrustes: rmse 7.551921e-06  max resid 5.292302e-05 
    ## ... Similar to previous best
    ## Run 2 stress 0.159298 
    ## ... Procrustes: rmse 0.0001800255  max resid 0.002468571 
    ## ... Similar to previous best
    ## Run 3 stress 0.159298 
    ## ... Procrustes: rmse 0.0001942133  max resid 0.002609492 
    ## ... Similar to previous best
    ## Run 4 stress 0.159298 
    ## ... Procrustes: rmse 0.000192946  max resid 0.002607318 
    ## ... Similar to previous best
    ## Run 5 stress 0.1592979 
    ## ... Procrustes: rmse 0.000176437  max resid 0.002553483 
    ## ... Similar to previous best
    ## Run 6 stress 0.173302 
    ## Run 7 stress 0.1592978 
    ## ... Procrustes: rmse 4.84758e-05  max resid 0.0004579921 
    ## ... Similar to previous best
    ## Run 8 stress 0.159298 
    ## ... Procrustes: rmse 0.0001916711  max resid 0.002599921 
    ## ... Similar to previous best
    ## Run 9 stress 0.1592978 
    ## ... New best solution
    ## ... Procrustes: rmse 3.074129e-06  max resid 2.601491e-05 
    ## ... Similar to previous best
    ## Run 10 stress 0.159298 
    ## ... Procrustes: rmse 0.0001863176  max resid 0.002559743 
    ## ... Similar to previous best
    ## Run 11 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001778677  max resid 0.002530041 
    ## ... Similar to previous best
    ## Run 12 stress 0.159298 
    ## ... Procrustes: rmse 0.0001897222  max resid 0.002575608 
    ## ... Similar to previous best
    ## Run 13 stress 0.159298 
    ## ... Procrustes: rmse 0.0001873852  max resid 0.002568287 
    ## ... Similar to previous best
    ## Run 14 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001742413  max resid 0.002513958 
    ## ... Similar to previous best
    ## Run 15 stress 0.1732947 
    ## Run 16 stress 0.159298 
    ## ... Procrustes: rmse 0.0001840396  max resid 0.002555721 
    ## ... Similar to previous best
    ## Run 17 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001811929  max resid 0.00254666 
    ## ... Similar to previous best
    ## Run 18 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001782825  max resid 0.002535227 
    ## ... Similar to previous best
    ## Run 19 stress 0.1592978 
    ## ... Procrustes: rmse 5.411506e-05  max resid 0.0004055244 
    ## ... Similar to previous best
    ## Run 20 stress 0.159298 
    ## ... Procrustes: rmse 0.0001910864  max resid 0.002578653 
    ## ... Similar to previous best
    ## Run 21 stress 0.173285 
    ## Run 22 stress 0.1732843 
    ## Run 23 stress 0.159298 
    ## ... Procrustes: rmse 0.0001881864  max resid 0.002569862 
    ## ... Similar to previous best
    ## Run 24 stress 0.4189868 
    ## Run 25 stress 0.1592978 
    ## ... New best solution
    ## ... Procrustes: rmse 4.706562e-05  max resid 0.0003575562 
    ## ... Similar to previous best
    ## Run 26 stress 0.1734078 
    ## Run 27 stress 0.1732956 
    ## Run 28 stress 0.159298 
    ## ... Procrustes: rmse 0.0001637931  max resid 0.002303683 
    ## ... Similar to previous best
    ## Run 29 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001538236  max resid 0.00225558 
    ## ... Similar to previous best
    ## Run 30 stress 0.1592978 
    ## ... New best solution
    ## ... Procrustes: rmse 3.723227e-05  max resid 0.0002737049 
    ## ... Similar to previous best
    ## Run 31 stress 0.1732933 
    ## Run 32 stress 0.159298 
    ## ... Procrustes: rmse 0.0001826249  max resid 0.002527583 
    ## ... Similar to previous best
    ## Run 33 stress 0.1592984 
    ## ... Procrustes: rmse 0.0001082946  max resid 0.001775236 
    ## ... Similar to previous best
    ## Run 34 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001684633  max resid 0.002465911 
    ## ... Similar to previous best
    ## Run 35 stress 0.159298 
    ## ... Procrustes: rmse 0.0001825345  max resid 0.002527201 
    ## ... Similar to previous best
    ## Run 36 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001761583  max resid 0.002504945 
    ## ... Similar to previous best
    ## Run 37 stress 0.159298 
    ## ... Procrustes: rmse 0.0001782641  max resid 0.002510155 
    ## ... Similar to previous best
    ## Run 38 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001769084  max resid 0.002510759 
    ## ... Similar to previous best
    ## Run 39 stress 0.1592978 
    ## ... New best solution
    ## ... Procrustes: rmse 1.846949e-05  max resid 0.0001653995 
    ## ... Similar to previous best
    ## Run 40 stress 0.159298 
    ## ... Procrustes: rmse 0.0001687765  max resid 0.002359936 
    ## ... Similar to previous best
    ## Run 41 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001575324  max resid 0.002324305 
    ## ... Similar to previous best
    ## Run 42 stress 0.1732835 
    ## Run 43 stress 0.159298 
    ## ... Procrustes: rmse 0.0001708777  max resid 0.002365934 
    ## ... Similar to previous best
    ## Run 44 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001571633  max resid 0.002319541 
    ## ... Similar to previous best
    ## Run 45 stress 0.1592978 
    ## ... Procrustes: rmse 1.798223e-05  max resid 0.0001684401 
    ## ... Similar to previous best
    ## Run 46 stress 0.159298 
    ## ... Procrustes: rmse 0.0001723134  max resid 0.002371885 
    ## ... Similar to previous best
    ## Run 47 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001652757  max resid 0.002352307 
    ## ... Similar to previous best
    ## Run 48 stress 0.159298 
    ## ... Procrustes: rmse 0.0001694262  max resid 0.002365387 
    ## ... Similar to previous best
    ## Run 49 stress 0.159298 
    ## ... Procrustes: rmse 0.000170241  max resid 0.002365286 
    ## ... Similar to previous best
    ## Run 50 stress 0.159298 
    ## ... Procrustes: rmse 0.0001702841  max resid 0.002363807 
    ## ... Similar to previous best
    ## Run 51 stress 0.159298 
    ## ... Procrustes: rmse 0.0001715177  max resid 0.002367923 
    ## ... Similar to previous best
    ## Run 52 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001667571  max resid 0.002351892 
    ## ... Similar to previous best
    ## Run 53 stress 0.159298 
    ## ... Procrustes: rmse 0.0001699057  max resid 0.00236564 
    ## ... Similar to previous best
    ## Run 54 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001663053  max resid 0.002348971 
    ## ... Similar to previous best
    ## Run 55 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001607658  max resid 0.002372779 
    ## ... Similar to previous best
    ## Run 56 stress 0.159298 
    ## ... Procrustes: rmse 0.0001704547  max resid 0.002365915 
    ## ... Similar to previous best
    ## Run 57 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001601755  max resid 0.002322454 
    ## ... Similar to previous best
    ## Run 58 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001564841  max resid 0.002307413 
    ## ... Similar to previous best
    ## Run 59 stress 0.1732954 
    ## Run 60 stress 0.159298 
    ## ... Procrustes: rmse 0.0001673804  max resid 0.002357523 
    ## ... Similar to previous best
    ## Run 61 stress 0.159298 
    ## ... Procrustes: rmse 0.0001668516  max resid 0.002356011 
    ## ... Similar to previous best
    ## Run 62 stress 0.1734242 
    ## Run 63 stress 0.1592984 
    ## ... Procrustes: rmse 0.0001028808  max resid 0.001753664 
    ## ... Similar to previous best
    ## Run 64 stress 0.159298 
    ## ... Procrustes: rmse 0.0001702217  max resid 0.002371519 
    ## ... Similar to previous best
    ## Run 65 stress 0.1592978 
    ## ... Procrustes: rmse 1.751806e-05  max resid 0.0001612634 
    ## ... Similar to previous best
    ## Run 66 stress 0.159298 
    ## ... Procrustes: rmse 0.0001698481  max resid 0.002363816 
    ## ... Similar to previous best
    ## Run 67 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001541675  max resid 0.002304338 
    ## ... Similar to previous best
    ## Run 68 stress 0.159298 
    ## ... Procrustes: rmse 0.0001687191  max resid 0.002360054 
    ## ... Similar to previous best
    ## Run 69 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001603081  max resid 0.002325551 
    ## ... Similar to previous best
    ## Run 70 stress 0.1592978 
    ## ... Procrustes: rmse 2.504735e-05  max resid 0.0002482711 
    ## ... Similar to previous best
    ## Run 71 stress 0.159298 
    ## ... Procrustes: rmse 0.0001530986  max resid 0.002279854 
    ## ... Similar to previous best
    ## Run 72 stress 0.159298 
    ## ... Procrustes: rmse 0.0001679558  max resid 0.002354842 
    ## ... Similar to previous best
    ## Run 73 stress 0.159298 
    ## ... Procrustes: rmse 0.0001691845  max resid 0.002363229 
    ## ... Similar to previous best
    ## Run 74 stress 0.1592978 
    ## ... Procrustes: rmse 3.191509e-05  max resid 0.0002341286 
    ## ... Similar to previous best
    ## Run 75 stress 0.159298 
    ## ... Procrustes: rmse 0.0001707753  max resid 0.002365066 
    ## ... Similar to previous best
    ## Run 76 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001641643  max resid 0.002342676 
    ## ... Similar to previous best
    ## Run 77 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001596832  max resid 0.002320403 
    ## ... Similar to previous best
    ## Run 78 stress 0.1732939 
    ## Run 79 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001673727  max resid 0.002354772 
    ## ... Similar to previous best
    ## Run 80 stress 0.159298 
    ## ... Procrustes: rmse 0.0001591718  max resid 0.002318047 
    ## ... Similar to previous best
    ## Run 81 stress 0.159298 
    ## ... Procrustes: rmse 0.0001691361  max resid 0.002360076 
    ## ... Similar to previous best
    ## Run 82 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001649742  max resid 0.002331985 
    ## ... Similar to previous best
    ## Run 83 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001662541  max resid 0.002348243 
    ## ... Similar to previous best
    ## Run 84 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001677812  max resid 0.002355366 
    ## ... Similar to previous best
    ## Run 85 stress 0.159298 
    ## ... Procrustes: rmse 0.0001713906  max resid 0.002368524 
    ## ... Similar to previous best
    ## Run 86 stress 0.159298 
    ## ... Procrustes: rmse 0.0001693362  max resid 0.002360799 
    ## ... Similar to previous best
    ## Run 87 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001578846  max resid 0.002311922 
    ## ... Similar to previous best
    ## Run 88 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001665726  max resid 0.002356299 
    ## ... Similar to previous best
    ## Run 89 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001656648  max resid 0.002349422 
    ## ... Similar to previous best
    ## Run 90 stress 0.159298 
    ## ... Procrustes: rmse 0.0001714476  max resid 0.002373879 
    ## ... Similar to previous best
    ## Run 91 stress 0.159298 
    ## ... Procrustes: rmse 0.0001680816  max resid 0.002362071 
    ## ... Similar to previous best
    ## Run 92 stress 0.159298 
    ## ... Procrustes: rmse 0.0001718693  max resid 0.002368652 
    ## ... Similar to previous best
    ## Run 93 stress 0.159298 
    ## ... Procrustes: rmse 0.0001719763  max resid 0.002373932 
    ## ... Similar to previous best
    ## Run 94 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001626132  max resid 0.002334844 
    ## ... Similar to previous best
    ## Run 95 stress 0.1732947 
    ## Run 96 stress 0.1592978 
    ## ... Procrustes: rmse 2.006637e-05  max resid 0.0001804203 
    ## ... Similar to previous best
    ## Run 97 stress 0.159298 
    ## ... Procrustes: rmse 0.0001671161  max resid 0.002352603 
    ## ... Similar to previous best
    ## Run 98 stress 0.159298 
    ## ... Procrustes: rmse 0.0001681026  max resid 0.002356779 
    ## ... Similar to previous best
    ## Run 99 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001548215  max resid 0.002298373 
    ## ... Similar to previous best
    ## Run 100 stress 0.1592979 
    ## ... Procrustes: rmse 0.0001641968  max resid 0.002346896 
    ## ... Similar to previous best
    ## *** Best solution repeated 57 times

``` r
nmds_comm_trajectories
```

    ## 
    ## Call:
    ## metaMDS(comm = comm_trajectories_relative, distance = "bray",      k = 2, try = 100, trymax = 100, engine = "monoMDS", autotransform = FALSE) 
    ## 
    ## global Multidimensional Scaling using monoMDS
    ## 
    ## Data:     comm_trajectories_relative 
    ## Distance: bray 
    ## 
    ## Dimensions: 2 
    ## Stress:     0.1592978 
    ## Stress type 1, weak ties
    ## Best solution was repeated 57 times in 100 tries
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
comm_AM1_nmds
```

    ##           MDS1        MDS2
    ## 1  -0.53257790 -0.33472191
    ## 2   0.05909854 -0.63761137
    ## 4  -0.47895339  0.04423619
    ## 5  -0.02494590 -0.92390787
    ## 6   0.14229480  0.40047447
    ## 7  -0.99943866  0.75492504
    ## 9   0.18249741 -0.37640119
    ## 10  0.77174316 -0.21161113
    ## 11  0.03885077 -0.53417163
    ## 13 -0.88750289  0.32365796
    ## 14  0.61716845  0.79863848
    ## 15  0.24362653  0.07145690
    ## 16 -0.34964218 -0.59375616
    ## 18 -0.23116413  0.36261064
    ## 19 -0.28434617  0.54353746
    ## 20 -0.06155094  0.06936839
    ## 22 -0.29239780 -0.07251778
    ## 23 -0.56009752  0.28288563
    ## 24  0.46413310  0.17622269
    ## 25  0.60791992 -0.13073093
    ## 26  1.47675372  0.25607174
    ## 27  0.53704901 -0.63904890
    ## 28 -0.74753453 -0.57292438
    ## 30  0.30901660  0.94331767

``` r
species <- nmds_comm_AM1_rm_total$species

sp_names <- rownames(species)
sp_names <- gsub("\\.", " ", sp_names)



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


rect(xleft = 0, ybottom = 0, xright = xmax*1.2, ytop = ymax*1.2 ,col = transparent("#56B4E9", trans.val = 0.75), border = "white", lwd = 4)
rect(xleft = xmin*1.2, ybottom = ymin*1.2, xright = 0, ytop = 0,col = transparent("#E69F00", trans.val = 0.75), border = "white", lwd = 4)
rect(xleft = 0, ybottom = ymin*1.2, xright = xmax*1.2, ytop = 0,col = transparent("#CC79A7", trans.val = 0.75), border = "white", lwd = 4)
rect(xleft = xmin*1.2, ybottom = 0, xright = 0, ytop = ymax*1.2,col = transparent("#009E73", trans.val = 0.75), border = "white", lwd = 4)


abline(h = 0, v = 0, lty = 2)

library(scales)
library(shape)


#pal <- col_numeric(palette = c("white", "black"), domain = urb, na.color = "grey50", alpha = FALSE, reverse = FALSE)
#col <-pal(urb)

points(comm_AM1_nmds[,1],comm_AM1_nmds[,2], col = "grey70", bg = "white", pch = 16, cex = 1.5)

Arrows(x0 <- rep(0, nrow(species)),
       y0 <- rep(0, nrow(species)),
       x1 <- species[,1],
       y1 <- species[,2], arr.type = "triangle", arr.length = 0.4, col = "#C49C94", lwd = 1.5)

text(x = species[,1]*1.2, y = species[,2]*1.2, labels = sp_names, cex = 0.8, col = "black", font = 2)

axis(1, cex.axis = 1)
axis(2, cex.axis = 1, las = 2)
title(xlab = "nMDS 1", cex.lab = 1.2, line = 3)
title(ylab = "nMDS 2", cex.lab = 1.2, line = 3)
#title(main = "Morphological traits", line = 0.5, adj = 0, cex.main = 1.5)
letters(x = 10, y = 97, "a)", cex = 1.5)


screen(2)

#set.seed(1);species_pred <- jitter(species_pred, amount = 0.25)

xmin <- min(c(points_high_MDS1_high_MDS2[,1], points_low_MDS1_high_MDS2[,1], points_high_MDS1_low_MDS2[,1], points_low_MDS1_low_MDS2[,1],species_pred[,1]))*1.1
xmax <- max(c(points_high_MDS1_high_MDS2[,1], points_low_MDS1_high_MDS2[,1], points_high_MDS1_low_MDS2[,1], points_low_MDS1_low_MDS2[,1],species_pred[,1]))*1.1

ymin <- min(c(points_high_MDS1_high_MDS2[,2], points_low_MDS1_high_MDS2[,2], points_high_MDS1_low_MDS2[,2], points_low_MDS1_low_MDS2[,2],species_pred[,2]))*1.1
ymax <- max(c(points_high_MDS1_high_MDS2[,2], points_low_MDS1_high_MDS2[,2], points_high_MDS1_low_MDS2[,2], points_low_MDS1_low_MDS2[,2],species_pred[,2]))*1.1

par(mar = c(4, 4, 1, 1))
  
plot(NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = "", ylab = "", xaxt = "n", yaxt = "n")

lines(x = points_high_MDS1_high_MDS2[,1], y = points_high_MDS1_high_MDS2[,2], lwd = 3, lty = 1, col = "#56B4E9") #controle

lines(x = points_low_MDS1_high_MDS2[,1], y = points_low_MDS1_high_MDS2[,2], lwd = 3, lty = 1, col = "#009E73") #atrasado

lines(x = points_high_MDS1_low_MDS2[,1], y = points_high_MDS1_low_MDS2[,2], lwd = 3, lty = 1, col = "#CC79A7") #fechado

lines(x = points_low_MDS1_low_MDS2[,1], y = points_low_MDS1_low_MDS2[,2], lwd = 3, lty = 1, col = "#E69F00") #fechado



#####################################################################



##################################################

points(x = points_high_MDS1_high_MDS2[find_day(day = 32, min = 32, max = 158),1],
       y = points_high_MDS1_high_MDS2[find_day(day = 32, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#56B4E9", lwd = 3, cex = 3.5)

text(x = points_high_MDS1_high_MDS2[find_day(day = 32, min = 32, max = 158),1],
       y = points_high_MDS1_high_MDS2[find_day(day = 32, min = 32, max = 158),2], labels = "32", cex = 0.7, col = "#56B4E9", font = 2)

points(x = points_high_MDS1_high_MDS2[find_day(day = 80, min = 32, max = 158),1],
       y = points_high_MDS1_high_MDS2[find_day(day = 80, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#56B4E9", lwd = 3, cex = 3.5)

text(x = points_high_MDS1_high_MDS2[find_day(day = 80, min = 32, max = 158),1],
       y = points_high_MDS1_high_MDS2[find_day(day = 80, min = 32, max = 158),2], labels = "80", cex = 0.7, col = "#56B4E9", font = 2)

points(x = points_high_MDS1_high_MDS2[find_day(day = 116, min = 32, max = 158),1],
       y = points_high_MDS1_high_MDS2[find_day(day = 116, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#56B4E9", lwd = 3, cex = 3.5)

text(x = points_high_MDS1_high_MDS2[find_day(day = 116, min = 32, max = 158),1],
       y = points_high_MDS1_high_MDS2[find_day(day = 116, min = 32, max = 158),2], labels = "116", cex = 0.7, col = "#56B4E9", font = 2)

points(x = points_high_MDS1_high_MDS2[find_day(day = 158, min = 32, max = 158),1],
       y = points_high_MDS1_high_MDS2[find_day(day = 158, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#56B4E9", lwd = 3, cex = 3.5)

text(x = points_high_MDS1_high_MDS2[find_day(day = 158, min = 32, max = 158),1],
       y = points_high_MDS1_high_MDS2[find_day(day = 158, min = 32, max = 158),2], labels = "158", cex = 0.7, col = "#56B4E9", font = 2)
#################################################


##################################################

points(x = points_low_MDS1_high_MDS2[find_day(day = 32, min = 32, max = 158),1],
       y = points_low_MDS1_high_MDS2[find_day(day = 32, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#009E73", lwd = 3, cex = 3.5)

text(x = points_low_MDS1_high_MDS2[find_day(day = 32, min = 32, max = 158),1],
       y = points_low_MDS1_high_MDS2[find_day(day = 32, min = 32, max = 158),2], labels = "32", cex = 0.7, font = 2, col = "#009E73")

points(x = points_low_MDS1_high_MDS2[find_day(day = 80, min = 32, max = 158),1],
       y = points_low_MDS1_high_MDS2[find_day(day = 80, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#009E73", lwd = 3, cex = 3.5)

text(x = points_low_MDS1_high_MDS2[find_day(day = 80, min = 32, max = 158),1],
       y = points_low_MDS1_high_MDS2[find_day(day = 80, min = 32, max = 158),2], labels = "80", cex = 0.7, font = 2, col = "#009E73")

points(x = points_low_MDS1_high_MDS2[find_day(day = 116, min = 32, max = 158),1],
       y = points_low_MDS1_high_MDS2[find_day(day = 116, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#009E73", lwd = 3, cex = 3.5)

text(x = points_low_MDS1_high_MDS2[find_day(day = 116, min = 32, max = 158),1],
       y = points_low_MDS1_high_MDS2[find_day(day = 116, min = 32, max = 158),2], labels = "116", cex = 0.7, font = 2, col = "#009E73")

points(x = points_low_MDS1_high_MDS2[find_day(day = 158, min = 32, max = 158),1],
       y = points_low_MDS1_high_MDS2[find_day(day = 158, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#009E73", lwd = 3, cex = 3.5)

text(x = points_low_MDS1_high_MDS2[find_day(day = 158, min = 32, max = 158),1],
       y = points_low_MDS1_high_MDS2[find_day(day = 158, min = 32, max = 158),2], labels = "158", cex = 0.7, font = 2, col = "#009E73")
#################################################



##################################################

points(x = points_high_MDS1_low_MDS2[find_day(day = 32, min = 32, max = 158),1],
       y = points_high_MDS1_low_MDS2[find_day(day = 32, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#CC79A7", lwd = 3, cex = 3.5)

text(x = points_high_MDS1_low_MDS2[find_day(day = 32, min = 32, max = 158),1],
       y = points_high_MDS1_low_MDS2[find_day(day = 32, min = 32, max = 158),2], labels = "32", cex = 0.7, font = 2, col = "#CC79A7")

points(x = points_high_MDS1_low_MDS2[find_day(day = 80, min = 32, max = 158),1],
       y = points_high_MDS1_low_MDS2[find_day(day = 80, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#CC79A7", lwd = 3, cex = 3.5)

text(x = points_high_MDS1_low_MDS2[find_day(day = 80, min = 32, max = 158),1],
       y = points_high_MDS1_low_MDS2[find_day(day = 80, min = 32, max = 158),2], labels = "80", cex = 0.7, font = 2, col = "#CC79A7")

points(x = points_high_MDS1_low_MDS2[find_day(day = 116, min = 32, max = 158),1],
       y = points_high_MDS1_low_MDS2[find_day(day = 116, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#CC79A7", lwd = 3, cex = 3.5)

text(x = points_high_MDS1_low_MDS2[find_day(day = 116, min = 32, max = 158),1],
       y = points_high_MDS1_low_MDS2[find_day(day = 116, min = 32, max = 158),2], labels = "116", cex = 0.7, font = 2, col = "#CC79A7")

points(x = points_high_MDS1_low_MDS2[find_day(day = 158, min = 32, max = 158),1],
       y = points_high_MDS1_low_MDS2[find_day(day = 158, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#CC79A7", lwd = 3, cex = 3.5)

text(x = points_high_MDS1_low_MDS2[find_day(day = 158, min = 32, max = 158),1],
       y = points_high_MDS1_low_MDS2[find_day(day = 158, min = 32, max = 158),2], labels = "158", cex = 0.7, font = 2, col = "#CC79A7")
#################################################





##################################################

points(x = points_low_MDS1_low_MDS2[find_day(day = 32, min = 32, max = 158),1],
       y = points_low_MDS1_low_MDS2[find_day(day = 32, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#E69F00", lwd = 3, cex = 3.5)

text(x = points_low_MDS1_low_MDS2[find_day(day = 32, min = 32, max = 158),1],
       y = points_low_MDS1_low_MDS2[find_day(day = 32, min = 32, max = 158),2], labels = "32", cex = 0.7, font = 2, col = "#E69F00")

points(x = points_low_MDS1_low_MDS2[find_day(day = 80, min = 32, max = 158),1],
       y = points_low_MDS1_low_MDS2[find_day(day = 80, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#E69F00", lwd = 3, cex = 3.5)

text(x = points_low_MDS1_low_MDS2[find_day(day = 80, min = 32, max = 158),1],
       y = points_low_MDS1_low_MDS2[find_day(day = 80, min = 32, max = 158),2], labels = "80", cex = 0.7, font = 2, col = "#E69F00")

points(x = points_low_MDS1_low_MDS2[find_day(day = 116, min = 32, max = 158),1],
       y = points_low_MDS1_low_MDS2[find_day(day = 116, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#E69F00", lwd = 3, cex = 3.5)

text(x = points_low_MDS1_low_MDS2[find_day(day = 116, min = 32, max = 158),1],
       y = points_low_MDS1_low_MDS2[find_day(day = 116, min = 32, max = 158),2], labels = "116", cex = 0.7, font = 2, col = "#E69F00")

points(x = points_low_MDS1_low_MDS2[find_day(day = 158, min = 32, max = 158),1],
       y = points_low_MDS1_low_MDS2[find_day(day = 158, min = 32, max = 158),2],
       pch = 21, bg = "white", col = "#E69F00", lwd = 3, cex = 3.5)

text(x = points_low_MDS1_low_MDS2[find_day(day = 158, min = 32, max = 158),1],
       y = points_low_MDS1_low_MDS2[find_day(day = 158, min = 32, max = 158),2], labels = "158", cex = 0.7, font = 2, col = "#E69F00")
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

letters(x = 10, y = 97, "b)", cex = 1.5)
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
    ##  Resampling run 300 finished. Time elapsed: 0.12 minutes...
    ##  Resampling run 400 finished. Time elapsed: 0.16 minutes...
    ##  Resampling run 500 finished. Time elapsed: 0.20 minutes...
    ##  Resampling run 600 finished. Time elapsed: 0.24 minutes...
    ##  Resampling run 700 finished. Time elapsed: 0.28 minutes...
    ##  Resampling run 800 finished. Time elapsed: 0.32 minutes...
    ##  Resampling run 900 finished. Time elapsed: 0.36 minutes...
    ## Time elapsed: 0 hr 0 min 23 sec

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
    ##  Resampling run 200 finished. Time elapsed: 0.07 minutes...
    ##  Resampling run 300 finished. Time elapsed: 0.10 minutes...
    ##  Resampling run 400 finished. Time elapsed: 0.14 minutes...
    ##  Resampling run 500 finished. Time elapsed: 0.17 minutes...
    ##  Resampling run 600 finished. Time elapsed: 0.20 minutes...
    ##  Resampling run 700 finished. Time elapsed: 0.23 minutes...
    ##  Resampling run 800 finished. Time elapsed: 0.27 minutes...
    ##  Resampling run 900 finished. Time elapsed: 0.30 minutes...
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
    ##  Resampling run 200 finished. Time elapsed: 0.07 minutes...
    ##  Resampling run 300 finished. Time elapsed: 0.10 minutes...
    ##  Resampling run 400 finished. Time elapsed: 0.13 minutes...
    ##  Resampling run 500 finished. Time elapsed: 0.17 minutes...
    ##  Resampling run 600 finished. Time elapsed: 0.20 minutes...
    ##  Resampling run 700 finished. Time elapsed: 0.24 minutes...
    ##  Resampling run 800 finished. Time elapsed: 0.27 minutes...
    ##  Resampling run 900 finished. Time elapsed: 0.30 minutes...
    ## Time elapsed: 0 hr 0 min 20 sec

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
    ##  Resampling run 100 finished. Time elapsed: 0.03 minutes...
    ##  Resampling run 200 finished. Time elapsed: 0.07 minutes...
    ##  Resampling run 300 finished. Time elapsed: 0.10 minutes...
    ##  Resampling run 400 finished. Time elapsed: 0.13 minutes...
    ##  Resampling run 500 finished. Time elapsed: 0.17 minutes...
    ##  Resampling run 600 finished. Time elapsed: 0.21 minutes...
    ##  Resampling run 700 finished. Time elapsed: 0.25 minutes...
    ##  Resampling run 800 finished. Time elapsed: 0.29 minutes...
    ##  Resampling run 900 finished. Time elapsed: 0.32 minutes...
    ## Time elapsed: 0 hr 0 min 21 sec

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
    ##  Resampling run 200 finished. Time elapsed: 0.09 minutes...
    ##  Resampling run 300 finished. Time elapsed: 0.13 minutes...
    ##  Resampling run 400 finished. Time elapsed: 0.20 minutes...
    ##  Resampling run 500 finished. Time elapsed: 0.25 minutes...
    ##  Resampling run 600 finished. Time elapsed: 0.30 minutes...
    ##  Resampling run 700 finished. Time elapsed: 0.35 minutes...
    ##  Resampling run 800 finished. Time elapsed: 0.39 minutes...
    ##  Resampling run 900 finished. Time elapsed: 0.44 minutes...
    ## Time elapsed: 0 hr 0 min 28 sec

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
    ##  Resampling run 100 finished. Time elapsed: 0.04 minutes...
    ##  Resampling run 200 finished. Time elapsed: 0.09 minutes...
    ##  Resampling run 300 finished. Time elapsed: 0.14 minutes...
    ##  Resampling run 400 finished. Time elapsed: 0.19 minutes...
    ##  Resampling run 500 finished. Time elapsed: 0.23 minutes...
    ##  Resampling run 600 finished. Time elapsed: 0.27 minutes...
    ##  Resampling run 700 finished. Time elapsed: 0.32 minutes...
    ##  Resampling run 800 finished. Time elapsed: 0.36 minutes...
    ##  Resampling run 900 finished. Time elapsed: 0.41 minutes...
    ## Time elapsed: 0 hr 0 min 27 sec

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
