library(vegan)
library(gllvm)

source("functions/pairwise_gllvm.R")
source("functions/my.anova.gllvm.R")
source("functions/My_coefplot.R")
source("functions/remove_sp.R")
source("functions/extract_mles.R")
source("functions/my_ordiplot.R")
source("functions/get_scaled_lvs.R")



comm_AM1 <- remove_sp(communities_AM1, 2)
comm_AM2 <- remove_sp(communities_AM2, 2)
comm_AM3 <- remove_sp(communities_AM3, 2)
comm_AM4 <- remove_sp(communities_AM4, 2)


#remove_AM1 <- colnames(comm_AM1) != "Scinax fuscovarius" & colnames(comm_AM1) != "Dendropsophus minutus"
#remove_AM2 <- colnames(comm_AM2) != "Scinax fuscovarius" & colnames(comm_AM2) != "Dendropsophus minutus"
#remove_AM3 <- colnames(comm_AM3) != "Scinax fuscovarius" & colnames(comm_AM3) != "Dendropsophus minutus"
#remove_AM4 <- colnames(comm_AM4) != "Scinax fuscovarius" & colnames(comm_AM4) != "Dendropsophus minutus"

#comm_AM1<-comm_AM1[,remove_AM1]
#comm_AM2<-comm_AM2[,remove_AM2]
#comm_AM3<-comm_AM3[,remove_AM3]
#comm_AM4<-comm_AM4[,remove_AM4]

ncol(comm_AM1)
ncol(comm_AM2)
ncol(comm_AM3)
ncol(comm_AM4)

################################### AM1 ####
rem_atradado <- Exp_design$treatments != "atrasado"

mod_block_gllvm_AM1 <- gllvm(y = comm_AM1[rem_atradado,], formula = ~ block2, X = Exp_design[rem_atradado,], family = "negative.binomial",
                             num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 2, starting.val = "res"), row.eff = "fixed", seed = 10)
mod_treat_gllvm_AM1 <- gllvm(y = comm_AM1[rem_atradado,], formula = ~ treatments, X = Exp_design[rem_atradado,], family = "negative.binomial",
                             num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 2, starting.val = "res"), row.eff = "fixed", seed = 10)
mod_block_treat_gllvm_AM1 <- gllvm(y = comm_AM1[rem_atradado,], formula = ~ block2 + treatments, X = Exp_design[rem_atradado,], family = "negative.binomial",
                                   num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 2, starting.val = "res"), row.eff = "fixed", seed = 10)

anova_block_AM1<-my.anova.gllvm(mod_treat_gllvm_AM1, mod_block_treat_gllvm_AM1)
anova_treatment_AM1<-my.anova.gllvm(mod_block_gllvm_AM1, mod_block_treat_gllvm_AM1)



pairwise_AM1_treatments  <- pairwise_gllvm(com = comm_AM1[rem_atradado,],
                                           pairwise = Exp_design[rem_atradado,]$treatments,
                                           covariable = Exp_design[rem_atradado,]$block2,
                                           family = "negative.binomial", remove_species =  1, num.lv = 1,
                                           seed = 10,
                                           n.init = 100, method = "VA", n.init.max = 25, starting.val = "res", row.eff = "fixed")
pairwise_AM1_treatments$results
########################################


################################### AM2 ####
mod_block_gllvm_AM2 <- gllvm(y = comm_AM2, formula = ~ block2, X = Exp_design, family = "negative.binomial", num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 2, starting.val = "res"), row.eff = "fixed", seed = 10)
mod_treat_gllvm_AM2 <- gllvm(y = comm_AM2, formula = ~ treatments, X = Exp_design, family = "negative.binomial", num.lv = 1, method = "VA", control.start    = list(n.init = 1000, n.init.max = 100, jitter.var = 2, starting.val = "res"), row.eff = "fixed", seed = 10)
mod_block_treat_gllvm_AM2 <- gllvm(y = comm_AM2, formula = ~ block2 + treatments, X = Exp_design, family = "negative.binomial", num.lv = 1, method = "VA", control.start    = list(n.init = 1000, n.init.max = 100, jitter.var = 2, starting.val = "res"), row.eff = "fixed", seed = 10)

anova_block_AM2<-my.anova.gllvm(mod_treat_gllvm_AM2, mod_block_treat_gllvm_AM2)
anova_treatment_AM2<-my.anova.gllvm(mod_block_gllvm_AM2, mod_block_treat_gllvm_AM2)


pairwise_AM2_treatments <- pairwise_gllvm(com = comm_AM2,
                                          pairwise = Exp_design$treatments,
                                          covariable = Exp_design$block2,
                                          family = "negative.binomial", remove_species =  1, num.lv = 1,
                                          seed = 1,
                                          n.init = 199, method = "VA", n.init.max = 25, starting.val = "res", row.eff = "fixed")

pairwise_AM2_treatments$results
########################################


################################### AM3 ####

mod_block_gllvm_AM3 <- gllvm(y = comm_AM3, formula = ~ block2, X = Exp_design, family = "negative.binomial", num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 2, starting.val = "res"), row.eff = "fixed", seed = 10)
mod_treat_gllvm_AM3 <- gllvm(y = comm_AM3, formula = ~ treatments, X = Exp_design, family = "negative.binomial", num.lv = 1, method = "VA", control.start    = list(n.init = 1000, n.init.max = 100, jitter.var = 2, starting.val = "res"), row.eff = "fixed", seed = 10)
mod_block_treat_gllvm_AM3 <- gllvm(y = comm_AM3, formula = ~ block2 + treatments, X = Exp_design, family = "negative.binomial", num.lv = 1, method = "VA", control.start    = list(n.init = 1000, n.init.max = 100, jitter.var = 2, starting.val = "res"), row.eff = "fixed", seed = 10)

anova_block_AM3<-my.anova.gllvm(mod_treat_gllvm_AM3, mod_block_treat_gllvm_AM3)
anova_treatment_AM3<-my.anova.gllvm(mod_block_gllvm_AM3, mod_block_treat_gllvm_AM3)


pairwise_AM3_treatments <- pairwise_gllvm(com = comm_AM3,
                                          pairwise = Exp_design$treatments,
                                          covariable = Exp_design$block2,
                                          family = "negative.binomial", remove_species =  1, num.lv = 1,
                                          seed = 1,
                                          n.init = 199, method = "VA", n.init.max = 25, starting.val = "res", row.eff = "fixed")

pairwise_AM3_treatments$results
########################################


################################### AM4 ####

mod_block_gllvm_AM4 <- gllvm_mc(y = comm_AM4,
                    formula = ~ block2,
                    X = Exp_design,
                    family = "negative.binomial",
                    num.lv = 1, method = "VA", starting.val = "res",
                    seed = 3, row.eff = "fixed",
                    export_list = c("gllvm", "Exp_design","rem_atradado"), n_cores = detectCores()-2, n_runs = 1000)
plot(mod_block_gllvm_AM4$AICc, type = "lines")
mod_block_gllvm_AM4$best
mod_block_gllvm_AM4$best$sd$Xcoef
mod_block_gllvm_AM4$best$sd$Xcoef



mod_treat_gllvm_AM4 <- gllvm_mc(y = comm_AM4,
                                formula = ~ treatments,
                                X = Exp_design,
                                family = "negative.binomial",
                                num.lv = 1, method = "VA", starting.val = "res",
                                seed = 3, row.eff = "fixed",
                                export_list = c("gllvm", "Exp_design","rem_atradado"), n_cores = detectCores()-2, n_runs = 1000)
plot(mod_treat_gllvm_AM4$AICc, type = "lines")
mod_treat_gllvm_AM4$best
mod_treat_gllvm_AM4$best$sd$Xcoef
mod_treat_gllvm_AM4$best$sd$Xcoef



mod_block_treat_gllvm_AM4 <- gllvm_mc(y = comm_AM4,
                                formula = ~ block2 + treatments,
                                X = Exp_design,
                                family = "negative.binomial",
                                num.lv = 1, method = "VA", starting.val = "res",
                                seed = 3, row.eff = "fixed",
                                export_list = c("gllvm", "Exp_design","rem_atradado"), n_cores = detectCores()-2, n_runs = 1000)
plot(mod_block_treat_gllvm_AM4$AICc, type = "lines")



#mod_block_gllvm_AM4 <- gllvm(y = comm_AM4,formula = ~ block2, X = Exp_design, family = "negative.binomial", num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 10, starting.val = "res"), row.eff = "fixed", seed = 10)
#mod_treat_gllvm_AM4 <- gllvm(y = comm_AM4, formula = ~ treatments, X = Exp_design, family = "negative.binomial", num.lv = 1, method = "VA", control.start    = list(n.init = 1000, n.init.max = 100, jitter.var = 10, starting.val = "res"), row.eff = "fixed", seed = 10)
#mod_block_treat_gllvm_AM4 <- gllvm(y = comm_AM4, formula = ~ block2 + treatments, X = Exp_design, family = "negative.binomial", num.lv = 1, method = "VA", control.start    = list(n.init = 1000, n.init.max = 100, jitter.var = 10, starting.val = "res"), row.eff = "fixed", seed = 10)

anova_block_AM4<-my.anova.gllvm(mod_treat_gllvm_AM4, mod_block_treat_gllvm_AM4)
anova_treatment_AM4<-my.anova.gllvm(mod_block_gllvm_AM4, mod_block_treat_gllvm_AM4)


pairwise_AM4_treatments <- pairwise_gllvm(com = comm_AM4,
                                          pairwise = Exp_design$treatments,
                                          covariable = Exp_design$block2,
                                          family = "negative.binomial", remove_species =  1, num.lv = 1,
                                          seed = 10,
                                          n.init = 199, method = "VA", n.init.max = 25, starting.val = "zero", row.eff = "fixed")

pairwise_AM4_treatments$results
########################################

###################################################
comm_all_drop_AM1 <- comm_all[Exp_design_all$AM != 1,]
Exp_design_drop_AM1 <- Exp_design_all[Exp_design_all$AM != 1,]

Exp_design_drop_AM1$sites_by_sample <- paste(Exp_design_drop_AM1$sites,Exp_design_drop_AM1$AM, sep ="_")
Exp_design_drop_AM1$AM <- as.numeric(Exp_design_drop_AM1$AM)
comm_all_drop_AM1 <- remove_sp(comm_all_drop_AM1,2)

comm_all_drop_atrasado <- comm_all[Exp_design_all$treatments != "atrasado",]
Exp_design_drop_atrasado <- Exp_design_all[Exp_design_all$treatments != "atrasado",]

Exp_design_drop_atrasado$AM <- as.numeric(Exp_design_drop_atrasado$AM)
comm_all_drop_atrasado <- remove_sp(comm_all_drop_atrasado,2)


Exp_design_drop_AM1$AM2  <- Exp_design_drop_AM1$AM^2
#####################################################


mod_block <- gllvm(y = comm_all_drop_AM1, formula = ~ block2, X = Exp_design_drop_AM1, family = "negative.binomial", num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 10, starting.val = "res"), row.eff = "fixed", seed = 10)
mod_block_AM <- gllvm(y = comm_all_drop_AM1, formula = ~ block2 + AM, X = Exp_design_drop_AM1, family = "negative.binomial", num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 10, starting.val = "res"), row.eff = "fixed", seed = 10)
mod_block_AM2 <- gllvm(y = comm_all_drop_AM1, formula = ~ block2 + AM + AM2, X = Exp_design_drop_AM1, family = "negative.binomial", num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 10, starting.val = "res"), row.eff = "fixed", seed = 10)
mod_block_treat <- gllvm(y = comm_all_drop_AM1, formula = ~ block2 + treatments, X = Exp_design_drop_AM1, family = "negative.binomial", num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 10, starting.val = "res"), row.eff = "fixed", seed = 10)
mod_block_treat_AM <- gllvm(y = comm_all_drop_AM1, formula = ~ block2 + treatments + AM, X = Exp_design_drop_AM1, family = "negative.binomial", num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 10, starting.val = "res"), row.eff = "fixed", seed = 10)
mod_block_treat_AM2 <- gllvm(y = comm_all_drop_AM1, formula = ~ block2 + treatments + AM + AM2, X = Exp_design_drop_AM1, family = "negative.binomial", num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 10, starting.val = "res"), row.eff = "fixed", seed = 10)
mod_block_interaction <- gllvm(y = comm_all_drop_AM1, formula = ~ block2 + treatments * AM, X = Exp_design_drop_AM1, family = "negative.binomial", num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 10, starting.val = "res"), row.eff = "fixed", seed = 10)
mod_block_interaction_quad <- gllvm(y = comm_all_drop_AM1, formula = ~ block2 + treatments * (AM + AM2), X = Exp_design_drop_AM1, family = "negative.binomial", num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 10, starting.val = "res"), row.eff = "fixed", seed = 10)

mod_block_interaction_quad$params$Xcoef

anova_AM<-my.anova.gllvm(mod_block_treat, mod_block_treat_AM)
anova_treat<-my.anova.gllvm(mod_block_AM, mod_block_treat_AM)
anova_interaction<-my.anova.gllvm(mod_block_treat_AM, mod_block_interaction)

anova_AM<-my.anova.gllvm(mod_block_treat, mod_block_treat_AM)
anova_AM2 <-my.anova.gllvm(mod_block_treat_AM, mod_block_treat_AM2)
anova_treat2 <-my.anova.gllvm(mod_block_AM2, mod_block_treat_AM2)
anova_interaction2 <-my.anova.gllvm(mod_block_treat_AM2, mod_block_interaction_quad)

###################################################


#####################################################

Exp_design_drop_atrasado$AM2  <- Exp_design_drop_atrasado$AM^2

mod_block <- gllvm(y = comm_all_drop_atrasado, formula = ~ block2, X = Exp_design_drop_atrasado, family = "negative.binomial", num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 10, starting.val = "res"), row.eff = "fixed", seed = 10)
mod_block_AM <- gllvm(y = comm_all_drop_atrasado, formula = ~ block2 + AM, X = Exp_design_drop_atrasado, family = "negative.binomial", num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 10, starting.val = "res"), row.eff = "fixed", seed = 10)
mod_block_AM2 <- gllvm(y = comm_all_drop_atrasado, formula = ~ block2 + AM + AM2, X = Exp_design_drop_atrasado, family = "negative.binomial", num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 10, starting.val = "res"), row.eff = "fixed", seed = 10)
mod_block_treat <- gllvm(y = comm_all_drop_atrasado, formula = ~ block2 + treatments, X = Exp_design_drop_atrasado, family = "negative.binomial", num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 10, starting.val = "res"), row.eff = "fixed", seed = 10)
mod_block_treat_AM <- gllvm(y = comm_all_drop_atrasado, formula = ~ block2 + treatments + AM, X = Exp_design_drop_atrasado, family = "negative.binomial", num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 10, starting.val = "res"), row.eff = "fixed", seed = 10)
mod_block_treat_AM2 <- gllvm(y = comm_all_drop_atrasado, formula = ~ block2 + treatments + AM + AM2, X = Exp_design_drop_atrasado, family = "negative.binomial", num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 10, starting.val = "res"), row.eff = "fixed", seed = 10)
mod_block_interaction <- gllvm(y = comm_all_drop_atrasado, formula = ~ block2 + treatments * AM, X = Exp_design_drop_atrasado, family = "negative.binomial", num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 10, starting.val = "res"), row.eff = "fixed", seed = 10)
mod_block_interaction_quad <- gllvm(y = comm_all_drop_atrasado, formula = ~ block2 + treatments * (AM + AM2), X = Exp_design_drop_atrasado, family = "negative.binomial", num.lv = 1, method = "VA", control.start  = list(n.init = 1000, n.init.max = 100, jitter.var = 10, starting.val = "res"), row.eff = "fixed", seed = 10)

mod_block_interaction_quad$params$Xcoef

anova_AM<-my.anova.gllvm(mod_block_treat, mod_block_treat_AM)
anova_treat<-my.anova.gllvm(mod_block_AM, mod_block_treat_AM)
anova_interaction<-my.anova.gllvm(mod_block_treat_AM, mod_block_interaction)

anova_AM<-my.anova.gllvm(mod_block_treat, mod_block_treat_AM)
anova_AM2 <-my.anova.gllvm(mod_block_treat_AM, mod_block_treat_AM2)
anova_treat2 <-my.anova.gllvm(mod_block_AM2, mod_block_treat_AM2)
anova_interaction2 <-my.anova.gllvm(mod_block_treat_AM2, mod_block_interaction_quad)



###################################################
