#Ms Cc repeated heat shock (developmental timing heat shock) preliminary analyses 


#load libraries
library(readr)
library(ggplot2)
library(tidyr)
library(Rmisc)
library(dplyr)
library(nlme)
library(lme4)
library(scales)


#load data
#load data file

rhs <- read_csv("data files/Ms+Cc_RHS_complete_clean.csv", 
                col_types = cols(shock.stage = col_factor(levels = c("control", 
                                                                     "early", "mid", "late"))))
View(rhs)


rhs.vlong <- read_csv("data files/Ms+Cc_RHS_comp_clean_vlong.csv", 
                      col_types = cols(shock.stage = col_factor(levels = c("control", 
                                                                           "early", "mid", "late"))))
View(rhs.vlong)

rhs.long<- read_csv("data files/Ms+Cc_RHS_comp_clean_long.csv", 
                    col_types = cols(shock.stage = col_factor(levels = c("control", 
                                                                         "early", "mid", "late"))))



#------------------------

#analysis of wasp survival by shock stage treatment

#remove early shock stage, since no wasps survived (complete separation)
rhs_ww <- subset(rhs, shock.stage!="early")

#remove individuals that died before the molt to 5th
rhs_ww <- subset(rhs_ww, died.bf5==0)

#make a column for total number of wasps that died for each host
rhs_ww$tot.died <- rhs_ww$tot.load - rhs_ww$num.ecl


#glm model with binomial distribution
#response variable = tot.died vs num.ecl
#predictors = shock.stage and load

waspsurv_mod <- glm(cbind(tot.died, num.ecl) ~ shock.stage * tot.load,
                    family = "quasibinomial",
                    data=rhs_ww,
                    na.action = na.omit)

anova(waspsurv_mod)
summary(waspsurv_mod)



#Very overdispersed, attempting a mixed effects model

#make id number a character
rhs_ww$id <- as.character(rhs_ww$id)

#rescale load
rhs_ww$resc.load <- rescale(rhs_ww$tot.load, to=c(0,1))

#glmm model
waspsurv_re_mod <- glmer(cbind(tot.died, num.ecl) ~ shock.stage * resc.load + (1|id),
                         family = "binomial",
                         data = rhs_ww,
                         na.action = na.omit)
anova(waspsurv_re_mod)
summary(waspsurv_re_mod)



#test fixed effects

waspsurv_re_mod_stage <- glmer(cbind(tot.died, num.ecl) ~ shock.stage + (1|id),
                               family = "binomial",
                               data = rhs_ww,
                               na.action = na.omit)

waspsurv_re_mod_load <- glmer(cbind(tot.died, num.ecl) ~ resc.load + (1|id),
                              family = "binomial",
                              data = rhs_ww,
                              na.action = na.omit)

waspsurv_re_mod_int <- glmer(cbind(tot.died, num.ecl) ~ shock.stage:resc.load + (1|id),
                             family = "binomial",
                             data = rhs_ww,
                             na.action = na.omit)

waspsurv_re_mod_null <- glmer(cbind(tot.died, num.ecl) ~ 1 + (1|id),
                              family = "binomial",
                              data = rhs_ww,
                              na.action = na.omit)

waspsurv_re_mod_add <- glmer(cbind(tot.died, num.ecl) ~ shock.stage + resc.load + (1|id),
                             family = "binomial",
                             data = rhs_ww,
                             na.action = na.omit)

#seems like shock stage and load have effects, but their interaction does not. 
#additive model does not seem significantly different than the full model

anova(waspsurv_re_mod_null, waspsurv_re_mod, waspsurv_re_mod_stage, waspsurv_re_mod_load,
      waspsurv_re_mod_int, waspsurv_re_mod_add)

AIC(waspsurv_re_mod_null, waspsurv_re_mod, waspsurv_re_mod_stage, waspsurv_re_mod_load,
    waspsurv_re_mod_int, waspsurv_re_mod_add)




#-------------------------

#analyze wasp adult mass

#make long dataframe of wasp mass and sex
rhs_wam <- rhs_ww %>% gather(sex, ad.mass, ind.fem.mass, ind.male.mass)
rhs_wam$sex<-gsub("ind.fem.mass", "female", rhs_wam$sex)
rhs_wam$sex<-gsub("ind.male.mass", "male", rhs_wam$sex)


#linear mixed effects model with mass as response, shock stage, sex and load as predictors
waspmss_re_mod <- lme(ad.mass ~ shock.stage * sex * tot.load,
                      random = ~1|id,
                      method = "ML",
                      data = rhs_wam,
                      na.action = na.omit)
anova(waspmss_re_mod)
summary(waspmss_re_mod)







