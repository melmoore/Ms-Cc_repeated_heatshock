#Ms Cc repeated heat shock (developmental timing heat shock) preliminary analyses 


#load libraries
library(scales)
library(readr)
library(ggplot2)
library(tidyr)
library(Rmisc)
library(dplyr)
library(nlme)
library(lme4)
library(mgcv)
library(MuMIn)


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

rhs_long<- read_csv("data files/Ms+Cc_RHS_comp_clean_long.csv", 
                    col_types = cols(shock.stage = col_factor(levels = c("control", 
                                                                         "early", "mid", "late"))))



#------------------------

#analysis of wasp survival by shock stage treatment

#remove early shock stage, since no wasps survived (complete separation)
rhs_ww <- subset(rhs, shock.stage!="early")

#remove individuals that died before the molt to 5th
rhs_ww <- subset(rhs_ww, died.bf5==0)

#drop rows with NAs in tot.load--1 that did not have num_em recorded, 2 that got lost in freezer and
#not dissected
rhs_ww <- drop_na(rhs_ww, tot.load)

#remove individuals with load greater than 300
rhs_ww <- subset(rhs_ww, tot.load<300)

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



#model selection using dredge() 

#dredge requires dataframe with no NAs--subsetting to only columns in the model
rhs_wwd <- select(rhs_ww, id, shock.stage, resc.load, tot.died, num.ecl)
rhs_wwd <- drop_na(rhs_wwd)

waspsurv_re_mod_fd <- glmer(cbind(tot.died, num.ecl) ~ shock.stage * resc.load + (1|id),
                         family = "binomial",
                         data = rhs_wwd,
                         na.action = na.fail)


ws_mod_dredge <- dredge(waspsurv_re_mod_fd)
ws_mod_dredge


#best model according to dredge
ws_mod_d <- glmer(cbind(tot.died, num.ecl) ~ shock.stage + (1|id),
                     family = "binomial",
                     data = rhs_ww,
                     na.action = na.omit) 



#test fixed effects

ws_re_mod_ss <- glmer(cbind(tot.died, num.ecl) ~ + (1|id),
                         family = "binomial",
                         data = rhs_ww,
                         na.action = na.omit)



anova(ws_mod_d_ld, ws_re_mod_ss)




#-------------------------------

#GLMM binomial analysis of wasp survival to emergence

wemsurv_re_mod <- glmer(cbind(num.em, tot.unem) ~ shock.stage * resc.load + (1|id),
                        family = "binomial",
                        data = rhs_ww,
                        na.action = na.omit)

anova(wemsurv_re_mod)
summary(wemsurv_re_mod)


#model selection using dredge() 

#dredge requires dataframe with no NAs--subsetting to only columns in the model
rhs_wwd_em <- select(rhs_ww, id, shock.stage, resc.load, num.em, tot.unem)
rhs_wwd_em <- drop_na(rhs_wwd_em)

wemsurv_re_mod_fd <- glmer(cbind(num.em, tot.unem) ~ shock.stage * resc.load + (1|id),
                            family = "binomial",
                            data = rhs_wwd_em,
                            na.action = na.fail)

wems_mod_dredge <- dredge(wemsurv_re_mod_fd)
wems_mod_dredge


#best model according to dredge
wems_mod_d <- glmer(cbind(num.em, tot.unem) ~ shock.stage + (1|id),
                    family = "binomial",
                    data = rhs_ww,
                    na.action = na.omit)



wems_mod_ss <- glmer(cbind(num.em, tot.unem) ~  (1|id),
                     family = "binomial",
                     data = rhs_ww,
                     na.action = na.omit)

anova(wems_mod_d, wems_mod_ss)


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



#model selection using dredge() 

#dredge requires dataframe with no NAs--subsetting to only columns in the model
rhs_wamd <- select(rhs_wam, id, shock.stage, sex, tot.load, ad.mass)
rhs_wamd <- drop_na(rhs_wamd)


wmss_mod_d <- lme(ad.mass ~ shock.stage * sex * tot.load,
                       random = ~1|id,
                       method = "ML",
                       data = rhs_wamd,
                       na.action = na.fail)

wmss_mod_dredge <- dredge(wmss_mod_d)
wmss_mod_dredge

#best model is full model



#raw data plot
awsm.plot<-ggplot(rhs_wam, aes(x=tot.load, y=ad.mass, group=interaction(shock.stage, sex), color=sex))
awsm.plot+geom_point(aes(shape=sex),
                     size=3
)+geom_smooth(method=lm,
              size=1.2
)+scale_color_manual(values=c("orange", "black"),
                     name="Sex",
                     label=c("Female", "Male")
)+facet_wrap(~shock.stage)



#test model fit

#create plotting data set by removing rows with NA in wasp adult mass
rhs_mod <- drop_na(rhs_wam, ad.mass)

rhs_mod$pred <- predict(waspmss_re_mod, level=0)
rhs_mod$resid <- residuals(waspmss_re_mod, level=0)

wspmss_pr_plot <- ggplot(rhs_mod, aes(x=pred, y=resid, color=shock.stage))
wspmss_pr_plot + geom_point(shape=1
)+geom_hline(aes(yintercept=0),
             color="black", linetype="dashed")


wspmss_rl_plot <- ggplot(rhs_mod, aes(x=resc.load, y=resid, color=shock.stage))
wspmss_rl_plot + geom_point(shape=1
)+geom_hline(aes(yintercept=0),
             color="black", linetype="dashed")


#pred lines not working for some reason
wspmss_pred_plot <- ggplot(rhs_mod, aes(x=resc.load, y=ad.mass, color=shock.stage))
wspmss_pred_plot + geom_point(size=5, shape=1
)+geom_line(data=rhs_mod, aes(x=resc.load, y=pred)
)+facet_wrap(~shock.stage)



#-----------------------

#investigating the effects of the lower late heat shock on wasp adult mass

#try removing low late heat shock, see what effects that has
rhs_wamh <- subset(rhs_wam, hs.cal!="low")


#linear mixed effects model with mass as response, shock stage, sex and load as predictors
#removing the low late heatshock does not seem to have an effect
waspmss_re_mod_hot <- lme(ad.mass ~ shock.stage * sex * tot.load,
                      random = ~1|id,
                      method = "ML",
                      data = rhs_wamh,
                      na.action = na.omit)

anova(waspmss_re_mod_hot)
summary(waspmss_re_mod_hot)



#creating summary of sex mass to directly compare shock treatments

awsm.sum<-summarySE(rhs_wamh, measurevar = "ad.mass",
                    groupvars = c("shock.stage", "sex"),
                    na.rm=TRUE)
awsm.sum


#plotting mean of wasp adult mass by shock treatment

awsm.sum.plot<-ggplot(awsm.sum, aes(x=shock.stage, y=ad.mass, group=sex, color=sex))
awsm.sum.plot+geom_point(aes(shape=sex),
                         size=5
)+geom_line(aes(linetype=sex),
            size=1.2
)+geom_errorbar(aes(ymin=ad.mass-se, ymax=ad.mass+se),
                width=.3, size=1
)+scale_color_manual(values=c("black","red"),
                     name="Sex",
                     label=c("Female", "Male")
)+scale_linetype_manual(values=c("solid", "dashed"),
                        name="Sex",
                        label=c("Female", "Male")
)+scale_shape_manual(values=c(16, 17),
                     name="Sex",
                     label=c("Female", "Male"))



#analyze mass for just females, without those exposed to low late shock
#shows same patterns--wasps are bigger than control at mid, and more affected by tot.load
#control and late shock don't seem to differ
rhs_wwh <- subset(rhs_ww, hs.cal!="low")


femmss_re_mod <- lme(ind.fem.mass ~ shock.stage * tot.load,
                     random = ~1|id,
                     data=rhs_wwh,
                     method = "ML",
                     na.action = na.omit)
anova(femmss_re_mod)
summary(femmss_re_mod)


#----------------------------------

#analyzing the effects of shock stage and load on host mass and development time 

#LINEAR MODEL

#creating a long dataframe that has an instar designation "end", which indicates mass and age
#at the end of dev, regradless of if they wandered, had emergence, or were culled

rhs_long$instar2 <- ifelse(rhs_long$instar=="wand" | rhs_long$instar=="em" |
                             rhs_long$instar=="cull", "end", rhs_long$instar)


#drop rows with NAs in mass (will remove individuals with missing data, and remove extra rows
#created by the wand, em, cull instar system I made originally)
rhs_long <- drop_na(rhs_long, mass)

#remove individuals that died before molting to 5th
rhs_long <- subset(rhs_long, died.bf5==0)

#remove individuals with load greater than 300
rhs_long <- subset(rhs_long, tot.load<300)


#make sure necessary columns are in the correct class format
rhs_long$mass <- as.numeric(rhs_long$mass)
rhs_long$id <- as.character(rhs_long$id)

#convert mass to log_mass
rhs_long$log_mass <- log(rhs_long$mass)


##log.mass is numeric, age is integer, shock stage is factor, load is numeric
##random effect is linear age by individual (character)
##syntax of (x+y+z)^2 removes the 4 way interactions, leaving only 3 way interactions with the different age terms

lms_mod1<-lme(log_mass~(age+I(age^2)):(shock.stage * tot.load),
              random=~age|id,
              data=rhs_long,
              na.action=na.omit,
              method="ML")

anova(lms_mod1)
summary(lms_mod1)


#testing model fit
rhs_modl <- rhs_long
rhs_modl <- drop_na(rhs_modl, age)

rhs_modl$pred <- predict(lms_mod1, level=0)
rhs_modl$resid <- residuals(lms_mod1, level=0)


#predicted values against actual values
lmss_lmfit <- ggplot(rhs_modl, aes(x=age, y=log_mass, group=shock.stage, color=shock.stage))
lmss_lmfit + geom_point(shape=1, size=3
)+geom_line(data=rhs_modl, aes(x=age, y=pred, group=shock.stage, color=shock.stage)
)+facet_wrap(~shock.stage)


#residuals against predicted values
lmss_pr_plot <- ggplot(rhs_modl, aes(x=pred, y=resid, group=shock.stage, color=shock.stage))
lmss_pr_plot+geom_point(shape=1, size=3
)+geom_hline(aes(yintercept=0),
             color="black", linetype="dashed", size=1.5)

#residuals against age
lmss_ra_plot <- ggplot(rhs_modl, aes(x=age, y=resid, group=shock.stage, color=shock.stage))
lmss_ra_plot+geom_point(shape=1, size=3
)+geom_hline(aes(yintercept=0),
             color="black", linetype="dashed", size=1.5)




#-------------------------

#analyzing the effects of shock stage and on host mass and development time--not including load, since it has no impact on early treatment group
#will do separate analysis on other shock groups to see the effect of load

#GAMM MODEL

#create log_mass column
rhs_long$log_mass <- log(rhs_long$mass)

#remove individuals that died
rhs_lngc <- subset(rhs_long, died.bf5==0)

#subset to only terms used in the model, drop NAs
rhs_mod <- select(rhs_long, id, shock.stage, instar, log_mass, age)
rhs_mod <- drop_na(rhs_mod)

#make id a factor so it will work as a random effect in the GAMM model
rhs_mod$id<-as.factor(rhs_mod$id)


#run a full GAMM model
gam_mass_mod<-gam(log_mass ~ s(age, by= shock.stage, k=30, bs="ts") + s(id, bs="re") + shock.stage,
                  method="ML", data=rhs_mod, na.action = na.omit)

anova(gam_mass_mod)
summary(gam_mass_mod)



#test model fit
rhs_mod$gam_pred <- predict(gam_mass_mod, level=0)
rhs_mod$gam_resid <- residuals(gam_mass_mod, level=0)



#gam model predict fit against actual data
gam_pred_fit <- ggplot(rhs_mod, aes(x=age, y=log_mass, group = id, color=shock.stage))
gam_pred_fit + geom_point(shape=1, size=3
)+geom_line(data=rhs_mod, aes(x=age, y=gam_pred)
)+facet_wrap(~shock.stage)


#gam model pred against resid
gam_pr_plot <- ggplot(rhs_mod, aes(x=gam_pred, y=gam_resid, color=shock.stage))
gam_pr_plot + geom_point(size=3, shape=1
)+geom_hline(aes(yintercept=0),
             size=1.5, color="black", linetype="dashed")



#gam model residuals against age
gam_ra_plot <- ggplot(rhs_mod, aes(x=age, y=gam_resid, color=shock.stage))
gam_ra_plot+geom_point(size=3, shape=1
)+geom_hline(aes(yintercept=0),
             size=1.5, color="black", linetype="dashed")



gam.check(gam_mass_mod)



#---------------------

#removing early shock stage so I can examine the effect of load on host growth
rhs_ne <- subset(rhs_lngc, shock.stage!="early")

#select columns only used in model, drop NAs
rhs_ne <- select(rhs_ne, id, shock.stage, tot.load, log_mass, age)
rhs_ne <- drop_na(rhs_ne)

#make id a factor so it works as a random effect in gamm model
rhs_ne$id <- as.factor(rhs_ne$id)


#run a full GAMM model--load and age in same smooth
gam_ml_mod<-gam(log_mass ~ s(age, tot.load,  by= shock.stage, bs="ts") + s(id, bs="re") + shock.stage,
                  method="ML", data=rhs_ne, na.action = na.omit)

anova(gam_ml_mod)
summary(gam_ml_mod)


#run full gamm model--load and age in separate smooths
gam_ml_noint_mod <- gam(log_mass ~ s(age, by=shock.stage, bs="ts")
                        + s(tot.load, by = shock.stage, bs="ts") + s(id, bs="re")
                        + shock.stage, method="ML", data=rhs_ne, na.action = na.omit)

anova(gam_ml_noint_mod)
summary(gam_ml_noint_mod)

anova(gam_ml_mod, gam_ml_noint_mod, test="Chisq")
AIC(gam_ml_mod, gam_ml_noint_mod) #no interaction model is much better


#test model residuals
rhs_ne$pred <- predict(gam_ml_noint_mod, level=0)
rhs_ne$resid <- residuals(gam_ml_noint_mod, level=0)



#gam model predict fit against actual data
gamml_pred_fit <- ggplot(rhs_ne, aes(x=age, y=log_mass, group = id, color=shock.stage))
gamml_pred_fit + geom_point(shape=1, size=3
)+geom_line(data=rhs_ne, aes(x=age, y=pred)
)+facet_wrap(~shock.stage)



#gam model pred against resid
gamml_pr_plot <- ggplot(rhs_ne, aes(x=pred, y=resid, color=shock.stage))
gamml_pr_plot + geom_point(size=3, shape=1
)+geom_hline(aes(yintercept=0),
             size=1.5, color="black", linetype="dashed")



#gam model residuals against age
gamml_ra_plot <- ggplot(rhs_ne, aes(x=age, y=resid, color=shock.stage))
gamml_ra_plot+geom_point(size=3, shape=1
)+geom_hline(aes(yintercept=0),
             size=1.5, color="black", linetype="dashed")


#gam model residuals against load
gamml_rl_plot <- ggplot(rhs_ne, aes(x=tot.load, y=resid, color=shock.stage))
gamml_rl_plot+geom_point(size=3, shape=1
)+geom_hline(aes(yintercept=0),
             size=1.5, color="black", linetype="dashed")






#adjusting k values for age smooth

#run full gamm model--load and age in separate smooths
gam_ml_nointk_mod <- gam(log_mass ~ s(age, by=shock.stage, k=20, bs="ts")
                        + s(tot.load, by = shock.stage, bs="ts") + s(id, bs="re")
                        + shock.stage, method="ML", data=rhs_ne, na.action = na.omit)

gam.check(gam_ml_nointk_mod)


#----------------------------

#Looking at knot specification and degrees of freedom for the GAMM model with age and load in 
#separate smooths


#GAMM model, k=10
gam_mass_mod_noint <- gam(log_mass ~ s(age, by=shock.stage, k=10, bs="ts")
                          + s(tot.load, by=shock.stage, k=10, bs="ts") + s(id, bs="re")
                          + shock.stage, method="ML", data=rhs_modl, na.action = na.omit)

gam.check(gam_mass_mod_noint)


#trying to follow the troubleshooting example in the ?choose.k help page
plot(gam_mass_mod_noint,pages=1,residuals=TRUE) 


#test residuals against each smooth to see if one needs to have the k adjusted
gam(gam_resid ~ s(age, by=shock.stage, k=20, bs="ts")+ s(id, bs="re")
    + shock.stage, method="ML", data=rhs_modl, na.action = na.omit)

gam(gam_resid ~ s(tot.load, by=shock.stage, k=20, bs="ts") + s(id, bs="re")
    + shock.stage, method="ML", data=rhs_modl, na.action = na.omit)



#GAMM model, k=20 for both smooths
#smooth with load looks ok, did not change. Age improved, try increasing the k again.
gammss_mod_ni_1 <- gam(log_mass ~ s(age, by=shock.stage, k=20, bs="ts")
                          + s(tot.load, by=shock.stage, k=20, bs="ts") + s(id, bs="re")
                          + shock.stage, method="ML", data=rhs_modl, na.action = na.omit)

gam.check(gammss_mod_ni_1)


#GAMM model, k=30 for age, k=10 for load
gammss_mod_ni_2 <- gam(log_mass ~ s(age, by=shock.stage, k=30, bs="ts")
                       + s(tot.load, by=shock.stage, k=10, bs="ts") + s(id, bs="re")
                       + shock.stage, method="ML", data=rhs_modl, na.action = na.omit)

gam.check(gammss_mod_ni_2)



#GAMM model, k=40 for age, k=10 for load--this one throws an error message
gammss_mod_ni_3 <- gam(log_mass ~ s(age, by=shock.stage, k=40, bs="ts")
                       + s(tot.load, by=shock.stage, k=10, bs="ts") + s(id, bs="re")
                       + shock.stage, method="ML", data=rhs_modl, na.action = na.omit)

gam.check(gammss_mod_ni_3)


#----------------

#looking at whether incidence of wanderers differs between treatments

#remove individuals that died before molt to 5th
rhs_cl <- subset(rhs, died.bf5==0)

#create "class" column indicating whether the individual wandered, had emergence, or was a WOWE
rhs_cl$date.em.j[is.na(rhs_cl$date.em.j)]<-0
rhs_cl$date.cull.j[is.na(rhs_cl$date.cull.j)]<-0
rhs_cl$date.wand.j[is.na(rhs_cl$date.wand.j)]<-0
rhs_cl$num.em[is.na(rhs_cl$num.em)]<-0


rhs_cl$class<-ifelse(rhs_cl$date.em.j>0 | rhs_cl$num.em>0, "em",
                     ifelse(rhs_cl$date.cull.j>0, "mongo",
                            ifelse(rhs_cl$date.wand.j>0, "wand", "unk")))

#remove individuals with class "unk"
rhs_cl <- subset(rhs_cl, class!="unk")

#create binary column where 1 indicates wandering and 0 indicates other outcome (em or WOWE)
rhs_cl$bin_wand <- ifelse(rhs_cl$class=="wand", 1, 0)


#glm analyzing incidence of wandering between shock stages--if I'm interpreting this correctly, there are
#fewer incidents of wanderers at the control than the shocked treatments
wand_mod <- glm(bin_wand ~ shock.stage,
                family = binomial,
                data = rhs_cl,
                na.action = na.omit)

summary(wand_mod)


#---------------

#removing wanderers for a quick and dirty analysis of outcomes for the other treatments
rhs_nw <- subset(rhs_cl, class!="wand")


#creating a binary class column--emergence = 1, WOWE = 0
rhs_nw$bin_class <- ifelse(rhs_nw$class=="em", 1, 0)


#glm analyzing occurence of emergence vs WOWE between shock stage treatments
class_mod <- glm(bin_class ~ shock.stage,
                 family = binomial,
                 data = rhs_nw,
                 na.action = na.omit)

summary(class_mod)














