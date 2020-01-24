#Ms+Cc RHS preliminary figures

#load libraries!

library(readr)
library(ggplot2)
library(tidyr)
library(Rmisc)
library(dplyr)
library(viridis)
library(cowplot)
library(extrafont)

#----------------------------

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



#-----------------------

#prelim plot of adult wasp mass

  ##need to think about how to deal with NAs/0s
  ##unlike survival, which is counts, doesn't make sense to include a 0 mass
  ##but will change the sample size between survival and mass for the means of the different sexes

fem.adwmass.plot<-ggplot(rhs, aes(x=shock.stage, y=ind.fem.mass, color=shock.stage))
fem.adwmass.plot+geom_jitter()

fem.adwmss.bp <- ggplot(rhs, aes(x=shock.stage, y=ind.fem.mass, fill=shock.stage))
fem.adwmss.bp + geom_boxplot()

male.adwmass.plot<-ggplot(rhs, aes(x=shock.stage, y=ind.male.mass, color=shock.stage))
male.adwmass.plot+geom_jitter()


#creating long data set for adult mass, so I can plot male and female together

rhs.wam <- rhs %>% gather(sex, ad.mass, ind.fem.mass, ind.male.mass)
rhs.wam$sex<-gsub("ind.fem.mass", "female", rhs.wam$sex)
rhs.wam$sex<-gsub("ind.male.mass", "male", rhs.wam$sex)

#Removing early treatment, since there are no wasps

rhs.wam.ne<-subset(rhs.wam, shock.stage!="early")



#plotting male and female mass together, by load
  ##color by shock stage

awsm.plot<-ggplot(rhs.wam.ne, aes(x=tot.load, y=ad.mass, group=interaction(shock.stage, sex), color=shock.stage))
awsm.plot+geom_point(aes(shape=sex),
                     size=3
)+geom_smooth(method=lm,
              size=1.2
)+scale_color_manual(values=c("#95D840", "#1F9F88", "#440D54"),
                     name="Shock Treatment",
                     label=c("Control", "Middle", "Late")
)+facet_wrap(~shock.stage)


#plotting male and female mass together, by load
##color by sex

awsm.plot2<-ggplot(rhs.wam, aes(x=tot.load, y=ad.mass, group=sex, color=sex))
awsm.plot2+geom_point(aes(shape=sex),
                      size=3
)+geom_smooth(method=lm,
              size=1.2
)+scale_color_manual(values=c("black", "red"),
                     name="Sex",
                     label=c("Female", "Male")
)+facet_wrap(~shock.stage)



#creating summary of sex mass to directly compare shock treatments

awsm.sum<-summarySE(rhs.wam, measurevar = "ad.mass",
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





#plotting adult mass, separating out late treatment that had a low heat shock

awsm.low.sum<-summarySE(rhs.wam, measurevar = "ad.mass",
                    groupvars = c("shock.stage", "hs.cal", "sex"),
                    na.rm=TRUE)
awsm.low.sum


#plotting mean of wasp adult mass by shock treatment

awsm.low.sum.plot<-ggplot(awsm.low.sum, aes(x=shock.stage, y=ad.mass, group=interaction(sex, hs.cal), color=sex))
awsm.low.sum.plot+geom_point(aes(shape=hs.cal),
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
)+scale_shape_manual(values=c(16, 17, 18),
                     name="heat shock temp")






#-----------------------

#prelim plots of wasp survival to eclosion

#Converting 0s to NAs except for early treatment

rhs$ps.ecl<-ifelse(rhs$shock.stage=="early", 0.000001, rhs$ps.ecl)

rhs$ps.ecl[rhs$ps.ecl==0]<-NA
rhs$ps.ecl[rhs$ps.ecl==0.000001000]<-0


#ps.ecl by load
psecl.plot<-ggplot(rhs, aes(x=tot.load, y=ps.ecl, group=shock.stage, color=shock.stage))
psecl.plot+geom_point(size=3
)+geom_smooth(method=lm, se=FALSE,
              size=1.2
)+scale_y_continuous(limits=c(0,1)
)+scale_color_manual(values=c("#95D840", "#1F9F88", "#3F87BC","#440D54"),
                     name="Shock Stage",
                     label=c("Control", "Early", "Middle", "Late"))


rhs$mass.in.hs<-as.numeric(rhs$mass.in.hs)

rhs.lm<-subset(rhs,shock.stage!="control")
rhs.lm<-subset(rhs.lm, shock.stage!="early")


#ps.ecl by size at heat shock
psecl.plot2<-ggplot(rhs.lm, aes(x=mass.in.hs, y=ps.ecl, group=shock.stage, color=shock.stage))
psecl.plot2+geom_point(
)+geom_smooth(method=lm
)+facet_wrap(~shock.stage)



#mean wasp survival to eclosion

psecl.sum<-summarySE(rhs, measurevar = "ps.ecl",
                     groupvars = "shock.stage",
                     na.rm = TRUE)
psecl.sum

psecl.sum$hs<-c("con", "test", "test", "test")


psecl.sum.plot<-ggplot(psecl.sum, aes(x=shock.stage, y=ps.ecl, color=hs, group=hs))
psecl.sum.plot+geom_point(aes(shape=hs),
                          size=7
)+geom_line(size=2
)+geom_errorbar(aes(ymin=ps.ecl-se, ymax=ps.ecl+se),
                width=.4, size=1.2
)+scale_color_manual(values=c("orange", "black"),
                     name="Treatment",
                     label=c("Control", "Heat Shock")
)+scale_shape_manual(values=c(16, 17),
                     name="Treatment",
                     label=c("Control", "Heat Shock")
)+scale_x_discrete(labels=c("Control", "Early", "Middle", "Late")
)+labs(x="Shock Stage", y="Survival to eclosion (%)"                     
)+theme(text = element_text(family=("Cambria"),face = "bold"),
        strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                        size = 3),
        strip.text = element_text(size=30),
        axis.line.x=element_line(colour = 'black', size = 2),
        axis.line.y=element_line(colour = 'black', size = 2),
        axis.ticks = element_line(colour = 'black', size = 2),
        axis.ticks.length = unit(2.5, "mm"),
        axis.text.x = element_text(size = 26, face = "bold"),
        axis.text.y = element_text(size = 26, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold",
                                    margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 30, face = "bold",
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.key.width=unit(15,"mm"),
        legend.key.height = unit(10,"mm"),
        legend.text=element_text(size=18, face = "bold"),
        legend.title=element_text(size=20, face = "bold"),
        legend.background = element_rect(color="black",linetype="solid",size=2))


#ps.ecl (for dissected hosts) boxplot of raw data

psecl.boxplot<-ggplot(rhs, aes(x=shock.stage, y=ps.ecl, fill=shock.stage))
psecl.boxplot+geom_boxplot()


#plotting ps.ecl by sex
  ##sex by total load

rhs$ps.ecl.totm<-rhs$male.ecl/rhs$tot.load
rhs$ps.ecl.totf<-rhs$fem.ecl/rhs$tot.load

rhs.tot<-gather(rhs, sex, surv, ps.ecl.totf, ps.ecl.totm)

rhs.tot$surv[is.na(rhs.tot$surv)]<-0

pseclsex.plot<-ggplot(rhs.tot, aes(x=shock.stage, y=surv, group=sex, color=sex))
pseclsex.plot+geom_jitter()


#mean ps.ecl by sex

psecl.sex.sum<-summarySE(rhs.tot, measurevar = "surv",
                         groupvars = c("shock.stage", "sex"),
                         na.rm=TRUE)
psecl.sex.sum


pseclsex.sum.plot<-ggplot(psecl.sex.sum, aes(x=shock.stage, y=surv, group=sex, color=sex))
pseclsex.sum.plot+geom_point(aes(shape=sex),
                             size=3
)+geom_line(aes(linetype=sex),
            size=1.2
)+geom_errorbar(aes(ymin=surv-se, ymax=surv+se),
                width=.3, size=1
)+scale_color_manual(values=c("black", "red"),
                     name="Sex",
                     label=c("Female", "Male")
)+scale_shape_manual(values=c(16,17),
                     name="Sex",
                     label=c("Female", "Male")
)+scale_linetype_manual(values=c("solid", "dashed"),
                        name="Sex",
                        label=c("Female", "Male"))



#------------------- 
#plotting mean number of wasps that eclosed, by sex

#subsetting out the low heat shock group (after 11/24/18)
rhs.hot<-subset(rhs, hs.cal!="low")

#subsetting out dead bf5 hosts
rhs.hot<-subset(rhs.hot, died.bf5=="0")

#subsetting hosts that have not had their eclosed wasps counted
rhs.hot$date.ecl.j[is.na(rhs.hot$date.ecl.j)]<-0

rhs.hot$keep.ecl<-ifelse(rhs.hot$shock.stage=="early", 1,
                      ifelse(rhs.hot$date.ecl.j>0 & rhs.hot$num.ecl>0, 1, 0))

rhs.hot<-subset(rhs.hot, keep.ecl=="1")


#Creating long data set for wasp num.ecl by sex
rhs.hot.lng<-gather(rhs.hot, sex, num.ecl2, fem.ecl, male.ecl)

rhs.hot.lng$num.ecl2[is.na(rhs.hot.lng$num.ecl2)]<-0

#Summary of num.ecl
numecl.sum<-summarySE(rhs.hot.lng, measurevar = "num.ecl2",
                      groupvars = c("shock.stage","sex"),
                      na.rm=TRUE)
numecl.sum

numecl.sum$group<-c("con", "con",
                    "eml", "eml",
                    "eml", "eml",
                    "eml", "eml")

#make early shock stage treatment == NA

numecl.sum$num.ecl2[numecl.sum$num.ecl2=="0"]<-NA
numecl.sum$se[numecl.sum$se=="0"]<-NA


#summary of total num of ecl wasps
totecl.sum<-summarySE(rhs.hot, measurevar = "num.ecl",
                      groupvars = "shock.stage",
                      na.rm=TRUE)

totecl.sum
colnames(totecl.sum)[colnames(totecl.sum)=="num.ecl"]<-"num.ecl2"

totecl.sum$group<-c("con", "eml", "eml", "eml")


totecl.sum.plot<-ggplot(totecl.sum, aes(x=shock.stage, y=num.ecl2, group=group, color=group))
totecl.sum.plot+geom_point(aes(shape=group),
                           size=5
)+geom_line(size=1.2
)+geom_errorbar(aes(ymin=num.ecl2-se, ymax=num.ecl2+se),
                width=.3, size=1
)+scale_color_manual(values=c("orange", "black"),
                     name="Treatment",
                     label=c("Control", "Heat Shock")
)+scale_shape_manual(values=c(16, 17),
                     name="Treatment",
                     label=c("Control", "Heat Shock")
)+labs(x="Shock Stage", y="Number of Eclosed Wasps"
)+theme(text = element_text(family=("Cambria")),
        strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                        size = 1),
        strip.text = element_text(size=30),
        axis.line.x=element_line(colour = 'black', size = 1),
        axis.line.y=element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = 'black', size = 1),
        axis.ticks.length = unit(2, "mm"),
        axis.text.x = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        legend.key.width=unit(15,"mm"),
        legend.key.height = unit(10,"mm"),
        legend.text=element_text(size=18),
        legend.title=element_text(size=20),
        legend.background = element_rect(color="black",linetype="solid",size=1))




numecl.plot<-ggplot(numecl.sum, aes(x=shock.stage, y=num.ecl2, group=interaction(group,sex), color=sex))
numecl.plot+geom_point(aes(shape=sex),
                       size=5
)+geom_line(aes(linetype=sex),
            size=1.2
)+geom_errorbar(aes(ymin=num.ecl2-se, ymax=num.ecl2+se),
                width=.3, size=1
)+scale_color_manual(values=c("black", "red"),
                     name="Sex",
                     label=c("Female", "Male")
)+scale_shape_manual(values=c(16,17),
                     name="Sex",
                     label=c("Female", "Male")
)+scale_linetype_manual(values=c("solid", "dashed"),
                        name="Sex",
                        label=c("Female", "Male"))

#----------------------

#plotting the mean proportions of wasp stages by total load in a bar plot (like class prop plot)

#creating column of wasps that emerged, but failed to eclose as adults
rhs$fail.ecl<-rhs$num.em-rhs$num.ecl

rhs.lng<-gather(rhs, wasp.stage, wasp.num, tot.unem, fail.ecl, num.ecl)

rhs.lng$wasp.stage<-factor(rhs.lng$wasp.stage, levels=c("tot.unem", "fail.ecl", "num.ecl"))


wstage.sum<-summarySE(rhs.lng, measurevar = "wasp.num",
                      groupvars = c("shock.stage", "wasp.stage"),
                      na.rm = TRUE)
wstage.sum


wstage.sum.plot<-ggplot(wstage.sum, aes(x=shock.stage, y=wasp.num, group=wasp.stage, color=wasp.stage))
wstage.sum.plot+geom_point(
)+geom_line()

#calculating mean load from mean wasp stage numbers
mncon.load<-wstage.sum[1,4]+wstage.sum[2,4]+wstage.sum[3,4]
mnmid.load<-wstage.sum[7,4]+wstage.sum[8,4]+wstage.sum[9,4]
mnlate.load<-wstage.sum[10,4]+wstage.sum[11,4]+wstage.sum[12,4]

#removing early data (from 1-2 hosts, not representative)
wstage.sum[4,4]<-0

#adding mean load to wstage.sum

wstage.sum$mn.load<-c(mncon.load,mncon.load,mncon.load,
                      0,0,0,
                      mnmid.load, mnmid.load, mnmid.load,
                      mnlate.load,mnlate.load,mnlate.load)

#calculating mean prop for each stage
wstage.sum$prop<-wstage.sum$wasp.num/wstage.sum$mn.load

wstage.sum$prop[is.nan(wstage.sum$prop)]<-0


#plotting mean prop in bar plot for each shock stage treatment

wstage.plot<-ggplot(wstage.sum,aes(x=shock.stage, y=prop, fill=wasp.stage))
wstage.plot+geom_bar(position="fill",stat="identity"
)+scale_fill_manual(values=c("#95D840", "#1F9F88", "#440D54"),
                    breaks=c("tot.unem", "fail.ecl", "num.ecl"),
                    labels=c("Unemerged", "Emerged", "Eclosed"),
                    name="Stage"
)+labs(x="Shock Stage", y="Proportion"
)+theme(text = element_text(family=("Cambria")),
        strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                        size = 1),
        strip.text = element_text(size=30),
        axis.line.x=element_line(colour = 'black', size = 1),
        axis.line.y=element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = 'black', size = 1),
        axis.ticks.length = unit(2, "mm"),
        axis.text.x = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        legend.key.width=unit(15,"mm"),
        legend.key.height = unit(10,"mm"),
        legend.text=element_text(size=18),
        legend.title=element_text(size=20),
        legend.background = element_rect(color="black",linetype="solid",size=1))





#---------------
#plotting sex ratio by shock stage--does not seem different between treatments

rhs$fem.ecl[is.na(rhs$fem.ecl)]<-0
rhs$male.ecl[is.na(rhs$male.ecl)]<-0

rhs$fem.rat<-rhs$fem.ecl/rhs$male.ecl

rhs$fem.rat[is.nan(rhs$fem.rat)]<-0

femrat.plot<-ggplot(rhs, aes(x=shock.stage, y=fem.rat, color=shock.stage))
femrat.plot+geom_jitter()

femrat.boxplot<-ggplot(rhs, aes(x=shock.stage, y=fem.rat, fill=shock.stage))
femrat.boxplot+geom_boxplot()


#-----------------------

#prelim plot of number of wasps emerged

rhs$num.em[rhs$num.em==0]<-NA
rhs<-subset(rhs, tot.load<280)

#remove ind that have not been dissected

rhs$tot.unem[is.na(rhs$tot.unem)]<-0
rhs.dis<-subset(rhs, tot.unem!=0)

#by load

nem.plot<-ggplot(rhs.dis, aes(x=tot.load, y=num.em, group=shock.stage, color=shock.stage))
nem.plot+geom_point(
)+geom_smooth(method=lm)


#by date of test

nem.plot2<-ggplot(rhs, aes(x=date.in.hs.j, y=num.em, group=shock.stage, color=shock.stage))
nem.plot2+geom_point(aes(shape=hs.cal)
)+geom_smooth(method=lm)

#subset out those shocked after 11/24/18, when chamber was recalibrated

rhs.2<-subset(rhs, hs.cal!="low")

nem.plot3<-ggplot(rhs.2, aes(x=date.in.hs.j, y=num.em, group=shock.stage, color=shock.stage))
nem.plot3+geom_point(
)+geom_smooth(method=lm)


#-------------------

#calculating the percent of mature and immature 2nds

rhs$p.unem.mat<-rhs$num.unem.mat/rhs$tot.unem
rhs$p.unem.im<-rhs$num.unem.im/rhs$tot.unem

unem.plot<-ggplot(rhs, aes(x=tot.load, y=p.unem.mat, group=shock.stage, color=shock.stage))
unem.plot+geom_point(
)+geom_smooth(method=lm)

unem.plot2<-ggplot(rhs, aes(x=tot.load, y=p.unem.im, group=shock.stage, color=shock.stage))
unem.plot2+geom_point(
)+geom_smooth(method=lm)


#Plotting distribution of % mature and immature 2nds

distmat.plot<-ggplot(rhs, aes(x=p.unem.mat, fill=shock.stage))
distmat.plot+geom_density(alpha=.5)

distim.plot<-ggplot(rhs, aes(x=p.unem.im, fill=shock.stage))
distim.plot+geom_density(alpha=.5)


#--------------

#plotting development time from molt to 5th

rhs$dt5<-rhs$date.em.j-rhs$date.5.j
rhs$dtcull<-rhs$date.cull.j-rhs$date.5.j

rhs$dt5<-coalesce(rhs$dt5, rhs$dtcull)


dt5.plot<-ggplot(rhs, aes(x=shock.stage, y=dt5, color=shock.stage))
dt5.plot+geom_jitter()

#creating summary of length of 5th instar

dt5.sum<-summarySE(rhs, measurevar = "dt5",
                   groupvars = c("shock.stage"),
                   na.rm=TRUE)
dt5.sum

#plot of mean length of the 5th instar

mn.dt5.plot<-ggplot(dt5.sum, aes(x=shock.stage, y=dt5, color=shock.stage))
mn.dt5.plot+geom_point(size=5
)+geom_errorbar(aes(ymin=dt5-se, ymax=dt5+se),
                width=.3, 
                size=1
)+labs(x="Shock Stage", y="Time to emergence from 5th")



#plotting time to emergence for host lifespan

ttemh.plot<-ggplot(rhs, aes(x=shock.stage, y=ttem.h, color=shock.stage))
ttemh.plot+geom_jitter()

ttemh.sum<-summarySE(rhs, measurevar = "ttem.h",
                     groupvars = "shock.stage",
                     na.rm = TRUE)
ttemh.sum

ttemhsum.plot<-ggplot(ttemh.sum, aes(x=shock.stage, y=ttem.h, color=shock.stage))
ttemhsum.plot+geom_point(size=3
)+geom_errorbar(aes(ymin=ttem.h-se, ymax=ttem.h+se)
)+labs(y="Time from Hatching to Emergence", x="Shock Stage")



#plotting time to emergence for wasp life span (ovp to emergence)

ttemw.plot<-ggplot(rhs, aes(x=shock.stage, y=ttem.w, color=shock.stage))
ttemw.plot+geom_jitter()


ttemw.sum<-summarySE(rhs, measurevar = "ttem.w",
                     groupvars = "shock.stage",
                     na.rm = TRUE)
ttemw.sum


ttemwsum.plot<-ggplot(ttemw.sum, aes(x=shock.stage, y=ttem.w, color=shock.stage))
ttemwsum.plot+geom_point(size=3
)+geom_errorbar(aes(ymin=ttem.w-se, ymax=ttem.w+se)
)+labs(y="Time from Oviposition to Emergence", x="Shock Stage")



#Boxplot of dt5

dt5.bxplot<-ggplot(rhs, aes(y=dt5, x=shock.stage, fill=shock.stage))
dt5.bxplot+geom_boxplot()


#plotting number emerged by dt5
ne.dt5.plot<-ggplot(rhs, aes(x=dt5, y=num.em, group=shock.stage, color=shock.stage))
ne.dt5.plot+geom_jitter(
)+geom_smooth(method=lm)

#plotting ps.em by dt5
psem.dt5.plot<-ggplot(rhs, aes(x=dt5, y=ps.em, group=shock.stage, color=shock.stage))
psem.dt5.plot+geom_jitter(
)+geom_smooth(method=lm)


#plotting ps.ecl by dt5
psecl.dt5.plot<-ggplot(rhs, aes(x=dt5, y=ps.ecl, group=shock.stage, color=shock.stage))
psecl.dt5.plot+geom_jitter(
)+geom_smooth(method=lm)


#------------------------------

#plotting mass by age for parasitized hosts

#make mass numeric
rhs.long$mass <- as.numeric(rhs.long$mass)

#removing an individual with a typo in mass

rhs.long<-subset(rhs.long, mass<20000)

#taking log of mass
rhs.long$log.mass<-log(rhs.long$mass)


#individual mass
am.plot<-ggplot(rhs.long, aes(x=age, y=log.mass, group=interaction(id, shock.stage), color=shock.stage))
am.plot+geom_point(aes(shape=instar)
)+geom_line(
)+facet_wrap(~shock.stage)


#constructing summarySE for age and mass to make mean plot

lm.sum<-summarySE(rhs.long, measurevar = "log.mass",
                  groupvars = c("shock.stage", "instar"),
                  na.rm=TRUE)
lm.sum


age.sum<-summarySE(rhs.long, measurevar = "age",
                   groupvars = c("shock.stage", "instar"),
                   na.rm = TRUE)
age.sum


lm.sum$age<-age.sum[,4]
lm.sum$age.se<-age.sum[,6]

lm.sum.em<-subset(lm.sum, N>10)

ammean.plot<-ggplot(lm.sum.em, aes(x=age, y=log.mass, group=shock.stage, color=shock.stage))
ammean.plot+geom_point(size=3
)+geom_line(size=1.2
)+geom_errorbar(aes(ymin=log.mass-se, ymax=log.mass+se),
                width=.3, size=1
)+geom_errorbarh(aes(xmin=age-age.se, xmax=age+age.se),
                 height=.3, size=1
)+scale_color_manual(values=c("#95D840", "#1F9F88", "#3F87BC","#440D54"),
                     name="Shock Stage",
                     label=c("Control", "Early", "Middle", "Late"))


#----------------------------

#Plotting % of em and mongo for each shock stage

#remove dead individuals
rhs.cl<-subset(rhs, died.bf5==0)

#removing control individual that molted poorly and died, but was recorded as cull instead of dead
rhs.cl<-subset(rhs.cl, id!="232")

#creating column "class"--em, mongo, wander
  ##some individuals are missing date.em info, try and locate to fill in missing info

rhs.cl$date.em.j[is.na(rhs.cl$date.em.j)]<-0
rhs.cl$date.cull.j[is.na(rhs.cl$date.cull.j)]<-0
rhs.cl$date.wand.j[is.na(rhs.cl$date.wand.j)]<-0
rhs.cl$num.em[is.na(rhs.cl$num.em)]<-0


rhs.cl$class<-ifelse(rhs.cl$date.em.j>0 | rhs.cl$num.em>0, "em",
                     ifelse(rhs.cl$date.cull.j>0, "mongo",
                            ifelse(rhs.cl$date.wand.j>0, "wand", "unk")))

check<-rhs.cl[,c("id","shock.stage", "date.em.j", "date.cull.j", "date.wand.j", "class")]
View(check)

#determining num of each class in each treatment, and the total number in each treatment

n.table<-table(rhs.cl$shock.stage, rhs.cl$class)
tot.table<-table(rhs.cl$shock.stage)

#converting tables to data frames

n.table<-data.frame(n.table)
tot.table<-data.frame(tot.table)

#renaming columns in tables
n.table<-dplyr::rename(n.table, shock.stage=Var1, class=Var2)
tot.table<-dplyr::rename(tot.table, shock.stage=Var1)

#creating a column in n.table that has the total numbers for each treatment from tot.table

n.table$tot<-ifelse(n.table$shock.stage=="control", tot.table[which(tot.table$shock.stage=="control"), 2],
                    ifelse(n.table$shock.stage=="early", tot.table[which(tot.table$shock.stage=="early"), 2],
                           ifelse(n.table$shock.stage=="mid", tot.table[which(tot.table$shock.stage=="mid"), 2],
                                  ifelse(n.table$shock.stage=="late", tot.table[which(tot.table$shock.stage=="late"), 2], 0))))


#calculating the proportion of each class in each shock stage treatment

n.table<-mutate(n.table, prop = Freq/tot)

#subsetting out the unknown classes--due to missing data

n.table<-subset(n.table, class!="unk")


#Plotting the prop of each class in each shock stage treatment

class.plot<-ggplot(n.table,aes(x=shock.stage, y=prop, fill=class))
class.plot+geom_bar(position="fill",stat="identity"
)+scale_fill_manual(values=c("#95D840", "#1F9F88", "#440D54"),
                    breaks=c("em", "mongo", "wand"),
                    labels=c("Emergence", "WOWE", "Wandering"),
                    name="Outcome"
)+labs(x="Shock Stage", y="Proportion"
)+scale_x_discrete(labels=c("Control", "Early", "Middle", "Late")                   
)+theme(text = element_text(family=("Cambria"),face = "bold"),
        strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                        size = 3),
        strip.text = element_text(size=30),
        axis.line.x=element_line(colour = 'black', size = 2),
        axis.line.y=element_line(colour = 'black', size = 2),
        axis.ticks = element_line(colour = 'black', size = 2),
        axis.ticks.length = unit(2.5, "mm"),
        axis.text.x = element_text(size = 26, face = "bold"),
        axis.text.y = element_text(size = 26, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold",
                                    margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 30, face = "bold",
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.key.width=unit(15,"mm"),
        legend.key.height = unit(10,"mm"),
        legend.text=element_text(size=18, face = "bold"),
        legend.title=element_text(size=20, face = "bold"),
        legend.background = element_rect(color="black",linetype="solid",size=2))

View(n.table)








