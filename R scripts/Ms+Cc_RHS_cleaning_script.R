#Ms+Cc RHS cleaning script

#Load libraries

library(readr)
library(ggplot2)
library(tidyr)
library(Rmisc)
library(dplyr)
library(viridis)
library(cowplot)
library(extrafont)

#-----------------------------


#load data files

rhs <- read_csv("data files/Ms+Cc_RHS_incomp_12-21-18.csv", 
                col_types = cols(shock.stage = col_factor(levels = c("control", 
                                                                     "early", "mid", "late"))))
View(rhs)

#---------------------


#removing empty rows at bottom of data sheet

rhs$id[is.na(rhs$id)]<-0

rhs<-subset(rhs, id>0)


#-----------------------

#Transform date columns into julian date

##Converts x into julian date
j.date<-function(x){
  strptime(x, "%m/%d")$yday+1
}


#Takes all columns that have "date." in the name, and converts contents to Julian day using j.date function. Renames columns (adds a 
##j to end of column name), and binds the out put julian day columns to the original data set

lapj.date<-function(df){
  date.j<-lapply(df[,grep("date.",colnames(df))],j.date)
  date.j<-as.data.frame(date.j)
  colnames(date.j)<-paste(colnames(date.j), "j", sep = ".")
  output.df<-cbind(df,date.j)
  output.df
}


rhs<-lapj.date(rhs)

#creating sorting column for those that died before the molt to 5th

rhs$date.died.j[is.na(rhs$date.died.j)]<-0
rhs$died.bf5<-ifelse(rhs$date.died.j==0, "0", "1")


#------------------------

#converting time columns to decimal time

#Function that turns turns time (x) into a character, splits it at the :, and adds it together to get decimal time
dec.time<-function(x) {
  x<-as.character(x)
  sapply(strsplit(x,":"),function(x){
    x <- as.numeric(x)
    y<-x[1]+x[2]/60
    
  })
}


#Function that applies the dec.time function to every column with "time." in the name, and adds decminal time columns to 
##dataframe
dec.time.col<-function(df){
  dct<-lapply(df[,grep("time.",colnames(df))],dec.time)
  dct<-as.data.frame(dct)
  colnames(dct)<-paste(colnames(dct), "dec", sep = ".")
  output.df<-cbind(df,dct)
  output.df
}

rhs<-dec.time.col(rhs)





#--------------------------

#Converting the blank out put of formula columns to NA

rhs[rhs=="#N/A"]<-NA


#---------------------

#For "lates" that had emergence in the 4th instar (and therefor before being shocked), changing treatment to "control"

rhs$instar.em[is.na(rhs$instar.em)]<-0

#ifelse returns numeric instead of factor, so have to specify that output is a factor with factor()
rhs$shock.stage<-factor(ifelse(rhs$shock.stage=="late" & rhs$instar.em==4, "control", rhs$shock.stage))

#rename factor levels to be the treatment names
rhs$shock.stage<-revalue(rhs$shock.stage, c("1"="control", "2"="early", "3"="mid", "4"="late"))


#----------------------

#Calculating development times

rhs$tt3<-rhs$date.3.j-rhs$date.hatch.j
rhs$tt4<-rhs$date.4.j-rhs$date.hatch.j
rhs$tt5<-rhs$date.5.j-rhs$date.hatch.j
rhs$ttp5.1<-rhs$date.p5.1.j-rhs$date.hatch.j
rhs$ttp5.2<-rhs$date.p5.2.j-rhs$date.hatch.j
rhs$ttp5.3<-rhs$date.p5.3.j-rhs$date.hatch.j
rhs$ttp5.4<-rhs$date.p5.4.j-rhs$date.hatch.j
rhs$ttwand<-rhs$date.wand.j-rhs$date.hatch.j
rhs$ttcull<-rhs$date.cull.j-rhs$date.hatch.j
rhs$ttem.h<-rhs$date.em.j-rhs$date.hatch.j
rhs$ttem.w<-rhs$date.em.j-rhs$date.ovp.j
rhs$ttecl<-rhs$date.ecl.j-rhs$date.ovp.j

#-------------------------

#Calculate wasp metrics



#Calculate percent survival

rhs$ps.em<-rhs$num.em/rhs$tot.load
rhs$ps.ecl<-rhs$num.ecl/rhs$tot.load

#Convert NaN to 0s

rhs$ps.em[is.nan(rhs$ps.em)]<-0
rhs$ps.ecl[is.nan(rhs$ps.ecl)]<-0


#calculate mass for adult male and female wasps

rhs$ind.male.mass<-rhs$male.mass/rhs$male.ecl
rhs$ind.fem.mass<-rhs$fem.mass/rhs$fem.ecl


#create column for total number unem, then calculate the perc that were mature or immature 2nds
rhs$tot.unem<-rhs$num.unem.mat+rhs$num.unem.im

rhs$p.unem.mat<-rhs$num.unem.mat/rhs$tot.unem
rhs$p.unem.im<-rhs$num.unem.im/rhs$tot.unem


#------------------------------------

#creating a sorting column for lates that were shocked at a lower temperature
  ##lates that were in the heat shock when it was recalibrated on 11/24/18 (j.date=328)
  ##lates that went in after the 11/22/18 would also have experienced this lower shock, so using this as cut off date (j.date=326)

rhs$date.in.hs.j[is.na(rhs$date.in.hs.j)]<-0

rhs$hs.cal<-ifelse(rhs$date.in.hs.j=="0", "none",
                   ifelse(rhs$date.in.hs.j<326, "high", "low"))



#-------------------------------------------

#creating a long dataset for age and mass of hosts

mss.lng<- rhs %>% gather(instar, mass, mass.3, mass.4, mass.5, mass.p5.1, mass.p5.2, mass.p5.3, mass.p5.4, mass.48em, mass.cull, mass.wand)
mss.lng$instar<-gsub("mass.", "", mss.lng$instar)
mss.lng$instar<-gsub("48", "", mss.lng$instar)

age.lng<-rhs %>% gather(instar, age, tt3, tt4, tt5, ttp5.1, ttp5.2, ttp5.3, ttp5.4, ttem.h, ttcull, ttwand)
age.lng$instar<-gsub("tt", "", age.lng$instar)
age.lng$instar<-gsub("em.h", "em", age.lng$instar)

age.lng<-age.lng %>% select(id, shock.stage, instar, age)

rhs.lng<-merge(mss.lng, age.lng, by=c("id", "shock.stage", "instar"))


#creating a long dataset for age and mass of hosts
  ##without p5 data, as it creates a very large dataframe

mss.lng2<- rhs %>% gather(instar, mass, mass.3, mass.4, mass.5, mass.48em, mass.cull, mass.wand)
mss.lng2$instar<-gsub("mass.", "", mss.lng2$instar)
mss.lng2$instar<-gsub("48", "", mss.lng2$instar)

age.lng2<-rhs %>% gather(instar, age, tt3, tt4, tt5, ttem.h, ttcull, ttwand)
age.lng2$instar<-gsub("tt", "", age.lng2$instar)
age.lng2$instar<-gsub("em.h", "em", age.lng2$instar)


age.lng2<-age.lng2 %>% select(id, shock.stage, instar, age)

rhs.lng2<-merge(mss.lng2, age.lng2, by=c("id", "shock.stage","instar"))


#--------------------------

#write to csv

write.csv(rhs, "Ms+Cc_RHS_incomplete_clean.csv",row.names = FALSE)
write.csv(rhs.lng, "Ms+Cc_RHS_incomp_clean_vlong.csv", row.names=FALSE)
write.csv(rhs.lng2,"Ms+Cc_RHS_incomp_clean_long.csv", row.names=FALSE)



