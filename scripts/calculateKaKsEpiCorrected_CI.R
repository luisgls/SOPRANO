rm(list=ls())
library(ggplot2)
library(reshape)
library(tidyr)
source("~/Rscripts/qqplot_check.R")


#input data sources
args <- commandArgs(trailingOnly = TRUE)
file.name <-  args[1]
#file.name <-"~/tools/local_SSBselection/tmp/test.data_epitopes"
file.name2 <-  args[2]  
#file.name2 <- "~/tools/local_SSBselection/tmp/test.epitope_NaNs.txt"
file.name3 <- args[3] 
#file.name3 <- "~/tools/local_SSBselection/tmp/test.nonepitope_NaNs.txt"

df<-read.csv(file.name,header=F,sep="\t")
colnames(df)<-c("EnsembleID","Total","Class")

##From long to wide format
df2<-spread(df,Class,Total,fill=0)

##Read calculation for total sites
df.sites.extra<-read.csv(file.name2,header=F,sep="\t")
df.sites.intra<-read.csv(file.name3,header=F,sep="\t")


na_epi<-sum(df2$extra_missense_variant)
ns_epi<-sum(df2$extra_synonymous_variant)

na_nonepi<-sum(df2$intra_missense_variant)
ns_nonepi<-sum(df2$intra_synonymous_variant)


ALL.extraKaKs<-(na_epi/df.sites.extra[1,1])/(ns_epi/df.sites.extra[1,2])
ALL.intraKaKs<-(na_nonepi/df.sites.intra[1,1])/(ns_nonepi/df.sites.intra[1,2])

#Mutations for the epitope part
m<-na_epi
s<-ns_epi
M<-df.sites.extra[1,1]
S<-df.sites.extra[1,2]
##Calculate conf interval using Katz method
p1<-(m/(M+1))
p2<-(s/(S+1))
globaldnds<-p1/p2
N1 <- M
N2 <- S
SE = sqrt( (1-p1)/(N1*p1) + (1-p2)/(N2*p2) )
finalLowCI = globaldnds * exp(-1.96*SE)
finalHighCI = globaldnds * exp(1.96*SE)
N<-m+s

#Mutations for the non epitope part
m<-na_nonepi
s<-ns_nonepi
M<-df.sites.intra[1,1]
S<-df.sites.intra[1,2]
##Calculate conf interval using Katz method
p1<-(m/(M+1))
p2<-(s/(S+1))
globaldnds<-p1/p2
N1 <- M
N2 <- S
SE = sqrt( (1-p1)/(N1*p1) + (1-p2)/(N2*p2) )
finalLowCI2 = globaldnds * exp(-1.96*SE)
finalHighCI2 = globaldnds * exp(1.96*SE)
N2<-m+s

#Print results
high = finalHighCI2
low = finalLowCI2
val = ALL.extraKaKs
SE = (log(high)-log(low) )/ (2 * 1.96)
EST = log(val)
z = EST/SE
PVAL1=exp(-0.717*z - 0.416*z^(2))
PVAL2=exp(-0.717*-z - 0.416*-z^(2))

if(PVAL2 > 0 & PVAL2 <= 1) {
  PVAL = PVAL2
} else {
  PVAL = PVAL1
}

cat(paste(paste(ALL.extraKaKs,finalLowCI,finalHighCI,N,ALL.intraKaKs,finalLowCI2,finalHighCI2,N2,PVAL,"\t"),"\n"))
