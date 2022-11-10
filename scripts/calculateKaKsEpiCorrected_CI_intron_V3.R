rm(list=ls())
suppressMessages(library("ggplot2"))
suppressMessages(library("reshape"))
suppressMessages(library("tidyr"))
suppressMessages(library("graphics"))


#input data sources
args <- commandArgs(trailingOnly = TRUE)
file.name <-  args[1]
#file.name <-"~/tools/local_SSBselection/tmp/test.data_epitopes"
file.name2 <-  args[2]  
#file.name2 <- "~/tools/local_SSBselection/tmp/test.epitope_NaNs.txt"
file.name3 <- args[3] 
#file.name3 <- "~/tools/local_SSBselection/tmp/test.nonepitope_NaNs.txt"
file.name4 <- args[4]
#file.name4 <- "~/tools/local_SSBselection/results/test.intronic.rate"

df<-read.csv(file.name,header=F,sep="\t")
colnames(df)<-c("EnsembleID","Total","Class")

##From long to wide format
#df2<-spread(df,Class,Total,fill=0)

#df2<-pivot_wider(df, names_from = Class, values_from = Total, values_fill = list(Total = 0), values_fn= list(Total = length))
df2<-spread(df,Class,Total,fill=0)
df2<-as.data.frame(df2)

##Read calculation for total sites
df.sites.extra<-read.csv(file.name2,header=F,sep="\t")
df.sites.intra<-read.csv(file.name3,header=F,sep="\t")

##Read intronic data
df.intron<-read.csv(file.name4,header=F,sep="\t")
colnames(df.intron)<-c("EnsemblID","intronrate","mutsintron","intronlength")

df3<-merge(df2,df.intron,by.x="EnsembleID",by.y="EnsemblID", all.x = T)
df3[is.na(df3)] <- 0

na_epi<-sum(df3$extra_missense_variant)
ns_epi<-sum(df3$extra_synonymous_variant)

na_nonepi<-sum(df3$intra_missense_variant)
ns_nonepi<-sum(df3$intra_synonymous_variant)

ns_intron<-sum(df3$mutsintron)
NS_intron<-sum(df3$intronlength)
NS2_intron<-ns_intron/( (ns_nonepi/df.sites.intra[1,2] + ns_epi/df.sites.extra[1,2])/2 )

ALL.extraKaKs<-(na_epi/df.sites.extra[1,1])/(ns_epi/df.sites.extra[1,2])
ALL.intraKaKs<-(na_nonepi/df.sites.intra[1,1])/(ns_nonepi/df.sites.intra[1,2])

ALL.intra_intronKaKs<-(na_nonepi/df.sites.intra[1,1]) / ( (ns_nonepi+ns_intron)/(df.sites.intra[1,2]+NS2_intron) )

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

#Estimate P-value from confidence interval of the non-target region
high = finalHighCI2
low = finalLowCI2
val = ALL.extraKaKs - ALL.intraKaKs
SE = ((high)-(low) )/ (2 * 1.96)
EST = (val)
z = EST/SE
PVAL1=exp(-0.717*z - 0.416*z^(2))
PVAL2=exp(-0.717*-z - 0.416*-z^(2))

if(PVAL2 > 0 & PVAL2 <= 1) {
  PVAL = PVAL2
} else {
  PVAL = PVAL1
}

if (PVAL < 0.0001) {
  PVAL = 0.0001
}

###PRint the results: TYPE,dNdS_target,lowCI_target,highCI_target,TotalMuts_Target,dNdS_nonTarget,lowCI_nontarget,highCI_nontarget,TotalMuts_nonTarget,Pvalue,
### nonsyn_target,nonsynsites_target,syn_target,synsites_target,nonsyn_nontarget,nonsynsites_nontarget,syn_nontarget,synsites_nontarget
cat(paste(paste("ExonicOnly",ALL.extraKaKs,finalLowCI,finalHighCI,N,
                ALL.intraKaKs,finalLowCI2,finalHighCI2,N2,
                PVAL,
                na_epi,df.sites.extra[1,1],ns_epi,df.sites.extra[1,2],
                na_nonepi,df.sites.intra[1,1],ns_nonepi,df.sites.intra[1,2],"\t"),"\n"))

#Estimate confidence interval using intron muts 

m<-na_nonepi
s<-ns_nonepi+ns_intron
M<-df.sites.intra[1,1]
S<-df.sites.intra[1,2]+NS2_intron

##Calculate conf interval using Katz method
p1<-(m/(M+1))
p2<-(s/(S+1))
globaldnds<-p1/p2
N1 <- M
N2 <- S
SE = sqrt( (1-p1)/(N1*p1) + (1-p2)/(N2*p2) )
finalLowCI3 = globaldnds * exp(-1.96*SE)
finalHighCI3 = globaldnds * exp(1.96*SE)
N3 <- m + s
##Pval using intronic
#Estimate P-value from confidence interval of the non-target region
high = finalHighCI3
low = finalLowCI3
val = ALL.extraKaKs - ALL.intra_intronKaKs
SE = (high-low)/(2 * 1.96)
EST = val
z = EST/SE
PVAL1intron=exp(-0.717*z - 0.416*z^(2))
PVAL2intron=exp(-0.717*-z - 0.416*-z^(2))

if(PVAL2intron > 0 & PVAL2intron <= 1) {
  PVALintron = PVAL2intron
} else {
  PVALintron = PVAL1intron
}

if (PVALintron < 0.0001) {
  PVALintron = 0.0001
}
#Print results

cat(paste(paste("ExonicIntronic",ALL.extraKaKs,finalLowCI,finalHighCI,N,
                ALL.intra_intronKaKs,finalLowCI3,finalHighCI3,N3,
                PVALintron,
                na_epi,df.sites.extra[1,1],ns_epi,df.sites.extra[1,2],
                na_nonepi,df.sites.intra[1,1],(ns_nonepi+ns_intron),(df.sites.intra[1,2]+NS2_intron),
                "\t"),"\n"))
