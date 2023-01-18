## remove (almost) everything in the working environment.
rm(list = ls())
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  args <- 1
}


setwd("/media/ssd2To/FishStrat/")
source("parIbasam.R")
pops <- c("Leff", "Trieux", "Jaudy", "Leguer", "Yar", "Douron", "Penze", "Elorn", "Aulne", "Goyen", "Odet", "Aven", "Laita", "Scorff", "Blavet")
nSIMUL <- 30
nYear <- 40
nInit <- 10 # Nb years at initialization
npop <- 15

#h = c(1, .95, .9, .85, .8, .75, .7) # Philopatry (homing) rates #edit al - 22/03/21 5% and 15%
# scenarioConnect=7 #scenario 1 for h=1.00, scenario 2 for h=0.95, scenario 3 pour h=0.80
# 
# # 0: control / 
# # 1: no fishing on sink populations 
# # 2: no fishing on source populations 
# scenarioFishing = 0

frates_vec =c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)
# scenarioFishingRate = 1

# Sc_name <- paste0("FishStrat_",scenarioConnect,scenarioFishing,scenarioFishingRate)


scn=c(
  101, 102, 103, 104, 105, 106, 107, 108, 109, 1010, 1011, 1012, 1013, NA, NA, NA, NA, NA, NA, NA, NA
  ,401, 402, 403, 404, 405, 406, 407, 408, 409, 4010, 4011, 4012, 4013, NA, NA, NA, NA, NA, NA, NA, NA
  ,401, 412, 413, 414, 415, 416, 417, 418, 419, 4110, 4111, 4112, 4113, 4114, 4115, 4116, 4117, 4118, 4119, 4120, 4121
  ,401, 422, 423, 424, 425, 426, 427, 428, 429, 4210, 4211, 4212, 4213, 4214, 4215, 4216, 4217, 4218, 4219, 4220, 4221
  ,701, 702, 703, 704, 705, 706, 707, 708, 709, 7010, 7011, 7012, 7013, NA, NA, NA, NA, NA, NA, NA, NA
  ,701, 712, 713, 714, 715, 716, 717, 718, 719, 7110, 7111, 7112, 7113, 7114, 7115, 7116, 7117, 7118, 7119, 7120, 7121
  ,701, 722, 723, 724, 725, 726, 727, 728, 729, 7210, 7211, 7212, 7213, 7214, 7215, 7216, 7217, 7218, 7219, 7220, 7221
)
SCN <-matrix(scn,length(frates_vec),7)

frates_vec =c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)
SCN <- SCN[1:length(frates_vec),]

##### METAPOPULATION ABUNDANCE

metapops.means=metapops.means.rel=metapops.exp.means=Pglobal.means=extinction.means.rel=array(,dim=c(nrow(SCN),nSIMUL,ncol(SCN)))
Metapop<-Exploitation<-Pexpl.final<-Extinction<-Prop.Immigrants<-Population<-RatioImExpl<-list()
pops.means=migs.means=migexpl.means=array(,dim=c(nrow(SCN),nSIMUL,npop,ncol(SCN)))

LFPARR<-LFSMOLT<-LFMSW<-LF1SW<-LFRETURNS<-RMSW<-list()
LFparr.means=LFsmolt.means=LFreturns.means=LFmsw.means=LF1sw.means=rmsw.means=array(,dim=c(nrow(SCN),nSIMUL,npop,ncol(SCN)))
for (i in 1:ncol(SCN)){
#  for (i in 1:1){

  LFparr.pop=LFsmolt.pop=LFreturns.pop=LFmsw.pop=LF1sw.pop=Rmsw.pop=list(list())
  j=0
  #for (iEXPE in EXPE){ # Loop over scenario 
  for (iEXPE in SCN[,i]){
    j=j+1
    if(is.na(iEXPE)) next;
    
    load(paste0("results/PHENOTYPE",iEXPE,".RData"))
    
    scn_id <- as.numeric(strsplit(as.character(iEXPE), "")[[1]])
    scenarioConnect = scn_id[1]
    scenarioFishing = scn_id[2]
    scenarioFishingRate = scn_id[3]
    
    
    #### POP
    for (pop in 1:npop){
      
      LFparr.tmp=LFsmolt.tmp=LFreturns.tmp=LFmsw.tmp=LF1sw.tmp=rmsw.tmp=NULL
      
      for (simul in 1:nSIMUL){
        
        tmp1=tmp2=tmp3=tmp4=tmp5=tmp6=NULL
        tmp1 <- LFParr[[paste0(iEXPE)]][[simul]][,pop]
        tmp2 <- LFSmolt[[paste0(iEXPE)]][[simul]][,pop]#/LFReturns[[paste0(iEXPE)]][[simul]][1,pop]
        tmp3 <- LFReturns[[paste0(iEXPE)]][[simul]][,pop]
        tmp4 <- LFmsw[[paste0(iEXPE)]][[simul]][,pop]
        tmp5 <- LF1sw[[paste0(iEXPE)]][[simul]][,pop]
        tmp6 <- rMSW[[paste0(iEXPE)]][[simul]][,pop]
        
        LFparr.means[j,simul,pop,i] <- mean(tail(tmp1,5))
        LFsmolt.means[j,simul,pop,i] <- mean(tail(tmp2,5))
        LFreturns.means[j,simul,pop,i] <- mean(tail(tmp3,5))
        LFmsw.means[j,simul,pop,i] <- mean(tail(tmp4,5))
        LF1sw.means[j,simul,pop,i] <- mean(tail(tmp5,5))
        rmsw.means[j,simul,pop,i] <- mean(tail(tmp6,5))
        
        LFparr.tmp <- cbind(LFparr.tmp, tmp1)
        LFsmolt.tmp <- cbind(LFsmolt.tmp, tmp2)
        LFreturns.tmp <- cbind(LFreturns.tmp, tmp3)
        LFmsw.tmp <- cbind(LFmsw.tmp, tmp4)
        LF1sw.tmp <- cbind(LF1sw.tmp, tmp5)
        rmsw.tmp <- cbind(rmsw.tmp, tmp6)
        
      }#end loop simul
      
      LFparr.pop[[paste0(iEXPE)]][[pop]]<-LFparr.tmp
      LFsmolt.pop[[paste0(iEXPE)]][[pop]]<-LFsmolt.tmp
      LFreturns.pop[[paste0(iEXPE)]][[pop]]<-LFreturns.tmp
      LFmsw.pop[[paste0(iEXPE)]][[pop]]<-LFmsw.tmp
      LF1sw.pop[[paste0(iEXPE)]][[pop]]<-LF1sw.tmp
      Rmsw.pop[[paste0(iEXPE)]][[pop]]<-rmsw.tmp
      
    }#end loop pop
    
    #### METAPOP
    load(paste0("results/DEMOGRAPHY",iEXPE,".RData"))
    metapops=exploit=NULL
    metapop_exp=NULL
    extinction=NULL
    for (simul in 1:nSIMUL){
      tmp=NULL
      tmp <- nReturns[[paste0(iEXPE)]][[simul]]
      
      metapops <- cbind(metapops, rowSums(tmp, na.rm=TRUE))
      metapops.means[j,simul,i] <- mean(tail(metapops[,simul],5))
      #Means[j,i] <- mean(tail(metapops,5))
      metapops.means.rel[j,simul,i] <- metapops.means[j,simul,i]/metapops.means[1,simul,i] #faire valeur median control toutes simul confondues ?
      
      metapop_exp <- cbind(metapop_exp, rowSums(nExploitation[[paste0(iEXPE)]][[simul]], na.rm=T)) #compute total number of fish exploited
      metapops.exp.means[j,simul,i] <- mean(tail(metapop_exp[,simul],5)) / metapops.means[j,simul,i] #exploitation rate
      
      tmp3 <- 1-(apply(tmp,1,function(x) length(which(!is.na(x))))/npop)
      extinction <- cbind(extinction, tmp3)
      extinction.means.rel[j,simul,i] <- mean(tail(tmp3,5))
      
      tmp2=NULL
      tmp2 <- nExploitation[[paste0(iEXPE)]][[simul]]
      exploit <- cbind(exploit, rowSums(tmp2, na.rm=TRUE))
      Pglobal <- exploit / metapops
      Pglobal.means[j,simul,i] <- mean(tail(Pglobal[,simul],5))
    }

  } # loop iEXPE
  
  LFPARR[[i]]<-LFparr.pop
  LFSMOLT[[i]]<-LFsmolt.pop
  LFRETURNS[[i]]<-LFreturns.pop
  LFMSW[[i]]<-LFmsw.pop
  LF1SW[[i]]<-LF1sw.pop
  RMSW[[i]]<-Rmsw.pop
  
  Pexpl.final[[i]]<- Pglobal.means
} # loop SCN






load(paste0("results/PHENOTYPE",101,".RData"))
init <- apply(tail(LFmsw[[paste0(101)]][[1]],10),2,mean,na.rm=TRUE)
tmp<-init/sum(init)
prop.sink <- sum(tmp[dat$Type!="source"],na.rm=TRUE)
prop.source <- sum(tmp[dat$Type!="sink"],na.rm=TRUE)


  






#### POPULATION DYNAMICS ####
#-------------------------------#

#note: always pass alpha on the 0-255 scale
makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

col.sources=c("dodgerblue2", "deepskyblue")
col.sinks=c("brown1", "lightcoral", "white")

# Create a list of gradients for each colour 2 to 10 over five steps from 
library(viridis)
#col <- viridis(13)
col <- (inferno(length(frates_vec)))
#col=c("#000000", "#1874CD", "#009ACD", "#00BFFF", "#228B22", "#66CD00", "#7FFF00")

pop.type<-as.numeric(factor(dat$Type))
col.type <- c("black",col.sinks[1],col.sources[1])







pdf(file="results/PHENOTYPE.pdf")


xlab="Time (years)"
ylab=" value" # Pre-Fisheries Abundance
lwd=2


#### METAPOP
par(mfcol=c(1,1))
x=2:50
main=paste0("Ratio MSW/1SW")
#tmp <- LFPARR
ylim=c(0,0.4)#c(0,max(unlist(tmp), na.rm=TRUE))
plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
#title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
rect(0,0,10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
for (pop in 1:npop){
  #tmp<-apply(RMSW[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
  tmp<-apply(RMSW[[1]][["101"]][[pop]],1,median,na.rm=TRUE)
  #tmp2.5<-apply(RMSW[[1]][["101"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  #tmp97.5<-apply(RMSW[[1]][["101"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  #polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[1],20))
  points(x,tmp[x],col=col[1],type='l',lwd=lwd)
}


par(mfcol=c(1,1))
x=2:50
main=paste0("Fork Length MSW")
#tmp <- LFPARR
ylim=c(700,800)#c(0,max(ylim))
plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
#title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
rect(0,0,10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
for (pop in 1:npop){
  #tmp<-apply(LFMSW[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
  tmp<-apply(LFMSW[[1]][["101"]][[pop]],1,median,na.rm=TRUE)
  #tmp2.5<-apply(LFMSW[[1]][["101"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  #tmp97.5<-apply(LFMSW[[1]][["101"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  #polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[1],20))
  points(x,tmp[x],col=col[1],type='l',lwd=lwd)
}










### PHENOTYPES

par(mfcol=c(3,2))
x=2:50
for (pop in 1:npop){
  
  main=paste0("Fork Length Parr 0+ (", pops[pop],")")
  #tmp <- LFPARR
  ylim=c(60,120)#c(0,max(unlist(tmp), na.rm=TRUE))
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main="")
  title(main, col.main = col.type[pop.type[pop]])
  rect(0,0,10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
  #tmp<-apply(LFPARR[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
  tmp<-apply(LFPARR[[1]][["101"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(LFPARR[[1]][["101"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(LFPARR[[1]][["101"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[1],20))
  points(x,tmp[x],col=col[1],type='l',lwd=lwd)
  
  tmp<-apply(LFPARR[[1]][["105"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(LFPARR[[1]][["105"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(LFPARR[[1]][["105"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[5],20))
  points(x,tmp[x],col=col[5],type='l',lwd=lwd)
  
  tmp<-apply(LFPARR[[1]][["1010"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(LFPARR[[1]][["1010"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(LFPARR[[1]][["1010"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[10],20))
  points(x,tmp[x],col=col[10],type='l',lwd=lwd)
  
  
  
  
  
  main=paste0("Fork Length Smolt 0+ (", pops[pop],")")
  #tmp <- apply(LFSMOLT[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
  tmp<-apply(LFSMOLT[[1]][["101"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(LFSMOLT[[1]][["101"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(LFSMOLT[[1]][["101"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  ylim=c(100,160)#c(0,max(ylim))
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  rect(0,0,10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[1],20))
  points(x,tmp[x],col=col[1],type='l',lwd=lwd)
  
  tmp<-apply(LFSMOLT[[1]][["105"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(LFSMOLT[[1]][["105"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(LFSMOLT[[1]][["105"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[5],20))
  points(x,tmp[x],col=col[5],type='l',lwd=lwd)
  
  tmp<-apply(LFSMOLT[[1]][["1010"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(LFSMOLT[[1]][["1010"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(LFSMOLT[[1]][["1010"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[10],20))
  points(x,tmp[x],col=col[10],type='l',lwd=lwd)
  
  
  
  
  
  # main=paste0("Parr females maturation threshold (", pops[pop],")")
  # tmp <- apply(GFMID2[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
  # ylim=c(550,650)#c(0,max(ylim))
  # plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  # #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  # rect(0,0,10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
  # points(x,tmp[x],col=1,type='l',lwd=lwd)
  
  main=paste0("Fork Length MSW (", pops[pop],")")
  #tmp <- apply(LFMSW[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
  tmp<-apply(LFMSW[[1]][["101"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(LFMSW[[1]][["101"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(LFMSW[[1]][["101"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  ylim=c(700,800)#c(0,max(ylim))
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  rect(0,0,10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[1],20))
  points(x,tmp[x],col=1,type='l',lwd=lwd)
  
  tmp<-apply(LFMSW[[1]][["105"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(LFMSW[[1]][["105"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(LFMSW[[1]][["105"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[5],20))
  points(x,tmp[x],col=col[5],type='l',lwd=lwd)
  
  tmp<-apply(LFMSW[[1]][["1010"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(LFMSW[[1]][["1010"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(LFMSW[[1]][["1010"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[10],20))
  points(x,tmp[x],col=col[10],type='l',lwd=lwd)
  
  
  
  
  
  main=paste0("Fork Length 1SW (", pops[pop],")")
  #tmp <- apply(LF1SW[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
  tmp<-apply(LF1SW[[1]][["101"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(LF1SW[[1]][["101"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(LF1SW[[1]][["101"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  ylim=c(550,650)#c(0,max(ylim))
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  rect(0,0,10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[1],20))
  points(x,tmp[x],col=1,type='l',lwd=lwd)
  
  tmp<-apply(LF1SW[[1]][["105"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(LF1SW[[1]][["105"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(LF1SW[[1]][["105"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[5],20))
  points(x,tmp[x],col=col[5],type='l',lwd=lwd)
  
  tmp<-apply(LF1SW[[1]][["1010"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(LF1SW[[1]][["1010"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(LF1SW[[1]][["1010"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[10],20))
  points(x,tmp[x],col=col[10],type='l',lwd=lwd)
  
  
  
  
  
  
  main=paste0("Ratio MSW/1SW (", pops[pop],")")
  #tmp <- apply(RMSW[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
  tmp<-apply(RMSW[[1]][["101"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(RMSW[[1]][["101"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(RMSW[[1]][["101"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  ylim=c(0,.4)#c(0,max(ylim))
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  rect(0,min(ylim),10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[1],20))
  points(x,tmp[x],col=col[1],type='l',lwd=lwd)
  
  tmp<-apply(RMSW[[1]][["105"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(RMSW[[1]][["105"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(RMSW[[1]][["105"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[7],20))
  points(x,tmp[x],col=col[5],type='l',lwd=lwd)
  
  tmp<-apply(RMSW[[1]][["1010"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(RMSW[[1]][["1010"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(RMSW[[1]][["1010"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[10],20))
  points(x,tmp[x],col=col[10],type='l',lwd=lwd)
  
  
  plot.new()
  #} # end loop pop
  
  
  
  
  
  
  
  
  
  par(mfrow=c(3,3))
  
  
  #for (pop in 1:npop){
  # pop=14
  
  xlab="local exploitation rate"
  ylab="Fork Length Parr 0+"
  
  xlim=c(0,max(frates_vec))
  ylim=c(80,120)#c(0,max(unlist(tmp), na.rm=TRUE))
  
  
  #### ALL
  main=paste0("h=1 / ALL")#paste0("Neutral gene (", pops[pop],")")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFparr.means[,,pop,1]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.85 / ALL")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFparr.means[,,pop,2]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / ALL")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFparr.means[,,pop,5]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  
  #### SOURCE
  plot.new()
  
  
  main=paste0("h=0.85 / SOURCE")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFparr.means[,,pop,3]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / SOURCE")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFparr.means[,,pop,6]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  
  #### SINK
  plot.new()
  
  
  main=paste0("h=0.85 / SINK")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFparr.means[,,pop,4]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / SINK")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFparr.means[,,pop,7]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  
  
  
  
  
  
  
  
  
  
  par(mfrow=c(3,3))
  
  
  #for (pop in 1:npop){
  #pop=14
  
  xlab="local exploitation rate"
  ylab=paste0("Fork Length Smolt (", pops[pop],")")
  
  #tmp <- LFPARR
  ylim=c(120,160)
  
  #### ALL
  main=paste0("h=1 / ALL")#paste0("Neutral gene (", pops[pop],")")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFsmolt.means[,,pop,1]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.85 / ALL")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFsmolt.means[,,pop,2]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / ALL")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFsmolt.means[,,pop,5]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  
  #### SOURCE
  plot.new()
  
  
  main=paste0("h=0.85 / SOURCE")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFsmolt.means[,,pop,3]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / SOURCE")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFsmolt.means[,,pop,6]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  
  #### SINK
  plot.new()
  
  
  main=paste0("h=0.85 / SINK")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFsmolt.means[,,pop,4]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / SINK")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFsmolt.means[,,pop,7]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  
  
  
  
  
  
  
  
  
  
  par(mfrow=c(3,3))
  
  
  #for (pop in 1:npop){
  #pop=14
  
  xlab="local exploitation rate"
  ylab=paste0("Fork Length MSW (", pops[pop],")")
  
  #tmp <- LFPARR
  ylim=c(700,800)
  
  #### ALL
  main=paste0("h=1 / ALL")#paste0("Neutral gene (", pops[pop],")")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFmsw.means[,,pop,1]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.85 / ALL")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFmsw.means[,,pop,2]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / ALL")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFmsw.means[,,pop,5]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  
  #### SOURCE
  plot.new()
  
  
  main=paste0("h=0.85 / SOURCE")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFmsw.means[,,pop,3]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / SOURCE")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFmsw.means[,,pop,6]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  
  #### SINK
  plot.new()
  
  
  main=paste0("h=0.85 / SINK")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFmsw.means[,,pop,4]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / SINK")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFmsw.means[,,pop,7]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  
  
  
  
  
  
  
  
  par(mfrow=c(3,3))
  
  
  #for (pop in 1:npop){
  #pop=14
  
  xlab="local exploitation rate"
  ylab=paste0("Fork Length 1SW (", pops[pop],")")
  
  #tmp <- LFPARR
  ylim=c(550,650)
  
  #### ALL
  main=paste0("h=1 / ALL")#paste0("Neutral gene (", pops[pop],")")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LF1sw.means[,,pop,1]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.85 / ALL")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LF1sw.means[,,pop,2]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / ALL")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LF1sw.means[,,pop,5]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  
  #### SOURCE
  plot.new()
  
  
  main=paste0("h=0.85 / SOURCE")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LF1sw.means[,,pop,3]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / SOURCE")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LF1sw.means[,,pop,6]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  
  #### SINK
  plot.new()
  
  
  main=paste0("h=0.85 / SINK")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LF1sw.means[,,pop,4]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / SINK")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LF1sw.means[,,pop,7]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  
  
  
  
  
  
  
  
  par(mfrow=c(3,3))
  
  
  #for (pop in 1:npop){
  #pop=13
  
  xlab="local exploitation rate"
  ylab=paste0("Ratio MSW/1SW (", pops[pop],")")
  
  #tmp <- LFPARR
  ylim=c(0,.4)
  
  #### ALL
  main=paste0("h=1 / ALL")#paste0("Neutral gene (", pops[pop],")")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-rmsw.means[,,pop,1]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.85 / ALL")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-rmsw.means[,,pop,2]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / ALL")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-rmsw.means[,,pop,5]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  
  #### SOURCE
  plot.new()
  
  
  main=paste0("h=0.85 / SOURCE")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-rmsw.means[,,pop,3]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / SOURCE")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-rmsw.means[,,pop,6]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  
  #### SINK
  plot.new()
  
  
  main=paste0("h=0.85 / SINK")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-rmsw.means[,,pop,4]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / SINK")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-rmsw.means[,,pop,7]
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(frates_vec,tmp2.5,frates_vec,tmp97.5,col=col[1:length(frates_vec)])
  points(frates_vec,tmp,col=col[1:length(frates_vec)],type='p',lwd=lwd,pch=20)
  
  
  
  
  
  
  
  
  
  
  
  
  ##### METAPOPULATION EXPLOITATION RATE
  
  
  par(mfrow=c(3,3))
  
  
  #for (pop in 1:npop){
  # pop=15
  
  xlab="Metapopulation exploitation rate"
  ylab=paste0("Fork Length Parr 0+ (", pops[pop],")")
  
  xlim=c(0,0.5)
  ylim=c(80,120)
  #### ALL
  main=paste0("h=1 / ALL")#paste0("Neutral gene (", pops[pop],")")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFparr.means[,,pop,1]
  var2<-Pglobal.means[,,1];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  
  main=paste0("h=0.85 / ALL")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFparr.means[,,pop,2]
  var2<-Pglobal.means[,,2];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / ALL")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFparr.means[,,pop,5]
  var2<-Pglobal.means[,,5];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  
  #### SOURCE
  plot.new()
  
  
  main=paste0("h=0.85 / SOURCE")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFparr.means[,,pop,3]
  var2<-Pglobal.means[,,3];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / SOURCE")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFparr.means[,,pop,6]
  var2<-Pglobal.means[,,6];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  
  #### SINK
  plot.new()
  
  
  main=paste0("h=0.85 / SINK")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFparr.means[,,pop,4]
  var2<-Pglobal.means[,,4];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / SINK")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFparr.means[,,pop,7]
  var2<-Pglobal.means[,,7];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  
  
  
  
  
  par(mfrow=c(3,3))
  
  
  #for (pop in 1:npop){
  # pop=15
  
  xlab="Metapopulation exploitation rate"
  ylab=paste0("Ratio MSW/1SW (", pops[pop],")")
  
  xlim=c(0,0.5)
  ylim=c(0,.4)
  
  #### ALL
  main=paste0("h=1 / ALL")#paste0("Neutral gene (", pops[pop],")")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-rmsw.means[,,pop,1]
  var2<-Pglobal.means[,,1];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  
  main=paste0("h=0.85 / ALL")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-rmsw.means[,,pop,2]
  var2<-Pglobal.means[,,2];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / ALL")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-rmsw.means[,,pop,5]
  var2<-Pglobal.means[,,5];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  
  #### SOURCE
  plot.new()
  
  
  main=paste0("h=0.85 / SOURCE")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-rmsw.means[,,pop,3]
  var2<-Pglobal.means[,,3];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / SOURCE")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-rmsw.means[,,pop,6]
  var2<-Pglobal.means[,,6];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  
  #### SINK
  plot.new()
  
  
  main=paste0("h=0.85 / SINK")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-rmsw.means[,,pop,4]
  var2<-Pglobal.means[,,4];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / SINK")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-rmsw.means[,,pop,7]
  var2<-Pglobal.means[,,7];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  
  
  
  
  
  
  
  
  
  
  par(mfrow=c(3,3))
  
  
  #for (pop in 1:npop){
  # pop=15
  
  xlab="Metapopulation exploitation rate"
  ylab=paste0("Fork Length MSW (", pops[pop],")")
  
  xlim=c(0,0.5)
  ylim=c(700,800)
  #### ALL
  main=paste0("h=1 / ALL")#paste0("Neutral gene (", pops[pop],")")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFmsw.means[,,pop,1]
  var2<-Pglobal.means[,,1];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  
  main=paste0("h=0.85 / ALL")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFmsw.means[,,pop,2]
  var2<-Pglobal.means[,,2];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / ALL")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFmsw.means[,,pop,5]
  var2<-Pglobal.means[,,5];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  
  #### SOURCE
  plot.new()
  
  
  main=paste0("h=0.85 / SOURCE")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFmsw.means[,,pop,3]
  var2<-Pglobal.means[,,3];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / SOURCE")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFmsw.means[,,pop,6]
  var2<-Pglobal.means[,,6];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  
  #### SINK
  plot.new()
  
  
  main=paste0("h=0.85 / SINK")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFmsw.means[,,pop,4]
  var2<-Pglobal.means[,,4];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / SINK")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LFmsw.means[,,pop,7]
  var2<-Pglobal.means[,,7];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  
  
  
  
  
  
  
  
  par(mfrow=c(3,3))
  
  
  #for (pop in 1:npop){
  # pop=15
  
  xlab="Metapopulation exploitation rate"
  ylab=paste0("Fork Length 1SW (", pops[pop],")")
  
  xlim=c(0,0.5)
  ylim=c(550,650)
  #### ALL
  main=paste0("h=1 / ALL")#paste0("Neutral gene (", pops[pop],")")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LF1sw.means[,,pop,1]
  var2<-Pglobal.means[,,1];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  
  main=paste0("h=0.85 / ALL")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LF1sw.means[,,pop,2]
  var2<-Pglobal.means[,,2];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / ALL")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LF1sw.means[,,pop,5]
  var2<-Pglobal.means[,,5];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  
  #### SOURCE
  plot.new()
  
  
  main=paste0("h=0.85 / SOURCE")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LF1sw.means[,,pop,3]
  var2<-Pglobal.means[,,3];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / SOURCE")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LF1sw.means[,,pop,6]
  var2<-Pglobal.means[,,6];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  
  #### SINK
  plot.new()
  
  
  main=paste0("h=0.85 / SINK")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LF1sw.means[,,pop,4]
  var2<-Pglobal.means[,,4];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  main=paste0("h=0.7 / SINK")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-LF1sw.means[,,pop,7]
  var2<-Pglobal.means[,,7];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  
  
  
  
} # end loop pop
dev.off()


