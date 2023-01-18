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

GNEUTRAL<-GFMID1<-GFMID2<-GFMID3<-GFMID4<-GGROWTH<-list()
gNEUTRAL.means=gFMID1.means=gFMID2.means=gFMID3.means=gFMID4.means=gGROWTH.means=array(,dim=c(nrow(SCN),nSIMUL,npop,ncol(SCN)))
for (i in 1:ncol(SCN)){
  #for (i in 1:5){

    gNEUTRAL.pop=gFMID1.pop=gFMID2.pop=gFMID3.pop=gFMID4.pop=gGROWTH.pop=list(list())
    
  j=0
  #for (iEXPE in EXPE){ # Loop over scenario 
  for (iEXPE in SCN[,i]){
    j=j+1
    if(is.na(iEXPE)) next;
    
    load(paste0("results/GENOTYPE",iEXPE,".RData"))
    
    scn_id <- as.numeric(strsplit(as.character(iEXPE), "")[[1]])
    scenarioConnect = scn_id[1]
    scenarioFishing = scn_id[2]
    scenarioFishingRate = scn_id[3]
    
    
    #### POP
    for (pop in 1:npop){
      
      gNEUTRAL.tmp=gFMID1.tmp=gFMID2.tmp=gFMID3.tmp=gFMID4.tmp=gGROWTH.tmp=NULL
      
      for (simul in 1:nSIMUL){
        
        tmp1=tmp2=tmp3=tmp4=tmp5=tmp6=NULL
        #gNEUTRAL,gFMID1,gFMID2,gFMID3,gFMID4,gGROWTH
        tmp1 <- gNEUTRAL[[paste0(iEXPE)]][[simul]][,pop]
        tmp2 <- gFMID1[[paste0(iEXPE)]][[simul]][,pop]#/LFReturns[[paste0(iEXPE)]][[simul]][1,pop]
        tmp3 <- gFMID2[[paste0(iEXPE)]][[simul]][,pop]
        tmp4 <- gFMID3[[paste0(iEXPE)]][[simul]][,pop]
        tmp5 <- gFMID4[[paste0(iEXPE)]][[simul]][,pop]
        tmp6 <- gGROWTH[[paste0(iEXPE)]][[simul]][,pop]
        
        gNEUTRAL.means[j,simul,pop,i] <- mean(tail(tmp1,5))
        gFMID1.means[j,simul,pop,i] <- mean(tail(tmp2,5))
        gFMID2.means[j,simul,pop,i] <- mean(tail(tmp3,5))
        gFMID3.means[j,simul,pop,i] <- mean(tail(tmp4,5))
        gFMID4.means[j,simul,pop,i] <- mean(tail(tmp5,5))
        gGROWTH.means[j,simul,pop,i] <- mean(tail(tmp6,5))
        
        gNEUTRAL.tmp <- cbind(gNEUTRAL.tmp, tmp1)
        gFMID1.tmp <- cbind(gFMID1.tmp, tmp2)
        gFMID2.tmp <- cbind(gFMID2.tmp, tmp3)
        gFMID3.tmp <- cbind(gFMID3.tmp, tmp4)
        gFMID4.tmp <- cbind(gFMID4.tmp, tmp5)
        gGROWTH.tmp <- cbind(gGROWTH.tmp, tmp6)
        
      }#end loop simul
      
      gNEUTRAL.pop[[paste0(iEXPE)]][[pop]]<-gNEUTRAL.tmp
      gFMID1.pop[[paste0(iEXPE)]][[pop]]<-gFMID1.tmp
      gFMID2.pop[[paste0(iEXPE)]][[pop]]<-gFMID2.tmp
      gFMID3.pop[[paste0(iEXPE)]][[pop]]<-gFMID3.tmp
      gFMID4.pop[[paste0(iEXPE)]][[pop]]<-gFMID4.tmp
      gGROWTH.pop[[paste0(iEXPE)]][[pop]]<-gGROWTH.tmp
      
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
  
  GNEUTRAL[[i]]<-gNEUTRAL.pop
  GFMID1[[i]]<-gFMID1.pop
  GFMID2[[i]]<-gFMID2.pop
  GFMID3[[i]]<-gFMID3.pop
  GFMID4[[i]]<-gFMID4.pop
  GGROWTH[[i]]<-gGROWTH.pop
  
  Pexpl.final[[i]]<- Pglobal.means
} # loop SCN






load(paste0("results/GENOTYPE",101,".RData"))
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







pdf(file="results/GENOTYPE.pdf")


xlab="Time (years)"
ylab="Genetic value" # Pre-Fisheries Abundance
lwd=2


#### METAPOP
par(mfcol=c(1,1))
x=2:50
main=paste0("Neutral gene")
#tmp <- LFPARR
ylim=c(0.4,0.6)#c(0,max(unlist(tmp), na.rm=TRUE))
plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
#title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
rect(0,0,10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
for (pop in 1:npop){
#tmp<-apply(GNEUTRAL[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
tmp<-apply(GNEUTRAL[[1]][["101"]][[pop]],1,median,na.rm=TRUE)
#tmp2.5<-apply(GNEUTRAL[[1]][["101"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
#tmp97.5<-apply(GNEUTRAL[[1]][["101"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
#polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[1],20))
points(x,tmp[x],col=col[1],type='l',lwd=lwd)
}


par(mfcol=c(1,1))
x=2:50
main=paste0("GFMID4")
#tmp <- LFPARR
ylim=c(50,150)#c(0,max(ylim))
plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
#title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
rect(0,0,10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
for (pop in 1:npop){
  #tmp<-apply(GNEUTRAL[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
  tmp<-apply(GFMID4[[1]][["101"]][[pop]],1,median,na.rm=TRUE)
  #tmp2.5<-apply(GNEUTRAL[[1]][["101"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  #tmp97.5<-apply(GNEUTRAL[[1]][["101"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  #polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[1],20))
  points(x,tmp[x],col=col[1],type='l',lwd=lwd)
}










### PHENOTYPES

par(mfcol=c(3,2))
x=2:50
for (pop in 1:npop){
  
  main=paste0("Neutral gene (", pops[pop],")")
  #tmp <- LFPARR
  ylim=c(0.4,0.6)#c(0,max(unlist(tmp), na.rm=TRUE))
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  rect(0,0,10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
  #tmp<-apply(GNEUTRAL[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
  tmp<-apply(GNEUTRAL[[1]][["101"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(GNEUTRAL[[1]][["101"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(GNEUTRAL[[1]][["101"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[1],20))
  points(x,tmp[x],col=col[1],type='l',lwd=lwd)
  
  tmp<-apply(GNEUTRAL[[1]][["105"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(GNEUTRAL[[1]][["105"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(GNEUTRAL[[1]][["105"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[5],20))
  points(x,tmp[x],col=col[5],type='l',lwd=lwd)
  
  tmp<-apply(GNEUTRAL[[1]][["1010"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(GNEUTRAL[[1]][["1010"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(GNEUTRAL[[1]][["1010"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[10],20))
  points(x,tmp[x],col=col[10],type='l',lwd=lwd)
  
  
  
  
  
  main=paste0("Parr males maturation threshold (", pops[pop],")")
  #tmp <- apply(GFMID1[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
  tmp<-apply(GFMID1[[1]][["101"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(GFMID1[[1]][["101"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(GFMID1[[1]][["101"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
   ylim=c(1,1.5)#c(0,max(ylim))
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  rect(0,0,10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[1],20))
  points(x,tmp[x],col=col[1],type='l',lwd=lwd)
  
  tmp<-apply(GFMID1[[1]][["105"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(GFMID1[[1]][["105"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(GFMID1[[1]][["105"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[5],20))
  points(x,tmp[x],col=col[5],type='l',lwd=lwd)
  
  tmp<-apply(GFMID1[[1]][["1010"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(GFMID1[[1]][["1010"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(GFMID1[[1]][["1010"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[10],20))
  points(x,tmp[x],col=col[10],type='l',lwd=lwd)
  
  
  
  
  
  # main=paste0("Parr females maturation threshold (", pops[pop],")")
  # tmp <- apply(GFMID2[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
  # ylim=c(550,650)#c(0,max(ylim))
  # plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  # #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  # rect(0,0,10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
  # points(x,tmp[x],col=1,type='l',lwd=lwd)
  
  main=paste0("Anadromous males maturation threshold (", pops[pop],")")
  #tmp <- apply(GFMID3[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
  tmp<-apply(GFMID3[[1]][["101"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(GFMID3[[1]][["101"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(GFMID3[[1]][["101"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  ylim=c(20,70)#c(0,max(ylim))
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  rect(0,0,10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[1],20))
  points(x,tmp[x],col=1,type='l',lwd=lwd)
  
  tmp<-apply(GFMID3[[1]][["105"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(GFMID3[[1]][["105"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(GFMID3[[1]][["105"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[5],20))
  points(x,tmp[x],col=col[5],type='l',lwd=lwd)
  
  tmp<-apply(GFMID3[[1]][["1010"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(GFMID3[[1]][["1010"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(GFMID3[[1]][["1010"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[10],20))
  points(x,tmp[x],col=col[10],type='l',lwd=lwd)
  
  
  

  
  main=paste0("Anadromous females maturation threshold (", pops[pop],")")
  #tmp <- apply(GFMID4[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
  tmp<-apply(GFMID4[[1]][["101"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(GFMID4[[1]][["101"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(GFMID4[[1]][["101"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
   ylim=c(50,150)#c(0,max(ylim))
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  rect(0,0,10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[1],20))
  points(x,tmp[x],col=1,type='l',lwd=lwd)
  
  tmp<-apply(GFMID4[[1]][["105"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(GFMID4[[1]][["105"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(GFMID4[[1]][["105"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[5],20))
  points(x,tmp[x],col=col[5],type='l',lwd=lwd)
  
  tmp<-apply(GFMID4[[1]][["1010"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(GFMID4[[1]][["1010"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(GFMID4[[1]][["1010"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[10],20))
  points(x,tmp[x],col=col[10],type='l',lwd=lwd)
  
  
  
  
  
  
  main=paste0("Anadromous growth potential (", pops[pop],")")
  #tmp <- apply(GGROWTH[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
  tmp<-apply(GGROWTH[[1]][["101"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(GGROWTH[[1]][["101"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(GGROWTH[[1]][["101"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  ylim=c(-0.1,.1)#c(0,max(ylim))
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  rect(0,min(ylim),10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[1],20))
  points(x,tmp[x],col=col[1],type='l',lwd=lwd)
  
  tmp<-apply(GGROWTH[[1]][["105"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(GGROWTH[[1]][["105"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(GGROWTH[[1]][["105"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[7],20))
  points(x,tmp[x],col=col[5],type='l',lwd=lwd)
  
  tmp<-apply(GGROWTH[[1]][["1010"]][[pop]],1,median,na.rm=TRUE)
  tmp2.5<-apply(GGROWTH[[1]][["1010"]][[pop]],1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(GGROWTH[[1]][["1010"]][[pop]],1,quantile,probs=0.975,na.rm=TRUE)
  polygon(c(x,rev(x)),c(tmp2.5[x],rev(tmp97.5[x])),border=NA, col=makeTransparent(col[10],20))
  points(x,tmp[x],col=col[10],type='l',lwd=lwd)
  
  
  plot.new()
#} # end loop pop









par(mfrow=c(3,3))


#for (pop in 1:npop){
# pop=14
 
 xlab="local exploitation rate"
 ylab="gNEUTRAL"

 xlim=c(0,max(frates_vec))
ylim=c(0.4,0.6)#c(0,max(unlist(tmp), na.rm=TRUE))

  
  #### ALL
  main=paste0("h=1 / ALL")#paste0("Neutral gene (", pops[pop],")")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-gNEUTRAL.means[,,pop,1]
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
  var<-gNEUTRAL.means[,,pop,2]
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
  var<-gNEUTRAL.means[,,pop,5]
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
  var<-gNEUTRAL.means[,,pop,3]
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
  var<-gNEUTRAL.means[,,pop,6]
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
  var<-gNEUTRAL.means[,,pop,4]
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
  var<-gNEUTRAL.means[,,pop,7]
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
  ylab=paste0("Precocious males maturation threshold (", pops[pop],")")
  
  #tmp <- LFPARR
  ylim=c(0.8,1.8)
  
  #### ALL
  main=paste0("h=1 / ALL")#paste0("Neutral gene (", pops[pop],")")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-gFMID1.means[,,pop,1]
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
  var<-gFMID1.means[,,pop,2]
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
  var<-gFMID1.means[,,pop,5]
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
  var<-gFMID1.means[,,pop,3]
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
  var<-gFMID1.means[,,pop,6]
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
  var<-gFMID1.means[,,pop,4]
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
  var<-gFMID1.means[,,pop,7]
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
  ylab=paste0("Anadromous males maturation threshold (", pops[pop],")")
  
  #tmp <- LFPARR
  ylim=c(0,80)
  
  #### ALL
  main=paste0("h=1 / ALL")#paste0("Neutral gene (", pops[pop],")")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-gFMID3.means[,,pop,1]
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
  var<-gFMID3.means[,,pop,2]
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
  var<-gFMID3.means[,,pop,5]
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
  var<-gFMID3.means[,,pop,3]
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
  var<-gFMID3.means[,,pop,6]
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
  var<-gFMID3.means[,,pop,4]
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
  var<-gFMID3.means[,,pop,7]
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
  ylab=paste0("Anadromous females maturation threshold (", pops[pop],")")
  
  #tmp <- LFPARR
  ylim=c(50,150)
  
  #### ALL
  main=paste0("h=1 / ALL")#paste0("Neutral gene (", pops[pop],")")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-gFMID4.means[,,pop,1]
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
  var<-gFMID4.means[,,pop,2]
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
  var<-gFMID4.means[,,pop,5]
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
  var<-gFMID4.means[,,pop,3]
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
  var<-gFMID4.means[,,pop,6]
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
  var<-gFMID4.means[,,pop,4]
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
  var<-gFMID4.means[,,pop,7]
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
  ylab=paste0("Anadromous growth potential (", pops[pop],")")
  
  #tmp <- LFPARR
  ylim=c(-0.1,.1)
  
  #### ALL
  main=paste0("h=1 / ALL")#paste0("Neutral gene (", pops[pop],")")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-gGROWTH.means[,,pop,1]
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
  var<-gGROWTH.means[,,pop,2]
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
  var<-gGROWTH.means[,,pop,5]
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
  var<-gGROWTH.means[,,pop,3]
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
  var<-gGROWTH.means[,,pop,6]
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
  var<-gGROWTH.means[,,pop,4]
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
  var<-gGROWTH.means[,,pop,7]
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
  ylab=paste0("Anadromous growth potential (", pops[pop],")")
  
  xlim=c(0,0.5)
  ylim=c(-0.1,.1)
  
  #### ALL
  main=paste0("h=1 / ALL")#paste0("Neutral gene (", pops[pop],")")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-gGROWTH.means[,,pop,1]
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
  var<-gGROWTH.means[,,pop,2]
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
  var<-gGROWTH.means[,,pop,5]
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
  var<-gGROWTH.means[,,pop,3]
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
  var<-gGROWTH.means[,,pop,6]
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
  var<-gGROWTH.means[,,pop,4]
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
  var<-gGROWTH.means[,,pop,7]
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
  ylab=paste0("Anadromous male threshold (", pops[pop],")")
  
  xlim=c(0,0.5)
  ylim=c(30,60)
  #### ALL
  main=paste0("h=1 / ALL")#paste0("Neutral gene (", pops[pop],")")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-gFMID3.means[,,pop,1]
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
  var<-gFMID3.means[,,pop,2]
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
  var<-gFMID3.means[,,pop,5]
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
  var<-gFMID3.means[,,pop,3]
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
  var<-gFMID3.means[,,pop,6]
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
  var<-gFMID3.means[,,pop,4]
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
  var<-gFMID3.means[,,pop,7]
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
  ylab=paste0("Anadromous female threshold (", pops[pop],")")
  
  xlim=c(0,0.5)
  ylim=c(50,150)
  #### ALL
  main=paste0("h=1 / ALL")#paste0("Neutral gene (", pops[pop],")")
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
  axis(1,frates_vec,frates_vec)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  var<-gFMID4.means[,,pop,1]
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
  var<-gFMID4.means[,,pop,2]
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
  var<-gFMID4.means[,,pop,5]
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
  var<-gFMID4.means[,,pop,3]
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
  var<-gFMID4.means[,,pop,6]
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
  var<-gFMID4.means[,,pop,4]
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
  var<-gFMID4.means[,,pop,7]
  var2<-Pglobal.means[,,7];ExpGl<-apply(var2,1,median,na.rm=TRUE)
  tmp<-apply(var,1,median,na.rm=TRUE)
  tmp2.5<-apply(var,1,quantile,probs=0.025,na.rm=TRUE)
  tmp97.5<-apply(var,1,quantile,probs=0.975,na.rm=TRUE)
  abline(h=tmp[1],lty=2,col="grey")
  segments(ExpGl,tmp2.5,ExpGl,tmp97.5,col=col[1:length(ExpGl)])
  points(ExpGl,tmp,col=col[1:length(ExpGl)],type='p',lwd=lwd,pch=20)
  
  
  
  
  
} # end loop pop
  dev.off()
  