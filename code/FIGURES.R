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


##### METAPOPULATION ABUNDANCE

#EXPE=c(101, 102, 103, 104, 105, 106)
#401 402 403 404 405 406 )

Means=MeansRel=Sds=matrix(NA,nrow(SCN),ncol(SCN))
metapops.means=metapops.means.rel=metapops.exp.means=Pglobal.means=extinction.means.rel=array(,dim=c(nrow(SCN),nSIMUL,ncol(SCN)))
Metapop<-Exploitation<-Pexpl.final<-Extinction<-Prop.Immigrants<-Population<-RatioImExpl<-list()
pops.means=migs.means=migexpl.means=array(,dim=c(nrow(SCN),nSIMUL,npop,ncol(SCN)))
for (i in 1:ncol(SCN)){
  res=res2=res3=list()
  res.pop=prop.mig=res.expl=list(list())
  j=0
  #for (iEXPE in EXPE){ # Loop over scenario 
  for (iEXPE in SCN[,i]){
    j=j+1
    if(is.na(iEXPE)) next;
    
    load(paste0("results/DEMOGRAPHY",iEXPE,".RData"))
    
    scn_id <- as.numeric(strsplit(as.character(iEXPE), "")[[1]])
    scenarioConnect = scn_id[1]
    scenarioFishing = scn_id[2]
    scenarioFishingRate = scn_id[3]
    
    
    #### POP
    for (pop in 1:npop){
      pop.tmp=mig.tmp=expl.tmp=NULL
      for (simul in 1:nSIMUL){
        tmp=tmp2=tmp3=NULL
        tmp <- nReturns[[paste0(iEXPE)]][[simul]][,pop]
        tmp2 <- nReturns[[paste0(iEXPE)]][[simul]][,pop]/nReturns[[paste0(iEXPE)]][[simul]][1,pop]
        tmp3 <- Mig[[paste0(iEXPE)]][[simul]]$NIm[,pop]/tmp
        tmp4 <- Mig[[paste0(iEXPE)]][[simul]]$NIm[,pop]/nExploitation[[paste0(iEXPE)]][[simul]][,pop]
        
        pops.means[j,simul,pop,i] <- mean(tail(tmp2,5))
        migs.means[j,simul,pop,i] <- mean(tail(tmp3,5))
        migexpl.means[j,simul,pop,i] <- mean(tail(tmp4,5))
        
        pop.tmp <- cbind(pop.tmp, tmp2)
        mig.tmp <- cbind(mig.tmp, tmp3)
        expl.tmp <- cbind(expl.tmp, tmp4)
        
      }#end loop simul
      res.pop[[paste0(iEXPE)]][[pop]]<-pop.tmp
      prop.mig[[paste0(iEXPE)]][[pop]]<-mig.tmp
      res.expl[[paste0(iEXPE)]][[pop]]<-expl.tmp
    }#end loop pop
    
    
    #### METAPOP
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
    
    #Means[j,i] <- mean(tail(metapops,5))
    #Sds[j,i] <- sd(tail(metapops,5))
    
    #MeansRel[j,i] <- Means[j,i]/Means[1,i]
    #res[[paste0(iEXPE)]]<- list(Means=Means,Sds=Sds)
    metapops.rel <- metapops/metapops[1,]
    res[[paste0(iEXPE)]]<-apply(metapops.rel,1,mean, na.rm=TRUE)
    res2[[paste0(iEXPE)]]<-apply(Pglobal,1,mean, na.rm=TRUE)
    res3[[paste0(iEXPE)]] <- apply(extinction,1,mean, na.rm=TRUE)
  } # loop iEXPE
  Population[[i]]<-res.pop
  Prop.Immigrants[[i]]<-prop.mig
  Metapop[[i]]<-res
  Exploitation[[i]]<-res2
  Pexpl.final[[i]]<- Pglobal.means
  Extinction[[i]]<-res3
  RatioImExpl[[i]]<-res.expl
} # loop SCN






load(paste0("results/DEMOGRAPHY",101,".RData"))
init <- apply(tail(nReturns[[paste0(101)]][[1]],10),2,mean,na.rm=TRUE)
tmp<-init/sum(init)
prop.sink <- sum(tmp[dat$Type!="source"],na.rm=TRUE)
prop.source <- sum(tmp[dat$Type!="sink"],na.rm=TRUE)




pdf(file="results/FishStrat_results.pdf")

#### METAPOPULATION DYNAMICS ####
#-------------------------------#



par(mfrow=c(3,3))

tmp <- Metapop




xlab="Time (years)"

#ylab="Metapopulation Abundance"
#ylab="PFA"
#ylab="Metapopulation Relative Abundance"
ylab="Relative PFA" # Pre-Fisheries Abundance

# Create a list of gradients for each colour 2 to 10 over five steps from 
library(viridis)
#col <- viridis(13)
col <- (inferno(length(frates_vec)))
#col=c("#000000", "#1874CD", "#009ACD", "#00BFFF", "#228B22", "#66CD00", "#7FFF00")

lwd=1

plot(NULL,xlim=c(1,nInit+nYear),ylim=c(0,max(unlist(tmp), na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry=1 / Fishing all")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[1]])){
  points(1:50,tmp[[1]][[j]],col=col[j],type='l',lwd=lwd)
}
#legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)

plot(NULL,xlim=c(1,nInit+nYear),ylim=c(0,max(unlist(tmp), na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry=0.85 / Fishing all")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[2]])){
  points(1:50,tmp[[2]][[j]],col=col[j],type='l',lwd=lwd)
}
#legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)


plot(NULL,xlim=c(1,nInit+nYear),ylim=c(0,max(unlist(tmp), na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry=0.7 / Fishing all")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[5]])){
  points(1:50,tmp[[5]][[j]],col=col[j],type='l',lwd=lwd)
}
#legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n")




#plot.new()
plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = "Local Exploitation rates")
#legend("bottomleft",paste0(frates_vec),col=col,lty=rep(1,length(frates_vec)),bty="n",title="Local Exploitation rates",cex=.6)
legend_image <- as.raster(matrix(rev(col), ncol=1))
text(x=0.5, y = 1.1, "Local Exploitation rates", cex=.6)
text(x=0.45, y = seq(0,1,l=length(frates_vec)), labels = paste0(frates_vec),adj=.2,cex=.5)
rasterImage(legend_image, 0.3, 0, .4,1)


plot(NULL,xlim=c(1,nInit+nYear),ylim=c(0,max(unlist(tmp), na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry=0.85 / Fishing sources")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[3]])){
  points(1:50,tmp[[3]][[j]],col=col[j],type='l',lwd=lwd)
}
abline(h=1-prop.source,lty=2)
#legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)

plot(NULL,xlim=c(1,nInit+nYear),ylim=c(0,max(unlist(tmp), na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry=0.7 / Fishing sources")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[6]])){
  points(1:50,tmp[[6]][[j]],col=col[j],type='l',lwd=lwd)
}
abline(h=1-prop.source,lty=2)
#legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n")


plot.new()
# legend_image <- as.raster(matrix(rev(col), ncol=1))
# plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
# text(x=0.1, y = seq(0,1,l=length(frates_vec)), labels = paste0(frates_vec),adj=.2,cex=.7)
# rasterImage(legend_image, 0, 0, .1,1)

plot(NULL,xlim=c(1,nInit+nYear),ylim=c(0,max(unlist(tmp), na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry=0.85 / Fishing sinks")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[4]])){
  points(1:50,tmp[[4]][[j]],col=col[j],type='l',lwd=lwd)
}
abline(h=1-prop.sink,lty=2)
#legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)

plot(NULL,xlim=c(1,nInit+nYear),ylim=c(0,max(unlist(tmp), na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry=0.7 / Fishing sinks")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[7]])){
  points(1:50,tmp[[7]][[j]],col=col[j],type='l',lwd=lwd)
}
abline(h=1-prop.sink,lty=2)
# #legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")

#legend("bottomleft",paste0(frates_vec),col=col[1:length(frates_vec)],lty=rep(1,length(frates_vec)),bty="n")






#### FISHING STRATEGIES / LOCAL EXPLOITATION RATES ####
#-------------------------------#



#tmp <- metapops.means
tmp <- metapops.means.rel


#col=c("#030303", "#454545", "#8C8C8C", "#CCCCCC", "#E3E3E3", "#F7F7F7")
#col=c("#000000", "#1874CD", "#009ACD", "#00BFFF", "#228B22", "#66CD00", "#7FFF00")
col.sources=c("dodgerblue2", "deepskyblue")
col.sinks=c("brown1", "lightcoral", "white")
par(mfcol=c(3,2))

lwd=3

xlab="Local Exploitation rates"

### FISHING ALL
plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Fishing all",xaxt='n')
axis(1,1:length(frates_vec), frates_vec)
#for (i in 1:ncol(SCN)){
points(1:length(frates_vec),apply(tmp[,,1],1,mean),col=col[1],pch=16,bg=col[1],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,2],1,mean),col="grey",pch=16,bg="grey",cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,5],1,mean),col="lightgrey",pch=16,bg="lightgrey",cex=2,type='l',lwd=lwd)
# for (j in 1:13){
#   #points(rep(j,nSIMUL),tmp[j,,1],col=col[1],pch=1,bg=col[1])
#   points(j,mean(tmp[j,,1]),col=col[1],pch=16,bg=col[1],cex=2,type='l')
#   #segments(j,mean(tmp[j,,1])-sd(tmp[j,,1]),j,mean(tmp[j,,1])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,1]
#   #lines(lowess(x,y,f = 1/3), col = col[1], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,2],col=col[2],pch=1,bg=col[2])
#   points(j,mean(tmp[j,,2]),col=col[2],pch=16,bg=col[2],cex=2)
#   #segments(j,mean(tmp[j,,2])-sd(tmp[j,,2]),j,mean(tmp[j,,2])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,2]
#   #lines(lowess(x,y,f = 1/3), col = col[2], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,5],col=col[5],pch=1,bg=col[5])
#   points(j,mean(tmp[j,,5]),col=col[5],pch=16,bg=col[5],cex=2)
#   #segments(j,mean(tmp[j,,5])-sd(tmp[j,,5]),j,mean(tmp[j,,5])-sd(tmp[j,,5]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,5]
#   #lines(lowess(x,y,f = 1/3), col = col[5], lty = 1)
# }
legend("bottomleft",c("h=1","h=0.85","h=0.7"),col=c(1,"grey","lightgrey"),lty=rep(1,3),lwd=rep(2,3),bty="n")


### FISHING SOURCES
plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Fishing source",xaxt='n')
axis(1,1:length(frates_vec), frates_vec)
points(1:length(frates_vec),apply(tmp[,,1],1,mean),col=col[1],pch=16,bg=col[1],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,3],1,mean),col=col.sources[1],pch=16,bg=col.sources[1],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,6],1,mean),col=col.sources[2],pch=16,bg=col.sources[2],cex=2,type='l',lwd=lwd)
#for (i in 1:ncol(SCN)){
# for (j in 1:13){
#   #points(rep(j,nSIMUL),tmp[j,,1],col=col[1],pch=16,bg=col[1])
#   points(j,mean(tmp[j,,1]),col=col[1],pch=16,bg=col[1],cex=2)
#   #segments(j,mean(tmp[j,,1])-sd(tmp[j,,1]),j,mean(tmp[j,,1])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,1]
#   #lines(lowess(x,y,f = 1/3), col = col[1], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,3],col=col[3],pch=16,bg=col[3])
#   points(j,mean(tmp[j,,3]),col=col[2],pch=16,bg=col[2],cex=2)
#   #segments(j,mean(tmp[j,,3])-sd(tmp[j,,3]),j,mean(tmp[j,,3])-sd(tmp[j,,3]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,3]
#   #lines(lowess(x,y,f = 1/3), col = col[2], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,5],col=col[5],pch=1,bg=col[5])
#   points(j,mean(tmp[j,,7]),col=col[5],pch=16,bg=col[5],cex=2)
#   #segments(j,mean(tmp[j,,5])-sd(tmp[j,,5]),j,mean(tmp[j,,5])-sd(tmp[j,,5]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,13]
#   #lines(lowess(x,y,f = 1/3), col = col[5], lty = 1)
# }
legend("bottomleft",c("h=1","h=0.85","h=0.7"),col=c(1,col.sources),lty=rep(1,3),bty="n")


### FISHING SINKS
plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Fishing sinks",xaxt='n')
axis(1,1:length(frates_vec), frates_vec)
points(1:length(frates_vec),apply(tmp[,,1],1,mean),col=col[1],pch=16,bg=col[1],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,4],1,mean),col=col.sinks[1],pch=16,bg=col.sinks[1],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,7],1,mean),col=col.sinks[2],pch=16,bg=col.sinks[2],cex=2,type='l',lwd=lwd)
#for (i in 1:ncol(SCN)){
# for (j in 1:13){
#   #points(rep(j,nSIMUL),tmp[j,,1],col=col[1],pch=16,bg=col[1])
#   points(j,mean(tmp[j,,1]),col=col[1],pch=16,bg=col[1],cex=2)
#   #segments(j,mean(tmp[j,,1])-sd(tmp[j,,1]),j,mean(tmp[j,,1])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,1]
#   #lines(lowess(x,y,f = 1/3), col = col[1], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,4],col=col[2],pch=1,bg=col[2])
#   points(j,mean(tmp[j,,4]),col=col[2],pch=16,bg=col[2],cex=2)
#   #segments(j,mean(tmp[j,,3])-sd(tmp[j,,3]),j,mean(tmp[j,,3])-sd(tmp[j,,3]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,4]
#   #lines(lowess(x,y,f = 1/3), col = col[2], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,5],col=col[5],pch=1,bg=col[5])
#   points(j,mean(tmp[j,,7]),col=col[5],pch=16,bg=col[5],cex=2)
#   #segments(j,mean(tmp[j,,5])-sd(tmp[j,,5]),j,mean(tmp[j,,5])-sd(tmp[j,,5]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,7]
#   #lines(lowess(x,y,f = 1/3), col = col[5], lty = 1)
# }
legend("bottomleft",c("h=1","h=0.85","h=0.7"),col=c(1,col.sinks),lty=rep(1,3),lwd=rep(2,3),bty="n")



# plot(NULL,xlim=c(1,13),ylim=c(0,max(MeansRel)),xlab=xlab,ylab=ylab)
# for (i in 1:ncol(SCN)){
#   points(1:13,MeansRel[,i],col=col[i], type='l')
#   points(1:6,MeansRel[,i],col=col[i], pch=paste0(i))
#   #segments(1:6,MeansRel[,i]-Sds[,i],1:6,Means[,i]+Sds[,i],col=col[i])
# }









# Npops <- array(,dim=c((nYears+nInit),nSIMUL,length(nReturns)))
# for (scn in 1:length(nReturns)){
#   #for (div in 1:length(DIV)){
#   for (simul in 1:nSIMUL){
#     Npops[,simul,scn] <- nReturns[[scn]][[1]]$NRet[[simul]][,13] #nb returns metapop per year, simu and scenario
#   }
#   #}
# }




#par(mfrow=c(2,1))

### DISPERSAL 1
plot.new()


### DISPERSAL 0.85 FUNCTION OF LOCAL EXPLOITATION RATE
plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry = 0.85",xaxt='n')
axis(1,1:length(frates_vec), frates_vec)
points(1:length(frates_vec),apply(tmp[,,2],1,mean),col="grey",pch=16,bg="grey",cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,3],1,mean),col=col.sources[1],pch=16,bg=col.sources[1],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,4],1,mean),col=col.sinks[1],pch=16,bg=col.sinks[1],cex=2,type='l',lwd=lwd)
#for (i in 1:ncol(SCN)){
# for (j in 1:13){
#   #points(rep(j,nSIMUL),tmp[j,,1],col=col[1],pch=1,bg=col[1])
#   points(j,mean(tmp[j,,2]),col=col[2],pch=16,bg=col[2],cex=2)
#   #segments(j,mean(tmp[j,,1])-sd(tmp[j,,1]),j,mean(tmp[j,,1])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:6,6,nSIMUL);y<-tmp[,,2]
#   #lines(lowess(x,y,f = 1/3), col = col[2], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,2],col=col[2],pch=1,bg=col[2])
#   points(j,mean(tmp[j,,3]),col=col[3],pch=16,bg=col[3],cex=2)
#   #segments(j,mean(tmp[j,,2])-sd(tmp[j,,2]),j,mean(tmp[j,,2])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:6,6,nSIMUL);y<-tmp[,,3]
#   #lines(lowess(x,y,f = 1/3), col = col[3], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,5],col=col[5],pch=1,bg=col[5])
#   points(j,mean(tmp[j,,4]),col=col[4],pch=16,bg=col[4],cex=2)
#   #segments(j,mean(tmp[j,,5])-sd(tmp[j,,5]),j,mean(tmp[j,,5])-sd(tmp[j,,5]),col=col[i])
#   #x<-matrix(1:6,6,nSIMUL);y<-tmp[,,4]
#   #lines(lowess(x,y,f = 1/3), col = col[4], lty = 1)
# }
legend("bottomleft",c("all","source","sink"),col=c("grey",col.sources[1],col.sinks[1]),lty=rep(1,3),lwd=rep(2,3),bty="n")


### DISPERSAL 0.7 FUNCTION OF LOCAL EXPLOITATION RATE
plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry = 0.7",xaxt='n')
axis(1,1:length(frates_vec), frates_vec)
points(1:length(frates_vec),apply(tmp[,,5],1,mean),col="lightgrey",pch=16,bg="grey",cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,6],1,mean),col=col.sources[2],pch=16,bg=col.sources[2],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,7],1,mean),col=col.sinks[2],pch=16,bg=col.sinks[2],cex=2,type='l',lwd=lwd)
#for (i in 1:ncol(SCN)){
# for (j in 1:13){
#   #points(rep(j,nSIMUL),tmp[j,,1],col=col[1],pch=1,bg=col[1])
#   points(j,mean(tmp[j,,5]),col=col[5],pch=16,bg=col[5],cex=2)
#   #segments(j,mean(tmp[j,,1])-sd(tmp[j,,1]),j,mean(tmp[j,,1])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:6,6,nSIMUL);y<-tmp[,,5]
#   #lines(lowess(x,y,f = 1/3), col = col[5], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,2],col=col[2],pch=1,bg=col[2])
#   points(j,mean(tmp[j,,6]),col=col[6],pch=16,bg=col[6],cex=2)
#   #segments(j,mean(tmp[j,,2])-sd(tmp[j,,2]),j,mean(tmp[j,,2])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:6,6,nSIMUL);y<-tmp[,,6]
#   #lines(lowess(x,y,f = 1/3), col = col[6], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,5],col=col[5],pch=1,bg=col[5])
#   points(j,mean(tmp[j,,7]),col=col[7],pch=16,bg=col[7],cex=2)
#   #segments(j,mean(tmp[j,,5])-sd(tmp[j,,5]),j,mean(tmp[j,,5])-sd(tmp[j,,5]),col=col[i])
#   #x<-matrix(1:6,6,nSIMUL);y<-tmp[,,7]
#   #lines(lowess(x,y,f = 1/3), col = col[7], lty = 1)
# }
legend("bottomleft",c("all","source","sink"),col=c("lightgrey",col.sources[2],col.sinks[2]),lty=rep(1,3),lwd=rep(2,3),bty="n")








#### EXPLOITATION ####
#---------------------------#


par(mfrow=c(3,3))
#col=c("#000000", "#1874CD", "#009ACD", "#00BFFF", "#228B22", "#66CD00", "#7FFF00")
tmp <- Exploitation

xlab="Time (years)"

ylab="Metapopulation Exploitation Rate"

lwd=1

plot(NULL,xlim=c(1,nInit+nYear),ylim=c(0,max(unlist(tmp), na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry=1 / Fishing all")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[1]])){
  points(1:50,tmp[[1]][[j]],col=col[j],type='l',lwd=lwd)
}
#legend("topleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("topleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)

plot(NULL,xlim=c(1,nInit+nYear),ylim=c(0,max(unlist(tmp), na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry=0.85 / Fishing all")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[2]])){
  points(1:50,tmp[[2]][[j]],col=col[j],type='l',lwd=lwd)
}
#legend("topleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("topleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)


plot(NULL,xlim=c(1,nInit+nYear),ylim=c(0,max(unlist(tmp), na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry=0.7 / Fishing all")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[5]])){
  points(1:50,tmp[[5]][[j]],col=col[j],type='l',lwd=lwd)
}
#legend("topleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("topleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n")




#plot.new()
#legend("topleft",paste0(frates_vec),col=col,lty=rep(1,length(frates_vec)),bty="n",title="Local Exploitation rates",cex=.6)
plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = "Local Exploitation rates")
#legend("bottomleft",paste0(frates_vec),col=col,lty=rep(1,length(frates_vec)),bty="n",title="Local Exploitation rates",cex=.6)
legend_image <- as.raster(matrix(rev(col), ncol=1))
text(x=0.5, y = 1.1, "Local Exploitation rates", cex=.6)
text(x=0.45, y = seq(0,1,l=length(frates_vec)), labels = paste0(frates_vec),adj=.2,cex=.5)
rasterImage(legend_image, 0.3, 0, .4,1)

plot(NULL,xlim=c(1,nInit+nYear),ylim=c(0,max(unlist(tmp), na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry=0.85 / Fishing sources")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[3]])){
  points(1:50,tmp[[3]][[j]],col=col[j],type='l',lwd=lwd)
}
#legend("topleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("topleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)

plot(NULL,xlim=c(1,nInit+nYear),ylim=c(0,max(unlist(tmp), na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry=0.7 / Fishing sources")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[6]])){
  points(1:50,tmp[[6]][[j]],col=col[j],type='l',lwd=lwd)
}
#legend("topleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("topleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n")


plot.new()

plot(NULL,xlim=c(1,nInit+nYear),ylim=c(0,max(unlist(tmp), na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry=0.85 / Fishing sinks")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[4]])){
  points(1:50,tmp[[4]][[j]],col=col[j],type='l',lwd=lwd)
}
#legend("topleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("topleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)

plot(NULL,xlim=c(1,nInit+nYear),ylim=c(0,max(unlist(tmp), na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry=0.7 / Fishing sinks")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[7]])){
  points(1:50,tmp[[7]][[j]],col=col[j],type='l',lwd=lwd)
}
# #legend("topleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
# legend("topleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n")










tmp <- Pglobal.means

xlab="Local Exploitation rates"
ylab="Metapopulation Exploitation Rate"

#col=c("#030303", "#454545", "#8C8C8C", "#CCCCCC", "#E3E3E3", "#F7F7F7")
#col=c("#000000", "#1874CD", "#009ACD", "#00BFFF", "#228B22", "#66CD00", "#7FFF00")
col.sources=c("dodgerblue2", "deepskyblue")
col.sinks=c("brown1", "lightcoral", "white")
par(mfcol=c(3,2))

lwd=3


### FISHING ALL
plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Fishing all",xaxt='n')
axis(1,1:length(frates_vec), frates_vec)
#for (i in 1:ncol(SCN)){
points(1:length(frates_vec),apply(tmp[,,1],1,mean),col=col[1],pch=16,bg=col[1],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,2],1,mean),col="grey",pch=16,bg="grey",cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,5],1,mean),col="lightgrey",pch=16,bg="lightgrey",cex=2,type='l',lwd=lwd)
# for (j in 1:13){
#   #points(rep(j,nSIMUL),tmp[j,,1],col=col[1],pch=1,bg=col[1])
#   points(j,mean(tmp[j,,1]),col=col[1],pch=16,bg=col[1],cex=2,type='l')
#   #segments(j,mean(tmp[j,,1])-sd(tmp[j,,1]),j,mean(tmp[j,,1])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,1]
#   #lines(lowess(x,y,f = 1/3), col = col[1], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,2],col=col[2],pch=1,bg=col[2])
#   points(j,mean(tmp[j,,2]),col=col[2],pch=16,bg=col[2],cex=2)
#   #segments(j,mean(tmp[j,,2])-sd(tmp[j,,2]),j,mean(tmp[j,,2])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,2]
#   #lines(lowess(x,y,f = 1/3), col = col[2], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,5],col=col[5],pch=1,bg=col[5])
#   points(j,mean(tmp[j,,5]),col=col[5],pch=16,bg=col[5],cex=2)
#   #segments(j,mean(tmp[j,,5])-sd(tmp[j,,5]),j,mean(tmp[j,,5])-sd(tmp[j,,5]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,5]
#   #lines(lowess(x,y,f = 1/3), col = col[5], lty = 1)
# }
legend("topleft",c("h=1","h=0.85","h=0.7"),col=c(1,"grey","lightgrey"),lty=rep(1,3),lwd=rep(2,3),bty="n")


### FISHING SOURCES
plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Fishing source",xaxt='n')
axis(1,1:length(frates_vec), frates_vec)
points(1:length(frates_vec),apply(tmp[,,1],1,mean),col=col[1],pch=16,bg=col[1],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,3],1,mean),col=col.sources[1],pch=16,bg=col.sources[1],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,6],1,mean),col=col.sources[2],pch=16,bg=col.sources[2],cex=2,type='l',lwd=lwd)
#for (i in 1:ncol(SCN)){
# for (j in 1:13){
#   #points(rep(j,nSIMUL),tmp[j,,1],col=col[1],pch=16,bg=col[1])
#   points(j,mean(tmp[j,,1]),col=col[1],pch=16,bg=col[1],cex=2)
#   #segments(j,mean(tmp[j,,1])-sd(tmp[j,,1]),j,mean(tmp[j,,1])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,1]
#   #lines(lowess(x,y,f = 1/3), col = col[1], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,3],col=col[3],pch=16,bg=col[3])
#   points(j,mean(tmp[j,,3]),col=col[2],pch=16,bg=col[2],cex=2)
#   #segments(j,mean(tmp[j,,3])-sd(tmp[j,,3]),j,mean(tmp[j,,3])-sd(tmp[j,,3]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,3]
#   #lines(lowess(x,y,f = 1/3), col = col[2], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,5],col=col[5],pch=1,bg=col[5])
#   points(j,mean(tmp[j,,7]),col=col[5],pch=16,bg=col[5],cex=2)
#   #segments(j,mean(tmp[j,,5])-sd(tmp[j,,5]),j,mean(tmp[j,,5])-sd(tmp[j,,5]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,13]
#   #lines(lowess(x,y,f = 1/3), col = col[5], lty = 1)
# }
legend("topleft",c("h=1","h=0.85","h=0.7"),col=c(1,col.sources),lty=rep(1,3),lwd=rep(2,3),bty="n")


### FISHING SINKS
plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Fishing sinks",xaxt='n')
axis(1,1:length(frates_vec), frates_vec)
points(1:length(frates_vec),apply(tmp[,,1],1,mean),col=col[1],pch=16,bg=col[1],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,4],1,mean),col=col.sinks[1],pch=16,bg=col.sinks[1],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,7],1,mean),col=col.sinks[2],pch=16,bg=col.sinks[2],cex=2,type='l',lwd=lwd)
#for (i in 1:ncol(SCN)){
# for (j in 1:13){
#   #points(rep(j,nSIMUL),tmp[j,,1],col=col[1],pch=16,bg=col[1])
#   points(j,mean(tmp[j,,1]),col=col[1],pch=16,bg=col[1],cex=2)
#   #segments(j,mean(tmp[j,,1])-sd(tmp[j,,1]),j,mean(tmp[j,,1])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,1]
#   #lines(lowess(x,y,f = 1/3), col = col[1], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,4],col=col[2],pch=1,bg=col[2])
#   points(j,mean(tmp[j,,4]),col=col[2],pch=16,bg=col[2],cex=2)
#   #segments(j,mean(tmp[j,,3])-sd(tmp[j,,3]),j,mean(tmp[j,,3])-sd(tmp[j,,3]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,4]
#   #lines(lowess(x,y,f = 1/3), col = col[2], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,5],col=col[5],pch=1,bg=col[5])
#   points(j,mean(tmp[j,,7]),col=col[5],pch=16,bg=col[5],cex=2)
#   #segments(j,mean(tmp[j,,5])-sd(tmp[j,,5]),j,mean(tmp[j,,5])-sd(tmp[j,,5]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,7]
#   #lines(lowess(x,y,f = 1/3), col = col[5], lty = 1)
# }
legend("topleft",c("h=1","h=0.85","h=0.7"),col=c(1,col.sinks),lty=rep(1,3),lwd=rep(2,3),bty="n")



# plot(NULL,xlim=c(1,13),ylim=c(0,max(MeansRel)),xlab=xlab,ylab=ylab)
# for (i in 1:ncol(SCN)){
#   points(1:13,MeansRel[,i],col=col[i], type='l')
#   points(1:6,MeansRel[,i],col=col[i], pch=paste0(i))
#   #segments(1:6,MeansRel[,i]-Sds[,i],1:6,Means[,i]+Sds[,i],col=col[i])
# }









# Npops <- array(,dim=c((nYears+nInit),nSIMUL,length(nReturns)))
# for (scn in 1:length(nReturns)){
#   #for (div in 1:length(DIV)){
#   for (simul in 1:nSIMUL){
#     Npops[,simul,scn] <- nReturns[[scn]][[1]]$NRet[[simul]][,length(frates_vec)] #nb returns metapop per year, simu and scenario
#   }
#   #}
# }




#par(mfrow=c(2,1))

### DISPERSAL 1
plot.new()


### DISPERSAL 0.85 FUNCTION OF LOCAL EXPLOITATION RATE
plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry = 0.85",xaxt='n')
axis(1,1:length(frates_vec), frates_vec)
points(1:length(frates_vec),apply(tmp[,,2],1,mean),col="grey",pch=16,bg="grey",cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,3],1,mean),col=col.sources[1],pch=16,bg=col.sources[1],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,4],1,mean),col=col.sinks[1],pch=16,bg=col.sinks[1],cex=2,type='l',lwd=lwd)
#for (i in 1:ncol(SCN)){
# for (j in 1:13){
#   #points(rep(j,nSIMUL),tmp[j,,1],col=col[1],pch=1,bg=col[1])
#   points(j,mean(tmp[j,,2]),col=col[2],pch=16,bg=col[2],cex=2)
#   #segments(j,mean(tmp[j,,1])-sd(tmp[j,,1]),j,mean(tmp[j,,1])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:6,6,nSIMUL);y<-tmp[,,2]
#   #lines(lowess(x,y,f = 1/3), col = col[2], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,2],col=col[2],pch=1,bg=col[2])
#   points(j,mean(tmp[j,,3]),col=col[3],pch=16,bg=col[3],cex=2)
#   #segments(j,mean(tmp[j,,2])-sd(tmp[j,,2]),j,mean(tmp[j,,2])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:6,6,nSIMUL);y<-tmp[,,3]
#   #lines(lowess(x,y,f = 1/3), col = col[3], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,5],col=col[5],pch=1,bg=col[5])
#   points(j,mean(tmp[j,,4]),col=col[4],pch=16,bg=col[4],cex=2)
#   #segments(j,mean(tmp[j,,5])-sd(tmp[j,,5]),j,mean(tmp[j,,5])-sd(tmp[j,,5]),col=col[i])
#   #x<-matrix(1:6,6,nSIMUL);y<-tmp[,,4]
#   #lines(lowess(x,y,f = 1/3), col = col[4], lty = 1)
# }
legend("topleft",c("all","source","sink"),col=c("grey",col.sources[1],col.sinks[1]),lty=rep(1,3),lwd=rep(2,3),bty="n")


### DISPERSAL 0.7 FUNCTION OF LOCAL EXPLOITATION RATE
plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main=" Philopatry = 0.7",xaxt='n')
axis(1,1:length(frates_vec), frates_vec)
points(1:length(frates_vec),apply(tmp[,,5],1,mean),col="lightgrey",pch=16,bg="grey",cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,6],1,mean),col=col.sources[2],pch=16,bg=col.sources[2],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,7],1,mean),col=col.sinks[2],pch=16,bg=col.sinks[2],cex=2,type='l',lwd=lwd)
#for (i in 1:ncol(SCN)){
# for (j in 1:13){
#   #points(rep(j,nSIMUL),tmp[j,,1],col=col[1],pch=1,bg=col[1])
#   points(j,mean(tmp[j,,5]),col=col[5],pch=16,bg=col[5],cex=2)
#   #segments(j,mean(tmp[j,,1])-sd(tmp[j,,1]),j,mean(tmp[j,,1])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:6,6,nSIMUL);y<-tmp[,,5]
#   #lines(lowess(x,y,f = 1/3), col = col[5], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,2],col=col[2],pch=1,bg=col[2])
#   points(j,mean(tmp[j,,6]),col=col[6],pch=16,bg=col[6],cex=2)
#   #segments(j,mean(tmp[j,,2])-sd(tmp[j,,2]),j,mean(tmp[j,,2])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:6,6,nSIMUL);y<-tmp[,,6]
#   #lines(lowess(x,y,f = 1/3), col = col[6], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,5],col=col[5],pch=1,bg=col[5])
#   points(j,mean(tmp[j,,7]),col=col[7],pch=16,bg=col[7],cex=2)
#   #segments(j,mean(tmp[j,,5])-sd(tmp[j,,5]),j,mean(tmp[j,,5])-sd(tmp[j,,5]),col=col[i])
#   #x<-matrix(1:6,6,nSIMUL);y<-tmp[,,7]
#   #lines(lowess(x,y,f = 1/3), col = col[7], lty = 1)
# }
legend("topleft",c("all","source","sink"),col=c("lightgrey",col.sources[2],col.sinks[2]),lty=rep(1,3),lwd=rep(2,3),bty="n")









#### GLOBAL EXPLOITATION ####
#---------------------------#

# getFishingScn <- function(scenarioFishing, scenarioFishingRate){
# dat <- read.csv2("/media/hdd4To/mbuoro/MetaIBASAM-Projects/FishStrat/data/dataPop.csv", header=TRUE, stringsAsFactors = FALSE, comment.char = "#")
# dat <- dat[-which(dat$Population == "Couesnon"),] # remove Couesnon
# 
# frates = c(0, frates_vec[scenarioFishingRate], 0, frates_vec[scenarioFishingRate]) # (1SW_init, 1SW, MSW_init, MSW)
# 
# fratesSource = frates # fishing rates by stage (1SW vs MSW) for source populations (13/04/2021 : ajout du taux d'exploitation initial et au bout de 10ans) 
# fratesSink = frates # fishing rates by stage (1SW vs MSW) for sink populations
# fratesNeutral = frates # fishing rates by stage (1SW vs MSW) for neutral populations
# 
# ## 3.1 Fishing scenarios
# # fishing rates (frates) are provided in parIbasam.R file
# tmp <- matrix(0, nrow = npop, ncol = length(fratesSink))
# if (scenarioFishing==0){ # fishing all
#   for (i in 1:npop) {
#     if (dat$Type[i]=="sink") tmp[i,] <- fratesSink
#     if (dat$Type[i]=="neutral") tmp[i,] <- fratesNeutral
#     if (dat$Type[i]=="source") tmp[i,] <- fratesSource
#   }}
# 
# if (scenarioFishing==1){ # fishing sources/neutral
#   for (i in 1:npop) {
#     if (dat$Type[i]=="sink") tmp[i,] <- 0
#     if (dat$Type[i]=="neutral") tmp[i,] <- fratesNeutral
#     if (dat$Type[i]=="source") tmp[i,] <- fratesSource
#   }}
# 
# if (scenarioFishing==2){ # fishing sinks/neutral
#   for (i in 1:npop) {
#     if (dat$Type[i]=="sink") tmp[i,] <- fratesSink
#     if (dat$Type[i]=="neutral") tmp[i,] <- fratesNeutral
#     if (dat$Type[i]=="source") tmp[i,] <- 0
#   }}
# fishing.rates <- tmp
# return(fishing.rates)
# }
# 
# getFishingScn(0,2)
# dat$Type



par(mfcol=c(1,1))

ylab="PFA relative"
xlab="Metapopulation Exploitation Rate"
#### DISPERSAL 0.85 FUNCTION OF GLOBAL EXPLOITATION RATE
plot(NULL,xlim=c(0,0.55),ylim=c(0,max(metapops.means.rel, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="")
points(apply(metapops.exp.means[,,2],1,mean),apply(metapops.means.rel[,,2],1,mean),col="grey",pch=16,bg="grey",cex=2,type='b',lwd=4)
points(apply(metapops.exp.means[,,3],1,mean),apply(metapops.means.rel[,,3],1,mean),col=col.sources[1],pch=16,bg=col.sources[1],cex=2,type='b',lwd=4)
points(apply(metapops.exp.means[,,4],1,mean),apply(metapops.means.rel[,,4],1,mean),col=col.sinks[1],pch=16,bg=col.sinks[1],cex=2,type='b',lwd=4)
# for (j in 1:nrow(SCN)){
#   points(mean(metapops.exp.means[j,,2]),mean(metapops.means.rel[j,,2]),col="grey",pch=16,bg="grey",cex=2,type="b")
#   
#   points(mean(metapops.exp.means[j,,3]),mean(metapops.means.rel[j,,3]),col=col.sources[1],pch=16,bg=col.sources[1],cex=2)
# 
#   points(mean(metapops.exp.means[j,,4]),mean(metapops.means.rel[j,,4]),col=col.sinks[1],pch=16,bg=col.sinks[1],cex=2)
# }
# #polynomial ? loess ?
# x_2<-y_2<-NULL
# for (j in 1:nrow(SCN)){
#   x_2 <- c(x_2, metapops.exp.means[j,,2])
#   y_2<- c(y_2, metapops.means.rel[j,,2])
#   
# }
# #lines(lowess(x_2,y_2,f = 1/3), col = col[2], lty = 1)
# intercept<-1
# k<-rep(intercept,length(y_2))
# poly_model <- lm(y_2 ~ -1 + x_2+ I(x_2^2)+offset(k))
# x.predict <- seq(0,0.5,length.out = 20)
# y.predict<-NULL
# for (i in 1:length(x.predict)) {
#   y.predict <- c(y.predict,1 + poly_model$coefficients[1] * x.predict[i] + poly_model$coefficients[2] * x.predict[i]^2)
# }
# #y.predict <- predict(poly_model, newdata = data.frame(x = x.predict))
# lines(x.predict, y.predict, col = col[2],lwd=2,xpd=FALSE)
# 
# x_3<-y_3<-NULL
# for (j in 1:nrow(SCN)){
#   x_3 <- c(x_3, metapops.exp.means[j,,3])
#   y_3<- c(y_3, metapops.means.rel[j,,3])
#   
# }
# #lines(lowess(x_3,y_3,f = 1/3), col = col[3], lty = 1)
# intercept<-1
# k<-rep(intercept,length(y_3))
# poly_model <- lm(y_3 ~ -1 + x_3+ I(x_3^2)+offset(k))
# x.predict <- seq(0,0.5,length.out = 20)
# y.predict<-NULL
# for (i in 1:length(x.predict)) {
#   y.predict <- c(y.predict,1 + poly_model$coefficients[1] * x.predict[i] + poly_model$coefficients[2] * x.predict[i]^2)
# }
# #y.predict <- predict(poly_model, newdata = data.frame(x = x.predict))
# #lines(x.predict, y.predict, col = col[3],lwd=2,xpd=FALSE)
# 
# x_4<-y_4<-NULL
# for (j in 1:nrow(SCN)){
#   x_4 <- c(x_4, metapops.exp.means[j,,4])
#   y_4<- c(y_4, metapops.means.rel[j,,4])
#   
# }#lines(lowess(x_4,y_4,f = 1/3), col = col[4], lty = 1)
# intercept<-1
# k<-rep(intercept,length(y_4))
# poly_model <- lm(y_4 ~ -1 + x_4+ I(x_4^2)+offset(k))
# x.predict <- seq(0,0.5,length.out = 20)
# y.predict<-NULL
# for (i in 1:length(x.predict)) {
#   y.predict <- c(y.predict,1 + poly_model$coefficients[1] * x.predict[i] + poly_model$coefficients[2] * x.predict[i]^2)
# }
# #y.predict <- predict(poly_model, newdata = data.frame(x = x.predict))
# #lines(x.predict, y.predict, col = col[4],lwd=2,xpd=FALSE)

#legend("topright",c("all","source","sink"),col=c("grey",col.sources[1],col.sinks[1]),lty=rep(1,3),bty="n")


#### DISPERSAL 0.7 FUNCTION OF GLOBAL EXPLOITATION RATE
#plot(NULL,xlim=c(0,0.55),ylim=c(0,max(metapops.means.rel, na.rm=TRUE)),xlab="Global Exploitation rates",ylab="Relative Mean Metapop Ab",main="h=0.7")
points(apply(metapops.exp.means[,,5],1,mean),apply(metapops.means.rel[,,5],1,mean),col="lightgrey",pch=16,bg="lightgrey",cex=2,type='b',lwd=4)
points(apply(metapops.exp.means[,,6],1,mean),apply(metapops.means.rel[,,6],1,mean),col=col.sources[2],pch=16,bg=col.sources[2],cex=2,type='b',lwd=4)
points(apply(metapops.exp.means[,,7],1,mean),apply(metapops.means.rel[,,7],1,mean),col=col.sinks[2],pch=16,bg=col.sinks[2],cex=2,type='b',lwd=4)
# for (j in 1:nrow(SCN)){
#   points(mean(metapops.exp.means[j,,5]),mean(metapops.means.rel[j,,5]),col="lightgrey",pch=16,bg="lightgrey",cex=2)
# 
#   points(mean(metapops.exp.means[j,,6]),mean(metapops.means.rel[j,,6]),col=col.sources[2],pch=16,bg=col.sources[2],cex=2)
# 
#   points(mean(metapops.exp.means[j,,7]),mean(metapops.means.rel[j,,7]),col=col.sinks[2],pch=16,bg=col.sinks[2],cex=2)
# }
# x_2<-y_2<-NULL
# for (j in 1:nrow(SCN)){
#   x_2 <- c(x_2, metapops.exp.means[j,,5])
#   y_2<- c(y_2, metapops.means.rel[j,,5])
# }
# lines(lowess(x_2,y_2,f = 1/3), col = col[5], lty = 1)
# 
# x_3<-y_3<-NULL
# for (j in 1:nrow(SCN)){
#   x_3 <- c(x_3, metapops.exp.means[j,,6])
#   y_3<- c(y_3, metapops.means.rel[j,,6])
# }
# lines(lowess(x_3,y_3,f = 1/3), col = col[6], lty = 1)
# 
# x_4<-y_4<-NULL
# for (j in 1:nrow(SCN)){
#   x_4 <- c(x_4, metapops.exp.means[j,,7])
#   y_4<- c(y_4, metapops.means.rel[j,,7])
# }
# lines(lowess(x_4,y_4,f = 1/3), col = col[7], lty = 1)
#legend("topright",c("all","source","sink"),col=c("lightgrey",col.sources[2],col.sinks[2]),lty=rep(1,3),bty="n")

legend("topright"
       ,c("fishing all","fishing source / h=0.85","fishing source / h=0.7","fishing sink / h=0.85","fishing sink / h=0.7")
       ,col=c("grey",col.sources,col.sinks)
       ,text.col=c("grey",col.sources,col.sinks)
       ,pch=rep(16,5)
       #,lty=c(1,1,2,1,2)
       ,bty="n")








#col=c("#030303", "#454545", "#8C8C8C", "#CCCCCC", "#E3E3E3", "#F7F7F7")
#col=c("#000000", "#1874CD", "#009ACD", "#00BFFF", "#228B22", "#66CD00", "#7FFF00")
par(mfcol=c(1,1))

tmp <- metapops.means.rel
tmp2 <- Pglobal.means


ylab="PFA relative"
xlab="Metapopulation Exploitation Rate"
cex=.5#+(frates_vec*2)
xlim=c(0,.5)

Points=FALSE
text=TRUE
lines=FALSE

### FISHING ALL
plot(NULL,xlim=xlim,ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="")#,xaxt='n')
#axis(1,1:6, frates_vec)
#for (i in 1:ncol(SCN)){
if(Points){
  points(apply(tmp2[,,1],1,mean),apply(tmp[,,1],1,mean),col="darkgrey",pch=16,bg="darkgrey",cex=cex,type='b',lwd=4)
  
  points(apply(tmp2[,,2],1,mean),apply(tmp[,,2],1,mean),col="grey",pch=16,bg="grey",cex=cex,type='b',lwd=4)
  points(apply(tmp2[,,3],1,mean),apply(tmp[,,3],1,mean),col=col.sources[1],pch=16,bg=col.sources[1],cex=cex,type='b',lwd=4)
  points(apply(tmp2[,,4],1,mean),apply(tmp[,,4],1,mean),col=col.sinks[1],pch=16,bg=col.sinks[1],cex=cex,type='b',lwd=4)
  
  #points(apply(tmp2[,,5],1,mean),apply(tmp[,,5],1,mean),col="lightgrey",pch=16,bg="lightgrey",cex=cex,type='b',lwd=4)
  points(apply(tmp2[,,6],1,mean),apply(tmp[,,6],1,mean),col=col.sources[2],pch=16,bg=col.sources[2],cex=cex,type='b',lwd=4)
  points(apply(tmp2[,,7],1,mean),apply(tmp[,,7],1,mean),col=col.sinks[2],pch=16,bg=col.sinks[2],cex=cex,type='b',lwd=4)
}

if(text){
  
  text(apply(tmp2[,,1],1,mean),apply(tmp[,,1],1,mean),frates_vec,col="darkgrey",cex=cex)
  
  #text(apply(tmp2[,,2],1,mean),apply(tmp[,,2],1,mean),frates_vec,col="grey",cex=cex)
  text(apply(tmp2[,,3],1,mean),apply(tmp[,,3],1,mean),frates_vec,col=col.sources[1],cex=cex)
  text(apply(tmp2[,,4],1,mean),apply(tmp[,,4],1,mean),frates_vec,col=col.sinks[1],cex=cex)
  
  #text(apply(tmp2[,,5],1,mean),apply(tmp[,,5],1,mean),frates_vec,col="lightgrey",cex=cex)
  text(apply(tmp2[,,6],1,mean),apply(tmp[,,6],1,mean),frates_vec,col=col.sources[2],cex=cex)
  text(apply(tmp2[,,7],1,mean),apply(tmp[,,7],1,mean),frates_vec,col=col.sinks[2],cex=cex)
}
# 
# for (j in 1:13){
#   #points(rep(j,nSIMUL),tmp[j,,1],col=col[1],pch=1,bg=col[1])
#   if(Points) points(mean(tmp2[j,,1]),mean(tmp[j,,1]),col=col[1],pch=16,bg=col[1],cex=cex[j],type='b',lwd=4)
#   if(text) text(mean(tmp2[j,,1]),mean(tmp[j,,1]),frates_vec[j],col=col[1],cex=cex[j])
#segments(j,mean(tmp[j,,1])-sd(tmp[j,,1]),j,mean(tmp[j,,1])-sd(tmp[j,,1]),col=col[i])
#x<-tmp2[,,1];y<-tmp[,,1]
#lines(lowess(x,y,f = 1/3), col = col[1], lty = 1)

#points(rep(j,nSIMUL),tmp[j,,2],col=col[2],pch=1,bg=col[2])
#points(mean(tmp2[j,,2]),mean(tmp[j,,2]),col=col[2],pch=16,bg=col[2],cex=2)
#segments(j,mean(tmp[j,,2])-sd(tmp[j,,2]),j,mean(tmp[j,,2])-sd(tmp[j,,1]),col=col[i])
#x<-tmp2[,,2];y<-tmp[,,2]
#lines(lowess(x,y,f = 1/3), col = col[2], lty = 1)

#points(rep(j,nSIMUL),tmp[j,,5],col=col[5],pch=1,bg=col[5])
#points(mean(tmp2[j,,5]),mean(tmp[j,,5]),col=col[5],pch=16,bg=col[5],cex=2)
#segments(j,mean(tmp[j,,5])-sd(tmp[j,,5]),j,mean(tmp[j,,5])-sd(tmp[j,,5]),col=col[i])
#x<-tmp2[,,5];y<-tmp[,,5]
#lines(lowess(x,y,f = 1/3), col = col[5], lty = 1)
#}
#legend("topleft",c("h=1","h=0.85","h=0.7"),col=c(1,1,1),lty=c(1,1,2),bty="n")

# if(lines){
# y=x=NULL
# for (j in 1:13){
#   y<-c(y,mean(tmp[j,,1]))
#   x<-c(x,mean(tmp2[j,,1]))
# }
# 
# intercept <- 1.0
# k<-rep(intercept,length(y))
# poly_model <- lm(y ~ -1 + x+ I(x^2)+offset(k))
# #x.predict <- seq(0,0.5,length.out = length(y))
# #y.predict <- predict(poly_model, newdata = data.frame(x = x.predict))
# x.predict <- seq(0,0.6,length.out = 20)
# y.predict=NULL
# for (i in 1:length(x.predict)){
#   y.predict[i] <- intercept + poly_model$coefficients[1]*x.predict[i] + poly_model$coefficients[2]*x.predict[i]*x.predict[i]
# }
# lines(x.predict, y.predict, col = col[1],lwd=2,xpd=FALSE)
# }



### FISHING SOURCES
# col.sources=c("dodgerblue2", "deepskyblue")
# #plot(NULL,xlim=xlim,ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Fishing source")#,xaxt='n')
# #axis(1,1:6, frates_vec)
# #for (i in 1:ncol(SCN)){
# for (j in 1:13){
#   #points(rep(j,nSIMUL),tmp[j,,1],col=col[1],pch=1,bg=col[1])
#   #points(mean(tmp2[j,,1]),mean(tmp[j,,1]),col=col[1],pch=16,bg=col[1],cex=2)
#   #segments(j,mean(tmp[j,,1])-sd(tmp[j,,1]),j,mean(tmp[j,,1])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:6,6,nSIMUL);y<-tmp[,,1]
#   #lines(lowess(x,y,f = 1/3), col = col[1], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,3],col=col[3],pch=1,bg=col[3])
#   if(Points)  points(mean(tmp2[j,,3]),mean(tmp[j,,3]),col=col.sources[1],pch=16,bg=col.sources[1],cex=cex[j])
#   if(text) text(mean(tmp2[j,,3]),mean(tmp[j,,3]),frates_vec[j],col=col.sources[1],cex=cex[j])
#   #segments(j,mean(tmp[j,,3])-sd(tmp[j,,3]),j,mean(tmp[j,,3])-sd(tmp[j,,3]),col=col[i])
#   #x<-tmp2[,,3];y<-tmp[,,3]
#   #lines(lowess(x,y,f = 1/3), col = col[2], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,5],col=col[5],pch=1,bg=col[5])
#   points(mean(tmp2[j,,6]),mean(tmp[j,,6]),col=col.sources[2],pch=16,bg=col.sources[2],cex=cex[j])
#   if(text) text(mean(tmp2[j,,6]),mean(tmp[j,,6]),frates_vec[j],col=col.sources[2],cex=cex[j])
#   #segments(j,mean(tmp[j,,5])-sd(tmp[j,,5]),j,mean(tmp[j,,5])-sd(tmp[j,,5]),col=col[i])
#   #x<-tmp2[,,6];y<-tmp[,,6]
#   #lines(lowess(x,y,f = 1/3), col = col[5], lty = 1)
# }
# #legend("topleft",c("h=1","h=0.85","h=0.7"),col=col[c(1,3,6)],lty=rep(1,3),bty="n")
# 
# if(lines){
# y=x=NULL
# for (j in 1:13){
#   y<-c(y,mean(tmp[j,,3]))
#   x<-c(x,mean(tmp2[j,,3]))
# }
# 
# intercept <- 1.0
# k<-rep(intercept,length(y))
# poly_model <- lm(y ~ -1 + x+ I(x^2)+offset(k))
# #x.predict <- seq(0,0.5,length.out = length(y))
# #y.predict <- predict(poly_model, newdata = data.frame(x = x.predict))
# x.predict <- seq(0,0.6,length.out = 20)
# y.predict=NULL
# for (i in 1:length(x.predict)){
#   y.predict[i] <- intercept + poly_model$coefficients[1]*x.predict[i] + poly_model$coefficients[2]*x.predict[i]*x.predict[i]
# }
# lines(x.predict, y.predict, col = col.sources[1],lwd=2,xpd=FALSE)
# 
# y=x=NULL
# for (j in 1:13){
#   y<-c(y,mean(tmp[j,,6]))
#   x<-c(x,mean(tmp2[j,,6]))
# }
# 
# intercept <- 1.0
# k<-rep(intercept,length(y))
# poly_model <- lm(y ~ -1 + x+ I(x^2)+offset(k))
# #x.predict <- seq(0,0.5,length.out = length(y))
# #y.predict <- predict(poly_model, newdata = data.frame(x = x.predict))
# x.predict <- seq(0,0.6,length.out = 20)
# y.predict=NULL
# for (i in 1:length(x.predict)){
#   y.predict[i] <- intercept + poly_model$coefficients[1]*x.predict[i] + poly_model$coefficients[2]*x.predict[i]*x.predict[i]
# }
# lines(x.predict, y.predict, col = col.sources[2],lwd=2,xpd=FALSE,lty=2)
# }
# 
# 
# 
# 
# ### FISHING SINKS
# col.sinks=c("brown1", "lightcoral", "white")
# #plot(NULL,xlim=xlim,ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Fishing sinks")#,xaxt='n')
# #axis(1,1:6, frates_vec)
# #for (i in 1:ncol(SCN)){
# for (j in 1:13){
#   #points(rep(j,nSIMUL),tmp[j,,1],col=col[1],pch=1,bg=col[1])
#   #points(j,mean(tmp[j,,1]),col=col[1],pch=16,bg=col[1],cex=2)
#   #segments(j,mean(tmp[j,,1])-sd(tmp[j,,1]),j,mean(tmp[j,,1])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:6,6,nSIMUL);y<-tmp[,,1]
#   #lines(lowess(x,y,f = 1/3), col = col[1], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,4],col=col[2],pch=1,bg=col[2])
#   if(Points) points(mean(tmp2[j,,4]),mean(tmp[j,,4]),col=col.sinks[1],pch=16,bg=col.sinks[1],cex=cex[j])
#   if(text) text(mean(tmp2[j,,4]),mean(tmp[j,,4]),frates_vec[j],col=col.sinks[1],cex=cex[j])
#   #segments(j,mean(tmp[j,,3])-sd(tmp[j,,3]),j,mean(tmp[j,,3])-sd(tmp[j,,3]),col=col[i])
#   #x<-tmp[,,4];y<-tmp[,,4]
#   #lines(lowess(x,y,f = 1/3), col = col[2], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,5],col=col[5],pch=1,bg=col[5])
#   if(Points) points(mean(tmp2[j,,7]),mean(tmp[j,,7]),col=col.sinks[2],pch=16,bg=col.sinks[2],cex=cex[j])
#   if(text) text(mean(tmp2[j,,7]),mean(tmp[j,,7]),frates_vec[j],col=col.sinks[2],cex=cex[j])
#   #segments(j,mean(tmp[j,,5])-sd(tmp[j,,5]),j,mean(tmp[j,,5])-sd(tmp[j,,5]),col=col[i])
#   #x<-tmp2[,,7];y<-tmp[,,7]
#   #lines(lowess(x,y,f = 1/3), col = col[5], lty = 1)
#   
#   
# }
# #legend("topleft",c("h=1","h=0.85","h=0.7"),col=col[c(1,3,6)],lty=rep(1,3),bty="n")
# if(lines){
# y=x=NULL
# for (j in 1:13){
#   y<-c(y,mean(tmp[j,,4]))
#   x<-c(x,mean(tmp2[j,,4]))
# }
# 
# intercept <- 1.0
# k<-rep(intercept,length(y))
# poly_model <- lm(y ~ -1 + x+ I(x^2)+offset(k))
# #x.predict <- seq(0,0.5,length.out = length(y))
# #y.predict <- predict(poly_model, newdata = data.frame(x = x.predict))
# x.predict <- seq(0,0.6,length.out = 20)
# y.predict=NULL
# for (i in 1:length(x.predict)){
#   y.predict[i] <- intercept + poly_model$coefficients[1]*x.predict[i] + poly_model$coefficients[2]*x.predict[i]*x.predict[i]
# }
# lines(x.predict, y.predict, col = col.sinks[1],lwd=2,xpd=FALSE)
# 
# 
# 
# 
# y=x=NULL
# for (j in 1:13){
#   y<-c(y,mean(tmp[j,,7]))
#   x<-c(x,mean(tmp2[j,,7]))
# }
# 
# intercept <- 1.0
# k<-rep(intercept,length(y))
# poly_model <- lm(y ~ -1 + x+ I(x^2)+offset(k))
# #x.predict <- seq(0,0.5,length.out = length(y))
# #y.predict <- predict(poly_model, newdata = data.frame(x = x.predict))
# x.predict <- seq(0,0.6,length.out = 20)
# y.predict=NULL
# for (i in 1:length(x.predict)){
#   y.predict[i] <- intercept + poly_model$coefficients[1]*x.predict[i] + poly_model$coefficients[2]*x.predict[i]*x.predict[i]
# }
# lines(x.predict, y.predict, col = col.sinks[2],lwd=2,xpd=FALSE,lty=2)
# }

legend("topright"
       ,c("fishing all","fishing source / h=0.85","fishing source / h=0.7","fishing sink / h=0.85","fishing sink / h=0.7")
       ,col=c("grey",col.sources,col.sinks)
       ,text.col=c("grey",col.sources,col.sinks)
       ,pch=rep(16,5)
       #,lty=c(1,1,2,1,2)
       ,bty="n")









#### EXTINCTION ####
#---------------------------#


par(mfrow=c(3,3))
#col=c("#000000", "#1874CD", "#009ACD", "#00BFFF", "#228B22", "#66CD00", "#7FFF00")
tmp <- Extinction

xlab="Time (years)"

ylab="Prop. population extinct"

lwd=1

ylim=c(0,max(unlist(tmp), na.rm=TRUE))

plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main="Philopatry=1 / Fishing all")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[1]])){
  points(1:50,tmp[[1]][[j]],col=col[j],type='l',lwd=lwd)
}
#legend("topleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("topleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)

plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main="Philopatry=0.85 / Fishing all")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[2]])){
  points(1:50,tmp[[2]][[j]],col=col[j],type='l',lwd=lwd)
}
#legend("topleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("topleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)


plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main="Philopatry=0.7 / Fishing all")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[5]])){
  points(1:50,tmp[[5]][[j]],col=col[j],type='l',lwd=lwd)
}
#legend("topleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("topleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n")




#plot.new()
#legend("topleft",paste0(frates_vec),col=col,lty=rep(1,length(frates_vec)),bty="n",title="Local Exploitation rates",cex=.6)
plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = "Local Exploitation rates")
#legend("bottomleft",paste0(frates_vec),col=col,lty=rep(1,length(frates_vec)),bty="n",title="Local Exploitation rates",cex=.6)
legend_image <- as.raster(matrix(rev(col), ncol=1))
text(x=0.5, y = 1.1, "Local Exploitation rates", cex=.6)
text(x=0.45, y = seq(0,1,l=length(frates_vec)), labels = paste0(frates_vec),adj=.2,cex=.5)
rasterImage(legend_image, 0.3, 0, .4,1)

plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main="Philopatry=0.85 / Fishing sources")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[3]])){
  points(1:50,tmp[[3]][[j]],col=col[j],type='l',lwd=lwd)
}
#legend("topleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("topleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)

plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main="Philopatry=0.7 / Fishing sources")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[6]])){
  points(1:50,tmp[[6]][[j]],col=col[j],type='l',lwd=lwd)
}
#legend("topleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("topleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n")


plot.new()

plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main="Philopatry=0.85 / Fishing sinks")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[4]])){
  points(1:50,tmp[[4]][[j]],col=col[j],type='l',lwd=lwd)
}
#legend("topleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("topleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)

plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main="Philopatry=0.7 / Fishing sinks")
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
for (j in 1:length(tmp[[7]])){
  points(1:50,tmp[[7]][[j]],col=col[j],type='l',lwd=lwd)
}
# #legend("topleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
# legend("topleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n")












tmp <- extinction.means.rel

xlab="Local Exploitation rates"
#ylab="Metapopulation Exploitation Rate"

#col=c("#030303", "#454545", "#8C8C8C", "#CCCCCC", "#E3E3E3", "#F7F7F7")
#col=c("#000000", "#1874CD", "#009ACD", "#00BFFF", "#228B22", "#66CD00", "#7FFF00")
col.sources=c("dodgerblue2", "deepskyblue")
col.sinks=c("brown1", "lightcoral", "white")
par(mfcol=c(3,2))

lwd=3

ylim=c(0,0.2)#max(tmp, na.rm=TRUE))

### FISHING ALL
plot(NULL,xlim=c(1,length(frates_vec)),ylim=ylim,xlab=xlab,ylab=ylab,main="Fishing all",xaxt='n')
axis(1,1:length(frates_vec), frates_vec)
#for (i in 1:ncol(SCN)){
points(1:length(frates_vec),apply(tmp[,,1],1,mean),col=col[1],pch=16,bg=col[1],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,2],1,mean),col="grey",pch=16,bg="grey",cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,5],1,mean),col="lightgrey",pch=16,bg="lightgrey",cex=2,type='l',lwd=lwd)
# for (j in 1:13){
#   #points(rep(j,nSIMUL),tmp[j,,1],col=col[1],pch=1,bg=col[1])
#   points(j,mean(tmp[j,,1]),col=col[1],pch=16,bg=col[1],cex=2,type='l')
#   #segments(j,mean(tmp[j,,1])-sd(tmp[j,,1]),j,mean(tmp[j,,1])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,1]
#   #lines(lowess(x,y,f = 1/3), col = col[1], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,2],col=col[2],pch=1,bg=col[2])
#   points(j,mean(tmp[j,,2]),col=col[2],pch=16,bg=col[2],cex=2)
#   #segments(j,mean(tmp[j,,2])-sd(tmp[j,,2]),j,mean(tmp[j,,2])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,2]
#   #lines(lowess(x,y,f = 1/3), col = col[2], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,5],col=col[5],pch=1,bg=col[5])
#   points(j,mean(tmp[j,,5]),col=col[5],pch=16,bg=col[5],cex=2)
#   #segments(j,mean(tmp[j,,5])-sd(tmp[j,,5]),j,mean(tmp[j,,5])-sd(tmp[j,,5]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,5]
#   #lines(lowess(x,y,f = 1/3), col = col[5], lty = 1)
# }
legend("topleft",c("h=1","h=0.85","h=0.7"),col=c(1,"grey","lightgrey"),lty=rep(1,3),lwd=rep(2,3),bty="n")


### FISHING SOURCES
plot(NULL,xlim=c(1,length(frates_vec)),ylim=ylim,xlab=xlab,ylab=ylab,main="Fishing source",xaxt='n')
axis(1,1:length(frates_vec), frates_vec)
points(1:length(frates_vec),apply(tmp[,,1],1,mean),col=col[1],pch=16,bg=col[1],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,3],1,mean),col=col.sources[1],pch=16,bg=col.sources[1],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,6],1,mean),col=col.sources[2],pch=16,bg=col.sources[2],cex=2,type='l',lwd=lwd)
#for (i in 1:ncol(SCN)){
# for (j in 1:13){
#   #points(rep(j,nSIMUL),tmp[j,,1],col=col[1],pch=16,bg=col[1])
#   points(j,mean(tmp[j,,1]),col=col[1],pch=16,bg=col[1],cex=2)
#   #segments(j,mean(tmp[j,,1])-sd(tmp[j,,1]),j,mean(tmp[j,,1])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,1]
#   #lines(lowess(x,y,f = 1/3), col = col[1], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,3],col=col[3],pch=16,bg=col[3])
#   points(j,mean(tmp[j,,3]),col=col[2],pch=16,bg=col[2],cex=2)
#   #segments(j,mean(tmp[j,,3])-sd(tmp[j,,3]),j,mean(tmp[j,,3])-sd(tmp[j,,3]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,3]
#   #lines(lowess(x,y,f = 1/3), col = col[2], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,5],col=col[5],pch=1,bg=col[5])
#   points(j,mean(tmp[j,,7]),col=col[5],pch=16,bg=col[5],cex=2)
#   #segments(j,mean(tmp[j,,5])-sd(tmp[j,,5]),j,mean(tmp[j,,5])-sd(tmp[j,,5]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,13]
#   #lines(lowess(x,y,f = 1/3), col = col[5], lty = 1)
# }
legend("topleft",c("h=1","h=0.85","h=0.7"),col=c(1,col.sources),lty=rep(1,3),lwd=rep(2,3),bty="n")


### FISHING SINKS
plot(NULL,xlim=c(1,length(frates_vec)),ylim=ylim,xlab=xlab,ylab=ylab,main="Fishing sinks",xaxt='n')
axis(1,1:length(frates_vec), frates_vec)
points(1:length(frates_vec),apply(tmp[,,1],1,mean),col=col[1],pch=16,bg=col[1],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,4],1,mean),col=col.sinks[1],pch=16,bg=col.sinks[1],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,7],1,mean),col=col.sinks[2],pch=16,bg=col.sinks[2],cex=2,type='l',lwd=lwd)
#for (i in 1:ncol(SCN)){
# for (j in 1:13){
#   #points(rep(j,nSIMUL),tmp[j,,1],col=col[1],pch=16,bg=col[1])
#   points(j,mean(tmp[j,,1]),col=col[1],pch=16,bg=col[1],cex=2)
#   #segments(j,mean(tmp[j,,1])-sd(tmp[j,,1]),j,mean(tmp[j,,1])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,1]
#   #lines(lowess(x,y,f = 1/3), col = col[1], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,4],col=col[2],pch=1,bg=col[2])
#   points(j,mean(tmp[j,,4]),col=col[2],pch=16,bg=col[2],cex=2)
#   #segments(j,mean(tmp[j,,3])-sd(tmp[j,,3]),j,mean(tmp[j,,3])-sd(tmp[j,,3]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,4]
#   #lines(lowess(x,y,f = 1/3), col = col[2], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,5],col=col[5],pch=1,bg=col[5])
#   points(j,mean(tmp[j,,7]),col=col[5],pch=16,bg=col[5],cex=2)
#   #segments(j,mean(tmp[j,,5])-sd(tmp[j,,5]),j,mean(tmp[j,,5])-sd(tmp[j,,5]),col=col[i])
#   #x<-matrix(1:13,13,nSIMUL);y<-tmp[,,7]
#   #lines(lowess(x,y,f = 1/3), col = col[5], lty = 1)
# }
legend("topleft",c("h=1","h=0.85","h=0.7"),col=c(1,col.sinks),lty=rep(1,3),lwd=rep(2,3),bty="n")



# plot(NULL,xlim=c(1,13),ylim=c(0,max(MeansRel)),xlab=xlab,ylab=ylab)
# for (i in 1:ncol(SCN)){
#   points(1:13,MeansRel[,i],col=col[i], type='l')
#   points(1:6,MeansRel[,i],col=col[i], pch=paste0(i))
#   #segments(1:6,MeansRel[,i]-Sds[,i],1:6,Means[,i]+Sds[,i],col=col[i])
# }









# Npops <- array(,dim=c((nYears+nInit),nSIMUL,length(nReturns)))
# for (scn in 1:length(nReturns)){
#   #for (div in 1:length(DIV)){
#   for (simul in 1:nSIMUL){
#     Npops[,simul,scn] <- nReturns[[scn]][[1]]$NRet[[simul]][,length(frates_vec)] #nb returns metapop per year, simu and scenario
#   }
#   #}
# }




#par(mfrow=c(2,1))

### DISPERSAL 1
plot.new()


### DISPERSAL 0.85 FUNCTION OF LOCAL EXPLOITATION RATE
plot(NULL,xlim=c(1,length(frates_vec)),ylim=ylim,xlab=xlab,ylab=ylab,main="Philopatry = 0.85",xaxt='n')
axis(1,1:length(frates_vec), frates_vec)
points(1:length(frates_vec),apply(tmp[,,2],1,mean),col="grey",pch=16,bg="grey",cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,3],1,mean),col=col.sources[1],pch=16,bg=col.sources[1],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,4],1,mean),col=col.sinks[1],pch=16,bg=col.sinks[1],cex=2,type='l',lwd=lwd)
#for (i in 1:ncol(SCN)){
# for (j in 1:13){
#   #points(rep(j,nSIMUL),tmp[j,,1],col=col[1],pch=1,bg=col[1])
#   points(j,mean(tmp[j,,2]),col=col[2],pch=16,bg=col[2],cex=2)
#   #segments(j,mean(tmp[j,,1])-sd(tmp[j,,1]),j,mean(tmp[j,,1])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:6,6,nSIMUL);y<-tmp[,,2]
#   #lines(lowess(x,y,f = 1/3), col = col[2], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,2],col=col[2],pch=1,bg=col[2])
#   points(j,mean(tmp[j,,3]),col=col[3],pch=16,bg=col[3],cex=2)
#   #segments(j,mean(tmp[j,,2])-sd(tmp[j,,2]),j,mean(tmp[j,,2])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:6,6,nSIMUL);y<-tmp[,,3]
#   #lines(lowess(x,y,f = 1/3), col = col[3], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,5],col=col[5],pch=1,bg=col[5])
#   points(j,mean(tmp[j,,4]),col=col[4],pch=16,bg=col[4],cex=2)
#   #segments(j,mean(tmp[j,,5])-sd(tmp[j,,5]),j,mean(tmp[j,,5])-sd(tmp[j,,5]),col=col[i])
#   #x<-matrix(1:6,6,nSIMUL);y<-tmp[,,4]
#   #lines(lowess(x,y,f = 1/3), col = col[4], lty = 1)
# }
legend("topleft",c("all","source","sink"),col=c("grey",col.sources[1],col.sinks[1]),lty=rep(1,3),lwd=rep(2,3),bty="n")


### DISPERSAL 0.7 FUNCTION OF LOCAL EXPLOITATION RATE
plot(NULL,xlim=c(1,length(frates_vec)),ylim=ylim,xlab=xlab,ylab=ylab,main=" Philopatry = 0.7",xaxt='n')
axis(1,1:length(frates_vec), frates_vec)
points(1:length(frates_vec),apply(tmp[,,5],1,mean),col="lightgrey",pch=16,bg="grey",cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,6],1,mean),col=col.sources[2],pch=16,bg=col.sources[2],cex=2,type='l',lwd=lwd)
points(1:length(frates_vec),apply(tmp[,,7],1,mean),col=col.sinks[2],pch=16,bg=col.sinks[2],cex=2,type='l',lwd=lwd)
#for (i in 1:ncol(SCN)){
# for (j in 1:13){
#   #points(rep(j,nSIMUL),tmp[j,,1],col=col[1],pch=1,bg=col[1])
#   points(j,mean(tmp[j,,5]),col=col[5],pch=16,bg=col[5],cex=2)
#   #segments(j,mean(tmp[j,,1])-sd(tmp[j,,1]),j,mean(tmp[j,,1])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:6,6,nSIMUL);y<-tmp[,,5]
#   #lines(lowess(x,y,f = 1/3), col = col[5], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,2],col=col[2],pch=1,bg=col[2])
#   points(j,mean(tmp[j,,6]),col=col[6],pch=16,bg=col[6],cex=2)
#   #segments(j,mean(tmp[j,,2])-sd(tmp[j,,2]),j,mean(tmp[j,,2])-sd(tmp[j,,1]),col=col[i])
#   #x<-matrix(1:6,6,nSIMUL);y<-tmp[,,6]
#   #lines(lowess(x,y,f = 1/3), col = col[6], lty = 1)
#   
#   #points(rep(j,nSIMUL),tmp[j,,5],col=col[5],pch=1,bg=col[5])
#   points(j,mean(tmp[j,,7]),col=col[7],pch=16,bg=col[7],cex=2)
#   #segments(j,mean(tmp[j,,5])-sd(tmp[j,,5]),j,mean(tmp[j,,5])-sd(tmp[j,,5]),col=col[i])
#   #x<-matrix(1:6,6,nSIMUL);y<-tmp[,,7]
#   #lines(lowess(x,y,f = 1/3), col = col[7], lty = 1)
# }
legend("topleft",c("all","source","sink"),col=c("lightgrey",col.sources[2],col.sinks[2]),lty=rep(1,3),lwd=rep(2,3),bty="n")














#col=c("#030303", "#454545", "#8C8C8C", "#CCCCCC", "#E3E3E3", "#F7F7F7")
#col=c("#000000", "#1874CD", "#009ACD", "#00BFFF", "#228B22", "#66CD00", "#7FFF00")
par(mfcol=c(1,1))

tmp <- extinction.means.rel
tmp2 <- Pglobal.means


#ylab="Extinction risk"
xlab="Metapopulation Exploitation Rate"
cex=.5#+(frates_vec*2)
xlim=c(0,.5)
ylim=c(0,.2)#max(tmp, na.rm=TRUE))

Points=TRUE
text=TRUE
lines=TRUE

### FISHING ALL
plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main="")#,xaxt='n')
#axis(1,1:6, frates_vec)
#for (i in 1:ncol(SCN)){
if(Points){
  points(apply(tmp2[,,1],1,mean),apply(tmp[,,1],1,mean),col="darkgrey",pch=16,bg="darkgrey",cex=cex,type='b',lwd=4)
  
  points(apply(tmp2[,,2],1,mean),apply(tmp[,,2],1,mean),col="grey",pch=16,bg="grey",cex=cex,type='b',lwd=4)
  points(apply(tmp2[,,3],1,mean),apply(tmp[,,3],1,mean),col=col.sources[1],pch=16,bg=col.sources[1],cex=cex,type='b',lwd=4)
  points(apply(tmp2[,,4],1,mean),apply(tmp[,,4],1,mean),col=col.sinks[1],pch=16,bg=col.sinks[1],cex=cex,type='b',lwd=4)
  
  points(apply(tmp2[,,5],1,mean),apply(tmp[,,5],1,mean),col="lightgrey",pch=16,bg="lightgrey",cex=cex,type='b',lwd=4)
  points(apply(tmp2[,,6],1,mean),apply(tmp[,,6],1,mean),col=col.sources[2],pch=16,bg=col.sources[2],cex=cex,type='b',lwd=4)
  points(apply(tmp2[,,7],1,mean),apply(tmp[,,7],1,mean),col=col.sinks[2],pch=16,bg=col.sinks[2],cex=cex,type='b',lwd=4)
}

if(text){
  
  text(apply(tmp2[,,1],1,mean),apply(tmp[,,1],1,mean),frates_vec,col="darkgrey",cex=cex)
  
  text(apply(tmp2[,,2],1,mean),apply(tmp[,,2],1,mean),frates_vec,col="grey",cex=cex)
  text(apply(tmp2[,,3],1,mean),apply(tmp[,,3],1,mean),frates_vec,col=col.sources[1],cex=cex)
  text(apply(tmp2[,,4],1,mean),apply(tmp[,,4],1,mean),frates_vec,col=col.sinks[1],cex=cex)
  
  text(apply(tmp2[,,5],1,mean),apply(tmp[,,5],1,mean),frates_vec,col="lightgrey",cex=cex)
  text(apply(tmp2[,,6],1,mean),apply(tmp[,,6],1,mean),frates_vec,col=col.sources[2],cex=cex)
  text(apply(tmp2[,,7],1,mean),apply(tmp[,,7],1,mean),frates_vec,col=col.sinks[2],cex=cex)
}


legend("topright"
       ,c("fishing all","fishing source / h=0.85","fishing source / h=0.7","fishing sink / h=0.85","fishing sink / h=0.7")
       ,col=c("grey",col.sources,col.sinks)
       ,text.col=c("grey",col.sources,col.sinks)
       ,pch=rep(16,5)
       #,lty=c(1,1,2,1,2)
       ,bty="n")













#### POPULATION DYNAMICS ####
#-------------------------------#


par(mfrow=c(3,3))

tmp <- Population

xlab="Time (years)"

#ylab="Metapopulation Abundance"
#ylab="PFA"
#ylab="Metapopulation Relative Abundance"
ylab="Relative PFA" # Pre-Fisheries Abundance

# Create a list of gradients for each colour 2 to 10 over five steps from 
library(viridis)
#col <- viridis(13)
col <- (inferno(length(frates_vec)))
#col=c("#000000", "#1874CD", "#009ACD", "#00BFFF", "#228B22", "#66CD00", "#7FFF00")

pop.type<-as.numeric(factor(dat$Type))
col.type <- c("black",col.sinks[1],col.sources[1])

lwd=1
ylim=c(0,3)#c(0,max(unlist(tmp), na.rm=TRUE))

#pop=1
for (pop in 1:npop){
  
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab)#,main="Philopatry=1 / Fishing all")
  title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
  j=0
  for (iEXPE in SCN[,1]){
    j=j+1
    if(is.na(iEXPE)) next;
    points(1:50,apply(tmp[[1]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
  }
  #legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
  #legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)
  
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab)#,main="Philopatry=0.85 / Fishing all")
  title("Philopatry=0.85 / Fishing all", col.main = col.type[pop.type[pop]])
  rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
  j=0
  for (iEXPE in SCN[,2]){
    j=j+1
    if(is.na(iEXPE)) next;
    points(1:50,apply(tmp[[2]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
  }
  #legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
  #legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)
  
  
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim, xlab=xlab,ylab=ylab)#,main="Philopatry=0.7 / Fishing all")
  title("Philopatry=0.7 / Fishing all", col.main = col.type[pop.type[pop]])
  rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
  j=0
  for (iEXPE in SCN[,5]){
    j=j+1
    if(is.na(iEXPE)) next;
    points(1:50,apply(tmp[[5]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
  }
  #legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
  #legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n")
  
  
  
  
  #plot.new()
  plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = "Local Exploitation rates")
  #legend("bottomleft",paste0(frates_vec),col=col,lty=rep(1,length(frates_vec)),bty="n",title="Local Exploitation rates",cex=.6)
  legend_image <- as.raster(matrix(rev(col), ncol=1))
  text(x=0.5, y = 1.1, "Local Exploitation rates", cex=.6)
  text(x=0.45, y = seq(0,1,l=length(frates_vec)), labels = paste0(frates_vec),adj=.2,cex=.5)
  rasterImage(legend_image, 0.3, 0, .4,1)
  
  
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim, xlab=xlab,ylab=ylab)#,main="Philopatry=0.85 / Fishing sources")
  title("Philopatry=0.85 / Fishing sources", col.main = col.type[pop.type[pop]])
  rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
  j=0
  for (iEXPE in SCN[,3]){
    j=j+1
    if(is.na(iEXPE)) next;
    points(1:50,apply(tmp[[3]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
  }
  abline(h=1-prop.source,lty=2)
  #legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
  #legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)
  
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim, xlab=xlab,ylab=ylab)#,main="Philopatry=0.7 / Fishing sources")
  title("Philopatry=0.7 / Fishing sources", col.main = col.type[pop.type[pop]])
  rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
  j=0
  for (iEXPE in SCN[,6]){
    j=j+1
    if(is.na(iEXPE)) next;
    points(1:50,apply(tmp[[6]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
  }
  abline(h=1-prop.source,lty=2)
  #legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
  #legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n")
  
  
  #plot.new()
  plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,2),xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  text(length(frates_vec)/2,1,paste0(dat$Population[pop]," / ", dat$Type[pop]),col=col.type[pop.type[pop]],cex=3)
  
  # legend_image <- as.raster(matrix(rev(col), ncol=1))
  # plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
  # text(x=0.1, y = seq(0,1,l=length(frates_vec)), labels = paste0(frates_vec),adj=.2,cex=.7)
  # rasterImage(legend_image, 0, 0, .1,1)
  
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim, xlab=xlab,ylab=ylab)#,main="Philopatry=0.85 / Fishing sinks")
  title("Philopatry=0.85 / Fishing sinks", col.main = col.type[pop.type[pop]])
  rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
  j=0
  for (iEXPE in SCN[,4]){
    j=j+1
    if(is.na(iEXPE)) next;
    points(1:50,apply(tmp[[4]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
  }
  abline(h=1-prop.sink,lty=2)
  #legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
  #legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)
  
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim, xlab=xlab,ylab=ylab)#,main="Philopatry=0.7 / Fishing sinks")
  title("Philopatry=0.7 / Fishing sinks", col.main = col.type[pop.type[pop]])
  rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
  j=0
  for (iEXPE in SCN[,7]){
    j=j+1
    if(is.na(iEXPE)) next;
    points(1:50,apply(tmp[[7]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
  }
  abline(h=1-prop.sink,lty=2)
  # #legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
  
  #legend("bottomleft",paste0(frates_vec),col=col[1:length(frates_vec)],lty=rep(1,length(frates_vec)),bty="n")
  
} # end loop pop







#### IMMIGRANTS DYNAMICS ####
#-------------------------------#


par(mfrow=c(3,3))

tmp <- Prop.Immigrants

xlab="Time (years)"

#ylab="Metapopulation Abundance"
#ylab="PFA"
#ylab="Metapopulation Relative Abundance"
ylab="Ratio Immigrants" # Pre-Fisheries Abundance

# Create a list of gradients for each colour 2 to 10 over five steps from 
library(viridis)
#col <- viridis(13)
col <- (inferno(length(frates_vec)))
#col=c("#000000", "#1874CD", "#009ACD", "#00BFFF", "#228B22", "#66CD00", "#7FFF00")

pop.type<-as.numeric(factor(dat$Type))
col.type <- c("black",col.sinks[1],col.sources[1])


lwd=1
ylim=c(0,1)#c(0,max(unlist(tmp), na.rm=TRUE))

#pop=1
for (pop in 1:npop){
  
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab)#,main="Philopatry=1 / Fishing all")
  title("Philopatry=1 / Fishing all", col.main =col.type[pop.type[pop]])
  rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
  j=0
  for (iEXPE in SCN[,1]){
    j=j+1
    if(is.na(iEXPE)) next;
    points(1:50,apply(tmp[[1]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
  }
  #legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
  #legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)
  
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab)#,main="Philopatry=0.85 / Fishing all")
  title("Philopatry=0.85 / Fishing all", col.main =col.type[pop.type[pop]])
  rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
  j=0
  for (iEXPE in SCN[,2]){
    j=j+1
    if(is.na(iEXPE)) next;
    points(1:50,apply(tmp[[2]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
  }
  #legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
  #legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)
  
  
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim, xlab=xlab,ylab=ylab)#,main="Philopatry=0.7 / Fishing all")
  title("Philopatry=0.7 / Fishing all", col.main =col.type[pop.type[pop]])
  rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
  j=0
  for (iEXPE in SCN[,5]){
    j=j+1
    if(is.na(iEXPE)) next;
    points(1:50,apply(tmp[[5]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
  }
  #legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
  #legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n")
  
  
  
  
  #plot.new()
  plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = "Local Exploitation rates")
  #legend("bottomleft",paste0(frates_vec),col=col,lty=rep(1,length(frates_vec)),bty="n",title="Local Exploitation rates",cex=.6)
  legend_image <- as.raster(matrix(rev(col), ncol=1))
  text(x=0.5, y = 1.1, "Local Exploitation rates", cex=.6)
  text(x=0.45, y = seq(0,1,l=length(frates_vec)), labels = paste0(frates_vec),adj=.2,cex=.5)
  rasterImage(legend_image, 0.3, 0, .4,1)
  
  
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim, xlab=xlab,ylab=ylab)#,main="Philopatry=0.85 / Fishing sources")
  title("Philopatry=0.85 / Fishing sources", col.main =col.type[pop.type[pop]])
  rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
  j=0
  for (iEXPE in SCN[,3]){
    j=j+1
    if(is.na(iEXPE)) next;
    points(1:50,apply(tmp[[3]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
  }
  abline(h=1-prop.source,lty=2)
  #legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
  #legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)
  
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim, xlab=xlab,ylab=ylab)#,main="Philopatry=0.7 / Fishing sources")
  title("Philopatry=0.7 / Fishing sources", col.main =col.type[pop.type[pop]])
  rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
  j=0
  for (iEXPE in SCN[,6]){
    j=j+1
    if(is.na(iEXPE)) next;
    points(1:50,apply(tmp[[6]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
  }
  abline(h=1-prop.source,lty=2)
  #legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
  #legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n")
  
  
  #plot.new()
  plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,2),xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  text(length(frates_vec)/2,1,paste0(dat$Population[pop]," / ", dat$Type[pop]),col=col.type[pop.type[pop]],cex=3)
  
  # legend_image <- as.raster(matrix(rev(col), ncol=1))
  # plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
  # text(x=0.1, y = seq(0,1,l=length(frates_vec)), labels = paste0(frates_vec),adj=.2,cex=.7)
  # rasterImage(legend_image, 0, 0, .1,1)
  
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim, xlab=xlab,ylab=ylab)#,main="Philopatry=0.85 / Fishing sinks")
  title("Philopatry=0.85 / Fishing sinks", col.main =col.type[pop.type[pop]])
  rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
  j=0
  for (iEXPE in SCN[,4]){
    j=j+1
    if(is.na(iEXPE)) next;
    points(1:50,apply(tmp[[4]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
  }
  abline(h=1-prop.sink,lty=2)
  #legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
  #legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)
  
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim, xlab=xlab,ylab=ylab)#,main="Philopatry=0.7 / Fishing sinks")
  title("Philopatry=0.7 / Fishing sinks", col.main =col.type[pop.type[pop]])
  rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
  j=0
  for (iEXPE in SCN[,7]){
    j=j+1
    if(is.na(iEXPE)) next;
    points(1:50,apply(tmp[[7]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
  }
  abline(h=1-prop.sink,lty=2)
  # #legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
  
  #legend("bottomleft",paste0(frates_vec),col=col[1:length(frates_vec)],lty=rep(1,length(frates_vec)),bty="n")
  
} # end loop pop
















#### FISHING STRATEGIES / LOCAL EXPLOITATION RATES ####
#-------------------------------#



#tmp <- metapops.means
tmp <- pops.means


#col=c("#030303", "#454545", "#8C8C8C", "#CCCCCC", "#E3E3E3", "#F7F7F7")
#col=c("#000000", "#1874CD", "#009ACD", "#00BFFF", "#228B22", "#66CD00", "#7FFF00")
col.sources=c("dodgerblue2", "deepskyblue")
col.sinks=c("brown1", "lightcoral", "white")
par(mfcol=c(3,2))

lwd=3

xlab="Local Exploitation rates"

#ylab="Metapopulation Abundance"
#ylab="PFA"
#ylab="Metapopulation Relative Abundance"
ylab="Relative PFA" # Pre-Fisheries Abundance

# Create a list of gradients for each colour 2 to 10 over five steps from 
library(viridis)
#col <- viridis(13)
col <- (inferno(length(frates_vec)))
#col=c("#000000", "#1874CD", "#009ACD", "#00BFFF", "#228B22", "#66CD00", "#7FFF00")
pop.type<-as.numeric(factor(dat$Type))
col.type <- c("black",col.sinks[1],col.sources[1])

ylim=c(0,max(tmp, na.rm=TRUE))
#pop=1
for (pop in 1:npop){
  
  
  ### FISHING ALL
  plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,xaxt='n')
  title("Fishing all", col.main = col.type[pop.type[pop]])
  axis(1,1:length(frates_vec), frates_vec)
  #for (i in 1:ncol(SCN)){
  points(1:length(frates_vec),apply(tmp[,,pop,1],1,mean,na.rm=TRUE),col=col[1],pch=16,bg=col[1],cex=2,type='l',lwd=lwd)
  points(1:length(frates_vec),apply(tmp[,,pop,2],1,mean,na.rm=TRUE),col="grey",pch=16,bg="grey",cex=2,type='l',lwd=lwd)
  points(1:length(frates_vec),apply(tmp[,,pop,5],1,mean,na.rm=TRUE),col="lightgrey",pch=16,bg="lightgrey",cex=2,type='l',lwd=lwd)
  
  legend("topright",c("h=1","h=0.85","h=0.7"),col=c(1,"grey","lightgrey"),lty=rep(1,3),lwd=rep(2,3),bty="n")
  
  
  ### FISHING SOURCES
  plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Fishing source",xaxt='n')
  axis(1,1:length(frates_vec), frates_vec)
  points(1:length(frates_vec),apply(tmp[,,pop,1],1,mean,na.rm=TRUE),col=col[1],pch=16,bg=col[1],cex=2,type='l',lwd=lwd)
  points(1:length(frates_vec),apply(tmp[,,pop,3],1,mean,na.rm=TRUE),col=col.sources[1],pch=16,bg=col.sources[1],cex=2,type='l',lwd=lwd)
  points(1:length(frates_vec),apply(tmp[,,pop,6],1,mean,na.rm=TRUE),col=col.sources[2],pch=16,bg=col.sources[2],cex=2,type='l',lwd=lwd)
  
  legend("topright",c("h=1","h=0.85","h=0.7"),col=c(1,col.sources),lty=rep(1,3),bty="n")
  
  
  ### FISHING SINKS
  plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Fishing sinks",xaxt='n')
  axis(1,1:length(frates_vec), frates_vec)
  points(1:length(frates_vec),apply(tmp[,,pop,1],1,mean,na.rm=TRUE),col=col[1],pch=16,bg=col[1],cex=2,type='l',lwd=lwd)
  points(1:length(frates_vec),apply(tmp[,,pop,4],1,mean,na.rm=TRUE),col=col.sinks[1],pch=16,bg=col.sinks[1],cex=2,type='l',lwd=lwd)
  points(1:length(frates_vec),apply(tmp[,,pop,7],1,mean,na.rm=TRUE),col=col.sinks[2],pch=16,bg=col.sinks[2],cex=2,type='l',lwd=lwd)
  
  legend("topright",c("h=1","h=0.85","h=0.7"),col=c(1,col.sinks),lty=rep(1,3),lwd=rep(2,3),bty="n")
  
  
  
  
  #par(mfrow=c(2,1))
  
  ### DISPERSAL 1
  #plot.new()
  plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,2),xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  text(length(frates_vec)/2,1,paste0(dat$Population[pop]," / ", dat$Type[pop]),col=col.type[pop.type[pop]],cex=3)
  
  
  ### DISPERSAL 0.85 FUNCTION OF LOCAL EXPLOITATION RATE
  plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry = 0.85",xaxt='n')
  axis(1,1:length(frates_vec), frates_vec)
  points(1:length(frates_vec),apply(tmp[,,pop,2],1,mean,na.rm=TRUE),col="grey",pch=16,bg="grey",cex=2,type='l',lwd=lwd)
  points(1:length(frates_vec),apply(tmp[,,pop,3],1,mean,na.rm=TRUE),col=col.sources[1],pch=16,bg=col.sources[1],cex=2,type='l',lwd=lwd)
  points(1:length(frates_vec),apply(tmp[,,pop,4],1,mean,na.rm=TRUE),col=col.sinks[1],pch=16,bg=col.sinks[1],cex=2,type='l',lwd=lwd)
  
  legend("topright",c("all","source","sink"),col=c("grey",col.sources[1],col.sinks[1]),lty=rep(1,3),lwd=rep(2,3),bty="n")
  
  
  ### DISPERSAL 0.7 FUNCTION OF LOCAL EXPLOITATION RATE
  plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry = 0.7",xaxt='n')
  axis(1,1:length(frates_vec), frates_vec)
  points(1:length(frates_vec),apply(tmp[,,pop,5],1,mean,na.rm=TRUE),col="lightgrey",pch=16,bg="grey",cex=2,type='l',lwd=lwd)
  points(1:length(frates_vec),apply(tmp[,,pop,6],1,mean,na.rm=TRUE),col=col.sources[2],pch=16,bg=col.sources[2],cex=2,type='l',lwd=lwd)
  points(1:length(frates_vec),apply(tmp[,,pop,7],1,mean,na.rm=TRUE),col=col.sinks[2],pch=16,bg=col.sinks[2],cex=2,type='l',lwd=lwd)
  
  legend("topright",c("all","source","sink"),col=c("lightgrey",col.sources[2],col.sinks[2]),lty=rep(1,3),lwd=rep(2,3),bty="n")
  
} # end loop pop







#### FISHING STRATEGIES / LOCAL EXPLOITATION RATES ####
#-------------------------------#



#tmp <- metapops.means
tmp <- migs.means


#col=c("#030303", "#454545", "#8C8C8C", "#CCCCCC", "#E3E3E3", "#F7F7F7")
#col=c("#000000", "#1874CD", "#009ACD", "#00BFFF", "#228B22", "#66CD00", "#7FFF00")
col.sources=c("dodgerblue2", "deepskyblue")
col.sinks=c("brown1", "lightcoral", "white")
par(mfcol=c(3,2))

lwd=3

xlab="Local Exploitation rates"

#ylab="Metapopulation Abundance"
#ylab="PFA"
#ylab="Metapopulation Relative Abundance"
ylab="Ratio Immigrants" # Pre-Fisheries Abundance

# Create a list of gradients for each colour 2 to 10 over five steps from 
library(viridis)
#col <- viridis(13)
col <- (inferno(length(frates_vec)))
#col=c("#000000", "#1874CD", "#009ACD", "#00BFFF", "#228B22", "#66CD00", "#7FFF00")
pop.type<-as.numeric(factor(dat$Type))
col.type <- c("black",col.sinks[1],col.sources[1])

ylim=c(0,max(tmp, na.rm=TRUE))


#pop=1
for (pop in 1:npop){
  
  
  ### FISHING ALL
  plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,xaxt='n')
  title("Fishing all", col.main = col.type[pop.type[pop]])
  axis(1,1:length(frates_vec), frates_vec)
  #for (i in 1:ncol(SCN)){
  points(1:length(frates_vec),apply(tmp[,,pop,1],1,mean,na.rm=TRUE),col=col[1],pch=16,bg=col[1],cex=2,type='l',lwd=lwd)
  points(1:length(frates_vec),apply(tmp[,,pop,2],1,mean,na.rm=TRUE),col="grey",pch=16,bg="grey",cex=2,type='l',lwd=lwd)
  points(1:length(frates_vec),apply(tmp[,,pop,5],1,mean,na.rm=TRUE),col="lightgrey",pch=16,bg="lightgrey",cex=2,type='l',lwd=lwd)
  
  legend("topright",c("h=1","h=0.85","h=0.7"),col=c(1,"grey","lightgrey"),lty=rep(1,3),lwd=rep(2,3),bty="n")
  
  
  ### FISHING SOURCES
  plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Fishing source",xaxt='n')
  axis(1,1:length(frates_vec), frates_vec)
  points(1:length(frates_vec),apply(tmp[,,pop,1],1,mean,na.rm=TRUE),col=col[1],pch=16,bg=col[1],cex=2,type='l',lwd=lwd)
  points(1:length(frates_vec),apply(tmp[,,pop,3],1,mean,na.rm=TRUE),col=col.sources[1],pch=16,bg=col.sources[1],cex=2,type='l',lwd=lwd)
  points(1:length(frates_vec),apply(tmp[,,pop,6],1,mean,na.rm=TRUE),col=col.sources[2],pch=16,bg=col.sources[2],cex=2,type='l',lwd=lwd)
  
  legend("topright",c("h=1","h=0.85","h=0.7"),col=c(1,col.sources),lty=rep(1,3),bty="n")
  
  
  ### FISHING SINKS
  plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Fishing sinks",xaxt='n')
  axis(1,1:length(frates_vec), frates_vec)
  points(1:length(frates_vec),apply(tmp[,,pop,1],1,mean,na.rm=TRUE),col=col[1],pch=16,bg=col[1],cex=2,type='l',lwd=lwd)
  points(1:length(frates_vec),apply(tmp[,,pop,4],1,mean,na.rm=TRUE),col=col.sinks[1],pch=16,bg=col.sinks[1],cex=2,type='l',lwd=lwd)
  points(1:length(frates_vec),apply(tmp[,,pop,7],1,mean,na.rm=TRUE),col=col.sinks[2],pch=16,bg=col.sinks[2],cex=2,type='l',lwd=lwd)
  
  legend("topright",c("h=1","h=0.85","h=0.7"),col=c(1,col.sinks),lty=rep(1,3),lwd=rep(2,3),bty="n")
  
  
  
  
  #par(mfrow=c(2,1))
  
  ### DISPERSAL 1
  #plot.new()
  plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,2),xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  text(length(frates_vec)/2,1,paste0(dat$Population[pop]," / ", dat$Type[pop]),col=col.type[pop.type[pop]],cex=3)
  
  
  ### DISPERSAL 0.85 FUNCTION OF LOCAL EXPLOITATION RATE
  plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry = 0.85",xaxt='n')
  axis(1,1:length(frates_vec), frates_vec)
  points(1:length(frates_vec),apply(tmp[,,pop,2],1,mean,na.rm=TRUE),col="grey",pch=16,bg="grey",cex=2,type='l',lwd=lwd)
  points(1:length(frates_vec),apply(tmp[,,pop,3],1,mean,na.rm=TRUE),col=col.sources[1],pch=16,bg=col.sources[1],cex=2,type='l',lwd=lwd)
  points(1:length(frates_vec),apply(tmp[,,pop,4],1,mean,na.rm=TRUE),col=col.sinks[1],pch=16,bg=col.sinks[1],cex=2,type='l',lwd=lwd)
  
  legend("topright",c("all","source","sink"),col=c("grey",col.sources[1],col.sinks[1]),lty=rep(1,3),lwd=rep(2,3),bty="n")
  
  
  ### DISPERSAL 0.7 FUNCTION OF LOCAL EXPLOITATION RATE
  plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,max(tmp, na.rm=TRUE)),xlab=xlab,ylab=ylab,main="Philopatry = 0.7",xaxt='n')
  axis(1,1:length(frates_vec), frates_vec)
  points(1:length(frates_vec),apply(tmp[,,pop,5],1,mean,na.rm=TRUE),col="lightgrey",pch=16,bg="grey",cex=2,type='l',lwd=lwd)
  points(1:length(frates_vec),apply(tmp[,,pop,6],1,mean,na.rm=TRUE),col=col.sources[2],pch=16,bg=col.sources[2],cex=2,type='l',lwd=lwd)
  points(1:length(frates_vec),apply(tmp[,,pop,7],1,mean,na.rm=TRUE),col=col.sinks[2],pch=16,bg=col.sinks[2],cex=2,type='l',lwd=lwd)
  
  legend("topright",c("all","source","sink"),col=c("lightgrey",col.sources[2],col.sinks[2]),lty=rep(1,3),lwd=rep(2,3),bty="n")
  
} # end loop pop












#### IMMIGRANTS DYNAMICS ####
#-------------------------------#


par(mfrow=c(3,3))

tmp <- RatioImExpl

xlab="Time (years)"

#ylab="Metapopulation Abundance"
#ylab="PFA"
#ylab="Metapopulation Relative Abundance"
ylab="Ratio Immigrants/Prelevement" # Pre-Fisheries Abundance

# Create a list of gradients for each colour 2 to 10 over five steps from 
library(viridis)
#col <- viridis(13)
col <- (inferno(length(frates_vec)))
#col=c("#000000", "#1874CD", "#009ACD", "#00BFFF", "#228B22", "#66CD00", "#7FFF00")

pop.type<-as.numeric(factor(dat$Type))
col.type <- c("black",col.sinks[1],col.sources[1])


lwd=1
ylim=c(0,3)#c(0,max(unlist(tmp), na.rm=TRUE))

#pop=14
for (pop in 1:npop){

plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab)#,main="Philopatry=1 / Fishing all")
title("Philopatry=1 / Fishing all", col.main =col.type[pop.type[pop]])
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
j=0
for (iEXPE in SCN[,1]){
  j=j+1
  if(is.na(iEXPE)) next;
  points(1:50,apply(tmp[[1]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
}
#legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)

plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab)#,main="Philopatry=0.85 / Fishing all")
title("Philopatry=0.85 / Fishing all", col.main =col.type[pop.type[pop]])
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
j=0
for (iEXPE in SCN[,2]){
  j=j+1
  if(is.na(iEXPE)) next;
  points(1:50,apply(tmp[[2]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
}
#legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)
abline(h=1,lty=2)

plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim, xlab=xlab,ylab=ylab)#,main="Philopatry=0.7 / Fishing all")
title("Philopatry=0.7 / Fishing all", col.main =col.type[pop.type[pop]])
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
j=0
for (iEXPE in SCN[,5]){
  j=j+1
  if(is.na(iEXPE)) next;
  points(1:50,apply(tmp[[5]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
}
#legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n")
abline(h=1,lty=2)



#plot.new()
plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = "Local Exploitation rates")
#legend("bottomleft",paste0(frates_vec),col=col,lty=rep(1,length(frates_vec)),bty="n",title="Local Exploitation rates",cex=.6)
legend_image <- as.raster(matrix(rev(col), ncol=1))
text(x=0.5, y = 1.1, "Local Exploitation rates", cex=.6)
text(x=0.45, y = seq(0,1,l=length(frates_vec)), labels = paste0(frates_vec),adj=.2,cex=.5)
rasterImage(legend_image, 0.3, 0, .4,1)


plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim, xlab=xlab,ylab=ylab)#,main="Philopatry=0.85 / Fishing sources")
title("Philopatry=0.85 / Fishing sources", col.main =col.type[pop.type[pop]])
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
j=0
for (iEXPE in SCN[,3]){
  j=j+1
  if(is.na(iEXPE)) next;
  points(1:50,apply(tmp[[3]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
}
abline(h=1,lty=2)
#legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)

plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim, xlab=xlab,ylab=ylab)#,main="Philopatry=0.7 / Fishing sources")
title("Philopatry=0.7 / Fishing sources", col.main =col.type[pop.type[pop]])
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
j=0
for (iEXPE in SCN[,6]){
  j=j+1
  if(is.na(iEXPE)) next;
  points(1:50,apply(tmp[[6]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
}
abline(h=1,lty=2)
#legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n")


#plot.new()
plot(NULL,xlim=c(1,length(frates_vec)),ylim=c(0,2),xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
text(length(frates_vec)/2,1,paste0(dat$Population[pop]," / ", dat$Type[pop]),col=col.type[pop.type[pop]],cex=3)

# legend_image <- as.raster(matrix(rev(col), ncol=1))
# plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
# text(x=0.1, y = seq(0,1,l=length(frates_vec)), labels = paste0(frates_vec),adj=.2,cex=.7)
# rasterImage(legend_image, 0, 0, .1,1)

plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim, xlab=xlab,ylab=ylab)#,main="Philopatry=0.85 / Fishing sinks")
title("Philopatry=0.85 / Fishing sinks", col.main =col.type[pop.type[pop]])
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
j=0
for (iEXPE in SCN[,4]){
  j=j+1
  if(is.na(iEXPE)) next;
  points(1:50,apply(tmp[[4]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
}
abline(h=1,lty=2)
#legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")
#legend("bottomleft",paste0(frates_vec),col=col[1:6],lty=rep(1,6),bty="n",title=xlab,cex=1)

plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim, xlab=xlab,ylab=ylab)#,main="Philopatry=0.7 / Fishing sinks")
title("Philopatry=0.7 / Fishing sinks", col.main =col.type[pop.type[pop]])
rect(0,0,10,max(unlist(tmp), na.rm=TRUE), col= "lightgrey",border ="lightgrey")
j=0
for (iEXPE in SCN[,7]){
  j=j+1
  if(is.na(iEXPE)) next;
  points(1:50,apply(tmp[[7]][[paste0(iEXPE)]][[pop]],1,mean,na.rm=TRUE),col=col[j],type='l',lwd=lwd)
}
abline(h=1,lty=2)
# #legend("bottomleft",names(tmp[[1]]),col=col[1:6],lty=rep(1,6),bty="n")

#legend("bottomleft",paste0(frates_vec),col=col[1:length(frates_vec)],lty=rep(1,length(frates_vec)),bty="n")

} # end loop pop






tmp <- pops.means
tmp2 <- migexpl.means


#col=c("#030303", "#454545", "#8C8C8C", "#CCCCCC", "#E3E3E3", "#F7F7F7")
#col=c("#000000", "#1874CD", "#009ACD", "#00BFFF", "#228B22", "#66CD00", "#7FFF00")
col.sources=c("dodgerblue2", "deepskyblue")
col.sinks=c("brown1", "lightcoral", "white")
par(mfcol=c(3,2))

lwd=3

xlab="Ratio Immigrants/Prelevement"

#ylab="Metapopulation Abundance"
#ylab="PFA"
ylab="Population Relative Abundance"
#ylab="Ratio Immigrants" # Pre-Fisheries Abundance

# Create a list of gradients for each colour 2 to 10 over five steps from 
library(viridis)
#col <- viridis(13)
col <- (inferno(length(frates_vec)))
#col=c("#000000", "#1874CD", "#009ACD", "#00BFFF", "#228B22", "#66CD00", "#7FFF00")
pop.type<-as.numeric(factor(dat$Type))
col.type <- c("black",col.sinks[1],col.sources[1])

ylim=c(0,2)
xlim=c(0,4)

#pop=3
for (pop in 1:npop){

par(mfcol=c(1,1))
### FISHING ALL
plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)#,xaxt='n')
title(paste0(dat$Population[pop]," / ", dat$Type[pop]), col.main = col.type[pop.type[pop]])
#axis(1,apply(tmp[,,pop,1],1,mean,na.rm=TRUE), frates_vec)
#for (i in 1:ncol(SCN)){
#points(apply(tmp2[,,pop,1],1,mean,na.rm=TRUE),apply(tmp[,,pop,1],1,mean,na.rm=TRUE),col=col[1],pch=16,bg=col[1],cex=1+frates_vec)#,type='l',lwd=lwd)
points(apply(tmp2[,,pop,2],1,mean,na.rm=TRUE),apply(tmp[,,pop,2],1,mean,na.rm=TRUE),col="grey",pch=16,bg="grey",cex=1+frates_vec)#,type='l',lwd=lwd)
points(apply(tmp2[,,pop,5],1,mean,na.rm=TRUE),apply(tmp[,,pop,5],1,mean,na.rm=TRUE),col="lightgrey",pch=16,bg="lightgrey",cex=1+frates_vec)#,type='l',lwd=lwd)
abline(h=1,lty=2);abline(v=1,lty=2)
legend("topleft",c("h=0.85","h=0.7"),col=c("grey","lightgrey"),pch=rep(16,2),bty="n")

#text(apply(tmp2[,,pop,2],1,mean,na.rm=TRUE),apply(tmp[,,pop,2],1,mean,na.rm=TRUE),frates_vec,col="grey",cex=1)
#text(apply(tmp2[,,pop,5],1,mean,na.rm=TRUE),apply(tmp[,,pop,5],1,mean,na.rm=TRUE),frates_vec,col="lightgrey",cex=1)

if(dat$Type[pop]=="source" | dat$Type[pop]=="neutral"){
  ### FISHING SOURCES
  #plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)#,xaxt='n')
  #title("Fishing Source", col.main = col.type[pop.type[pop]])
  #axis(1,apply(tmp[,,pop,1],1,mean,na.rm=TRUE), frates_vec)
  #points(apply(tmp2[,,pop,1],1,mean,na.rm=TRUE),apply(tmp[,,pop,1],1,mean,na.rm=TRUE),col=col[1],pch=16,bg=col[1],cex=1+frates_vec)#,type='l',lwd=lwd)
  points(apply(tmp2[,,pop,3],1,mean,na.rm=TRUE),apply(tmp[,,pop,3],1,mean,na.rm=TRUE),col=col.sources[1],pch=16,bg=col.sources[1],cex=1+frates_vec)#,type='l',lwd=lwd)
  points(apply(tmp2[,,pop,6],1,mean,na.rm=TRUE),apply(tmp[,,pop,6],1,mean,na.rm=TRUE),col=col.sources[2],pch=16,bg=col.sources[2],cex=1+frates_vec)#,type='l',lwd=lwd)
  #abline(h=1,lty=2);abline(v=1,lty=2)
  legend("topright",c("h=0.85","h=0.7"),col=c(col.sources),pch=rep(16,2),bty="n")
  
  #text(apply(tmp2[,,pop,3],1,mean,na.rm=TRUE),apply(tmp[,,pop,3],1,mean,na.rm=TRUE),frates_vec,col=col.sinks[1],cex=1)
  #text(apply(tmp2[,,pop,6],1,mean,na.rm=TRUE),apply(tmp[,,pop,6],1,mean,na.rm=TRUE),frates_vec,col=col.sinks[2],cex=1)
  
}

if(dat$Type[pop]=="sink" | dat$Type[pop]=="neutral"){
  ### FISHING SINKS
  #plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main="Fishing sinks")#,xaxt='n')
  #axis(1,apply(tmp2[,,pop,1],1,mean,na.rm=TRUE), frates_vec)
  #points(apply(tmp2[,,pop,1],1,mean,na.rm=TRUE),apply(tmp[,,pop,1],1,mean,na.rm=TRUE),col=col[1],pch=16,bg=col[1],cex=1+frates_vec)#,type='l',lwd=lwd)
  points(apply(tmp2[,,pop,4],1,mean,na.rm=TRUE),apply(tmp[,,pop,4],1,mean,na.rm=TRUE),col=col.sinks[1],pch=16,bg=col.sinks[1],cex=1+frates_vec)#,type='l',lwd=lwd)
  points(apply(tmp2[,,pop,7],1,mean,na.rm=TRUE),apply(tmp[,,pop,7],1,mean,na.rm=TRUE),col=col.sinks[2],pch=16,bg=col.sinks[2],cex=1+frates_vec)#,type='l',lwd=lwd)
  abline(h=1,lty=2);abline(v=1,lty=2)
  legend("topright",c("h=0.85","h=0.7"),col=c(col.sinks),pch=rep(16,2),bty="n")
  
  #text(apply(tmp2[,,pop,4],1,mean,na.rm=TRUE),apply(tmp[,,pop,4],1,mean,na.rm=TRUE),frates_vec,col=col.sinks[1],cex=1)
  #text(apply(tmp2[,,pop,7],1,mean,na.rm=TRUE),apply(tmp[,,pop,7],1,mean,na.rm=TRUE),frates_vec,col=col.sinks[2],cex=1)
}



} # end loop pop


dev.off()


