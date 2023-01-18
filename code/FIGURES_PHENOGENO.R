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
#SCN<-SCN[,1]

##### METAPOPULATION ABUNDANCE

#EXPE=c(101, 102, 103, 104, 105, 106)
#401 402 403 404 405 406 )

LFPARR<-LFSMOLT<-LFMSW<-LF1SW<-LFRETURNS<-RMSW<-list()
LFparr.means=LFsmolt.means=LFreturns.means=LFmsw.means=LF1sw.means=rmsw.means=array(,dim=c(nrow(SCN),nSIMUL,npop,ncol(SCN)))
#for (i in 1:ncol(SCN)){
  for (i in 1:1){

  LFparr.pop=LFsmolt.pop=LFreturns.pop=LFmsw.pop=LF1sw.pop=Rmsw.pop=list(list())
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
    

  } # loop iEXPE
  
  LFPARR[[i]]<-LFparr.pop
  LFSMOLT[[i]]<-LFsmolt.pop
  LFRETURNS[[i]]<-LFreturns.pop
  LFMSW[[i]]<-LFmsw.pop
  LF1SW[[i]]<-LF1sw.pop
  RMSW[[i]]<-Rmsw.pop
} # loop SCN






load(paste0("results/DEMOGRAPHY",101,".RData"))
init <- apply(tail(LFmsw[[paste0(101)]][[1]],10),2,mean,na.rm=TRUE)
tmp<-init/sum(init)
prop.sink <- sum(tmp[dat$Type!="source"],na.rm=TRUE)
prop.source <- sum(tmp[dat$Type!="sink"],na.rm=TRUE)


  





#pdf(file="results/FishStrat_results_phenotype.pdf")






#### POPULATION DYNAMICS ####
#-------------------------------#
col.sources=c("dodgerblue2", "deepskyblue")
col.sinks=c("brown1", "lightcoral", "white")


xlab="Time (years)"

#ylab="Metapopulation Abundance"
#ylab="PFA"
#ylab="Metapopulation Relative Abundance"
ylab="Fork length (mm)" # Pre-Fisheries Abundance

# Create a list of gradients for each colour 2 to 10 over five steps from 
library(viridis)
#col <- viridis(13)
col <- (inferno(length(frates_vec)))
#col=c("#000000", "#1874CD", "#009ACD", "#00BFFF", "#228B22", "#66CD00", "#7FFF00")

pop.type<-as.numeric(factor(dat$Type))
col.type <- c("black",col.sinks[1],col.sources[1])

lwd=1








par(mfrow=c(4,5))
for (pop in 1:npop){
  
  main=paste0("Parr (", pops[pop],")")
  #tmp <- LFPARR
  tmp<-apply(LFPARR[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
  ylim=c(80,100)#c(0,max(unlist(tmp), na.rm=TRUE))
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  rect(0,0,10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
  points(1:50,tmp,col=1,type='l',lwd=lwd)
  
  main=paste0("Smolt (", pops[pop],")")
  tmp <- apply(LFSMOLT[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
  ylim=c(120,150)#c(0,max(ylim))
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  rect(0,0,10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
  points(1:50,tmp,col=1,type='l',lwd=lwd)
  
  main=paste0("1SW (", pops[pop],")")
  tmp <- apply(LF1SW[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
  ylim=c(550,650)#c(0,max(ylim))
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  rect(0,0,10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
  points(1:50,tmp,col=1,type='l',lwd=lwd)
  
  main=paste0("MSW (", pops[pop],")")
  tmp <- apply(LFMSW[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
  ylim=c(740,800)#c(0,max(ylim))
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  rect(0,0,10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
  points(1:50,tmp,col=1,type='l',lwd=lwd)
  
  
  main=paste0("Ratio MSW/1SW (", pops[pop],")")
  tmp <- apply(RMSW[[1]][["101"]][[pop]],1,mean,na.rm=TRUE)
  ylim=c(0,.5)#c(0,max(ylim))
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  rect(0,0,10,max(ylim)+10, col= "lightgrey",border ="lightgrey")
  points(1:50,tmp,col=1,type='l',lwd=lwd)
}









lwd=1
ylim=c(0.9,1.1)#c(0,max(unlist(tmp), na.rm=TRUE))

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



dev.off()



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


