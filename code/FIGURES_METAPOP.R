## remove (almost) everything in the working environment.
rm(list = ls())
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  args <- 1
}


#setwd("/results/METAPOP/")

source("R/Rfunctions.R")

source("parIbasam.R")
h = c(1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7)
pops <- c("Leff", "Trieux", "Jaudy", "Leguer", "Yar", "Douron", "Penze", "Elorn", "Aulne", "Goyen", "Odet", "Aven", "Laita", "Scorff", "Blavet")
nSIMUL <- 30
nYear <- 40
nInit <- 10 # Nb years at initialization
npop <- 15

#h = c(1, .95, .9, .85, .8, .75, .7) # Philopatry (homing) rates #edit al - 22/03/21 5% and 15%
# scenarioConnect=7 #scenario 1 for Disp. 0%.00, scenario 2 for h=0.95, scenario 3 pour h=0.80
# 
# # 0: control / 
# # 1: no fishing on sink populations 
# # 2: no fishing on source populations 
# scenarioFishing = 0

frates_vec =c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)
# scenarioFishingRate = 1

# Sc_name <- paste0("FishStrat_",scenarioConnect,scenarioFishing,scenarioFishingRate)


scn=c(
  101, 102, 103, 104, 105, 106, 107, 108, 109, 1010, 1011, 1012, NA, NA, NA, NA, NA, NA, NA, NA, NA
  #,201, 212, 213, 214, 215, 216, 217, 218, 219, 2121, 2111, 2112, NA, NA, NA, NA, NA, NA, NA, NA, NA
  #,201, 222, 223, 224, 225, 226, 227, 228, 229, 2222, 2211, 2212, NA, NA, NA, NA, NA, NA, NA, NA, NA
  ,401, 402, 403, 404, 405, 406, 407, 408, 409, 4010, 4011, 4012, NA, NA, NA, NA, NA, NA, NA, NA, NA
  #,401, 412, 413, 414, 415, 416, 417, 418, 419, 4110, 4111, 4112, 4113, 4114, 4115, 4116, 4117, 4118, 4119, 4120, 4121
  ,401, 422, 423, 424, 425, 426, 427, 428, 429, 4210, 4211, 4212, 4213, 4214, 4215, 4216, 4217, 4218, 4219, 4220, 4221
  ,401, 432, 433, 434, 435, 436, 437, 438, 439, 4310, 4311, 4312, 4313, 4314, 4315, 4316, 4317, 4318, 4319, 4320, 4321
  
  ,701, 702, 703, 704, 705, 706, 707, 708, 709, 7010, 7011, 7012, NA, NA, NA, NA, NA, NA, NA, NA, NA
  #,701, 712, 713, 714, 715, 716, 717, 718, 719, 7110, 7111, 7112, 7113, 7114, 7115, 7116, 7117, 7118, 7119, 7120, 7121
  ,701, 722, 723, 724, 725, 726, 727, 728, 729, 7210, 7211, 7212, 7213, 7214, 7215, 7216, 7217, 7218, 7219, 7220, 7221
  ,701, 732, 733, 734, 735, 736, 737, 738, 739, 7310, 7311, 7312, 7313, 7314, 7315, 7316, 7317, 7318, 7319, 7320, 7321
  
)
SCN <-matrix(scn,length(frates_vec),7)

frates_vec =c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)#, 0.55)
SCN <- SCN[1:length(frates_vec),]

##### METAPOPULATION ABUNDANCE

RatioMigrants <- list()
RatioMigrantsVsExpl <- list()
## Demography
nParr <- nParrMature <- list() # Nb returns
nSmolt <- list() # Nb returns
nReturns <- list() # Nb returns
Mig <- list() # Nb immigrnats
nMSW<- list() # Nb immigrnats
n1SW<- list() # Nb immigrnats
#NRet.res <- list() # Nb returns
#PE <- list()
gExpl <- dExpl <- list()
PE <- list()
Mig <- list()

## Phenotypes
LFParr <- list()
LFSmolt <- list()
LFReturns <- list()
LF1sw <- list()
LFmsw <- list()
rMSW<- list() # Nb immigrnats

## Genotypes
gNEUTRAL <- list() # Nb returns
gFMID1 <- list() # Nb returns
gFMID2 <- list() # Nb returns
gFMID3 <- list() # Nb immigrnats
gFMID4<- list() # Nb immigrnats
gGROWTH<- list() # Nb immigrnats


CVReturns<-CVParr<-CVSmolt<-CVmsw<-CV1sw<-CVgNEUTRAL<-CVgFMID1<-CVgFMID2<-CVgFMID3<-CVgFMID4<-CVgGROWTH<-list()


for (i in 1:ncol(SCN)){
  #for (i in 1:3){
  j=0
  #for (iEXPE in EXPE){ # Loop over scenario 
  for (iEXPE in SCN[,i]){
    j=j+1
    if(is.na(iEXPE)) next;
    
    load(paste0("results/METAPOP/METAPOP",iEXPE,".RData"))
    #load(paste0("results/EXPLOIT/EXPLOIT",iEXPE,".RData"))
    load(paste0("results/MIGRANTS/MIGRANTS",iEXPE,".RData"))
    
    
    RatioMigrants[[paste0(iEXPE)]] <- ratioMigrants
    RatioMigrantsVsExpl[[paste0(iEXPE)]] <- ratioMigrantsVsExpl
    
    nParr[[paste0(iEXPE)]] <- Nparr0
    nParrMature[[paste0(iEXPE)]] <- Nparr0.mature
    nSmolt[[paste0(iEXPE)]] <- Nsmolt0
    #Mig[[paste0(iEXPE)]] <- list(Nhom=Nhomers, NIm=NIm)
    nReturns[[paste0(iEXPE)]] <- NRet
    nMSW[[paste0(iEXPE)]] <- NMSW
    n1SW[[paste0(iEXPE)]] <- N1SW
    
    Mig[[paste0(iEXPE)]] <- RETURNS
    
    gExpl[[paste0(iEXPE)]] <- Pexpl.global
    #dExpl[[paste0(iEXPE)]] <- tmp1
    PE[[paste0(iEXPE)]] <- PE_mv
    
    ## Phenotypes
    LFParr[[paste0(iEXPE)]] <- LFparr0
    LFSmolt[[paste0(iEXPE)]] <- LFsmolt0
    LFReturns[[paste0(iEXPE)]] <- LFreturns
    LFmsw[[paste0(iEXPE)]] <- LFMSW
    LF1sw[[paste0(iEXPE)]] <- LF1SW
    rMSW[[paste0(iEXPE)]] <- RatioMSW
    
    CVParr[[paste0(iEXPE)]] <- CVparr0
    CVSmolt[[paste0(iEXPE)]] <- CVsmolt0
    CVReturns[[paste0(iEXPE)]] <- CVreturns
    CVmsw[[paste0(iEXPE)]] <- CVMSW
    CV1sw[[paste0(iEXPE)]] <- CV1SW
    
    
    ## Genotypes
    gNEUTRAL[[paste0(iEXPE)]] <- gNeutral
    gFMID1[[paste0(iEXPE)]] <- gFmid1
    gFMID2[[paste0(iEXPE)]] <- gFmid2
    gFMID3[[paste0(iEXPE)]] <- gFmid3
    gFMID4[[paste0(iEXPE)]] <- gFmid4
    gGROWTH[[paste0(iEXPE)]] <- gG
    
    CVgNEUTRAL[[paste0(iEXPE)]] <- CVgNeutral
    CVgFMID1[[paste0(iEXPE)]] <- CVgFmid1
    CVgFMID2[[paste0(iEXPE)]] <- CVgFmid2
    CVgFMID3[[paste0(iEXPE)]] <- CVgFmid3
    CVgFMID4[[paste0(iEXPE)]] <- CVgFmid4
    CVgGROWTH[[paste0(iEXPE)]] <- CVgG
    
  } # loop iEXPE
  
} # loop SCN

# contrib=matrix(,npop,50)
# nRetTot <- aggregate(count~year,tmp1,sum)
# for (y in 1:max(tmp1$year)){
#   contrib[1:npop,y] <- (tmp1$count[tmp1$year==y]/nRetTot[y,"count"])*100
# }




# load(paste0("results/PHENOTYPE",101,".RData"))
# init <- apply(tail(LFmsw[[paste0(101)]][[1]],10),2,mean,na.rm=TRUE)
# tmp<-init/sum(init)
# prop.sink <- sum(tmp[dat$Type!="source"],na.rm=TRUE)
# prop.source <- sum(tmp[dat$Type!="sink"],na.rm=TRUE)
# 



col.sources=c("dodgerblue2", "deepskyblue")
col.sinks=c("brown1", "lightcoral", "white")


# Create a list of gradients for each colour 2 to 10 over five steps from 
library(viridis)
#col <- viridis(13)
col <- (inferno(length(frates_vec)))
#col=c("#000000", "#1874CD", "#009ACD", "#00BFFF", "#228B22", "#66CD00", "#7FFF00")

pop.type<-as.numeric(factor(dat$Type))
col.type <- c("black",col.sinks[1],col.sources[1])













pdf(file="results/METAPOP.pdf")

#### PORTFOLIO ####
#-------------------------------#
# 0: control /
# 1: no fishing on sink populations = fishing source and neutral
# 2: no fishing on source populations = fishing sink and neutral
# 3: no fishing on sink populations & stronger fishing pressure on sources

#col.type <- c("#FF6347", "#4F94CD", "#48D1CC") # tomato/blue/green
col.type <- c( "#48D1CC",  "#4F94CD","#FF6347") # tomato/blue/green

lgd<-function(){
  legend("topright"
         #,c("fishing all","fishing sink+neutral","fishing source+neutral","fishing source only")
         #,c("fishing all","fishing sink+neutral","fishing source only")
         ,c("Exploitation (all pop.)","Conservation of source pop.","Exploitation of source pop.")
         ,col=c("black",col.type[c(1,3)])
         ,text.col=c("black",col.type[c(1,3)])
         #,pch=rep(16,5)
         ,lty=rep(1,4)
         ,bty="n"
         ,cex=0.5)
}
distPlot <- function(nscn, X, ylab,synchrony=FALSE,plot=TRUE){
  
  xlab="Metapopulation exploitation rate"
  xlim=c(0,max(frates_vec))
  lwd=2
  
  
  EXPE <- SCN[,nscn]
  scn_id <- as.numeric(strsplit(as.character(EXPE[2]), "")[[1]])
  scenarioConnect = scn_id[1]
  scenarioFishing = scn_id[2]
  scenarioFishingRate = scn_id[3]
  
  if (scenarioFishing==0){ fishingStrat <- "Expl. All"}
  if (scenarioFishing==1){ fishingStrat <- "Expl. Source+neutral"}
  if (scenarioFishing==2){ fishingStrat <- "Expl. Sink+neutral"}
  if (scenarioFishing==3){ fishingStrat <- "Expl. Source only"}
  
  #ylab="# Parr 0+"
  #main=paste0("h=",h[scenarioConnect],"/",fishingStrat)#paste0("Neutral gene (", pops[pop],")")
  main=""
  if(plot==TRUE){
    #ylim=c(0,1.1)#c(min(sapply(X,min,na.rm=TRUE)),max(sapply(X,max,na.rm=TRUE)))#c(0,max(unlist(tmp), na.rm=TRUE))
    xlim=c(0,.3)
    plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
    axis(1,frates_vec,frates_vec)
    #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
    lgd()
  }
  
  
  tmp=tmp2.5=tmp97.5=pGlob=NULL
  j=0
  for (i in EXPE){
    j=j+1
    if(is.na(i)) next;
    
    #var<-X[[paste0(i)]]
    varTail=NULL
    
    if(synchrony) {
      for (sim in 1:nSIMUL) {varTail[sim]<-X[[paste0(i)]][[sim]]$Synchrony}
    } else {
      for (sim in 1:nSIMUL) {varTail[sim]<-X[[paste0(i)]][[sim]]$pe}
    }
    #varTail <- apply(var, 2, function(x) mean(tail(x , n = 5),na.rm=TRUE))
    
    # tmp[j]<-quantile(varTail,probs=0.5,na.rm=TRUE)
    # tmp2.5[j]<-quantile(varTail,probs=0.025,na.rm=TRUE)
    # tmp97.5[j]<-quantile(varTail,probs=0.975,na.rm=TRUE)
    
    tmp[j]<-mean(varTail,na.rm=TRUE)
    tmp2.5[j]<-tmp[j]-se(varTail)
    tmp97.5[j]<-tmp[j]+se(varTail)
    
    pGlob[j] <- mean(gExpl[[paste0(i)]],na.rm=TRUE)
    
    # if(i==EXPE[1] & j==1) ref <- tmp[1]
  }
  
  #colfunc <- colorRampPalette(c("darkgrey", "black"))
  if (nscn==1) mycol <- "lightgrey"
  
  if (nscn==2) mycol <- "darkgrey"
  #if (nscn==3) mycol <- col.type[2] # source+neutral
  if (nscn==3) mycol <- makeTransparent(col.type[1], 80) #col.type[1] # sink+neutral
  if (nscn==4) mycol <- makeTransparent(col.type[3], 80) #col.type[3] # source only
  
  if (nscn==5) mycol <- "black"
  #if (nscn==7) mycol <- makeTransparent(col.type[2], 80) # source+neutral
  if (nscn==6) mycol <- col.type[1] #makeTransparent(col.type[1], 80) # sink+neutral
  if (nscn==7) mycol <- col.type[3] #makeTransparent(col.type[3], 80) # source only
  
  #colfunc <- colorRampPalette(c(col.sinks[1], col.sinks[2]))}
  #  mycol <- col.sinks[1]} else { mycol <- col.sinks[2]}
  x=pGlob
  polygon(c(x,rev(x)),c(tmp2.5,rev(tmp97.5)),border=NA, col=makeTransparent(mycol,20))
  
  
  #scatter.smooth(tmp ~ pGlob, span = 2/3, degree = 2,col=mycol)
  lines(loess.smooth(pGlob,tmp,span = 2/3, degree = 1), col = mycol,lwd=3)
  
  #if(text){
  text(pGlob,tmp,frates_vec, col=mycol,cex=.6)
  #}
  #points(pGlob,tmp,col=mycol,type='p',lwd=1,pcDisp. 0%,cex=.6)
  
  #abline(h=1,lty=2,col="grey")
  #abline(h=ref,lty=2,col="grey")
  #segments(pGlob,tmp2.5,pGlob,tmp97.5,col=mycol)
  #points(pGlob,tmp,col=mycol,type='b',lwd=lwd,pch=20)
  
}

distPred <- function(nscn, X, pglob, ylab,synchrony=FALSE,plot=TRUE){
  
  xlab=""
  #xlim=c(0,7)#max(frates_vec))
  lwd=2
  
  
  EXPE <- SCN[,nscn]
  scn_id <- as.numeric(strsplit(as.character(EXPE[2]), "")[[1]])
  scenarioConnect = scn_id[1]
  scenarioFishing = scn_id[2]
  scenarioFishingRate = scn_id[3]
  
  if (scenarioFishing==0){ fishingStrat <- "Expl. All"}
  if (scenarioFishing==1){ fishingStrat <- "Expl. Source+neutral"}
  if (scenarioFishing==2){ fishingStrat <- "Expl. Sink+neutral"}
  if (scenarioFishing==3){ fishingStrat <- "Expl. Source only"}
  
  #ylab="# Parr 0+"
  #main=paste0("h=",h[scenarioConnect],"/",fishingStrat)#paste0("Neutral gene (", pops[pop],")")
  main=""
  if(plot==TRUE){
    #ylim=c(0,1.1)#c(min(sapply(X,min,na.rm=TRUE)),max(sapply(X,max,na.rm=TRUE)))#c(0,max(unlist(tmp), na.rm=TRUE))
    xlim=c(0, ncol(SCN)+0.5)
    plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
    xtext<-c("Control","All"
             #,"Source & neutral"
             ,"Sink & neutral","Source only",
             "All"
             #,"Source & neutral"
             ,"Sink & neutral","Source only"
    )
    # xtext<-c("Control","No conservation (fishing all)"
    #          #,"Source & neutral"
    #          ,"Conserv. Sources","Conserv. Sinks+Neutrals",
    #          "No conservation (fishing all)"
    #          #,"Source & neutral"
    #          ,"Conserv. Sources","Conserv. Sinks+Neutrals"
    # )
    axis(1,labels=F,at=1:length(xtext),las=1,cex=0.4)
    #axis(1,line = 3, labels="Management strategies",at=4,las=1,cex=0.6)
    # text(1:length(xtext), min(ylim) - 0.1*min(ylim), xtext, srt = 45, xpd = T,cex=0.7
    #      ,adj = c(1,0),pos=2
    #      , col=c("black","darkgrey",col.type[c(1,3)],"lightgrey",makeTransparent(col.type[c(1,3)], 95))
    # )
    # mtext(xtext,side=1,line=1,outer=FALSE,at=1:length(xtext),adj=1,las=3,cex=.5
    #       ,col=c("black","darkgrey",col.type[c(1,3)],"lightgrey",makeTransparent(col.type[c(1,3)], 95))
    # )
    
    # xtext2<-c("Control","Exploitation"
    #           #,"Source & neutral"
    #           ,"Conservation","Exploitation",
    #           "Exploitation"
    #           #,"Source & neutral"
    #           ,"Conservation","Exploitation"
    # )
    # # mtext(xtext2,
    # #       side=1,line=0.5,outer=FALSE,at=1:length(xtext)-0.3,adj=0,las=1,cex=.28
    # #       ,col=c("black","black",col.type[c(1,3)],"black",col.type[c(1,3)])
    # # )
    # text(1:length(xtext2), min(ylim) - 0.1*min(ylim), xtext2, srt = 45, xpd = T,cex=0.7
    #      ,adj = 1,pos=2
    #      ,col=c("black","black",col.type[c(1,3)],"black",col.type[c(1,3)])
    # )
    # 
    # xtext2<-c("","all pop."
    #           #,"Source & neutral"
    #           ,"Source pop.","Source pop.",
    #           "all pop."
    #           #,"Source & neutral"
    #           ,"Source pop.","Source pop."
    # )
    # # mtext(xtext2,
    # #       side=1,line=1,outer=FALSE,at=1:length(xtext)-0.3,adj=0,las=1,cex=.28
    # #       ,col=c("black","black",col.type[c(1,3)],"black",col.type[c(1,3)])
    # # )
    # # 
    # text(1:length(xtext2)+0.3, min(ylim) - 0.1*min(ylim), xtext2, srt = 45, xpd = T,cex=0.7
    #      ,adj = 1,pos=2
    #      ,col=c("black","black",col.type[c(1,3)],"black",col.type[c(1,3)])
    # )
    
    #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
    #lgd()
    abline(v=1.5,lty=3)
    abline(v=4.5,lty=3)
    
    #text(c(0.5,3, 6),rep(max(ylim)*0.95,3),c("Disp. 0%","Disp. 15%","Disp. 30%"),cex=1,col="darkgrey")
    text(c(0.5,3, 6),rep(max(ylim)*0.95,3),c("Dispersal: 0%","Dispersal: 15%","Dispersal: 30%"),cex=0.5,col="darkgrey")
  }
  
  
  
  
  tmp=tmp2.5=tmp97.5=pGlob=NULL
  j=0
  for (i in EXPE){
    j=j+1
    if(is.na(i)) next;
    
    #var<-X[[paste0(i)]]
    varTail=NULL
    
    if(synchrony) {
      for (sim in 1:nSIMUL) {varTail[sim]<-X[[paste0(i)]][[sim]]$Synchrony}
    } else {
      for (sim in 1:nSIMUL) {varTail[sim]<-X[[paste0(i)]][[sim]]$pe}
    }
    #varTail <- apply(var, 2, function(x) mean(tail(x , n = 5),na.rm=TRUE))
    
    # tmp[j]<-quantile(varTail,probs=0.5,na.rm=TRUE)
    # tmp2.5[j]<-quantile(varTail,probs=0.025,na.rm=TRUE)
    # tmp97.5[j]<-quantile(varTail,probs=0.975,na.rm=TRUE)
    
    tmp[j]<-mean(varTail,na.rm=TRUE)
    tmp2.5[j]<-tmp[j]-se(varTail)
    tmp97.5[j]<-tmp[j]+se(varTail)
    
    pGlob[j] <- mean(gExpl[[paste0(i)]],na.rm=TRUE)
    
    # if(i==EXPE[1] & j==1) ref <- tmp[1]
  }
  
  #colfunc <- colorRampPalette(c("darkgrey", "black"))
  if (nscn==1) mycol <- "lightgrey"
  
  if (nscn==2) mycol <- "darkgrey"
  #if (nscn==3) mycol <- col.type[2] # source+neutral
  if (nscn==3) mycol <- col.type[1] # sink+neutral
  if (nscn==4) mycol <- col.type[3] # source only
  
  if (nscn==5) mycol <- "black"
  #if (nscn==7) mycol <- makeTransparent(col.type[2], 80) # source+neutral
  if (nscn==6) mycol <- col.type[1]#makeTransparent(col.type[1], 80) # sink+neutral
  if (nscn==7) mycol <- col.type[3]#makeTransparent(col.type[3], 80) # source only
  
  #colfunc <- colorRampPalette(c(col.sinks[1], col.sinks[2]))}
  #  mycol <- col.sinks[1]} else { mycol <- col.sinks[2]}
  
  
  #scatter.smooth(tmp ~ pGlob, span = 2/3, degree = 2,col=mycol)
  #lines(loess.smooth(pGlob,tmp,span = 2/3, degree = 1), col = mycol,lwd=3)
  lo <- loess(tmp ~ pGlob)
  pred <- predict(lo, data.frame(pGlob = pglob), se = TRUE)
  #predict(lo, data.frame(pGlob = seq(5, 30, 1)), se = TRUE)
  
  #if(text){
  #text(pGlob,tmp,frates_vec, col=mycol,cex=.6)
  #}
  #points(pGlob,tmp,col=mycol,type='p',lwd=1,pcDisp. 0%,cex=.6)
  
  #abline(h=1,lty=2,col="grey")
  #abline(h=ref,lty=2,col="grey")
  segments(nscn,pred$fit-pred$se.fit,nscn,pred$fit+pred$se.fit,col=mycol)
  points(nscn,pred$fit,col=mycol,pch=20)
  
}


Xlab<-function(y,cex){
  xtext2<-c("Control","Exploitation"
            #,"Source & neutral"
            ,"Conservation","Exploitation",
            "Exploitation"
            #,"Source & neutral"
            ,"Conservation","Exploitation"
  )
  # mtext(xtext2,
  #       side=1,line=0.5,outer=FALSE,at=1:length(xtext)-0.3,adj=0,las=1,cex=.28
  #       ,col=c("black","black",col.type[c(1,3)],"black",col.type[c(1,3)])
  # )
  text(1:length(xtext2), y, xtext2, srt = 45, xpd = T,cex=cex
       ,adj = 1,pos=2
       ,col=c("black","black",col.type[c(1,3)],"black",col.type[c(1,3)])
  )
  
  xtext2<-c("","all pop."
            #,"Source & neutral"
            ,"Source pop.","Source pop.",
            "all pop."
            #,"Source & neutral"
            ,"Source pop.","Source pop."
  )
  # mtext(xtext2,
  #       side=1,line=1,outer=FALSE,at=1:length(xtext)-0.3,adj=0,las=1,cex=.28
  #       ,col=c("black","black",col.type[c(1,3)],"black",col.type[c(1,3)])
  # )
  # 
  text(1:length(xtext2)+0.3, y, xtext2, srt = 45, xpd = T,cex=cex
       ,adj = 1,pos=2
       ,col=c("black","black",col.type[c(1,3)],"black",col.type[c(1,3)])
  )
  
}



# jpeg(filename = "article/images/PE.jpeg",
#      width = 580, height = 780, 
#      quality = 100)
# ggsave("article/images/PE.png", device = "png",
#        width = 8, height = 6, units = "in", dpi = 600, type = "cairo-png")

par(mfrow = c(2,2), mar=c(4, 4, 1, 1) + 0.1)#it goes c(bottom, left, top, right)
#layout(matrix(c(1,2,3,3), nrow = 2, ncol = 2, byrow = TRUE))

x<-PE
ylab <- "PORTFOLIO EFFECT"
ylim=c(1.5,2.5)

#### ALL
distPlot(1,x,ylab,FALSE)
#abline(h=1,lty=2,col="grey")
text(0,max(ylim),labels="A",cex=1.5)

#### Disp. 15%
#plot.new()
distPlot(2,x,ylab,FALSE,FALSE)
##distPlot(3,x,ylab,FALSE,FALSE) # source+neutral
distPlot(3,x,ylab,FALSE,FALSE) # sink+neutral
distPlot(4,x,ylab,FALSE,FALSE) # source only

#### Disp. 30%
#plot.new()
distPlot(5,x,ylab,FALSE,FALSE)
##distPlot(7,x,ylab,FALSE,FALSE) # source+neutral
distPlot(6,x,ylab,FALSE,FALSE) # sink+neutral
distPlot(7,x,ylab,FALSE,FALSE) # source only

#### PREDICTION
pglob=0.15
abline(v=pglob,lty=2,col="grey")
ylab <- ""#paste0("Portfolio effect (pGlob=",pglob,")")
ylim=c(1.5,2.5)
####ALL
distPred(1,x, pglob,ylab,FALSE)
#### Disp. 15%
#plot.new()
distPred(2,x, pglob,ylab,FALSE,FALSE)
##distPred(3,x, pglob,ylab,FALSE,FALSE) # source+neutral
distPred(3,x, pglob,ylab,FALSE,FALSE) # sink+neutral
distPred(4,x, pglob,ylab,FALSE,FALSE) # source only

#### Disp. 30%
#plot.new()
distPred(5,x, pglob,ylab,FALSE,FALSE)
##distPred(7,x, pglob,ylab,FALSE,FALSE) # source+neutral
distPred(6,x, pglob,ylab,FALSE,FALSE) # sink+neutral
distPred(7,x, pglob,ylab,FALSE,FALSE) # source only

Xlab(min(ylim)-.1, cex=0.7)
text(0,max(ylim),labels="B",cex=1.5)



x<-PE
ylab <- "SYNCHRONY"
ylim=c(0,.8)

#### ALL
distPlot(1,x,ylab,TRUE)
text(0,max(ylim),labels="C",cex=1.5)
#### Disp. 15%
#plot.new()
distPlot(2,x,ylab,TRUE,FALSE)
##distPlot(3,x,ylab,TRUE,FALSE) # source+neutral
distPlot(3,x,ylab,TRUE,FALSE) # sink+neutral
distPlot(4,x,ylab,TRUE,FALSE) # source only

#### Disp. 30%
#plot.new()
distPlot(5,x,ylab,TRUE,FALSE)
##distPlot(7,x,ylab,TRUE,FALSE) # source+neutral
distPlot(6,x,ylab,TRUE,FALSE) # sink+neutral
distPlot(7,x,ylab,TRUE,FALSE) # source only

#### PREDICTION
#pglob=0.1
abline(v=pglob,lty=2,col="grey")
ylab <- ""#paste0("Synchrony (pGlob=",pglob,")")
ylim=c(0.1,0.45)
distPred(1,x, pglob,ylab,TRUE)
#### Disp. 15%
#plot.new()
distPred(2,x, pglob,ylab,TRUE,FALSE)
#distPred(3,x, pglob,ylab,TRUE,FALSE) # source+neutral
distPred(3,x, pglob,ylab,TRUE,FALSE) # sink+neutral
distPred(4,x, pglob,ylab,TRUE,FALSE) # source only

#### Disp. 30%
#plot.new()
distPred(5,x, pglob,ylab,TRUE,FALSE)
#distPred(7,x, pglob,ylab,TRUE,FALSE) # source+neutral
distPred(6,x, pglob,ylab,TRUE,FALSE) # sink+neutral
distPred(7,x, pglob,ylab,TRUE,FALSE) # source only

Xlab(min(ylim)-.03, cex=0.7)
text(0,max(ylim),labels="D",cex=1.5)














lgd<-function(cex){
  legend("topright"
         #,c("fishing all","fishing sink+neutral","fishing source+neutral","fishing source only")
         #,c("fishing all","fishing sink+neutral","fishing source only")
         ,c("Exploitation (all pop.)","Conservation of source pop.","Exploitation of source pop.")
         ,col=c("black",col.type[c(1,3)])
         ,text.col=c("black",col.type[c(1,3)])
         #,pch=rep(16,5)
         ,lty=rep(1,4)
         ,bty="n"
         ,cex=cex)
}

distPlot <- function(nscn, X, ylab,plot=TRUE,relative=TRUE){
  
  xlab="Metapopulation exploitation rate"
  xlim=c(0,max(frates_vec))
  lwd=2
  
  
  EXPE <- SCN[,nscn]
  scn_id <- as.numeric(strsplit(as.character(EXPE[2]), "")[[1]])
  scenarioConnect = scn_id[1]
  scenarioFishing = scn_id[2]
  scenarioFishingRate = scn_id[3]
  
  if (scenarioFishing==0){ fishingStrat <- "Expl. All"}
  if (scenarioFishing==1){ fishingStrat <- "Expl. Source+neutral"}
  if (scenarioFishing==2){ fishingStrat <- "Expl. Sink+neutral"}
  if (scenarioFishing==3){ fishingStrat <- "Expl. Source only"}
  
  #ylab="# Parr 0+"
  #main=paste0("h=",h[scenarioConnect],"/",fishingStrat)#paste0("Neutral gene (", pops[pop],")")
  main=""
  if(plot==TRUE){
    #ylim=c(0,1.1)#c(min(sapply(X,min,na.rm=TRUE)),max(sapply(X,max,na.rm=TRUE)))#c(0,max(unlist(tmp), na.rm=TRUE))
    xlim=c(0, .3)
    plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
    axis(1,frates_vec,frates_vec)
    #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
   # lgd()
  }
  
  
  tmp=tmp2.5=tmp97.5=pGlob=NULL
  j=0
  for (i in EXPE){
    j=j+1
    if(is.na(i)) next;
    
    var<-X[[paste0(i)]]
    varTail <- apply(var, 2, function(x) mean(tail(x , n = 5),na.rm=TRUE))
    
    # tmp[j]<-quantile(varTail,probs=0.5,na.rm=TRUE)
    # tmp2.5[j]<-quantile(varTail,probs=0.025,na.rm=TRUE)
    # tmp97.5[j]<-quantile(varTail,probs=0.975,na.rm=TRUE)
    
    tmp[j]<-mean(varTail,na.rm=TRUE)
    tmp2.5[j]<-tmp[j]-se(varTail)
    tmp97.5[j]<-tmp[j]+se(varTail)
    
    pGlob[j] <- mean(gExpl[[paste0(i)]],na.rm=TRUE)
    
    if (relative){
      if(i==EXPE[1] & j==1) ref <- tmp[1]
    } else {ref <-1}
  }
  
  #colfunc <- colorRampPalette(c("darkgrey", "black"))
  if (nscn==1) mycol <- "lightgrey"
  
  if (nscn==2) mycol <- "darkgrey"
  #if (nscn==3) mycol <- col.type[2] # source+neutral
  if (nscn==3) mycol <- makeTransparent(col.type[1], 80) # sink+neutral col.type[1] # sink+neutral
  if (nscn==4) mycol <- makeTransparent(col.type[3], 80) # source only col.type[3] # source only
  
  if (nscn==5) mycol <- "black"
  #if (nscn==7) mycol <- makeTransparent(col.type[2], 80) # source+neutral
  if (nscn==6) mycol <- col.type[1] # sink+neutral
  if (nscn==7) mycol <- col.type[3] # source only
  #colfunc <- colorRampPalette(c(col.sinks[1], col.sinks[2]))}
  #  mycol <- col.sinks[1]} else { mycol <- col.sinks[2]}
  
  x=pGlob
  polygon(c(x,rev(x)),c(tmp2.5/ref,rev(tmp97.5/ref)),border=NA, col=makeTransparent(mycol,20))
  
  #scatter.smooth(tmp/ref ~ pGlob, span = 2/3, degree = 2,col=mycol)
  lines(loess.smooth(pGlob,tmp/ref,span = 2/3, degree = 1), col = mycol,lwd=3)
  
  #if(text){
  text(pGlob,tmp/ref, frates_vec, col=mycol,cex=.6)
  #}
  #points(pGlob,tmp/ref,col=mycol,type='p',lwd=1,pcDisp. 0%,cex=.6)
  
  abline(h=1,lty=2,col="grey")
  #abline(h=ref,lty=2,col="grey")
  #segments(pGlob,tmp2.5/ref,pGlob,tmp97.5/ref,col=col[1:length(frates_vec)])
  #points(pGlob,tmp/ref,col=col[1:length(frates_vec)],type='b',lwd=lwd,pch=20)
  
}
distPred <- function(nscn, X, pglob, ylab,plot=TRUE,relative=TRUE){
  
  xlab=""
  xlim=c(0,max(frates_vec))
  lwd=2
  
  
  EXPE <- SCN[,nscn]
  scn_id <- as.numeric(strsplit(as.character(EXPE[2]), "")[[1]])
  scenarioConnect = scn_id[1]
  scenarioFishing = scn_id[2]
  scenarioFishingRate = scn_id[3]
  
  if (scenarioFishing==0){ fishingStrat <- "Expl. All"}
  if (scenarioFishing==1){ fishingStrat <- "Expl. Source+neutral"}
  if (scenarioFishing==2){ fishingStrat <- "Expl. Sink+neutral"}
  if (scenarioFishing==3){ fishingStrat <- "Expl. Source only"}
  
  #ylab="# Parr 0+"
  #main=paste0("h=",h[scenarioConnect],"/",fishingStrat)#paste0("Neutral gene (", pops[pop],")")
  main=""
  if(plot==TRUE){
    #ylim=c(0,1.1)#c(min(sapply(X,min,na.rm=TRUE)),max(sapply(X,max,na.rm=TRUE)))#c(0,max(unlist(tmp), na.rm=TRUE))
    xlim=c(0, ncol(SCN)+0.5)
    plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
    xtext<-c("Control","All"
             #,"Source & neutral"
             ,"Sink & neutral","Source only",
             "All"
             #,"Source & neutral"
             ,"Sink & neutral","Source only"
    )
    # xtext<-c("Control","No conservation (fishing all)"
    #          #,"Source & neutral"
    #          ,"Conserv. Sources","Conserv. Sinks+Neutrals",
    #          "No conservation (fishing all)"
    #          #,"Source & neutral"
    #          ,"Conserv. Sources","Conserv. Sinks+Neutrals"
    # )
    axis(1,labels=F,at=1:length(xtext),las=1,cex=0.4)
    #axis(1,line = 2, labels="Management strategies",at=4,las=1,cex=0.6)
    # text(1:length(xtext), min(ylim) - 0.1*min(ylim), xtext, srt = 45, xpd = T,cex=0.7
    #      ,adj = c(1,0),pos=2
    #      , col=c("black","darkgrey",col.type[c(1,3)],"lightgrey",makeTransparent(col.type[c(1,3)], 95))
    # )
    # mtext(xtext,side=1,line=1,outer=FALSE,at=1:length(xtext),adj=1,las=3,cex=.5
    #       ,col=c("black","darkgrey",col.type[c(1,3)],"lightgrey",makeTransparent(col.type[c(1,3)], 95))
    # )
    
    # xtext2<-c("Control","Exploitation"
    #           #,"Source & neutral"
    #           ,"Conservation","Exploitation",
    #           "Exploitation"
    #           #,"Source & neutral"
    #           ,"Conservation","Exploitation"
    # )
    # mtext(xtext2,
    #       side=1,line=.5,outer=FALSE,at=1:length(xtext)-0.3,adj=0,las=1,cex=.25
    #       ,col=c("black","black",col.type[c(1,3)],"black",col.type[c(1,3)])
    # )
    # 
    # xtext2<-c("","all pop."
    #           #,"Source & neutral"
    #           ,"Source pop.","Source pop.",
    #           "all pop."
    #           #,"Source & neutral"
    #           ,"Source pop.","Source pop."
    # )
    # mtext(xtext2,
    #       side=1,line=1,outer=FALSE,at=1:length(xtext)-0.3,adj=0,las=1,cex=.25
    #       ,col=c("black","black",col.type[c(1,3)],"black",col.type[c(1,3)])
    # )
    
    
    #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
    #lgd()
    abline(v=1.5,lty=3)
    abline(v=4.5,lty=3)
    
    #text(c(0.5,3, 6),rep(max(ylim)*0.95,3),c("Disp. 0%","Disp. 15%","Disp. 30%"),cex=1,col="darkgrey")
    #text(c(0.5,3, 6),rep(max(ylim)*0.95,3),c("Dispersal: 0%","Dispersal: 15%","Dispersal: 30%"),cex=1,col="darkgrey")
  }
  
  
  tmp=tmp2.5=tmp97.5=pGlob=NULL
  j=0
  for (i in EXPE){
    j=j+1
    if(is.na(i)) next;
    
    var<-X[[paste0(i)]]
    varTail <- apply(var, 2, function(x) mean(tail(x , n = 5),na.rm=TRUE))
    
    # tmp[j]<-quantile(varTail,probs=0.5,na.rm=TRUE)
    # tmp2.5[j]<-quantile(varTail,probs=0.025,na.rm=TRUE)
    # tmp97.5[j]<-quantile(varTail,probs=0.975,na.rm=TRUE)
    
    tmp[j]<-mean(varTail,na.rm=TRUE)
    tmp2.5[j]<-tmp[j]-se(varTail)
    tmp97.5[j]<-tmp[j]+se(varTail)
    
    pGlob[j] <- mean(gExpl[[paste0(i)]],na.rm=TRUE)
    
    if (relative){
      if(i==EXPE[1] & j==1) ref <- tmp[1]
    } else {ref <-1}
    
  }
  
  #colfunc <- colorRampPalette(c("darkgrey", "black"))
  if (nscn==1) mycol <- "lightgrey"
  
  if (nscn==2) mycol <- "darkgrey"
  #if (nscn==3) mycol <- col.type[2] # source+neutral
  if (nscn==3) mycol <- col.type[1] # sink+neutral
  if (nscn==4) mycol <- col.type[3] # source only
  
  if (nscn==5) mycol <- "black"
  #if (nscn==7) mycol <- makeTransparent(col.type[2], 80) # source+neutral
  if (nscn==6) mycol <- col.type[1] #makeTransparent(col.type[1], 80) # sink+neutral
  if (nscn==7) mycol <- col.type[3] #makeTransparent(col.type[3], 80) # source only
  
  #colfunc <- colorRampPalette(c(col.sinks[1], col.sinks[2]))}
  #  mycol <- col.sinks[1]} else { mycol <- col.sinks[2]}
  
  
  #scatter.smooth(tmp ~ pGlob, span = 2/3, degree = 2,col=mycol)
  #lines(loess.smooth(pGlob,tmp,span = 2/3, degree = 1), col = mycol,lwd=3)
  lo <- loess(tmp/ref ~ pGlob)
  pred <- predict(lo, data.frame(pGlob = pglob), se = TRUE)
  #predict(lo, data.frame(pGlob = seq(5, 30, 1)), se = TRUE)
  
  #if(text){
  #text(pGlob,tmp,frates_vec, col=mycol,cex=.6)
  #}
  #points(pGlob,tmp,col=mycol,type='p',lwd=1,pcDisp. 0%,cex=.6)
  
  #abline(h=1,lty=2,col="grey")
  #abline(h=ref,lty=2,col="grey")
  segments(nscn,pred$fit-pred$se.fit,nscn,pred$fit+pred$se.fit,col=mycol)
  points(nscn,pred$fit,col=mycol,pch=20)
  
}

#par(mfrow = c(1,2), mar=c(4, 4, 1, 1) + 0.1)#it goes c(bottom, left, top, right)
par(mfrow = c(2,2), mar=c(4, 4, 1, 1) + 0.1)#it goes c(bottom, left, top, right)


x<-nReturns
ylab <- "Relative Abundance (PFA)"
ylim=c(0.5,1.1)
pglob=0.15
#### ALL
distPlot(1,x,ylab,TRUE)
# distPlot(2,x,ylab,FALSE)
# distPlot(5,x,ylab,FALSE)
#abline(v=0.1,lty=2)

#### Disp. 15%
#plot.new()
distPlot(2,x,ylab,FALSE)
#distPlot(3,x,ylab,FALSE) # source
distPlot(3,x,ylab,FALSE) # sink
distPlot(4,x,ylab,FALSE) # source only

#### Disp. 30%
#plot.new()
distPlot(5,x,ylab, FALSE)
#distPlot(7,x,ylab,FALSE) # source
distPlot(6,x,ylab,FALSE) # sink
distPlot(7,x,ylab,FALSE) # source only

lgd(0.5)
text(0,max(ylim),labels="A",cex=1.5)
#legend("topleft",legend="A",bty="n",xjust=0,cex=1.5)

#### PREDICTION
#pglob=0.1
abline(v=pglob,lty=2,col="grey")
ylab <- ""#paste0("Relative returns (pGlob=",pglob,")")
ylim=c(0.5,1.1)
distPred(1,x, pglob,ylab)
#### Disp. 15%
#plot.new()
distPred(2,x, pglob,ylab,FALSE)
#distPred(3,x, pglob,ylab,FALSE) # source+neutral
distPred(3,x, pglob,ylab,FALSE) # sink+neutral
distPred(4,x, pglob,ylab,FALSE) # source only

#### Disp. 30%
#plot.new()
distPred(5,x, pglob,ylab,FALSE)
#distPred(7,x, pglob,ylab,FALSE) # source+neutral
distPred(6,x, pglob,ylab,FALSE) # sink+neutral
distPred(7,x, pglob,ylab,FALSE) # source only

Xlab(min(ylim)-.05, cex=0.7)
text(0,max(ylim),labels="B",cex=1.5)
text(c(0.55,3, 6),rep(max(ylim)*0.95,3),c("Dispersal: 0%","Dispersal: 15%","Dispersal: 30%"),cex=0.5,col="darkgrey")

plot.new()
plot.new()




par(mfrow = c(2,2), mar=c(4, 4, 1, 1) + 0.1)#it goes c(bottom, left, top, right)

#x<-lapply(rMSW,log)
x<-rMSW
ylab <- "Ratio MSW/1SW"
#ylim=c(-1.5,-1.3)
ylim=c(0.21,0.26)
#### ALL
distPlot(1,x,ylab,TRUE,FALSE)
# distPlot(2,x,ylab,FALSE,FALSE)
# distPlot(5,x,ylab,FALSE,FALSE)
#abline(v=0.1,lty=2)
for (sc in 2:7 ){
  distPlot(sc,x,ylab,FALSE,FALSE)
}
lgd(0.5)
text(0,max(ylim),labels="A",cex=1.5)

#### PREDICTION
#pglob=0.15
abline(v=pglob,lty=2,col="grey")
ylab <- ""#paste0("Relative ratio MSW/1SW (pGlob=",pglob,")")
#ylim=c(0.21,0.26)
distPred(1,x, pglob,ylab,TRUE,FALSE)
for (sc in 2:7 ){
  distPred(sc,x,pglob,ylab,FALSE,FALSE)
}
text(0,max(ylim),labels="B",cex=1.5)
text(c(0.5,3, 6),rep(max(ylim)*0.98,3),c("Dispersal: 0%","Dispersal: 15%","Dispersal: 30%"),cex=.5,col="darkgrey")
Xlab(min(ylim)-.005, cex=0.7)




#par(mfrow=c(1,2))
rParrMature<-list()
for (l in names(nParr)){
  rParrMature[[paste0(l)]]<-nParrMature[[l]]/nParr[[l]]
}

x<-rParrMature
ylab <- "Proportion mature parr 0+"
ylim=c(0.02,0.04)
#### ALL
distPlot(1,x,ylab,TRUE,FALSE)
for (sc in 2:7 ){
  distPlot(sc,x,ylab,FALSE,FALSE)
}
lgd(0.5)
text(0,max(ylim),labels="A",cex=1.5)

#### PREDICTION
#pglob=0.1
abline(v=pglob,lty=2,col="grey")
ylab <- ""# paste0("Relative ratio mature parr 0+ (pGlob=",pglob,")")
ylim=c(0.025,0.035)
distPred(1,x, pglob,ylab,TRUE,FALSE)
for (sc in 2:7 ){
  distPred(sc,x,pglob,ylab,FALSE,FALSE)
}
text(0,max(ylim),labels="B",cex=1.5)
text(c(0.5,3, 6),rep(max(ylim)*0.98,3),c("Dispersal: 0%","Dispersal: 15%","Dispersal: 30%"),cex=.5,col="darkgrey")
Xlab(min(ylim)-.001, cex=0.7)







par(mfrow = c(4,2), mar=c(4, 4, 1, 1) + 0.1)#it goes c(bottom, left, top, right)

x<-LFParr
ylab <- "Relative size Parr 0+"
ylim=c(1,1.1)
#### ALL
distPlot(1,x,ylab)
for (sc in 2:7 ){
  distPlot(sc,x,ylab,FALSE,TRUE)
}
lgd(0.5)
text(0,max(ylim),pos=1,labels="A",cex=1.5)

#### PREDICTION
#pglob=0.15
abline(v=pglob,lty=2,col="grey")
ylab <-""# paste0("Relative size Parr 0+ (pGlob=",pglob,")")
ylim=c(1.03,1.05)
distPred(1,x, pglob,ylab)
for (sc in 2:7 ){
  distPred(sc,x,pglob,ylab,FALSE,TRUE)
}
text(0,max(ylim),pos=1,labels="B",cex=1.5)
text(c(0.5,3, 6),rep(max(ylim)*0.995,3),c("Dispersal: 0%","Dispersal: 15%","Dispersal: 30%"),cex=.5,col="darkgrey")
Xlab(min(ylim)-.004, cex=0.7)




#### Disp. 15%
x<-LFSmolt
ylab <- "Relative size Smolt"
ylim=c(1,1.1)
distPlot(1,x,ylab)
for (sc in 2:7 ){
  distPlot(sc,x,ylab,FALSE,TRUE)
}
lgd(0.5)
text(0,max(ylim),pos=1,labels="C",cex=1.5)

#### PREDICTION
#pglob=0.15
abline(v=pglob,lty=2,col="grey")
ylab <-""# paste0("Relative size Smolt 0+ (pGlob=",pglob,")")
ylim=c(1.015,1.025)
distPred(1,x, pglob,ylab)
for (sc in 2:7 ){
  distPred(sc,x,pglob,ylab,FALSE,TRUE)
}
text(0,max(ylim),pos=1,labels="D",cex=1.5)
text(c(0.5,3, 6),rep(max(ylim)*0.995,3),c("Dispersal: 0%","Dispersal: 15%","Dispersal: 30%"),cex=.5,col="darkgrey")
Xlab(min(ylim)-.002, cex=0.7)




#### Disp. 15%
x<-LF1sw
ylab <- "Relative size 1SW"
ylim=c(1,1.01)
distPlot(1,x,ylab)
for (sc in 2:7 ){
  distPlot(sc,x,ylab,FALSE,TRUE)
}
lgd(0.5)
text(0,max(ylim),pos=1,labels="E",cex=1.5)
#### PREDICTION
pglob=0.15
abline(v=pglob,lty=2,col="grey")
ylab <-""# paste0("Relative size 1SW (pGlob=",pglob,")")
ylim=c(1.001,1.003)
distPred(1,x, pglob,ylab)
for (sc in 2:7 ){
  distPred(sc,x,pglob,ylab,FALSE,TRUE)
}
text(0,max(ylim),pos=1,labels="F",cex=1.5)
text(c(0.5,3, 6),rep(max(ylim)*0.995,3),c("Dispersal: 0%","Dispersal: 15%","Dispersal: 30%"),cex=.5,col="darkgrey")
Xlab(min(ylim)-.0004, cex=0.7)



#### Disp. 15%
x<-LFmsw
ylab <- "Relative size MSW"
ylim=c(0.995,1.01)
distPlot(1,x,ylab)
for (sc in 2:7 ){
  distPlot(sc,x,ylab,FALSE,TRUE)
}
lgd(0.5)
text(0,max(ylim),pos=1,labels="G",cex=1.5)
#### PREDICTION
#pglob=0.15
abline(v=pglob,lty=2,col="grey")
ylab <-""# paste0("Relative size MSW (pGlob=",pglob,")")
ylim=c(0.998,1.001)
distPred(1,x, pglob,ylab)
for (sc in 2:7 ){
  distPred(sc,x,pglob,ylab,FALSE,TRUE)
}
text(0,max(ylim),pos=1,labels="H",cex=1.5)
text(c(0.5,3, 6),rep(max(ylim)*0.995,3),c("Dispersal: 0%","Dispersal: 15%","Dispersal: 30%"),cex=.5,col="darkgrey")
Xlab(min(ylim)-.0005, cex=0.7)







## Genotypes
# gNEUTRAL[[paste0(iEXPE)]] <- gNeutral
# gFMID1[[paste0(iEXPE)]] <- gFmid1
# gFMID2[[paste0(iEXPE)]] <- gFmid2
# gFMID3[[paste0(iEXPE)]] <- gFmid3
# gFMID4[[paste0(iEXPE)]] <- gFmid4
# gGROWTH[[paste0(iEXPE)]] <- gG

par(mfrow = c(2,2), mar=c(4, 4, 1, 1) + 0.1)#it goes c(bottom, left, top, right)


x<-gNEUTRAL
ylab <- "Neutral gene"
ylim=c(.95,1.1)
#### ALL
distPlot(1,x,ylab)
for (sc in 2:7 ){
  distPlot(sc,x,ylab,FALSE,TRUE)
}
lgd(0.5)
text(0,max(ylim),pos=1,labels="A",cex=1.5)
#### PREDICTION
pglob=0.15
abline(v=pglob,lty=2,col="grey")
ylab <-""# paste0("Neutral gene (pGlob=",pglob,")")
ylim=c(0.99,1.01)
distPred(1,x, pglob,ylab)
for (sc in 2:7 ){
  distPred(sc,x,pglob,ylab,FALSE,TRUE)
}
text(0,max(ylim),pos=1,labels="B",cex=1.5)
text(c(0.5,3, 6),rep(max(ylim)*0.995,3),c("Dispersal: 0%","Dispersal: 15%","Dispersal: 30%"),cex=.5,col="darkgrey")
Xlab(min(ylim)-.002, cex=0.7)



x<-gGROWTH
ylab <- "Growth potential (anadromous)"
ylim=c(0,1.1)
#### ALL
distPlot(1,x,ylab)
for (sc in 2:7 ){
  distPlot(sc,x,ylab,FALSE,TRUE)
}
lgd(0.5)
text(0,max(ylim),pos=1,labels="C",cex=1.5)
#### PREDICTION
#pglob=0.15
abline(v=pglob,lty=2,col="grey")
ylab <-""# paste0("Growth potential (anadromous) (pGlob=",pglob,")")
ylim=c(0.5,1.01)
distPred(1,x, pglob,ylab)
for (sc in 2:7 ){
  distPred(sc,x,pglob,ylab,FALSE,TRUE)
}
text(0,max(ylim),pos=1,labels="D",cex=1.5)
text(c(0.5,3, 6),rep(max(ylim)*0.95,3),c("Dispersal: 0%","Dispersal: 15%","Dispersal: 30%"),cex=.5,col="darkgrey")
Xlab(min(ylim)-.05, cex=0.7)






par(mfrow = c(3,2), mar=c(4, 4, 1, 1) + 0.1)#it goes c(bottom, left, top, right)


x<-gFMID1
ylab <- "Precocious mat. threshold (male)"
ylim=c(.95,1.01)
#### ALL
distPlot(1,x,ylab)
for (sc in 2:7 ){
  distPlot(sc,x,ylab,FALSE,TRUE)
}
lgd(0.5)
text(0,max(ylim),pos=1,labels="A",cex=2)
#### PREDICTION
#pglob=0.15
abline(v=pglob,lty=2,col="grey")
ylab <-""# paste0("Anadromous maturation threshold (male) (pGlob=",pglob,")")
ylim=c(0.95,1.01)
distPred(1,x, pglob,ylab)
for (sc in 2:7 ){
  distPred(sc,x,pglob,ylab,FALSE,TRUE)
}
text(0,max(ylim),pos=1,labels="B",cex=1.5)
text(c(0.5,3, 6),rep(max(ylim)*0.95,3),c("Dispersal: 0%","Dispersal: 15%","Dispersal: 30%"),cex=.5,col="darkgrey")
Xlab(min(ylim)-.01, cex=0.7)

x<-gFMID3
ylab <- "Anadromous mat. threshold (male)"
ylim=c(.95,1.02)
#### ALL
distPlot(1,x,ylab)
for (sc in 2:7 ){
  distPlot(sc,x,ylab,FALSE,TRUE)
}
lgd(0.5)
text(0,max(ylim),pos=1,labels="C",cex=2)
#### PREDICTION
#pglob=0.15
abline(v=pglob,lty=2,col="grey")
ylab <-""# paste0("Anadromous maturation threshold (female) (pGlob=",pglob,")")
ylim=c(0.95,1.02)
distPred(1,x, pglob,ylab)
for (sc in 2:7 ){
  distPred(sc,x,pglob,ylab,FALSE,TRUE)
}
text(0,max(ylim),pos=1,labels="D",cex=1.5)
text(c(0.5,3, 6),rep(max(ylim)*0.95,3),c("Dispersal: 0%","Dispersal: 15%","Dispersal: 30%"),cex=.5,col="darkgrey")
Xlab(min(ylim)-.01, cex=0.7)


x<-gFMID4
ylab <- "Anadromous mat. threshold (female)"
ylim=c(.95,1.02)
#### ALL
distPlot(1,x,ylab)
for (sc in 2:7 ){
  distPlot(sc,x,ylab,FALSE,TRUE)
}
lgd(0.5)
text(0,max(ylim),pos=1,labels="E",cex=2)
#### PREDICTION
#pglob=0.15
abline(v=pglob,lty=2,col="grey")
ylab <-""# paste0("Anadromous maturation threshold (female) (pGlob=",pglob,")")
ylim=c(0.98,1.02)
distPred(1,x, pglob,ylab)
for (sc in 2:7 ){
  distPred(sc,x,pglob,ylab,FALSE,TRUE)
}
text(0,max(ylim),pos=1,labels="F",cex=1.5)
text(c(0.5,3, 6),rep(max(ylim)*0.95,3),c("Dispersal: 0%","Dispersal: 15%","Dispersal: 30%"),cex=.5,col="darkgrey")
Xlab(min(ylim)-.005, cex=0.7)



dev.off()

























#### METAPOPULATION DYNAMICS ####
#-------------------------------#

lgd<-function() {legend(
  "topleft",
  #x = 0, y = max(ylim),
  legend = paste0(frates_vec),
  fill = col,
  border = NA,
  y.intersp = 0.5,
  bty="n",
  cex = .8, text.font = .1, title="Exploitation rate")
}


dynPlot <- function(VAR,main,ylab,ylim=NULL,relative=TRUE){
  #main=paste0("h=",h[scenarioConnect],"/ # Parr 0+")
  #ylab="Abundance (PFA)" # Pre-Fisheries Abundance
  if(is.null(ylim)){
    ylim=c(min(sapply(VAR,min,na.rm=TRUE)),max(sapply(VAR,max,na.rm=TRUE)))#c(0,max(unlist(tmp), na.rm=TRUE))
  }
  plot(NULL,xlim=c(1,nInit+nYear),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
  rect(0,min(ylim),10,max(ylim), col= "lightgrey",border ="lightgrey")
  rect(45,min(ylim),50,max(ylim), col= NULL,border ="darkgrey",lty=2)
  
  lgd()
  j=0
  for (i in SCN[,strat]){
    j=j+1
    # tmp<-apply(VAR[[paste0(i)]],1,median,na.rm=TRUE)
    # tmp2.5<-apply(VAR[[paste0(i)]],1,quantile,probs=0.025,na.rm=TRUE)
    # tmp97.5<-apply(VAR[[paste0(i)]],1,quantile,probs=0.975,na.rm=TRUE)
    
    tmp<-apply(VAR[[paste0(i)]],1,mean,na.rm=TRUE)
    tmp2.5<-tmp-apply(VAR[[paste0(i)]],1,se)
    tmp97.5<-tmp+apply(VAR[[paste0(i)]],1,se)
    
    if(i==SCN[1,strat] & j==1) ref <- tmp[11]
    if (relative){
      x=2:50
      polygon(c(x,rev(x)),c(tmp2.5[x]/ref,rev(tmp97.5[x]/ref)),border=NA, col=makeTransparent(col[j],20))
      points(x,tmp[x]/ref,col=col[j],type='l',lwd=lwd)
    } else{
      ref=1
      x=2:50
      polygon(c(x,rev(x)),c(tmp2.5[x]/ref,rev(tmp97.5[x]/ref)),border=NA, col=makeTransparent(col[j],20))
      points(x,tmp[x]/ref,col=col[j],type='l',lwd=lwd)
      
    }
    
  }
}


# jpeg(filename = "results/Figure_explrates.jpeg",
#      width = 580, height = 780, 
#      quality = 75)

#xlab="Time (years)"
xlab=""
lwd=2

### CONTROL
strat=1
EXPE <- SCN[,strat]
scn_id <- as.numeric(strsplit(as.character(EXPE[1]), "")[[1]])
scenarioConnect = scn_id[1]
scenarioFishing = scn_id[2]
scenarioFishingRate = scn_id[3]

par(mfrow=c(4,2), mai = c(.5, .7, .1, .7))
#dynPlot(nReturns, main=paste0("Philopatry strict - No dispersal"), ylab="Numbers of returns",ylim=c(0,1.5))
#dynPlot(rMSW, main=paste0("h=",h[1],"/ # MSW"), ylab="Ratio MSW/1SW",ylim=c(0,0.4))

#par(mfrow=c(2,2))
#dynPlot(nParr, main=paste0("h=",h[1],"/ # Parr 0+"), ylab="Abundance (mean±se)",ylim=c(0,1.5))
#dynPlot(nSmolt, main=paste0("h=",h[1],"/ # Smolt 0+"), ylab="Abundance (mean±se)",ylim=c(0,1.5))
#dynPlot(nMSW, main=paste0("h=",h[1],"/ # MSW"), ylab="Abundance (mean±se)")
#dynPlot(n1SW, main=paste0("h=",h[1],"/ # 1SW"), ylab="Abundance (mean±se)")
#dynPlot(nReturns, main=paste0("h=",h[1],"/ # Returns"), ylab="Abundance (mean±se)",ylim=c(0,1.5))
dynPlot(rMSW, main=paste0("Philopatry strict - No dispersal"), ylab="Ratio of MSW/1SW",ylim=c(0,0.4),relative=FALSE)

# rParrMature<-list()
# for (l in names(nParr)){
#   rParrMature[[paste0(l)]]<-nParrMature[[l]]/nParr[[l]]
# }
ylab <- "Proportion mature parr 0+"
dynPlot(rParrMature, main="", ylab,ylim=c(0,0.05),relative=FALSE)


#par(mfrow=c(2,2))
#dynPlot(LFParr, main=paste0("h=",h[1],"/ Fork Length Parr 0+"), ylab="Fork Length (mm)",ylim=c(.9,1.1))
#dynPlot(LFSmolt, main=paste0("h=",h[1],"/ Fork Length Smolt 0+"), ylab="Fork Length (mm)",ylim=c(0.9,1.1))
dynPlot(LF1sw, main="", ylab="Body size of 1SW",ylim=c(0.98,1.02))
dynPlot(LFmsw, main="", ylab="Body size of MSW",ylim=c(0.98,1.02))

## GENETIC
#par(mfrow=c(3,2))
#dynPlot(gNEUTRAL, main=paste0("h=",h[1],"/ gNEUTRAL"), ylab="Neutral gene",ylim=c(0.45,0.55),relative=FALSE)
dynPlot(gGROWTH, main="", ylab="Growth potential (anadromous)",ylim=c(-0.01,.025),relative=FALSE)

dynPlot(gFMID1, main="", ylab="Precocious male mat. Threshold",ylim=c(.9,1.15))
#dynPlot(gFMID2, main=paste0("h=",h[1],"/ gFMID2"), ylab="Precocious female maturation Threshold",ylim=c(.9,1.15))
xlab="Time (years)"
dynPlot(gFMID3, main="", ylab="Anadromous male maturation Threshold",ylim=c(.9,1.25))
dynPlot(gFMID4, main="", ylab="Anadromous female maturation Threshold",ylim=c(.9,1.25))


#dev.off()


# for (strat in 1:ncol(SCN)){
# #for (strat in 1:3){
# #strat=2
# 
#   EXPE <- SCN[,strat]
#   
#   scn_id <- as.numeric(strsplit(as.character(EXPE[1]), "")[[1]])
#   scenarioConnect = scn_id[1]
#   scenarioFishing = scn_id[2]
#   scenarioFishingRate = scn_id[3]
#   
#   
#   
# dynPlot(nReturns, main=paste0("h=",h[scenarioConnect],"/ # Returns"), ylab="Abundance (PFA)")
# 
#   
# par(mfrow=c(2,2))
# 
# dynPlot(nParr, main=paste0("h=",h[scenarioConnect],"/ # Parr 0+"), ylab="Abundance (PFA)")
# dynPlot(LFParr, main=paste0("h=",h[scenarioConnect],"/ Fork Length Parr 0+"), ylab="Fork Length (mm)")
# 
# dynPlot(nSmolt, main=paste0("h=",h[scenarioConnect],"/ # Smolt 0+"), ylab="Abundance (PFA)")
# dynPlot(LFSmolt, main=paste0("h=",h[scenarioConnect],"/ Fork Length Smolt 0+"), ylab="Fork Length (mm)")
# 
# par(mfrow=c(3,2))
# dynPlot(nMSW, main=paste0("h=",h[scenarioConnect],"/ # MSW"), ylab="Abundance (PFA)")
# dynPlot(LFmsw, main=paste0("h=",h[scenarioConnect],"/ Fork Length MSW"), ylab="Fork Length (mm)")
# dynPlot(n1SW, main=paste0("h=",h[scenarioConnect],"/ # 1SW"), ylab="Abundance (PFA)")
# dynPlot(LF1sw, main=paste0("h=",h[scenarioConnect],"/ Fork Length 1SW"), ylab="Fork Length (mm)")
# dynPlot(rMSW, main=paste0("h=",h[scenarioConnect],"/ Ratio MSW/1SW"), ylab="Ratio")
# plot.new()
# 
# ## GENETIC
# par(mfrow=c(3,2))
# dynPlot(gNEUTRAL, main=paste0("h=",h[scenarioConnect],"/ gNEUTRAL"), ylab="Neutral gene")
# dynPlot(gFMID1, main=paste0("h=",h[scenarioConnect],"/ gFMID1"), ylab="Precocious male maturation Threshold")
# #dynPlot(gFMID2, main=paste0("h=",h[scenarioConnect],"/ gFMID2"), ylab="Precocious female maturation Threshold")
# dynPlot(gFMID3, main=paste0("h=",h[scenarioConnect],"/ gFMID3"), ylab="Anadromous male maturation Threshold")
# dynPlot(gFMID4, main=paste0("h=",h[scenarioConnect],"/ gFMID4"), ylab="Anadromous female maturation Threshold")
# dynPlot(gGROWTH, main=paste0("h=",h[scenarioConnect],"/ gGROWTH"), ylab="Growth potential (anadromous)")
# 
# 
# ### Coefficient variation CV
# # par(mfrow=c(3,2))
# # dynPlot(CVgNEUTRAL, main=paste0("h=",h[scenarioConnect],"/ CVgNEUTRAL"), ylab="Neutral gene")
# # dynPlot(CVgFMID1, main=paste0("h=",h[scenarioConnect],"/ CVgFMID1"), ylab="Precocious male maturation Threshold")
# # #dynPlot(CVgFMID2, main=paste0("h=",h[scenarioConnect],"/ CVgFMID2"), ylab="Precocious female maturation Threshold")
# # dynPlot(CVgFMID3, main=paste0("h=",h[scenarioConnect],"/ CVgFMID3"), ylab="Anadromous male maturation Threshold")
# # dynPlot(CVgFMID4, main=paste0("h=",h[scenarioConnect],"/ CVgFMID4"), ylab="Anadromous female maturation Threshold")
# # dynPlot(CVgGROWTH, main=paste0("h=",h[scenarioConnect],"/ CVgGROWTH"), ylab="Growth potential (anadromous)")
# 
# 
# 
# } # end loop strat





#dev.off()





#### EXTINCTION RISK ####
Extinction<-list()
for (l in names(PE)){
  extinction=NULL
  for (sim in 1:nSIMUL){ 
    tmp <- PE[[paste0(l)]][[sim]]$Returns
    #tmp3 <- (apply(tmp,1,function(x) length(which(x<5)))/npop)
    tmp[tmp==0]<-NA
    tmp3 <- (apply(tmp,1,function(x) length(which(is.na(x))))/npop)
    #tmp3 <- 1-(apply(tmp,1,function(x) length(which(!is.na(x))))/npop)
    extinction <- cbind(extinction, tmp3)
  }
  Extinction[[paste0(l)]]<-extinction
}

par(mfrow = c(1,1))#it goes c(bottom, left, top, right)
distPlot <- function(nscn, X, ylab,plot=TRUE,relative=FALSE){
  
  xlab="Metapopulation exploitation rate"
  xlim=c(0,max(frates_vec))
  lwd=2
  
  
  EXPE <- SCN[,nscn]
  scn_id <- as.numeric(strsplit(as.character(EXPE[2]), "")[[1]])
  scenarioConnect = scn_id[1]
  scenarioFishing = scn_id[2]
  scenarioFishingRate = scn_id[3]
  
  if (scenarioFishing==0){ fishingStrat <- "Expl. All"}
  if (scenarioFishing==1){ fishingStrat <- "Expl. Source+neutral"}
  if (scenarioFishing==2){ fishingStrat <- "Expl. Sink+neutral"}
  if (scenarioFishing==3){ fishingStrat <- "Expl. Source only"}
  
  #ylab="# Parr 0+"
  #main=paste0("h=",h[scenarioConnect],"/",fishingStrat)#paste0("Neutral gene (", pops[pop],")")
  main=""
  if(plot==TRUE){
    #ylim=c(0,1.1)#c(min(sapply(X,min,na.rm=TRUE)),max(sapply(X,max,na.rm=TRUE)))#c(0,max(unlist(tmp), na.rm=TRUE))
    xlim=c(0, .3)
    plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
    axis(1,frates_vec,frates_vec)
    #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
    lgd()
  }
  
  
  tmp=tmp2.5=tmp97.5=pGlob=NULL
  j=0
  for (i in EXPE){
    j=j+1
    if(is.na(i)) next;
    
    var<-X[[paste0(i)]]
    varTail <- apply(var, 2, function(x) mean(tail(x , n = 5),na.rm=TRUE))
    
    # tmp[j]<-quantile(varTail,probs=0.5,na.rm=TRUE)
    # tmp2.5[j]<-quantile(varTail,probs=0.025,na.rm=TRUE)
    # tmp97.5[j]<-quantile(varTail,probs=0.975,na.rm=TRUE)
    
    tmp[j]<-mean(varTail,na.rm=TRUE)
    tmp2.5[j]<-tmp[j]-se(varTail)
    tmp97.5[j]<-tmp[j]+se(varTail)
    
    pGlob[j] <- mean(gExpl[[paste0(i)]],na.rm=TRUE)
    
    if (relative){
      if(i==EXPE[1] & j==1) ref <- tmp[1]
    } else {ref <-1}
  }
  
  #colfunc <- colorRampPalette(c("darkgrey", "black"))
  if (nscn==1) mycol <- "black"
  
  if (nscn==2) mycol <- "darkgrey"
  #if (nscn==3) mycol <- col.type[2] # source+neutral
  if (nscn==3) mycol <- col.type[1] # sink+neutral
  if (nscn==4) mycol <- col.type[3] # source only
  
  if (nscn==5) mycol <- "lightgrey"
  #if (nscn==7) mycol <- makeTransparent(col.type[2], 80) # source+neutral
  if (nscn==6) mycol <- makeTransparent(col.type[1], 80) # sink+neutral
  if (nscn==7) mycol <- makeTransparent(col.type[3], 80) # source only
  #colfunc <- colorRampPalette(c(col.sinks[1], col.sinks[2]))}
  #  mycol <- col.sinks[1]} else { mycol <- col.sinks[2]}
  
  x=pGlob
  #polygon(c(x,rev(x)),c(tmp2.5/ref,rev(tmp97.5/ref)),border=NA, col=makeTransparent(mycol,20))
  
  #scatter.smooth(tmp/ref ~ pGlob, span = 2/3, degree = 2,col=mycol)
  #lines(loess.smooth(pGlob,tmp/ref,span = 2/3, degree = 1), col = mycol,lwd=3)
  
  #if(text){
  text(pGlob,tmp/ref, frates_vec, col=mycol,cex=.6)
  #}
  #points(pGlob,tmp/ref,col=mycol,type='p',lwd=1,pcDisp. 0%,cex=.6)
  
  abline(h=1,lty=2,col="grey")
  #abline(h=ref,lty=2,col="grey")
  #segments(pGlob,tmp2.5/ref,pGlob,tmp97.5/ref,col=col[1:length(frates_vec)])
  #points(pGlob,tmp/ref,col=col[1:length(frates_vec)],type='b',lwd=lwd,pch=20)
  
}
x<-Extinction
ylab <- "Extinction risk (nReturns=0)"
ylim=c(0,.2)
#### ALL
distPlot(1,x,ylab)
for (sc in 2:7 ){
  distPlot(sc,x,ylab,FALSE,FALSE)
}







id.all<-1:npop
id.sinks <- which(dat$Type=="sink" | dat$Type=="neutral")
id.sources <- which(dat$Type=="source")

Immigrants.all<-Immigrants.sinks<-Immigrants.sources<-list()
for (scn in names(Mig)){
  
  mig=capt=NULL
  for (sim in 1:nSIMUL){ # loop over simulations
    tmp <- rowSums(Mig[[scn]][[sim]]$Immigrants[,id.all])
    mig<-cbind(mig,tmp)
    
    tmp <- rowSums(Mig[[scn]][[sim]]$Captured[,id.all])
    capt<-cbind(capt,tmp)
  }
  mig=mig+1 # to avoid 0
  capt=ceiling(capt+1) # to avoid 0
  ratio.all <- mig/capt
  Immigrants.all[[scn]] <-ratio.all
  #captured.sinks[[scn]] <-capt
  #Immigrants.sinks[[scn]] <-mig
  
  
  mig=capt=NULL
  for (sim in 1:nSIMUL){ # loop over simulations
    tmp <- rowSums(Mig[[scn]][[sim]]$Immigrants[,id.sinks])
    mig<-cbind(mig,tmp)
    
    tmp <- rowSums(Mig[[scn]][[sim]]$Captured[,id.sinks])
    capt<-cbind(capt,tmp)
  }
  mig=mig+1 # to avoid 0
  capt=ceiling(capt+1) # to avoid 0
  ratio.sinks <- mig/capt
  Immigrants.sinks[[scn]] <-ratio.sinks
  #captured.sinks[[scn]] <-capt
  #Immigrants.sinks[[scn]] <-mig
  
  
  mig=capt=NULL
  for (sim in 1:nSIMUL){ # loop over simulations
    tmp <- rowSums(Mig[[scn]][[sim]]$Immigrants[,id.sources])
    mig<-cbind(mig,tmp)
    
    tmp <- rowSums(Mig[[scn]][[sim]]$Captured[,id.sources])
    capt<-cbind(capt,tmp)
  }
  mig=mig+1 # to avoid 0
  capt=ceiling(capt+1) # to avoid 0
  ratio.sources <- mig/capt
  Immigrants.sources[[scn]] <-ratio.sources
}


#### Ratio Immigrants / Exploitation
# lgd<-function(){
#   legend("topright"
#          #,c("fishing all","fishing sink+neutral","fishing source+neutral","fishing source only")
#          ,c("fishing all","fishing sink+neutral","fishing source only")
#          ,col=c("black",col.type[c(1,3)])
#          ,text.col=c("black",col.type[c(1,3)])
#          #,pch=rep(16,5)
#          ,lty=rep(1,4)
#          ,bty="n"
#          ,cex=0.5)
# }
distPlot <- function(nscn, X,type,ylab,plot=TRUE,relative=FALSE){
  
  xlab="Metapopulation exploitation rate"
  xlim=c(0,max(frates_vec))
  lwd=2
  
  
  EXPE <- SCN[,nscn]
  scn_id <- as.numeric(strsplit(as.character(EXPE[2]), "")[[1]])
  scenarioConnect = scn_id[1]
  scenarioFishing = scn_id[2]
  scenarioFishingRate = scn_id[3]
  
  if (scenarioFishing==0){ fishingStrat <- "Expl. All"}
  if (scenarioFishing==1){ fishingStrat <- "Expl. Source+neutral"}
  if (scenarioFishing==2){ fishingStrat <- "Expl. Sink+neutral"}
  if (scenarioFishing==3){ fishingStrat <- "Expl. Source only"}
  
  #ylab="# Parr 0+"
  #main=paste0("h=",h[scenarioConnect],"/",fishingStrat)#paste0("Neutral gene (", pops[pop],")")
  main=type
  if(plot==TRUE){
    #ylim=c(0,1.1)#c(min(sapply(X,min,na.rm=TRUE)),max(sapply(X,max,na.rm=TRUE)))#c(0,max(unlist(tmp), na.rm=TRUE))
    xlim=c(0,.3)
    plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
    axis(1,frates_vec,frates_vec)
    #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
    #lgd()
    # rect(-1, 1, 10, 10, col = makeTransparent(col.type[1], 30), border = "transparent")
    # rect(-1, 0, 10, 1, col = makeTransparent(col.type[3], 30), border = "transparent")
  }
  
  
  tmp=tmp2.5=tmp97.5=pGlob=NULL
  j=0
  for (i in EXPE){
    j=j+1
    if(is.na(i)) next;
    
    var<-X[[paste0(i)]]
    varTail <- apply(var, 2, function(x) mean(tail(x , n = 5),na.rm=TRUE))
    #var<-X[[paste0(i)]]
    #varTail <- mean(tail(var , n = 5),na.rm=TRUE)
    
    # tmp[j]<-quantile(varTail,probs=0.5,na.rm=TRUE)
    # tmp2.5[j]<-quantile(varTail,probs=0.025,na.rm=TRUE)
    # tmp97.5[j]<-quantile(varTail,probs=0.975,na.rm=TRUE)
    
    tmp[j]<-mean(varTail,na.rm=TRUE)
    tmp2.5[j]<-tmp[j]-se(varTail)
    tmp97.5[j]<-tmp[j]+se(varTail)
    
    pGlob[j] <- mean(gExpl[[paste0(i)]],na.rm=TRUE)
    
    if (relative){
      if(i==EXPE[1] & j==1) ref <- tmp[1]
    } else {ref <-1}
  }
  
  #colfunc <- colorRampPalette(c("darkgrey", "black"))
  if (nscn==1) mycol <- "lightgrey"
  
  if (nscn==2) mycol <- "darkgrey"
  #if (nscn==3) mycol <- col.type[2] # source+neutral
  if (nscn==3) mycol <- makeTransparent(col.type[1], 80) #col.type[1] # sink+neutral
  if (nscn==4) mycol <- makeTransparent(col.type[3], 80) #col.type[3] # source only
  
  if (nscn==5) mycol <- "black"
  #if (nscn==7) mycol <- makeTransparent(col.type[2], 80) # source+neutral
  if (nscn==6) mycol <- col.type[1] #makeTransparent(col.type[1], 80) # sink+neutral
  if (nscn==7) mycol <- col.type[3] #makeTransparent(col.type[3], 80) # source only
  #colfunc <- colorRampPalette(c(col.sinks[1], col.sinks[2]))}
  #  mycol <- col.sinks[1]} else { mycol <- col.sinks[2]}
  
  x=pGlob
  #polygon(c(x,rev(x)),c(tmp2.5/ref,rev(tmp97.5/ref)),border=NA, col=makeTransparent(mycol,20))
  
  #scatter.smooth(tmp/ref ~ pGlob, span = 2/3, degree = 2,col=mycol)
  #lines(loess.smooth(pGlob,tmp/ref,span = 1/3, degree = 2), col = mycol,lwd=3)
  
  #if(text){
  text(pGlob,tmp/ref, frates_vec, col=mycol,cex=1)
  #}
  #points(pGlob,tmp/ref,col=mycol,type='p',lwd=1,pcDisp. 0%,cex=.6)
  
  #abline(h=1,lty=2,col="black")
  #abline(h=ref,lty=2,col="grey")
  #segments(pGlob,tmp2.5/ref,pGlob,tmp97.5/ref,col=col[1:length(frates_vec)])
  #points(pGlob,tmp/ref,col=col[1:length(frates_vec)],type='b',lwd=lwd,pch=20)
  
}

lgd<-function(cex){
  legend("topright"
         #,c("fishing all","fishing sink+neutral","fishing source+neutral","fishing source only")
         #,c("fishing all","fishing sink+neutral","fishing source only")
         ,c("Exploitation (all pop.)","Conservation of source pop.","Exploitation of source pop.")
         ,col=c("black",col.type[c(1,3)])
         ,text.col=c("black",col.type[c(1,3)])
         #,pch=rep(16,5)
         ,lty=rep(1,4)
         ,bty="n"
         ,cex=cex)
}



pdf(file="results/RATIO.pdf")
#par(mfrow = c(1,1))
par(mfrow = c(1,1), mar=c(4, 4, 1, 1) + 0.1)#it goes c(bottom, left, top, right)

#x<-Immigrants.all
x<-lapply(Immigrants.all,log)
ylab <- "Log-Ratio Immigrants/Harvest"
ylim=c(-1,10)

#### ALL
#pop=12
distPlot(1,x,"All",ylab)
for (sc in c(1,2,5) ){
  distPlot(sc,x,"All",ylab,FALSE,FALSE)
}
lgd(1)
abline(h=0,lty=2,col="black")

#x<-Immigrants.sinks
x<-lapply(Immigrants.sinks,log)
#ylab <- "Ratio Immigrants/Harvest"
#ylim=c(0,10)

#### ALL
#pop=12
distPlot(1,x,"sinks&neutral",ylab,FALSE,FALSE)
for (sc in c(3,6) ){
  #for (sc in 2:7 ){
  distPlot(sc,x,"sinks&neutral",ylab,FALSE,FALSE)
}


#x<-Immigrants.sources
x<-lapply(Immigrants.sources,log)
#ylab <- "Ratio Immigrants/Harvest"
#ylim=c(0,10)
#### ALL
#pop=12
distPlot(1,x,"Sources",ylab,FALSE,FALSE)
for (sc in c(4,7) ){
  #for (sc in 2:7 ){
  distPlot(sc,x,"Sources",ylab,FALSE,FALSE)
}



dev.off()






##### RATIO MIGRANTS VS EXPLOITATION
pdf(file="results/MIGRANTS.pdf")

distPlot <- function(nscn, X, pop,ylab,plot=TRUE,relative=FALSE){
  
  xlab="Metapopulation exploitation rate"
  xlim=c(0,max(frates_vec))
  lwd=2
  
  
  EXPE <- SCN[,nscn]
  scn_id <- as.numeric(strsplit(as.character(EXPE[2]), "")[[1]])
  scenarioConnect = scn_id[1]
  scenarioFishing = scn_id[2]
  scenarioFishingRate = scn_id[3]
  
  if (scenarioFishing==0){ fishingStrat <- "Expl. All"}
  if (scenarioFishing==1){ fishingStrat <- "Expl. Source+neutral"}
  if (scenarioFishing==2){ fishingStrat <- "Expl. Sink+neutral"}
  if (scenarioFishing==3){ fishingStrat <- "Expl. Source only"}
  
  #ylab="# Parr 0+"
  #main=paste0("h=",h[scenarioConnect],"/",fishingStrat)#paste0("Neutral gene (", pops[pop],")")
  main=pops[pop]
  if(plot==TRUE){
    #ylim=c(0,1.1)#c(min(sapply(X,min,na.rm=TRUE)),max(sapply(X,max,na.rm=TRUE)))#c(0,max(unlist(tmp), na.rm=TRUE))
    xlim=c(0,.3)
    plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,xaxt = "n")
    axis(1,frates_vec,frates_vec)
    #title("Philopatry=1 / Fishing all", col.main = col.type[pop.type[pop]])
    #lgd()
    rect(-1, 1, 10, 10, col = makeTransparent(col.type[1], 30), border = "transparent")
    rect(-1, 0, 10, 1, col = makeTransparent(col.type[3], 30), border = "transparent")
  }
  
  
  tmp=tmp2.5=tmp97.5=pGlob=NULL
  j=0
  for (i in EXPE){
    j=j+1
    if(is.na(i)) next;
    
    # var<-X[[paste0(i)]][[1]][,14]
    # varTail <- apply(var, 2, function(x) mean(tail(x , n = 5),na.rm=TRUE))
    var<-X[[paste0(i)]][[1]][,pop]
    varTail <- mean(tail(var , n = 5),na.rm=TRUE)
    
    # tmp[j]<-quantile(varTail,probs=0.5,na.rm=TRUE)
    # tmp2.5[j]<-quantile(varTail,probs=0.025,na.rm=TRUE)
    # tmp97.5[j]<-quantile(varTail,probs=0.975,na.rm=TRUE)
    
    tmp[j]<-mean(varTail,na.rm=TRUE)
    tmp2.5[j]<-tmp[j]-se(varTail)
    tmp97.5[j]<-tmp[j]+se(varTail)
    
    pGlob[j] <- mean(gExpl[[paste0(i)]],na.rm=TRUE)
    
    if (relative){
      if(i==EXPE[1] & j==1) ref <- tmp[1]
    } else {ref <-1}
  }
  
  #colfunc <- colorRampPalette(c("darkgrey", "black"))
  if (nscn==1) mycol <- "lightgrey"
  
  if (nscn==2) mycol <- "darkgrey"
  #if (nscn==3) mycol <- col.type[2] # source+neutral
  if (nscn==3) mycol <- makeTransparent(col.type[1], 80) #col.type[1] # sink+neutral
  if (nscn==4) mycol <- makeTransparent(col.type[3], 80) #col.type[3] # source only
  
  if (nscn==5) mycol <- "black"
  #if (nscn==7) mycol <- makeTransparent(col.type[2], 80) # source+neutral
  if (nscn==6) mycol <- col.type[1] #makeTransparent(col.type[1], 80) # sink+neutral
  if (nscn==7) mycol <- col.type[3] #makeTransparent(col.type[3], 80) # source only
  #colfunc <- colorRampPalette(c(col.sinks[1], col.sinks[2]))}
  #  mycol <- col.sinks[1]} else { mycol <- col.sinks[2]}
  
  x=pGlob
  #polygon(c(x,rev(x)),c(tmp2.5/ref,rev(tmp97.5/ref)),border=NA, col=makeTransparent(mycol,20))
  
  #scatter.smooth(tmp/ref ~ pGlob, span = 2/3, degree = 2,col=mycol)
  #lines(loess.smooth(pGlob,tmp/ref,span = 1/3, degree = 2), col = mycol,lwd=3)
  
  #if(text){
  text(pGlob,tmp/ref, frates_vec, col=mycol,cex=.6)
  #}
  #points(pGlob,tmp/ref,col=mycol,type='p',lwd=1,pcDisp. 0%,cex=.6)
  
  abline(h=1,lty=2,col="black")
  #abline(h=ref,lty=2,col="grey")
  #segments(pGlob,tmp2.5/ref,pGlob,tmp97.5/ref,col=col[1:length(frates_vec)])
  #points(pGlob,tmp/ref,col=col[1:length(frates_vec)],type='b',lwd=lwd,pch=20)
  
}

#par(mfrow = c(1,1))
par(mfrow = c(2,2), mar=c(4, 4, 1, 1) + 0.1)#it goes c(bottom, left, top, right)
x<-RatioMigrants
#x<-lapply(RatioMigrants,log)
ylab <- "Ratio Immigrants/Emigrants"
ylim=c(0.8,1.8)

#### ALL
#pop=12
for (pop in 1:15){
  distPlot(1,x,pop,ylab)
  lgd(1)
  abline(h=1,lty=2,col="black")
  for (sc in 2:7 ){
    distPlot(sc,x,pop,ylab,FALSE,FALSE)
  }
}



#par(mfrow = c(1,1))
par(mfrow = c(2,2), mar=c(4, 4, 1, 1) + 0.1)#it goes c(bottom, left, top, right)
x<-RatioMigrantsVsExpl
#x<-lapply(RatioMigrantsVsExpl,log)
ylab <- "Ratio Immigrants/Harvest"
ylim=c(0,10)

#### ALL
#pop=12
for (pop in 1:15){
  distPlot(1,x,pop,ylab)
  lgd(1)
  abline(h=1,lty=2,col="black")
  for (sc in 2:7 ){
    distPlot(sc,x,pop,ylab,FALSE,FALSE)
  }
}
dev.off()



