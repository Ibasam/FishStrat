## remove (almost) everything in the working environment.
rm(list = ls())
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  args <- 1
}


source("R/Rfunctions.R")

source("parIbasam.R")
pops <- c("Leff", "Trieux", "Jaudy", "Leguer", "Yar", "Douron", "Penze", "Elorn", "Aulne", "Goyen", "Odet", "Aven", "Laita", "Scorff", "Blavet")
nSIMUL <- 50
nYear <- 40
nInit <- 10 # Nb years at initialization
npop <- 15

frates_vec =c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)

scn=c(

  ### 5% dispersal ###
  #exploitation all populations
  2011, NA, 2013, NA, 2015, NA, 2017, NA, 2019, NA, 20111, NA, 20113, NA, 20115, NA, 20117, NA, 20119, NA, 20121 #selective MSW
  ,2021, 2022, 2023, 2024, 2025, 2026, 2027, 2028, 2029, 20210, 20211, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA #non selective
  #source conservation
  ,2211, NA, 2213, NA, 2215, NA, 2217, NA, 2219, NA, 22111, NA, 22113, NA, 22115, NA, 22117, NA, 22119, NA, 22121 #selective MSW
  ,2021, NA, 2223, NA, 2225, NA, 2227, NA, 2229, NA, 22211, NA, 22213, NA, 22215, NA, 22217, NA, 22219, NA, 22221 #non selective
  #source exploitation
  ,2311, NA, 2313, NA, 2315, NA, 2317, NA, 2319, NA, 23111, NA, 23113, NA, 23115, NA, 23117, NA, 23119, NA, 23121 #selective MSW
  ,2021, NA, 2323, NA, 2325, NA, 2327, NA, 2329, NA, 23211, NA, 23213, NA, 23215, NA, 23217, NA, 23219, NA, 23221 #non selective
  
  ### 15% dispersal ###
  #exploitation all populations
  ,4011, NA, 4013, NA, 4015, NA, 4017, NA, 4019, NA, 40111, NA, 40113, NA, 40115, NA, 40117, NA, 40119, NA, 40121 #selective MSW
  ,4021, 4022, 4023, 4024, 4025, 4026, 4027, 4028, 4029, 40210, 40211, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA #non selective
  #source conservation
  ,4211, NA, 4213, NA, 4215, NA, 4217, NA, 4219, NA, 42111, NA, 42113, NA, 42115, NA, 42117, NA, 42119, NA, 42121 #selective MSW
  ,4021, NA, 4223, NA, 4225, NA, 4227, NA, 4229, NA, 42211, NA, 42213, NA, 42215, NA, 42217, NA, 42219, NA, 42221 #non selective
  #source exploitation
  ,4311, NA, 4313, NA, 4315, NA, 4317, NA, 4319, NA, 43111, NA, 43113, NA, 43115, NA, 43117, NA, 43119, NA, 43121 #selective MSW
  ,4021, NA, 4323, NA, 4325, NA, 4327, NA, 4329, NA, 43211, NA, 43213, NA, 43215, NA, 43217, NA, 43219, NA, 43221 #non selective

  #### 30% dispersal ###
  #exploitation all populations
  ,7011, NA, 7013, NA, 7015, NA, 7017, NA, 7019, NA, 70111, NA, 70113, NA, 70115, NA, 70117, NA, 70119, NA, 70121 #selective MSW
  ,7021, 7022, 7023, 7024, 7025, 7026, 7027, 7028, 7029, 70210, 70211, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA #non selective
  #source conservation
  ,7211, NA, 7213, NA, 7215, NA, 7217, NA, 7219, NA, 72111, NA, 72113, NA, 72115, NA, 72117, NA, 72119, NA, 72121 #selective MSW
  ,7021, NA, 7223, NA, 7225, NA, 7227, NA, 7229, NA, 72211, NA, 72213, NA, 72215, NA, 72217, NA, 72219, NA, 72221 #non selective
  #source exploitation
  ,7311, NA, 7313, NA, 7315, NA, 7317, NA, 7319, NA, 73111, NA, 73113, NA, 73115, NA, 73117, NA, 73119, NA, 73121 #selective MSW
  ,7021, NA, 7323, NA, 7325, NA, 7327, NA, 7329, NA, 73211, NA, 73213, NA, 73215, NA, 73217, NA, 73219, NA, 73221 #non selective
  

)
SCN <-matrix(scn,length(frates_vec),18)


#####------------------------ GET EXTRACTED RESULTS ----------------------######

## Demography
nReturns <- list() # Nb returns
Mig <- list() # Nb immigrnats
gExpl <- gCapt <- gCaptMSW <- gCapt1SW <- gCapt.5 <- gCaptMSW.5 <- gCapt1SW.5 <- list()
PE <- list()
PQEXT <- list

## Phenotypes
LFReturns <- list()
rMSW<- list() # Nb immigrnats

## Genotypes
gNEUTRAL <- list() # Nb returns
gFMID1 <- list() # Nb returns
gFMID3 <- list() # Nb immigrnats
gFMID4<- list() # Nb immigrnats
gGROWTH<- list() # Nb immigrnats


CVReturns<-CVParr<-CVSmolt<-CVmsw<-CV1sw<-CVgNEUTRAL<-CVgFMID1<-CVgFMID3<-CVgFMID4<-CVgGROWTH<-list()
SDgNEUTRAL<-SDgFMID1<-SDgFMID3<-SDgFMID4<-SDgGROWTH<-list()



for (i in 1:ncol(SCN)){
  j=0
  for (iEXPE in SCN[,i]){
    j=j+1
    if(is.na(iEXPE)) next;
    
    load(paste0("results/METAPOP",iEXPE,".RData"))

    
    nReturns[[paste0(iEXPE)]] <- NRet

    Mig[[paste0(iEXPE)]] <- RETURNS
    
    gExpl[[paste0(iEXPE)]] <- Pexpl.global
    
    gCapt[[paste0(iEXPE)]] <- Ncapt.global
    gCaptMSW[[paste0(iEXPE)]] <- Ncapt.MSW
    gCapt1SW[[paste0(iEXPE)]] <- Ncapt.1SW
    gCapt.5[[paste0(iEXPE)]] <- Ncapt.global.5
    gCaptMSW.5[[paste0(iEXPE)]] <- Ncapt.MSW.5
    gCapt1SW.5[[paste0(iEXPE)]] <- Ncapt.1SW.5

    PE[[paste0(iEXPE)]] <- PE_mv
    
    PQEXT[[paste0(iEXPE)]] <- PQext
    
    ## Phenotypes
    LFReturns[[paste0(iEXPE)]] <- LFreturns
    rMSW[[paste0(iEXPE)]] <- RatioMSW

    CVParr[[paste0(iEXPE)]] <- CVparr0
    CVSmolt[[paste0(iEXPE)]] <- CVsmolt0
    CVReturns[[paste0(iEXPE)]] <- CVreturns
    CVmsw[[paste0(iEXPE)]] <- CVMSW
    CV1sw[[paste0(iEXPE)]] <- CV1SW
    
    
    ## Genotypes
    gNEUTRAL[[paste0(iEXPE)]] <- gNeutral
    gFMID1[[paste0(iEXPE)]] <- gFmid1
    gFMID3[[paste0(iEXPE)]] <- gFmid3
    gFMID4[[paste0(iEXPE)]] <- gFmid4
    gGROWTH[[paste0(iEXPE)]] <- gG
    
    CVgNEUTRAL[[paste0(iEXPE)]] <- CVgNeutral
    CVgFMID1[[paste0(iEXPE)]] <- CVgFmid1
    CVgFMID3[[paste0(iEXPE)]] <- CVgFmid3
    CVgFMID4[[paste0(iEXPE)]] <- CVgFmid4
    CVgGROWTH[[paste0(iEXPE)]] <- CVgG
    
    SDgNEUTRAL[[paste0(iEXPE)]] <- SDgNeutral
    SDgFMID1[[paste0(iEXPE)]] <- SDgFmid1
    SDgFMID3[[paste0(iEXPE)]] <- SDgFmid3
    SDgFMID4[[paste0(iEXPE)]] <- SDgFmid4
    SDgGROWTH[[paste0(iEXPE)]] <- SDgG

    
  } # loop iEXPE
  
} # loop SCN


#####------------------------ COLORS AND PLOT FUNCTIONS ----------------------######

library("NatParksPalettes")

color_ns <- natparks.pals("Yellowstone")[4] #3
color_s <- "tomato2" #9

colors <- c("royalblue3", natparks.pals("Yellowstone")[4],"orange")

lgd_spatial<-function(cex, position){
  legend(position
         ,c("Exploitation source","Exploitation sink + neutral", "Exploitation all pop.")
         ,col=colors[c(3,2,1)]
         ,text.col=colors[c(3,2,1)]
         #,pch=rep(16,5)
         ,bty="n",
         cex=cex)
}

distPlot_PE <- function(nscn, X,xaxis, ylab,synchrony=FALSE,plot=TRUE, smooth=TRUE){
  
  if (xaxis == "MER") {
    xlab="Metapopulation exploitation rate"
  }
  if (xaxis == "Captures") {
    xlab="Captures"
  }
  if (xaxis == "LER") {
    xlab="Local exploitation rate"
  }
  xlim=c(0,max(frates_vec))
  
  EXPE <- SCN[,nscn]
  scn_id <- as.numeric(strsplit(as.character(EXPE[3]), "")[[1]])
  scenarioConnect = scn_id[1]
  scenarioFishing = scn_id[2]
  scenarioFishingRate = scn_id[4]
  
  if (scenarioFishing==0){ fishingStrat <- "Expl. All"}
  if (scenarioFishing==2){ fishingStrat <- "Expl. Sink+neutral"}
  if (scenarioFishing==3){ fishingStrat <- "Expl. Source only"}
  
  main=""
  if(plot==TRUE){
    if (xaxis == "MER") {
      xlim=c(0, .25)
    }
    if (xaxis == "Captures") {
      xlim=c(0,500)
    }
    if (xaxis == "LER") {
      xlim=c(0,0.8)
    }
    
    plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main, cex.lab=1.3, cex.axis=1.1)
    
  }
  
  tmp=tmp2.5=tmp97.5=pGlob=nCapt=NULL
  j=0
  for (i in EXPE){
    j=j+1
    if(is.na(i)) next;
    
    varTail=NULL
    
    if(synchrony) {
      for (sim in 1:nSIMUL) {varTail[sim]<-X[[paste0(i)]][[sim]]$Synchrony}
    } else {
      for (sim in 1:nSIMUL) {varTail[sim]<-X[[paste0(i)]][[sim]]$pe}
    }
    
    tmp[j]<-mean(varTail,na.rm=TRUE)
    tmp2.5[j]<-tmp[j]-se(varTail)
    tmp97.5[j]<-tmp[j]+se(varTail)
    
    pGlob[j] <- mean(gExpl[[paste0(i)]],na.rm=TRUE)
    nCapt[j] <- mean(gCapt[[paste0(i)]],na.rm=TRUE)
    
  }
  
  lty=1
  
  #5%-15-30 dispersal
  if (nscn %in% c(1,2,7,8,13,14)) mycol <- colors[1] #exploitation all populations
  if (nscn %in% c(3,4,9,10,15,16)) {mycol <- colors[2]} # exploitation sink and neutral
  if (nscn %in% c(5,6,11,12,17,18)) {mycol <- colors[3]} # exploitation source only
  
  if (xaxis == "MER") {
    xx=pGlob
  }
  if (xaxis == "Captures") {
    xx=nCapt
  }
  if (xaxis == "LER") {
    xx=frates_vec
  
  }
  
  if(smooth == TRUE) {
    lines(loess.smooth(xx,tmp,span = 0.6, degree = 1), col = mycol,lwd=5)
  } else {
    lines(na.omit(xx),na.omit(tmp), col = mycol,lwd=5)
  }
  
  
  if(nscn==2 | nscn==8 | nscn==14){
    text(xx[c(1,5,9)],tmp[c(1,5,9)]+0.01*mean(tmp, na.rm=T), frates_vec[c(1,5,9)], col=mycol,cex=.6)
  }else{
    text(xx[c(1,5,9,13,17,21)],tmp[c(1,5,9,13,17,21)]+0.01*mean(tmp, na.rm=T), frates_vec[c(1,5,9,13,17,21)], col=mycol,cex=.6)
  }
  
  
}
plot_PE_5_NS <- function (scnToPlot, horizLines, x) {
  distPlot_PE(scnToPlot[1],x,"MER",ylab, synchrony = FALSE, plot = TRUE, smooth = FALSE)
  abline(h=horizLines, lty=2, col="lightgrey")
  distPlot_PE(scnToPlot[2],x,"MER",ylab, synchrony = FALSE, plot = FALSE, smooth = FALSE) 
  distPlot_PE(scnToPlot[3],x,"MER",ylab, synchrony = FALSE, plot = FALSE, smooth = FALSE) 
  lgd_spatial(0.8,"bottomleft")
}
plot_PE_all <- function (scnToPlot, horizLines, x) {
  distPlot_PE(scnToPlot[1],x,"MER",ylab, synchrony = FALSE, plot = TRUE)
  abline(h=horizLines, lty=2, col="lightgrey")
  distPlot_PE(scnToPlot[2],x,"MER",ylab, synchrony = FALSE, plot = FALSE) # sink+neutral
  distPlot_PE(scnToPlot[3],x,"MER",ylab, synchrony = FALSE, plot = FALSE) # source only
  lgd_spatial(0.8,"bottomleft")
}

distPlot <- function(nscn, X,xaxis, ylab, plot=TRUE, relative=TRUE, smooth=TRUE){
  
  if (xaxis == "MER") {
    xlab="Metapopulation exploitation rate"
  }
  if (xaxis == "Captures") {
    xlab="Captures"
  }
  if (xaxis == "LER") {
    xlab="Local exploitation rate"
  }
  xlim=c(0,max(frates_vec))
  
  EXPE <- SCN[,nscn]
  scn_id <- as.numeric(strsplit(as.character(EXPE[3]), "")[[1]])
  scenarioConnect = scn_id[1]
  scenarioFishing = scn_id[2]
  scenarioFishingRate = scn_id[4]
  
  if (scenarioFishing==0){ fishingStrat <- "Expl. All"}
  if (scenarioFishing==2){ fishingStrat <- "Expl. Sink+neutral"}
  if (scenarioFishing==3){ fishingStrat <- "Expl. Source only"}
  
  main=""
  if(plot==TRUE){
    if (xaxis == "MER") {
      xlim=c(0, .25)
    }
    if (xaxis == "Captures") {
      xlim=c(0,500)
    }
    if (xaxis == "LER") {
      xlim=c(0,0.8)
    }
    
    plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main, cex.lab=1.3, cex.axis=1.1)
    
  }
  
  
  tmp=tmp2.5=tmp97.5=pGlob=nCapt=NULL
  j=0
  for (i in EXPE){
    j=j+1
    if(is.na(i)) next;
    
    var<-X[[paste0(i)]]
    if(length(var)==nSIMUL) {
      varTail <- var
    }else{
      if (var[49,18] > 1.7e+16) {var[49,18]<-var[48,18]} #bug for this simul 18
      if (var[46,41] > 1.5e+16) {var[46,41]<-var[45,41]} #bug for this simul 18
      if (var[48,40] > 1.5e+14) {var[48,40]<-var[47,40]} #bug for this simul 18
      if (any(var[,45] > 10000, na.rm=T)) {var[,45]<-var[,44]} #bug for this simul 18
      if (any(var[,9] > 10000, na.rm=T)) {var[,9]<-var[,8]} #bug for this simul 18
      
      varTail <- apply(var, 2, function(x) mean(tail(x , n = 5),na.rm=TRUE))
    }
    
    tmp[j]<-mean(varTail,na.rm=TRUE)
    tmp2.5[j]<-tmp[j]-se(varTail)
    tmp97.5[j]<-tmp[j]+se(varTail)
    
    pGlob[j] <- mean(gExpl[[paste0(i)]],na.rm=TRUE)
    nCapt[j] <- mean(gCapt[[paste0(i)]],na.rm=TRUE)
    
    if (relative){
      if(i==EXPE[1] & j==1) ref <- tmp[1]
    } else {ref <-1}
  }
  
  lty=1
  
  #5%-15-30 dispersal
  if (nscn %in% c(1,2,7,8,13,14)) mycol <- colors[1] #exploitation all populations
  if (nscn %in% c(3,4,9,10,15,16)) {mycol <- colors[2]} # exploitation sink and neutral
  if (nscn %in% c(5,6,11,12,17,18)) {mycol <- colors[3]} # exploitation source only
  
  if (xaxis == "MER") {
    xx=pGlob
  }
  if (xaxis == "Captures") {
    xx=nCapt
  }
  if (xaxis == "LER") {
    xx=frates_vec
  }
  
  if(smooth == TRUE) {
    lines(loess.smooth(xx,tmp,span = 0.55, degree = 1), col = mycol,lwd=5)
  } else {
    lines(na.omit(xx),na.omit(tmp), col = mycol,lwd=5)
  }
  
  if(nscn==2 | nscn==8 | nscn==14){
    text(xx[c(1,5,9)],(tmp/ref)[c(1,5,9)] + 0.01*mean(tmp/ref, na.rm=T), frates_vec[c(1,5,9)], col=mycol,cex=.6)
  }else{
    text(xx[c(1,5,9,13,17,21)],(tmp/ref)[c(1,5,9,13,17,21)]+ 0.01*mean(tmp/ref, na.rm=T), frates_vec[c(1,5,9,13,17,21)], col=mycol,cex=.6)
  }
  
}
plot_5_NS <- function (scnToPlot, horizLines, x) {
  distPlot(scnToPlot[1],x,"MER",ylab, plot = TRUE, relative = FALSE, smooth = FALSE)
  abline(h=horizLines, lty=2, col="lightgrey")
  distPlot(scnToPlot[2],x,"MER",ylab, plot = FALSE, relative = FALSE, smooth = FALSE) # sink+neutral
  distPlot(scnToPlot[3],x,"MER",ylab, plot = FALSE, relative = FALSE, smooth = FALSE) # source only
  lgd_spatial(0.8,"bottomleft")

}
plot_all <- function (scnToPlot, horizLines, x) {
  distPlot(scnToPlot[1],x,"MER",ylab, plot = TRUE, relative = FALSE)
  abline(h=horizLines, lty=2, col="lightgrey")
  distPlot(scnToPlot[2],x,"MER",ylab, plot = FALSE, relative = FALSE) # sink+neutral
  distPlot(scnToPlot[3],x,"MER",ylab, plot = FALSE, relative = FALSE) # source only
  lgd_spatial(0.8,"bottomleft")
}

distPlot_Extinction <- function(nscn, X,xaxis, ylab, plot=TRUE, relative=TRUE, pop="all", weighted=FALSE, smooth=TRUE){
  
  if (xaxis == "MER") {
    xlab="Metapopulation exploitation rate"
  }
  if (xaxis == "Captures") {
    xlab="Captures"
  }
  if (xaxis == "LER") {
    xlab="Local exploitation rate"
  }
  xlim=c(0,max(frates_vec))
  
  EXPE <- SCN[,nscn]
  scn_id <- as.numeric(strsplit(as.character(EXPE[3]), "")[[1]])
  scenarioConnect = scn_id[1]
  scenarioFishing = scn_id[2]
  scenarioFishingRate = scn_id[4]
  
  if (scenarioFishing==0){ fishingStrat <- "Expl. All"}
  if (scenarioFishing==2){ fishingStrat <- "Expl. Sink+neutral"}
  if (scenarioFishing==3){ fishingStrat <- "Expl. Source only"}
  
  main=""
  if(plot==TRUE){
    if (xaxis == "MER") {
      xlim=c(0, .25)
    }
    if (xaxis == "Captures") {
      xlim=c(0,500)
    }
    if (xaxis == "LER") {
      xlim=c(0,0.8)
    }
    
    plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main, cex.lab=1.3, cex.axis=1.1)
    
  }
  
  
  tmp=tmp2.5=tmp97.5=pGlob=nCapt=NULL
  j=0
  for (i in EXPE){
    j=j+1
    if(is.na(i)) next;
    
    var<-X[[paste0(i)]]
    
    if(pop=="sink") {
      pops_id <- c(1,3,5,10,12,14)
      if(weighted == TRUE) {
        tmp[j]<- weighted.mean(var[pops_id], Area[pops_id])
      } else {
        tmp[j]<-mean(var[pops_id],na.rm=TRUE)
        #tmp2.5[j]<-tmp[j]-se(var)
        #tmp97.5[j]<-tmp[j]+se(var)        
      }
    }
    
    if(pop=="source") {
      pops_id <- c(2,4,9,11,13)
      if(weighted == TRUE) {
        tmp[j]<- weighted.mean(var[pops_id], Area[pops_id])
      } else {
        tmp[j]<-mean(var[pops_id],na.rm=TRUE)
        #tmp2.5[j]<-tmp[j]-se(var)
        #tmp97.5[j]<-tmp[j]+se(var) 
      }
    }
    
    if(pop=="all") {
      if(weighted == TRUE) {
        tmp[j]<- weighted.mean(var, Area)
      } else {
        tmp[j]<-mean(var,na.rm=TRUE)
        #tmp2.5[j]<-tmp[j]-se(var)
        #tmp97.5[j]<-tmp[j]+se(var) 
      }
    }
    
    pGlob[j] <- mean(gExpl[[paste0(i)]],na.rm=TRUE)
    nCapt[j] <- mean(gCapt[[paste0(i)]],na.rm=TRUE)
    
    if (relative){
      if(i==EXPE[1] & j==1) ref <- tmp[1]
    } else {ref <-1}
  }
  
  lty=1
  
  #5%-15-30 dispersal
  if (nscn %in% c(1,2,7,8,13,14)) mycol <- colors[1] #exploitation all populations
  if (nscn %in% c(3,4,9,10,15,16)) {mycol <- colors[2]} # exploitation sink and neutral
  if (nscn %in% c(5,6,11,12,17,18)) {mycol <- colors[3]} # exploitation source only
  
  if (xaxis == "MER") {
    xx=pGlob
  }
  if (xaxis == "Captures") {
    xx=nCapt
  }
  if (xaxis == "LER") {
    xx=frates_vec
  }
  
  if(smooth == TRUE) {
    lines(loess.smooth(xx,tmp,span = 0.55, degree = 1), col = mycol,lwd=5)
  } else {
    lines(na.omit(xx),na.omit(tmp), col = mycol,lwd=5)
  }
  
  if(nscn==2 | nscn==8 | nscn==14){
    text(xx[c(1,5,9)],(tmp/ref)[c(1,5,9)] + 0.01*mean(tmp/ref, na.rm=T), frates_vec[c(1,5,9)], col=mycol,cex=.6)
  }else{
    text(xx[c(1,5,9,13,17,21)],(tmp/ref)[c(1,5,9,13,17,21)]+ 0.01*mean(tmp/ref, na.rm=T), frates_vec[c(1,5,9,13,17,21)], col=mycol,cex=.6)
  }
  
}
plot_extinction_5_NS <- function (scnToPlot, horizLines, x, pop.type) {
  distPlot_Extinction(scnToPlot[1],x,"MER",ylab, plot=TRUE, relative=FALSE, pop=pop.type, smooth=FALSE)
  abline(h=horizLines, lty=2, col="lightgrey")
  distPlot_Extinction(scnToPlot[2],x,"MER",ylab, plot=FALSE, relative=FALSE, pop=pop.type, smooth=FALSE) # sink+neutral
  distPlot_Extinction(scnToPlot[3],x,"MER",ylab, plot=FALSE, relative=FALSE, pop=pop.type, smooth=FALSE) # source only
  lgd_spatial(0.8,"topleft")
}
plot_extinction_all <- function (scnToPlot, horizLines, x, pop.type) {
  distPlot_Extinction(scnToPlot[1],x,"MER",ylab, plot=TRUE, relative=FALSE, pop=pop.type)
  abline(h=horizLines, lty=2, col="lightgrey")
  distPlot_Extinction(scnToPlot[2],x,"MER",ylab, plot=FALSE, relative=FALSE, pop=pop.type) # sink+neutral
  distPlot_Extinction(scnToPlot[3],x,"MER",ylab, plot=FALSE, relative=FALSE, pop=pop.type) # source only
  lgd_spatial(0.8,"topleft")
}

#####------------------------ PORTFOLIO EFFECT ----------------------######

x<-PE
ylab <- "Portfolio effect"
ylim=c(0.9,2.5)
horizLines <- c(1,1.5,2)

#### Disp. 5%
scnToPlot <- c(2,4,6)
plot_PE_5_NS(scnToPlot, horizLines, x)
title("5% dispersal - Not selective")
#### Disp. 15%
scnToPlot <- c(8,10,12)
plot_PE_all(scnToPlot, horizLines, x)
title("15% dispersal - Not selective")
#### Disp. 30%
scnToPlot <- c(14,16,18)
plot_PE_all(scnToPlot, horizLines, x)
title("30% dispersal - Not selective")

##Selective fisheries
#### Disp. 5%
scnToPlot <- c(1,3,5)
plot_PE_all(scnToPlot, horizLines, x)
title("5% dispersal - Selective MSW")
#### Disp. 15%
scnToPlot <- c(7,9,11)
plot_PE_all(scnToPlot, horizLines, x)
title("15% dispersal - Selective MSW")
#### Disp. 30%
scnToPlot <- c(13,15,17)
plot_PE_all(scnToPlot, horizLines, x)
title("30% dispersal - Selective MSW")



#####------------------------ PRE FISHERIES ABUNDANCE ----------------------######

x<-nReturns  
#ylab <- "Relative Abundance (PFA)"
ylab <- "Abundance (PFA)"
#ylim=c(0.2,1.1) #Relative
ylim=c(600,2000)
horizLines <- c(800,1200,1600)

#### Disp. 5%
scnToPlot <- c(2,4,6)
plot_5_NS(scnToPlot, horizLines, x)
title("5% dispersal - Not selective")
#### Disp. 15%
scnToPlot <- c(8,10,12)
plot_all(scnToPlot, horizLines, x)
title("15% dispersal - Not selective")
##30% dispersal
scnToPlot <- c(14,16,18)
plot_all(scnToPlot, horizLines, x)
title("30% dispersal - Not selective")

#Selective fisheries
#5% dispersal
scnToPlot <- c(1,3,5)
plot_all(scnToPlot, horizLines, x)
title("5% dispersal - Selective MSW")
#15% dispersal
scnToPlot <- c(7,9,11)
plot_all(scnToPlot, horizLines, x)
title("15% dispersal - Selective MSW")
##30% dispersal
scnToPlot <- c(13,15,17)
plot_all(scnToPlot, horizLines, x)
title("30% dispersal - Selective MSW")



#####------------------------ EXPLOITATION ----------------------######

### ----------- TOTAL 
#x<-gCapt  
x<-gCapt.5
#ylab <- "Relative Abundance (PFA)"
ylab <- "Total catch"
#ylim=c(0.2,1.1) #Relative
ylim=c(0,350)
horizLines <- c(100,200,300)

#### Disp. 5%
scnToPlot <- c(2,4,6)
plot_5_NS(scnToPlot, horizLines, x)
title("5% dispersal - Not selective")
#### Disp. 15%
scnToPlot <- c(8,10,12)
plot_all(scnToPlot, horizLines, x)
title("15% dispersal - Not selective")
##30% dispersal
scnToPlot <- c(14,16,18)
plot_all(scnToPlot, horizLines, x)
title("30% dispersal - Not selective")

#Selective fisheries
#5% dispersal
scnToPlot <- c(1,3,5)
plot_all(scnToPlot, horizLines, x)
title("5% dispersal - Selective MSW")
#15% dispersal
scnToPlot <- c(7,9,11)
plot_all(scnToPlot, horizLines, x)
title("15% dispersal - Selective MSW")
##30% dispersal
scnToPlot <- c(13,15,17)
plot_all(scnToPlot, horizLines, x)
title("30% dispersal - Selective MSW")


###----------- MSW ONLY
#x<-gCaptMSW  
x<-gCaptMSW.5  
#ylab <- "Relative Abundance (PFA)"
ylab <- "MSW catch"
#ylim=c(0.2,1.1) #Relative
ylim=c(0,100)
horizLines <- c(40,80)

#### Disp. 5%
scnToPlot <- c(2,4,6)
plot_5_NS(scnToPlot, horizLines, x)
title("5% dispersal - Not selective")
#### Disp. 15%
scnToPlot <- c(8,10,12)
plot_all(scnToPlot, horizLines, x)
title("15% dispersal - Not selective")
##30% dispersal
scnToPlot <- c(14,16,18)
plot_all(scnToPlot, horizLines, x)
title("30% dispersal - Not selective")

#Selective fisheries
#5% dispersal
scnToPlot <- c(1,3,5)
plot_all(scnToPlot, horizLines, x)
title("5% dispersal - Selective MSW")
#15% dispersal
scnToPlot <- c(7,9,11)
plot_all(scnToPlot, horizLines, x)
title("15% dispersal - Selective MSW")
##30% dispersal
scnToPlot <- c(13,15,17)
plot_all(scnToPlot, horizLines, x)
title("30% dispersal - Selective MSW")


###---------------- OSW ONLY
#x<-gCapt1SW  
x<-gCapt1SW.5  
#ylab <- "Relative Abundance (PFA)"
ylab <- "1SW catch"
#ylim=c(0.2,1.1) #Relative
ylim=c(0,200)
horizLines <- c(50,150)

#### Disp. 5%
scnToPlot <- c(2,4,6)
plot_5_NS(scnToPlot, horizLines, x)
title("5% dispersal - Not selective")
#### Disp. 15%
scnToPlot <- c(8,10,12)
plot_all(scnToPlot, horizLines, x)
title("15% dispersal - Not selective")
##30% dispersal
scnToPlot <- c(14,16,18)
plot_all(scnToPlot, horizLines, x)
title("30% dispersal - Not selective")

#Selective fisheries
#5% dispersal
scnToPlot <- c(1,3,5)
plot_all(scnToPlot, horizLines, x)
title("5% dispersal - Selective MSW")
#15% dispersal
scnToPlot <- c(7,9,11)
plot_all(scnToPlot, horizLines, x)
title("15% dispersal - Selective MSW")
##30% dispersal
scnToPlot <- c(13,15,17)
plot_all(scnToPlot, horizLines, x)
title("30% dispersal - Selective MSW")



#####------------------------ MATURATION STRATEGY ----------------------######

x<-rMSW
ylab <- "Ratio MSW/1SW"
#ylim=c(-1.5,-1.3) #Relative
ylim=c(0,0.3)
horizLines <- c(0.05,0.15,0.25)

#### Disp. 5%
scnToPlot <- c(2,4,6)
plot_all(scnToPlot, horizLines, x)
title("5% dispersal - Not selective")
#### Disp. 15%
scnToPlot <- c(8,10,12)
plot_all(scnToPlot, horizLines, x)
title("15% dispersal - Not selective")
##30% dispersal
scnToPlot <- c(14,16,18)
plot_all(scnToPlot, horizLines, x)
title("30% dispersal - Not selective")

#Selective fisheries
#5% dispersal
scnToPlot <- c(1,3,5)
plot_all(scnToPlot, horizLines, x)
title("5% dispersal - Selective MSW")
#15% dispersal
scnToPlot <- c(7,9,11)
plot_all(scnToPlot, horizLines, x)
title("15% dispersal - Selective MSW")
##30% dispersal
scnToPlot <- c(13,15,17)
plot_all(scnToPlot, horizLines, x)
title("30% dispersal - Selective MSW")


#####------------------------ GROWTH POTENTIAL ----------------------######

x<-gGROWTH
ylab <- "Mean growth potential"
ylim=c(0.003,0.022)
horizLines <- c(0.01,0.02)

x<-SDgGROWTH
ylab <- "SD Growth potential"
ylim=c(0.07,0.075)
horizLines <- c(0.071,0.074)

#x<-CVgGROWTH
#ylab <- "CV Growth potential"
#ylim=c(0,10)
#ylim=c(0,0.022)

#### Disp. 5%
scnToPlot <- c(2,4,6)
plot_5_NS(scnToPlot, horizLines, x)
title("5% dispersal - Not selective")
#### Disp. 15%
scnToPlot <- c(8,10,12)
plot_all(scnToPlot, horizLines, x)
title("15% dispersal - Not selective")
##30% dispersal
scnToPlot <- c(14,16,18)
plot_all(scnToPlot, horizLines, x)
title("30% dispersal - Not selective")

#Selective fisheries
#5% dispersal
scnToPlot <- c(1,3,5)
plot_all(scnToPlot, horizLines, x)
title("5% dispersal - Selective MSW")
#15% dispersal
scnToPlot <- c(7,9,11)
plot_all(scnToPlot, horizLines, x)
title("15% dispersal - Selective MSW")
##30% dispersal
scnToPlot <- c(13,15,17)
plot_all(scnToPlot, horizLines, x)
title("30% dispersal - Selective MSW")


#####------------------------ MATURATION THRESHOLD (MALE PARR) ----------------------######

x<-gFMID1
ylab <- "Mean maturation threshold (male parr)"
ylim=c(1.2,1.35)
horizLines <- c(1.25,1.3)

#### Disp. 5%
scnToPlot <- c(2,4,6)
plot_5_NS(scnToPlot, horizLines, x)
title("5% dispersal - Not selective")
#### Disp. 15%
scnToPlot <- c(8,10,12)
plot_all(scnToPlot, horizLines, x)
title("15% dispersal - Not selective")
##30% dispersal
scnToPlot <- c(14,16,18)
plot_all(scnToPlot, horizLines, x)
title("30% dispersal - Not selective")

#Selective fisheries
#5% dispersal
scnToPlot <- c(1,3,5)
plot_all(scnToPlot, horizLines, x)
title("5% dispersal - Selective MSW")
#15% dispersal
scnToPlot <- c(7,9,11)
plot_all(scnToPlot, horizLines, x)
title("15% dispersal - Selective MSW")
##30% dispersal
scnToPlot <- c(13,15,17)
plot_all(scnToPlot, horizLines, x)
title("30% dispersal - Selective MSW")



#####------------------------ MATURATION THRESHOLD (FEMALE ANADROMOUS) ----------------------######

x<-gFMID4
ylab <- "Mean maturation threshold (fem.)"
ylim=c(30,120)
horizLines <- c(55,50)
# x<-CVgFMID4
# ylab <- "CV Anadromous mat. threshold (female)"
# ylim=c(0,1)
x<-SDgFMID4
ylab <- "SD maturation threshold (fem.)"
ylim=c(45,60)
horizLines <- c(55,50)

#### Disp. 5%
scnToPlot <- c(2,4,6)
plot_all(scnToPlot, horizLines, x)
title("5% dispersal - Not selective")
#### Disp. 15%
scnToPlot <- c(8,10,12)
plot_all(scnToPlot, horizLines, x)
title("15% dispersal - Not selective")
##30% dispersal
scnToPlot <- c(14,16,18)
plot_all(scnToPlot, horizLines, x)
title("30% dispersal - Not selective")

#Selective fisheries
#5% dispersal
scnToPlot <- c(1,3,5)
plot_all(scnToPlot, horizLines, x)
title("5% dispersal - Selective MSW")
#15% dispersal
scnToPlot <- c(7,9,11)
plot_all(scnToPlot, horizLines, x)
title("15% dispersal - Selective MSW")
##30% dispersal
scnToPlot <- c(13,15,17)
plot_all(scnToPlot, horizLines, x)
title("30% dispersal - Selective MSW")


#####------------------------ EXTINCTION RISK ----------------------######
x<-PQEXT
ylim=c(0,1)
horizLines <- c(0.2,0.6)

###------ ALL POP
ylab <- "Quasi-extinction probability (all pop.)"
#### Disp. 5%
scnToPlot <- c(2,4,6)
plot_extinction_5_NS(scnToPlot, horizLines, x, pop.type = "all")
title("5% dispersal - All pop (not weighted) - 20%Rmax - Not selective")
#### Disp. 15%
scnToPlot <- c(8,10,12)
plot_extinction_all(scnToPlot, horizLines, x, pop.type = "all")
title("15% dispersal - All pop (not weighted) - 20%Rmax - Not selective")
##30% dispersal
scnToPlot <- c(14,16,18)
plot_extinction_all(scnToPlot, horizLines, x, pop.type = "all")
title("30% dispersal - All pop (not weighted) - 20%Rmax - Not selective")

#Selective fisheries
#5% dispersal
scnToPlot <- c(1,3,5)
plot_extinction_all(scnToPlot, horizLines, x, pop.type = "all")
title("5% dispersal - All pop (not weighted) - 20%Rmax - Selective MSW")
#15% dispersal
scnToPlot <- c(7,9,11)
plot_extinction_all(scnToPlot, horizLines, x, pop.type = "all")
title("15% dispersal - All pop (not weighted) - 20%Rmax - Selective MSW")
##30% dispersal
scnToPlot <- c(13,15,17)
plot_extinction_all(scnToPlot, horizLines, x, pop.type = "all")
title("30% dispersal - All pop (not weighted) - 20%Rmax - Selective MSW")


###------ SINK POP ONLY
ylab <- "Quasi-extinction probability (sink pop.)"
#### Disp. 5%
scnToPlot <- c(2,4,6)
plot_extinction_5_NS(scnToPlot, horizLines, x, pop.type = "sink")
title("5% dispersal - sink pop (not weighted) - 20%Rmax - Not selective")
#### Disp. 15%
scnToPlot <- c(8,10,12)
plot_extinction_all(scnToPlot, horizLines, x, pop.type = "sink")
title("15% dispersal - sink pop (not weighted) - 20%Rmax - Not selective")
##30% dispersal
scnToPlot <- c(14,16,18)
plot_extinction_all(scnToPlot, horizLines, x, pop.type = "sink")
title("30% dispersal - sink pop (not weighted) - 20%Rmax - Not selective")

#Selective fisheries
#5% dispersal
scnToPlot <- c(1,3,5)
plot_extinction_all(scnToPlot, horizLines, x, pop.type = "sink")
title("5% dispersal - sink pop (not weighted) - 20%Rmax - Selective MSW")
#15% dispersal
scnToPlot <- c(7,9,11)
plot_extinction_all(scnToPlot, horizLines, x, pop.type = "sink")
title("15% dispersal - sink pop (not weighted) - 20%Rmax - Selective MSW")
##30% dispersal
scnToPlot <- c(13,15,17)
plot_extinction_all(scnToPlot, horizLines, x, pop.type = "sink")
title("30% dispersal - sink pop (not weighted) - 20%Rmax - Selective MSW")


###------ SOURCE POP ONLY
ylab <- "Quasi-extinction probability (source pop.)"
#### Disp. 5%
scnToPlot <- c(2,4,6)
plot_extinction_5_NS(scnToPlot, horizLines, x, pop.type = "source")
title("5% dispersal - source pop (not weighted) - 20%Rmax - Not selective")
#### Disp. 15%
scnToPlot <- c(8,10,12)
plot_extinction_all(scnToPlot, horizLines, x, pop.type = "source")
title("15% dispersal - source pop (not weighted) - 20%Rmax - Not selective")
##30% dispersal
scnToPlot <- c(14,16,18)
plot_extinction_all(scnToPlot, horizLines, x, pop.type = "source")
title("30% dispersal - source pop (not weighted) - 20%Rmax - Not selective")

#Selective fisheries
#5% dispersal
scnToPlot <- c(1,3,5)
plot_extinction_all(scnToPlot, horizLines, x, pop.type = "source")
title("5% dispersal - source pop (not weighted) - 20%Rmax - Selective MSW")
#15% dispersal
scnToPlot <- c(7,9,11)
plot_extinction_all(scnToPlot, horizLines, x, pop.type = "source")
title("15% dispersal - source pop (not weighted) - 20%Rmax - Selective MSW")
##30% dispersal
scnToPlot <- c(13,15,17)
plot_extinction_all(scnToPlot, horizLines, x, pop.type = "source")
title("30% dispersal - source pop (not weighted) - 20%Rmax - Selective MSW")



###########################################################
### Fig.1 Populations network for 15% dispersal rate ###
###########################################################

EXPE <- 4021 #15% dispersal, no exploitation
load("results/METAPOP4021.RData")
Mig[[paste0(EXPE)]] <- RETURNS

#Migrants flows for arrows
a <- array(,dim=c(npop,npop,nSIMUL))
  for (simul in 1:nSIMUL) {
    for (pop in 1:npop) {
      a[pop,,simul]<-apply(Mig[[paste0(EXPE)]][[simul]]$Im[[pop]][46:50,],2,mean)
      for (pop2 in 1:npop) {
        if (a[pop,pop2,simul]=="NaN" || is.na(a[pop,pop2,simul])) {
          a[pop,pop2,simul]=0
        }
      }
    }
  }

  arr <- array( unlist(a[,,]) , c(15,15,50) )
  mflow<-apply( arr , 1:2 , mean )#mean simulations
  colnames(mflow)<-pops
  rownames(mflow)<-pops


library(reshape2)
mflow2<-melt(mflow) #for absolute values 
mflow3<-mflow2[c(2,1,3)] #for absolute values

mflow3[,1]<-as.factor(mflow3[,1])
mflow3[,2]<-as.factor(mflow3[,2])

#populations size for nodes
Parpop<-array(dim=c(nSIMUL, npop))
  for (simul in 1:nSIMUL) {
    for (pop in 1:npop) {
      Parpop[simul, ] <- colMeans(Mig[[paste0(EXPE)]][[simul]]$Returns[46:50,], na.rm=T) #5last years
    }
  }
  Parpop2<-colMeans(Parpop)

dat <- read.csv2("data/dataPop.csv", header=TRUE, stringsAsFactors = FALSE, comment.char = "#")
Type=dat$Type
#network
nodes<-data.frame(id=c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12","p13","p14","p15"),
                  pop=pops,
                  #type=c("1","2","1","2","1","3","2","3","2","1","2","1","2","1","3"), #control
                  type=c(1,2,1,2,1,3,2,3,2,1,2,1,2,1,3), #control
                  type.label=Type,
                  size=Parpop2#choice of scenarios
)
mflow3$Var2=as.character(mflow3$Var2)
mflow3$Var1=as.character(mflow3$Var1)

mflow3$Var2[which(mflow3$Var2=="Leff")]<-"p1" ;mflow3$Var1[which(mflow3$Var1=="Leff")]<-"p1"
mflow3$Var2[which(mflow3$Var2=="Trieux")]<-"p2" ;mflow3$Var1[which(mflow3$Var1=="Trieux")]<-"p2"
mflow3$Var2[which(mflow3$Var2=="Jaudy")]<-"p3" ;mflow3$Var1[which(mflow3$Var1=="Jaudy")]<-"p3"
mflow3$Var2[which(mflow3$Var2=="Leguer")]<-"p4" ;mflow3$Var1[which(mflow3$Var1=="Leguer")]<-"p4"
mflow3$Var2[which(mflow3$Var2=="Yar")]<-"p5" ;mflow3$Var1[which(mflow3$Var1=="Yar")]<-"p5"
mflow3$Var2[which(mflow3$Var2=="Douron")]<-"p6" ;mflow3$Var1[which(mflow3$Var1=="Douron")]<-"p6"
mflow3$Var2[which(mflow3$Var2=="Penze")]<-"p7" ;mflow3$Var1[which(mflow3$Var1=="Penze")]<-"p7"
mflow3$Var2[which(mflow3$Var2=="Elorn")]<-"p8" ;mflow3$Var1[which(mflow3$Var1=="Elorn")]<-"p8"
mflow3$Var2[which(mflow3$Var2=="Aulne")]<-"p9" ;mflow3$Var1[which(mflow3$Var1=="Aulne")]<-"p9"
mflow3$Var2[which(mflow3$Var2=="Goyen")]<-"p10" ;mflow3$Var1[which(mflow3$Var1=="Goyen")]<-"p10"
mflow3$Var2[which(mflow3$Var2=="Odet")]<-"p11" ;mflow3$Var1[which(mflow3$Var1=="Odet")]<-"p11"
mflow3$Var2[which(mflow3$Var2=="Aven")]<-"p12" ;mflow3$Var1[which(mflow3$Var1=="Aven")]<-"p12"
mflow3$Var2[which(mflow3$Var2=="Laita")]<-"p13" ;mflow3$Var1[which(mflow3$Var1=="Laita")]<-"p13"
mflow3$Var2[which(mflow3$Var2=="Scorff")]<-"p14" ;mflow3$Var1[which(mflow3$Var1=="Scorff")]<-"p14"
mflow3$Var2[which(mflow3$Var2=="Blavet")]<-"p15" ;mflow3$Var1[which(mflow3$Var1=="Blavet")]<-"p15"

links<-mflow3
colnames(links)<-c("from","to","weigth")
links<-links[-which(links$weigth==0),]
vis.nodes<-nodes
vis.links<-links

vis.nodes$shape  <- "dot"  
vis.nodes$shadow <- TRUE # Nodes will drop shadow
#vis.nodes$title  <- nodes$pop # Text on click
vis.nodes$label  <- nodes$pop # Node label
vis.nodes$size   <- nodes$size/5 # Node size
vis.nodes$borderWidth <- 2 # Node border width
vis.nodes$color.background <- c("mediumseagreen", "orange", "grey")[nodes$type]
vis.nodes$color.border <- c("mediumseagreen", "orange", "grey")[nodes$type]
vis.nodes$color.highlight.background <- "red"
vis.nodes$color.highlight.border <- "gold"
tryagain<-links$weigth*1.2 #+1 for absolute values #*10 for proportion
vis.links$width <- tryagain # line width
vis.links$arrows <- "to" # arrows: 'from', 'to', or 'middle'
vis.links$smooth <- TRUE    # should the edges be curved?
vis.links$shadow <- FALSE    # edge shadow

#nodes (pop) geographical coordinates
lat<-c(48.701791,48.677231, 48.714110, 48.651108, 48.646443, 48.636448, 48.600963, 48.479118, 48.204225, 48.040924, 48.002719, 47.869874, 47.870947, 47.859305, 47.834972)
lon<-c(-3.056085, -3.157408 , -3.262510, -3.420931, -3.577629, -3.659678, -3.937004, -4.191071, -4.050434,-4.478585, -4.111829, -3.726740,  -3.545196, -3.400563, -3.207688)
#plot(lon, lat)
vis.nodes$x<- lon*1000
vis.nodes$y <- -lat*1000

library(visNetwork)
a<-visNetwork(vis.nodes, vis.links)
a<-visEdges(a, arrows=list(to=list(enable=T, scaleFactor=1.5)),color = list(color = "black", highlight = "red"), smooth = list(enabled = TRUE, type = "diagonalCross"))
a<-visNodes(a,fixed = TRUE,physics=T, font=list(color="black", size=0))
a<-visOptions(a, highlightNearest = TRUE, selectedBy = "type.label")
a<-visLegend(a,addNodes = data.frame(label = c("Sink","Source","Neutral"), shape = c("dot"), 
                                     size = c(10,10,10), 
                                     color = c("mediumseagreen","orange","grey"),font.size=20),
             addEdges = data.frame(label = "Emigration",color="black"), useGroups = FALSE)
a

min(Parpop2)
max(Parpop2)
min(mflow3$value)
max(mflow3$value)
