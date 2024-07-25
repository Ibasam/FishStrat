## remove (almost) everything in the working environment.
rm(list = ls())
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  args <- 1
}


options(tidyverse.quiet = TRUE)

dir_results <- "results/"

#_________________ PACKAGES _________________#
library(fst)
library(tidyverse,warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

source("R/pe_mv.R")

CV <- function(x) (sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE))

#_________________ PARAMETERS _________________#
nSIMUL <- 50
nYear <- 40
nInit <- 10 # Nb years at initialization
npop <- 15

Rmax=rep(10,npop)
Rmax5=Rmax*5/100
Rmax20=Rmax*20/100
theta=Rmax5/100
theta20=Rmax20/100
dat <- read.csv2("data/dataPop.csv", header=TRUE, stringsAsFactors = FALSE, comment.char = "#")
Area=dat$Area
rPROP=.25

#_________________ SCENARIOS _________________#
iEXPE <- args

scn_id <- as.numeric(strsplit(as.character(iEXPE), "")[[1]])
scenarioConnect = scn_id[1]
scenarioFishing = scn_id[2]
scenarioFishingRate = scn_id[4]


#_________________ DATA _________________#

cat("Composition for ", iEXPE, "\n")

Nparr0 <-Nparr0.mature<- Nsmolt0 <- NRet <- NMSW <- N1SW <- matrix(NA,nYear+nInit,nSIMUL)
LFparr0 <- LFsmolt0 <- LFreturns <-LFMSW <-LF1SW<-RatioMSW<- matrix(NA,nYear+nInit,nSIMUL)
gNeutral <- gFmid1 <- gFmid2 <-gFmid3 <-gFmid4<-gG<- matrix(NA,nYear+nInit,nSIMUL)

CVparr0 <- CVsmolt0 <- CVreturns <-CVMSW <-CV1SW<- matrix(NA,nYear+nInit,nSIMUL)
CVgNeutral <- CVgFmid1 <- CVgFmid2 <-CVgFmid3 <-CVgFmid4<-CVgG<- matrix(NA,nYear+nInit,nSIMUL)
SDgNeutral <- SDgFmid1 <- SDgFmid2 <-SDgFmid3 <-SDgFmid4<-SDgG<- matrix(NA,nYear+nInit,nSIMUL)

results=NULL

Pexpl.global=Ncapt.global=Ncapt.MSW=Ncapt.1SW=Ncapt.global.5=Ncapt.MSW.5=Ncapt.1SW.5=NULL
PE_mv=list()
RETURNS=list()

PQext <- list() 
Ext <- array(,dim=c(nSIMUL,npop))

for (sim in 1:nSIMUL){ # loop over simulations
  
  cat(iEXPE,"-Simulation :",sim," / ")
  perc <- (sim/nSIMUL)*100
  if(perc %in% (seq(0,90,10))) cat(perc,"% ~ ","\n")
  if(perc == 100) cat(perc,"%","\n")
  
  tmp=tmp1=tmp2=tmp3=NULL
  NIm <- Exploitation <- matrix(0,nYear+nInit,npop)
  
  df=NULL
  df <- read.fst(paste0(dir_results,"results_FishStrat_",iEXPE,"/Metapop_Sc",iEXPE,"_Sim",sim,".fst"))
  load(paste0(dir_results,"results_FishStrat_",iEXPE,"/Fishing_rates_Sc",iEXPE,"_Sim",sim,".RData"))

  nyears=nYear+nInit
  df <- subset(df,year<=nyears)
  
  
  ##------- PARR 0+ -------##
  tmp<-subset(df,Parr==1 & date==273 & AgeRiver<1)
  tmp1 <- as.vector(by(tmp$Lf, tmp$year, length)); 
  tmp2 <- as.vector(by(tmp$Lf, tmp$year, mean ,na.rm=TRUE))
  tmp3 <- as.vector(by(tmp$Lf, tmp$year, CV))
  
  Nparr0[1:length(tmp1),sim] <- tmp1
  LFparr0[1:length(tmp2),sim] <- tmp2
  CVparr0[1:length(tmp3),sim] <- tmp3
  
  tmp<-subset(tmp,Mature==1)
  tmp4<-as.vector(by(tmp$Lf, tmp$year, length));
  Nparr0.mature[1:length(tmp4),sim] <- tmp4
  
  ## PARR MALE
  tmp <-subset(df,Parr==1 & date==273 & AgeRiver<1 & Female==0 )
  tmp2 <- as.vector(by(tmp$gFmid1, tmp$year, mean ,na.rm=TRUE))
  tmp3 <- as.vector(by(tmp$gFmid1, tmp$year, CV))
  tmp4 <- as.vector(by(tmp$gFmid1, tmp$year, sd))
  gFmid1[1:length(tmp2),sim] <- tmp2
  CVgFmid1[1:length(tmp3),sim] <- tmp3
  SDgFmid1[1:length(tmp4),sim] <- tmp4
  
  ## PARR FEMALE
  tmp <-subset(df,Parr==1 & date==273 & AgeRiver<1 & Female==1 )
  tmp2 <- as.vector(by(tmp$gFmid2, tmp$year, mean ,na.rm=TRUE))
  tmp3 <- as.vector(by(tmp$gFmid2, tmp$year, CV))
  tmp4 <- as.vector(by(tmp$gFmid2, tmp$year, sd))
  gFmid2[1:length(tmp2),sim] <- tmp2
  CVgFmid2[1:length(tmp3),sim] <- tmp3
  SDgFmid2[1:length(tmp4),sim] <- tmp4
  
  
  ##--------- SMOLT 0+ ----------##
  tmp<-subset(df,Smolt==1 & date==90 & AgeRiver==1)
  tmp1 <- as.vector(by(tmp$Lf, tmp$year, length)); 
  tmp2 <- as.vector(by(tmp$Lf, tmp$year, mean ,na.rm=TRUE))
  tmp3 <- as.vector(by(tmp$Lf, tmp$year, CV))
  Nsmolt0[1:length(tmp1),sim] <- tmp1
  LFsmolt0[1:length(tmp2),sim] <- tmp2
  CVsmolt0[1:length(tmp3),sim] <- tmp3
  
  
  ##----------- ANADROMOUS (ALL RETURNS) ----------##
  tmp<-subset(df,Returns>0 & Atsea==0 & Mature==1 & date==273)
  tmp1 <- as.vector(by(tmp$Lf, tmp$year, length)); 
  NRet[1:length(tmp1),sim] <- tmp1
  
  tmp1=NULL
  tmp1 <-
    tmp %>%
    group_by(Pop, year) %>%
    summarize(count = n()) %>%
    ungroup() %>%
    complete(Pop, year, fill = list(count = 0))
  
  tmp1 <-tmp1[order(tmp1$year, tmp1$Pop),]

  Nreturns.matrix <- matrix(tmp1$count, max(tmp1$year), max(tmp1$Pop),byrow=TRUE)
  PE_mv[[sim]]<-pe_mv(Nreturns.matrix, type = "loess_detrended")
  PE_mv[[sim]]$Synchrony<-synchrony(Nreturns.matrix)
  PE_mv[[sim]]$Returns <- Nreturns.matrix
  
  
  tmp2<-apply(Nreturns.matrix,1,CV)
  CVreturns[1:length(tmp2),sim] <- tmp2
  
  ## MSW
  tmp<-subset(df,Returns>0 & Atsea==0 & Mature==1 & date==273 & AgeSea >=2)
  tmp1 <- as.vector(by(tmp$Lf, tmp$year, length)); 
  tmp2 <- as.vector(by(tmp$Lf, tmp$year, mean ,na.rm=TRUE))
  tmp3 <- as.vector(by(tmp$Lf, tmp$year, CV))
  NMSW[1:length(tmp1),sim] <- tmp1
  LFMSW[1:length(tmp2),sim] <- tmp2
  CVMSW[1:length(tmp3),sim] <- tmp3
  
  tmp1_MSW=NULL
  tmp1_MSW <-
    tmp %>%
    group_by(Pop, year) %>%
    summarize(count = n()) %>%
    ungroup() %>%
    complete(Pop, year, fill = list(count = 0))
  
  tmp1_MSW <-tmp1_MSW[order(tmp1_MSW$year, tmp1_MSW$Pop),]
  if (any (tmp1_MSW$year == 1)) {
  } else {
    year1 <- as.data.frame(cbind(c(1:15),rep(1,15),rep(0,15)))
    colnames(year1) <- c("Pop","year","count")
    tmp1_MSW <- rbind(year1, tmp1_MSW)
  }
  tmp1_MSW$ExplRate <- c(rep(fishing.rates[,4],nyears)) 
  tmp1_MSW$nFish <- ceiling(tmp1_MSW$count*tmp1_MSW$ExplRate)
  
  
  ## 1SW
  tmp<-subset(df,Returns>0 & Atsea==0 & Mature==1 & date==273 & AgeSea <2)
  tmp1 <- as.vector(by(tmp$Lf, tmp$year, length)); 
  tmp2 <- as.vector(by(tmp$Lf, tmp$year, mean ,na.rm=TRUE))
  tmp3 <- as.vector(by(tmp$Lf, tmp$year, CV))
  N1SW[1:length(tmp1),sim] <- tmp1
  LF1SW[1:length(tmp2),sim] <- tmp2
  CV1SW[1:length(tmp3),sim] <- tmp3
  
  tmp1_1SW=NULL
  tmp1_1SW <-
    tmp %>%
    group_by(Pop, year) %>%
    summarize(count = n()) %>%
    ungroup() %>%
    complete(Pop, year, fill = list(count = 0))
  
  tmp1_1SW <-tmp1_1SW[order(tmp1_1SW$year, tmp1_1SW$Pop),]
  tmp1_1SW$ExplRate <- c(rep(fishing.rates[,2],nyears))
  tmp1_1SW$nFish <- ceiling(tmp1_1SW$count*tmp1_1SW$ExplRate)
  
  
  tot_MSW <-aggregate(tmp1_MSW$count, by = list(year=tmp1_MSW$year), FUN = sum)
  tot_1SW <-aggregate(tmp1_1SW$count, by = list(year=tmp1_1SW$year), FUN = sum)
  tot <- tot_MSW + tot_1SW
  capt_MSW <- aggregate(tmp1_MSW$nFish, by = list(year=tmp1_MSW$year), FUN = sum)
  capt_1SW <- aggregate(tmp1_1SW$nFish, by = list(year=tmp1_1SW$year), FUN = sum)
  capt <- capt_MSW + capt_1SW
  Pexpl.global[sim] <- mean(capt$x/tot$x,na.rm=TRUE)
  #nCaptures
  Ncapt.global[sim] <- mean(capt$x,na.rm=TRUE)
  Ncapt.MSW[sim] <- mean(capt_MSW$x,na.rm=TRUE)
  Ncapt.1SW[sim] <- mean(capt_1SW$x,na.rm=TRUE)
  #nCaptures last 5 years
  Ncapt.global.5[sim] <- mean(capt$x[45:50],na.rm=TRUE)
  Ncapt.MSW.5[sim] <- mean(capt_MSW$x[45:50],na.rm=TRUE)
  Ncapt.1SW.5[sim] <- mean(capt_1SW$x[45:50],na.rm=TRUE)
  
  
  ## PropMSW
  RatioMSW[,sim] <- NMSW[,sim] / N1SW[,sim]
  
  
  ##-------------- ANADROMOUS (HOMERS) ------------##
  tmp<-subset(df,Returns>0 & Atsea==0 & Mature==1 & date==273 & CollecID == Pop)
  
  tmp2 <- as.vector(by(tmp$Lf, tmp$year, mean ,na.rm=TRUE))
  tmp3 <- as.vector(by(tmp$Lf, tmp$year, CV))
  LFreturns[1:length(tmp2),sim] <- tmp2
  #CVreturns[1:length(tmp3),sim] <- tmp3
  
  # Neutral
  tmp2 <- as.vector(by(tmp$gNeutral, tmp$year, mean ,na.rm=TRUE))
  tmp3 <- as.vector(by(tmp$gNeutral, tmp$year, CV))
  tmp4 <- as.vector(by(tmp$gNeutral, tmp$year, sd))
  gNeutral[1:length(tmp2),sim] <- tmp2
  CVgNeutral[1:length(tmp3),sim] <- tmp3
  SDgNeutral[1:length(tmp4),sim] <- tmp4
  
  # Growth potential
  tmp2 <- as.vector(by(tmp$gG, tmp$year, mean ,na.rm=TRUE))
  tmp3 <- as.vector(by(tmp$gG, tmp$year, CV))
  tmp4 <- as.vector(by(tmp$gG, tmp$year, sd))
  gG[1:length(tmp2),sim] <- tmp2
  CVgG[1:length(tmp3),sim] <- tmp3
  SDgG[1:length(tmp4),sim] <- tmp4
  
  ## ANADROMOUS MALE (HOMERS)
  tmp <-subset(df,Returns>0 & Atsea==0 & Mature==1 & date==273 & Female==0 & CollecID == Pop)
  tmp2 <- as.vector(by(tmp$gFmid3, tmp$year, mean ,na.rm=TRUE))
  tmp3 <- as.vector(by(tmp$gFmid3, tmp$year, CV))
  tmp4 <- as.vector(by(tmp$gFmid3, tmp$year, sd))
  gFmid3[1:length(tmp2),sim] <- tmp2
  CVgFmid3[1:length(tmp3),sim] <- tmp3
  SDgFmid3[1:length(tmp4),sim] <- tmp4
  
  ## ANADROMOUS FEMALE (HOMERS)
  tmp <-subset(df,Returns>0 & Atsea==0 & Mature==1 & date==273 & Female==1 & CollecID == Pop)
  tmp2 <- as.vector(by(tmp$gFmid4, tmp$year, mean ,na.rm=TRUE))
  tmp3 <- as.vector(by(tmp$gFmid4, tmp$year, CV))
  tmp4 <- as.vector(by(tmp$gFmid4, tmp$year, sd))
  gFmid4[1:length(tmp2),sim] <- tmp2
  CVgFmid4[1:length(tmp3),sim] <- tmp3
  SDgFmid4[1:length(tmp4),sim] <- tmp4

  
  ##----------- ANADROMOUS (IMMIGRANTS) ---------##
  tmp_mig <- NULL
  Im.table <- array(,dim=c(nyears,npop))

for (pop in 1:npop){ # loop over popualtions
  for (i in 1:nyears){
      tmp5 <-subset(tmp,Pop==pop & CollecID!=pop & year==i)
      NIm[i,pop] <- nrow(tmp5)
      
      # Nb immigrants by population of origin
      for (pop2 in 1:npop){
        if (pop2 != pop){
          Im.table[i,pop2] <- nrow(subset(tmp, Pop==pop & CollecID==pop2 & year==i))
        } else {Im.table[i,pop2] <- 0}
      }
      
      ## FISHING
      if (i<=nInit){
        Exploitation[i,pop] <- Nreturns.matrix[i,pop]*fishing.rates[pop,1]
      }else{
        Exploitation[i,pop] <- Nreturns.matrix[i,pop]*fishing.rates[pop,2]
      }
    } # end loop years
  tmp_mig[[pop]] <- Im.table
  
  } # end loop pop
  RETURNS[[sim]] <- list()
  RETURNS[[sim]]$Returns <- Nreturns.matrix
  RETURNS[[sim]]$Immigrants <- NIm
  RETURNS[[sim]]$Im <- tmp_mig
  RETURNS[[sim]]$Captured <- Exploitation
  
  
  ##------- POPULATIONS QUASI-EXTINCTION RISK -------##
  
  for (pop in 1:npop){ # loop over populations
    
    demo <- NULL
    demo <- subset(df, Pop==pop)
    nyears <- max(demo$year,na.rm=TRUE)-1

    ## PARR 0+
    parr0<-subset(demo,Parr==1 & date==273 & AgeRiver<1)
    nparr0 <- as.vector(table(factor(parr0$year, levels=1:50)))
    density <- nparr0 / (Area[pop]*rPROP)
    
    c=0
    for (i in 1:(nyears-1)){
      
      if (density[i] < theta20[pop] && density[i+1] < theta20[pop]) {
        c=c+1
      }
    }
    if (c >= 1) {
      Ext[repmod,pop]<-1
    }
    else {
      Ext[repmod,pop]<-0
    }
    
  } # end loop population
  
  
  
} # end loop sim

PQext <- colSums(Ext[,])/nSIMUL

### Save results

save(RETURNS,PE_mv,Pexpl.global,Ncapt.global,Ncapt.MSW,Ncapt.1SW,Ncapt.global.5,Ncapt.MSW.5,Ncapt.1SW.5,
     Nparr0,Nparr0.mature,Nsmolt0,NRet, NMSW, N1SW
     ,CVreturns
     ,RatioMSW,LFparr0,LFsmolt0,LFreturns,LFMSW,LF1SW
     ,CVparr0,CVsmolt0,CVMSW,CV1SW
     
     ,gNeutral,gFmid1,gFmid2,gFmid3,gFmid4,gG
     ,CVgNeutral,CVgFmid1,CVgFmid2,CVgFmid3,CVgFmid4,CVgG
     ,SDgNeutral,SDgFmid1,SDgFmid2,SDgFmid3,SDgFmid4,SDgG
     
     ,PQext
      ,file=paste0("results/METAPOP",iEXPE,".RData"))

q('no')
