## remove (almost) everything in the working environment.
rm(list = ls())
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  args <- 1
}


options(tidyverse.quiet = TRUE)

#dir_results <- "/media/ssd2To/FishStrat/results/"
dir_results <- "/media/hdd8To/FishStrat/results/"

#_________________ PACKAGES _________________#
library(fst)
library(tidyverse,warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

source("R/pe_mv.R")


# 'data.frame':	1569975 obs. of  45 variables:
#   $ Lf         : num  724 808 132 144 138 ...
# $ W          : num  4049 5177.1 21.8 28.4 25.2 ...
# $ Fat        : num  627.24 425.21 1.02 1.52 1.37 ...
# $ AgeRiver   : num  1 1 1 1 1 1 1 1 1 1 ...
# $ AgeSea     : num  2.37 2.37 0 0 0 ...
# $ Female     : num  1 1 1 0 1 1 1 1 1 1 ...
# $ Parr       : num  0 0 0 0 0 0 0 0 0 0 ...
# $ Mature     : num  0 0 0 0 0 0 0 0 0 0 ...
# $ Atsea      : num  1 1 0 0 0 0 0 0 0 0 ...
# $ Smolt      : num  0 0 1 1 1 1 1 1 1 1 ...
# $ ID         : num  20960513 20956213 817413 823413 835713 ...
# $ DW         : num  10.32 7.98 6.66 9.04 7.53 ...
# $ Returns    : num  1 1 0 0 0 0 0 0 0 0 ...
# $ Psurvival  : num  0.848 0.875 0.694 0.614 0.684 ...
# $ date       : num  90 90 90 90 90 90 90 90 90 90 ...
# $ year       : num  1 1 1 1 1 1 1 1 1 1 ...
# $ CollecID   : num  13 13 13 13 13 13 13 13 13 13 ...
# $ Wini       : num  0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 ...
# $ SpecificGR : num  104 102 104 104 108 ...
# $ motherID   : num  0 0 0 0 0 0 0 0 0 0 ...
# $ fatherID   : num  0 0 0 0 0 0 0 0 0 0 ...
# $ gG         : num  -0.071 0.0473 0.0237 0.0237 -0.0473 ...
# $ gG_sea     : num  -0.1183 -0.1893 0 -0.0473 0.071 ...
# $ gPercF     : num  0.122 0.132 0.137 0.118 0.113 ...
# $ pG         : num  -0.161103 0.000222 -0.149366 -0.020968 -0.130997 ...
# $ pG_sea     : num  0.1001 0.0043 0.119 0.1898 -0.1325 ...
# $ pPercF     : num  0.122 0.121 0.103 0.121 0.117 ...
# $ gSLmid     : num  89 89 89 89 89 ...
# $ galphaS    : num  0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 ...
# $ pSLmid     : num  89 89 89 89 89 ...
# $ palphaS    : num  0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 ...
# $ gFmid1     : num  1.25 0.73 1.51 1.77 0.99 1.12 1.25 1.38 1.77 2.03 ...
# $ pFmid1     : num  1.425 1.305 0.894 1.694 1.307 ...
# $ gFmid2     : num  42 40 36 43 38 49 46 37 40 40 ...
# $ pFmid2     : num  41.7 46.7 36.3 41.6 44.6 ...
# $ gFmid3     : num  64 48 40 104 48 32 -8 80 32 16 ...
# $ pFmid3     : num  33.42 51 2.43 130.57 43.62 ...
# $ gFmid4     : num  -50 145 85 10 145 160 85 100 100 40 ...
# $ pFmid4     : num  -120.8 240.6 46.6 -37.7 209.4 ...
# $ gNeutral   : num  0.475 0.4 0.45 0.475 0.4 0.475 0.525 0.725 0.475 0.525 ...
# $ motherStrat: num  0 0 0 0 0 0 0 0 0 0 ...
# $ fatherStrat: num  0 0 0 0 0 0 0 0 0 0 ...
# $ NreproYear : num  1 1 0 0 0 0 0 0 0 0 ...
# $ Nrepro     : num  1 1 0 0 0 0 0 0 0 0 ...
# $ Pop        : int  13 13 13 13 13 13 13 13 13 13 ...



#_________________ FUNCTIONS_________________#
# proportions.population <- function (population, window = TRUE, plotting = TRUE, titles = "") {
#   
#   parr0 <- (population$Parr==1 & population$AgeRiver<1)
#   parr1 <- (population$Parr==1 & population$AgeRiver>1)
#   parr0.mature <- (population$Parr==1 & population$AgeRiver<1 & population$Mature==1)
#   parr1.mature <- (population$Parr==1 & population$AgeRiver>1 & population$Mature==1)
#   grilses <- (population$Returns > 0 & population$AgeSea < 2)
#   MSW <- (population$Returns > 0 & population$AgeSea >= 2)
#   
#   parr1ratio <- sum(parr0)/sum(parr1)
#   OnevsMSWratio <- sum(grilses)/sum(MSW)
#   
#   sexratioParr <- sum(parr1 & population$Female == 0)/sum(parr1 & population$Female == 1)
#   sexratioGrilses <- sum(grilses & population$Female == 0)/sum(grilses & population$Female == 1)
#   sexratioMSW <- sum(MSW & population$Female == 0)/sum(MSW &  population$Female == 1)
#   
#   Nmale.mature <- sum(parr0.mature & population$Female == 0) + sum(parr1.mature & population$Female == 0) +  sum(grilses & population$Female == 0) + sum(MSW & population$Female == 0)
#   NFemale.mature <- sum(parr0.mature & population$Female == 1) + sum(parr1.mature & population$Female == 1) +  sum(grilses & population$Female == 1) + sum(MSW & population$Female == 1)
#   
#   return(list(
#     Ovs1parrratio = parr1ratio
#     , sexratioParr = sexratioParr
#     , OnevsMSWratio = OnevsMSWratio
#     , sexratioGrilses = sexratioGrilses
#     , sexratioMSW = sexratioMSW
#     , nReturns = sum(grilses) +  sum(MSW)
#     , OSR=Nmale.mature / NFemale.mature))
# }
# #source("code/pe_mv.R")

CV <- function(x) (sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE))

#_________________ PARAMETERS _________________#
#source("/media/hdd4To/mbuoro/MetaIBASAM-Projects/FishStrat/dataIbasam.R")
#source("/media/hdd4To/mbuoro/MetaIBASAM-Projects/FishStrat/parIbasam.R")

nSIMUL <- 50
nYear <- 40
nInit <- 10 # Nb years at initialization
npop <- 15
#popo=15


# #h = c(1, .95, .9, .85, .8, .75, .7) # Philopatry (homing) rates #edit al - 22/03/21 5% and 15%
# scenarioConnect=7 #scenario 1 for h=1.00, scenario 2 for h=0.95, scenario 3 pour h=0.80
# 
# # 0: control / 
# # 1: no fishing on sink populations 
# # 2: no fishing on source populations 
# scenarioFishing = 0
# 
# #frates_vec =c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
# scenarioFishingRate = 1
# 
# 
# Sc_name <- paste0("FishStrat_",scenarioConnect,scenarioFishing,scenarioFishingRate)



# SCENARIOS
iEXPE <- args
#EXPE <- c(110, 210,220)#,120,123,124, 310, 320, 323, 324)
#EXPE <- c(120, 121, 122, 123, 124, 320, 321, 322, 323, 324)

scn_id <- as.numeric(strsplit(as.character(iEXPE), "")[[1]])
scenarioConnect = scn_id[1]
scenarioFishing = scn_id[2]
scenarioFishingRate = scn_id[3]

#for (iEXPE in EXPE){ # Loop over scenario 

  
  #_________________ DATA _________________#
  #fishing.rates = scn = NULL
  #scn <- as.numeric(strsplit(as.character(iEXPE), "")[[1]])
  #scenarioConnect = scn[1]
  #scenarioEnvi = scn[2]
  #scenarioFishing = scn[3]
  

  cat("Composition for ", iEXPE, "\n")
  
  #pe=petr=sync=CV_est=CV_obs=CV_esttr=CV_obstr=NULL
  Nparr0 <-Nparr0.mature<- Nsmolt0 <- NRet <- NMSW <- N1SW <- matrix(NA,nYear+nInit,nSIMUL)
  LFparr0 <- LFsmolt0 <- LFreturns <-LFMSW <-LF1SW<-RatioMSW<- matrix(NA,nYear+nInit,nSIMUL)
  gNeutral <- gFmid1 <- gFmid2 <-gFmid3 <-gFmid4<-gG<- matrix(NA,nYear+nInit,nSIMUL)
  
  CVreturns<-CVparr0 <- CVsmolt0 <- CVreturns <-CVMSW <-CV1SW<- matrix(NA,nYear+nInit,nSIMUL)
  CVgNeutral <- CVgFmid1 <- CVgFmid2 <-CVgFmid3 <-CVgFmid4<-CVgG<- matrix(NA,nYear+nInit,nSIMUL)
  
  fishing_rates=NULL
  results=NULL
  
  Pexpl.global=Sync=NULL
  PE_mv=list()
  RETURNS=list()
  
  for (sim in 1:nSIMUL){ # loop over simulations
    #cat(indice,"/",sim,"- ")
    cat(iEXPE,"-Simulation :",sim," / ")
    perc <- (sim/nSIMUL)*100
    if(perc %in% (seq(0,90,10))) cat(perc,"% ~ ","\n")
    if(perc == 100) cat(perc,"%","\n")
    
    tmp=tmp1=tmp2=tmp3=NULL
    NIm <- Exploitation <- matrix(0,nYear+nInit,npop)
    
    #load(paste0(dir_results,"results_FishStrat_",iEXPE,"/Metapop_Sc",iEXPE,"_Sim",sim,".RData"))
    df=NULL
    df <- read.fst(paste0(dir_results,"results_FishStrat_",iEXPE,"/Metapop_Sc",iEXPE,"_Sim",sim,".fst"))
    load(paste0(dir_results,"results_FishStrat_",iEXPE,"/Fishing_rates_Sc",iEXPE,"_Sim",sim,".RData"))
    
    #nyears <- max(df$year,na.rm=TRUE)-1
    nyears=nYear+nInit
    df <- subset(df,year<=nyears)
    
    
    ## PARR 0+
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
    
    ## SMOLT 0+
    tmp<-subset(df,Smolt==1 & date==90 & AgeRiver==1)
    tmp1 <- as.vector(by(tmp$Lf, tmp$year, length)); 
    tmp2 <- as.vector(by(tmp$Lf, tmp$year, mean ,na.rm=TRUE))
    tmp3 <- as.vector(by(tmp$Lf, tmp$year, CV))
    Nsmolt0[1:length(tmp1),sim] <- tmp1
    LFsmolt0[1:length(tmp2),sim] <- tmp2
    CVsmolt0[1:length(tmp3),sim] <- tmp3
    
    
    

    ## ANADROMOUS
    tmp<-subset(df,Returns>0 & Atsea==0 & Mature==1 & date==273)
    tmp1 <- as.vector(by(tmp$Lf, tmp$year, length)); 
    tmp2 <- as.vector(by(tmp$Lf, tmp$year, mean ,na.rm=TRUE))
    tmp3 <- as.vector(by(tmp$Lf, tmp$year, CV))
    NRet[1:length(tmp1),sim] <- tmp1
    LFreturns[1:length(tmp2),sim] <- tmp2
    #CVreturns[1:length(tmp3),sim] <- tmp3
    

    
    tmp1=NULL
    tmp1 <-
      tmp %>%
      group_by(Pop, year) %>%
      summarize(count = n()) %>%
      ungroup() %>%
      complete(Pop, year, fill = list(count = 0))

    tmp1 <-tmp1[order(tmp1$year, tmp1$Pop),]
    tmp1$ExplRate <- c(rep(fishing.rates[,2],nyears))
    tmp1$nFish <- ceiling(tmp1$count*tmp1$ExplRate)

    Nreturns.matrix <- matrix(tmp1$count, max(tmp1$year), max(tmp1$Pop),byrow=TRUE)
    PE_mv[[sim]]<-pe_mv(Nreturns.matrix, type = "loess_detrended")
    PE_mv[[sim]]$Synchrony<-synchrony(Nreturns.matrix)
    PE_mv[[sim]]$Returns <- Nreturns.matrix

    
    tmp2<-apply(Nreturns.matrix,1,CV)
    CVreturns[1:length(tmp2),sim] <- tmp2
      
    tot <-aggregate(tmp1$count, by = list(year=tmp1$year), FUN = sum)
    capt <- aggregate(tmp1$nFish, by = list(year=tmp1$year), FUN = sum)
    Pexpl.global[sim] <- mean(capt$x/tot$x,na.rm=TRUE)
    
    
    ## HOMERS
    #tmp<-subset(demo,Returns>0 & Atsea==0 & Mature==1 & date==273 & CollecID == pop & year==i)
    #Nhomers[i,pop] <- nrow(tmp)
    
    ## IMMIGRANTS
    for (i in 1:nyears){
      for (pop in 1:npop){ # loop over popualtions
        tmp5 <-subset(tmp,Pop==pop & CollecID!=pop & year==i)
        NIm[i,pop] <- nrow(tmp5)
        
        
        ## FISHING
        if (i<=nInit){
          Exploitation[i,pop] <- Nreturns.matrix[i,pop]*fishing.rates[pop,1]
        }else{
          Exploitation[i,pop] <- Nreturns.matrix[i,pop]*fishing.rates[pop,2]
        }
      } # end loop pop
    } # end loop years
    RETURNS[[sim]] <- list()
    RETURNS[[sim]]$Returns <- Nreturns.matrix
    RETURNS[[sim]]$Immigrants <- NIm
    RETURNS[[sim]]$Captured <- Exploitation
    
    ## Neutral
    #tmp<-subset(df,Returns>0 & Atsea==0 & Mature==1 & date==273)
    tmp2 <- as.vector(by(tmp$gNeutral, tmp$year, mean ,na.rm=TRUE))
    tmp3 <- as.vector(by(tmp$gNeutral, tmp$year, CV))
    gNeutral[1:length(tmp2),sim] <- tmp2
    CVgNeutral[1:length(tmp3),sim] <- tmp3
    
    ## Growth potential
    #tmp<-subset(df,Returns>0 & Atsea==0 & Mature==1 & date==273)
    tmp2 <- as.vector(by(tmp$gG, tmp$year, mean ,na.rm=TRUE))
    tmp3 <- as.vector(by(tmp$gG, tmp$year, CV))
    gG[1:length(tmp2),sim] <- tmp2
    CVgG[1:length(tmp3),sim] <- tmp3
    
    
    
    
    
    
    ## MSW
    tmp<-subset(df,Returns>0 & Atsea==0 & Mature==1 & date==273 & AgeSea >=2)
    tmp1 <- as.vector(by(tmp$Lf, tmp$year, length)); 
    tmp2 <- as.vector(by(tmp$Lf, tmp$year, mean ,na.rm=TRUE))
    tmp3 <- as.vector(by(tmp$Lf, tmp$year, CV))
    NMSW[1:length(tmp1),sim] <- tmp1
    LFMSW[1:length(tmp2),sim] <- tmp2
    CVMSW[1:length(tmp3),sim] <- tmp3
    
    ## 1SW
    tmp<-subset(df,Returns>0 & Atsea==0 & Mature==1 & date==273 & AgeSea <2)
    tmp1 <- as.vector(by(tmp$Lf, tmp$year, length)); 
    tmp2 <- as.vector(by(tmp$Lf, tmp$year, mean ,na.rm=TRUE))
    tmp3 <- as.vector(by(tmp$Lf, tmp$year, CV))
    N1SW[1:length(tmp1),sim] <- tmp1
    LF1SW[1:length(tmp2),sim] <- tmp2
    CV1SW[1:length(tmp3),sim] <- tmp3
    
    ## PropMSW
    RatioMSW[,sim] <- NMSW[,sim] / N1SW[,sim]
    
    
      ## PARR MALE
      tmp <-subset(df,Parr==1 & date==273 & AgeRiver<1 & Female==0 )
      tmp2 <- as.vector(by(tmp$gFmid1, tmp$year, mean ,na.rm=TRUE))
      tmp3 <- as.vector(by(tmp$gFmid1, tmp$year, CV))
      gFmid1[1:length(tmp2),sim] <- tmp2
      CVgFmid1[1:length(tmp3),sim] <- tmp3

      ## PARR FEMALE
      tmp <-subset(df,Parr==1 & date==273 & AgeRiver<1 & Female==1 )
      tmp2 <- as.vector(by(tmp$gFmid2, tmp$year, mean ,na.rm=TRUE))
      tmp3 <- as.vector(by(tmp$gFmid2, tmp$year, CV))
      gFmid2[1:length(tmp2),sim] <- tmp2
      CVgFmid2[1:length(tmp3),sim] <- tmp3

      ## ANADROMOUS MALE
      tmp <-subset(df,Returns>0 & Atsea==0 & Mature==1 & date==273 & Female==0 )
      tmp2 <- as.vector(by(tmp$gFmid3, tmp$year, mean ,na.rm=TRUE))
      tmp3 <- as.vector(by(tmp$gFmid3, tmp$year, CV))
      gFmid3[1:length(tmp2),sim] <- tmp2
      CVgFmid3[1:length(tmp3),sim] <- tmp3

      ## ANADROMOUS FEMALE
      tmp <-subset(df,Returns>0 & Atsea==0 & Mature==1 & date==273 & Female==1 )
      tmp2 <- as.vector(by(tmp$gFmid4, tmp$year, mean ,na.rm=TRUE))
      tmp3 <- as.vector(by(tmp$gFmid4, tmp$year, CV))
      gFmid4[1:length(tmp2),sim] <- tmp2
      CVgFmid4[1:length(tmp3),sim] <- tmp3
    
  } # end loop sim
    
### Save results

save(RETURNS,PE_mv,Pexpl.global,
     Nparr0,Nparr0.mature,Nsmolt0,NRet, NMSW, N1SW
     ,CVreturns
     ,RatioMSW,LFparr0,LFsmolt0,LFreturns,LFMSW,LF1SW
     ,CVparr0,CVsmolt0,CVMSW,CV1SW
     
     ,gNeutral,gFmid1,gFmid2,gFmid3,gFmid4,gG
     ,CVgNeutral,CVgFmid1,CVgFmid2,CVgFmid3,CVgFmid4,CVgG
     
     ,file=paste0("results/METAPOP/METAPOP",iEXPE,".RData"))

q('no')
