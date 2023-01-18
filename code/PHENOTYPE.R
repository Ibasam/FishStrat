## remove (almost) everything in the working environment.
rm(list = ls())
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  args <- 1
}


#dir_results <- "/media/ssd2To/FishStrat/results/"
dir_results <- "/media/hdd8To/FishStrat/results/"

#_________________ PACKAGES _________________#
library(fst)
#devtools::install_github("ecofolio", username="seananderson")
#library(ecofolio)


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

#CV <- function(x) (sd(x/mean(x)))*100

#_________________ PARAMETERS _________________#
#source("/media/hdd4To/mbuoro/MetaIBASAM-Projects/FishStrat/dataIbasam.R")
#source("/media/hdd4To/mbuoro/MetaIBASAM-Projects/FishStrat/parIbasam.R")

nSIMUL <- 30
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

nParr <- list() # Nb returns
nSmolt <- list() # Nb returns
nReturns <- list() # Nb returns
Mig <- list() # Nb immigrnats
nMSW<- list() # Nb immigrnats
n1SW<- list() # Nb immigrnats
rMSW<- list() # Nb immigrnats
#NRet.res <- list() # Nb returns
#PE <- list()
nExploitation <- list()
LFParr <- list()
LFSmolt <- list()
LFReturns <- list()
LF1sw <- list()
LFmsw <- list()
Fishing_rates <- list()

# SCENARIOS
iEXPE <- args
#EXPE <- c(110, 210,220)#,120,123,124, 310, 320, 323, 324)
#EXPE <- c(120, 121, 122, 123, 124, 320, 321, 322, 323, 324)

#for (iEXPE in EXPE){ # Loop over scenario 

 
  tmp1 <- tmp2 <- tmp3 <- tmp4 <- tmp5 <- tmp6 <- tmp7 <- tmp8 <- tmp9<-tmp10<-tmp11<-tmp12<-tmp13<-list()
  
  #_________________ DATA _________________#
  #fishing.rates = scn = NULL
  #scn <- as.numeric(strsplit(as.character(iEXPE), "")[[1]])
  #scenarioConnect = scn[1]
  #scenarioEnvi = scn[2]
  #scenarioFishing = scn[3]
  

  cat("Composition for ", iEXPE, "\n")
  
  #pe=petr=sync=CV_est=CV_obs=CV_esttr=CV_obstr=NULL
  
  for (repmod in 1:nSIMUL){ # loop over simulations
    #cat(indice,"/",repmod,"- ")
    cat("Simulation :",repmod," / ")
    perc <- (repmod/nSIMUL)*100
    if(perc %in% (seq(0,90,10))) cat(perc,"% ~ ","\n")
    if(perc == 100) cat(perc,"%","\n")
    
    # Loading data files
    #NParr.table=NSmolt.table=NRet.table = NIm.table = Nhomers.table = NULL
    #Npeche.table = NULL
    #Ppeche.table = NULL
    
    LFparr0 <- LFsmolt0 <- LFreturns <-LFMSW <-LF1SW<-NMSW<-N1SW<-RatioMSW<- matrix(NA,nYear+nInit,npop)
    fishing_rates=NULL
    results=NULL
    
    
    #load(paste0(dir_results,"results_FishStrat_",iEXPE,"/Metapop_Sc",iEXPE,"_Sim",repmod,".RData"))
    df=NULL
    df <- read.fst(paste0(dir_results,"results_FishStrat_",iEXPE,"/Metapop_Sc",iEXPE,"_Sim",repmod,".fst"))
    load(paste0(dir_results,"results_FishStrat_",iEXPE,"/Fishing_rates_Sc",iEXPE,"_Sim",repmod,".RData"))
    
    for (pop in 1:npop){ # loop over popualtions
      
      demo <- NULL
      demo <- subset(df, Pop==pop)
      #demo <- results[[pop]]$pop
      #demo <- demo[demo$year>nInit,]
      #demo$year	<- demo$year - nInit
      nyears <- max(demo$year,na.rm=TRUE)-1
      #redds <- NULL
      #redds <- results[[pop]]$redds
      
      #fishing_tmp<-NULL
      #fishing_tmp <- results[[pop]]$fishing_rate
      #fishing_rates <- rbind(fishing_rates, fishing_tmp)
      
      #-------------------------------#
      #  Composition des populations  #
      #-------------------------------#  

      
      
      for (i in 1:nyears){
        
        ## PARR 0+
        tmp<-subset(demo,Parr==1 & date==273 & AgeRiver<1 & year==i)
        #Nparr0[i,pop] <- nrow(tmp)
        LFparr0[i,pop] <- mean(tmp$Lf,na.rm=TRUE)
        
        ## SMOLT 0+
        tmp<-subset(demo,Smolt==1 & date==90 & AgeRiver==1 & year==i)
        #Nsmolt0[i,pop] <- nrow(tmp)
        LFsmolt0[i,pop] <- mean(tmp$Lf,na.rm=TRUE)
        
        ## ANADROMOUS
        tmp<-subset(demo,Returns>0 & Atsea==0 & Mature==1 & date==273 & year==i)
        #NRet[i,pop] <- nrow(tmp)
        LFreturns[i,pop] <- mean(tmp$Lf,na.rm=TRUE)
        
        ## MSW
        tmp<-subset(demo,Returns>0 & Atsea==0 & Mature==1 & date==273 & AgeSea >= 2 & year==i)
        NMSW[i,pop] <- nrow(tmp)
        LFMSW[i,pop] <- mean(tmp$Lf,na.rm=TRUE)
        
        ## 1SW
        tmp<-subset(demo,Returns>0 & Atsea==0 & Mature==1 & date==273 & AgeSea < 2 & year==i)
        N1SW[i,pop] <- nrow(tmp)
        LF1SW[i,pop] <- mean(tmp$Lf,na.rm=TRUE)
        
        ## PropMSW
        RatioMSW[i,pop] <- NMSW[i,pop] / N1SW[i,pop]
        
      } # end loop years
  
      
    } # end loop population
    
    
    #tmp1[[repmod]] <- Nparr0
    #tmp2[[repmod]] <- Nsmolt0
    #tmp3[[repmod]] <- list(Nhom=Nhomers, NIm=NIm)
    #tmp4[[repmod]] <- NRet
    tmp5[[repmod]] <- NMSW
    tmp6[[repmod]] <- N1SW
    
    tmp7[[repmod]] <- RatioMSW
    
    tmp8[[repmod]] <- LFparr0
    tmp9[[repmod]] <- LFsmolt0
    tmp10[[repmod]] <- LFreturns
    tmp11[[repmod]] <- LFMSW
    tmp12[[repmod]] <- LF1SW
    
    
    
  } # end loop simul
  
  
  rMSW[[paste0(iEXPE)]] <- tmp7
  
  LFParr[[paste0(iEXPE)]] <- tmp8
  LFSmolt[[paste0(iEXPE)]] <- tmp9
  LFReturns[[paste0(iEXPE)]] <- tmp10
  LFmsw[[paste0(iEXPE)]] <- tmp11
  LF1sw[[paste0(iEXPE)]] <- tmp12
  
  #nExploitation[[paste0(iEXPE)]] <- tmp13
  
  #Fishing_rates[[paste0(iEXPE)]] <- fishing_rates


  #} # end scenario

### Save results

save(rMSW,LFParr,LFSmolt,LFReturns,LFmsw,LF1sw,file=paste0("results/PHENOTYPE",iEXPE,".RData"))

q('no')
