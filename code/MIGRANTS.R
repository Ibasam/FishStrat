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


nSIMUL <- 50
nYear <- 40
nInit <- 10 # Nb years at initialization
npop <- 15
#popo=15
nyears=nYear+nInit

# SCENARIOS
iEXPE <- args
#EXPE <- c(110, 210,220)#,120,123,124, 310, 320, 323, 324)
#EXPE <- c(120, 121, 122, 123, 124, 320, 321, 322, 323, 324)

#for (iEXPE in EXPE){ # Loop over scenario 

  ratioMigrants <- ratioMigrantsVsExpl <- list() # Nb returns
  tmp1 <- list()
  
  #_________________ DATA _________________#
  #fishing.rates = scn = NULL
  #scn <- as.numeric(strsplit(as.character(iEXPE), "")[[1]])
  #scenarioConnect = scn[1]
  #scenarioEnvi = scn[2]
  #scenarioFishing = scn[3]
  

  cat("Composition for ", iEXPE, "\n")
  
  #pe=petr=sync=CV_est=CV_obs=CV_esttr=CV_obstr=NULL
  ratioMig=ratioMigExpl=array(0,dim=c(50,nSIMUL,npop))
  
  for (repmod in 1:nSIMUL){ # loop over simulations
    #cat(indice,"/",repmod,"- ")
    cat(iEXPE,"-Simulation :",repmod," (")
    perc <- (repmod/nSIMUL)*100
    if(perc %in% (seq(0,90,10))) cat(perc,"%)","\n")
    if(perc == 100) cat(perc,"%)","\n")
    

    
    #load(paste0(dir_results,"results_FishStrat_",iEXPE,"/Metapop_Sc",iEXPE,"_Sim",repmod,".RData"))
    df=NULL
    df <- read.fst(paste0(dir_results,"results_FishStrat_",iEXPE,"/Metapop_Sc",iEXPE,"_Sim",repmod,".fst"))
    load(paste0(dir_results,"results_FishStrat_",iEXPE,"/Fishing_rates_Sc",iEXPE,"_Sim",repmod,".RData"))
    
    
    nMig=array(0,dim=c(npop,npop))
    Exploitation <- matrix(NA,nYear+nInit,npop)

      #-------------------------------#
      #  Composition des populations  #
      #-------------------------------#  
      
    returns <- subset(df,Returns>0 & Atsea==0 & Mature==1 & date==273)
      
      for (i in 1:nyears){
        for (pop in 1:npop){ # loop over popualtions
        tmp <-subset(returns,Pop==pop & year==i)
        migs <- table(tmp$CollecID)
        nMig[pop,as.numeric(names(migs))] <- migs   
        
        ## FISHING
        if (i<=nInit){
          Exploitation[i,pop] <- sum(migs)*fishing.rates[pop,1]
        }else{
          Exploitation[i,pop] <- sum(migs)*fishing.rates[pop,2]
        }
        
        } # end loop population
        
        for (pop in 1:npop){
        ratioMig[i,repmod,pop] <- sum(nMig[pop,],na.rm=TRUE)/sum(nMig[,pop],na.rm=TRUE) # Immigrants/Emigrants
        
        ratioMigExpl[i,repmod,pop] <- (sum(nMig[pop,],na.rm=TRUE) - nMig[pop,pop] +1) / (Exploitation[i,pop]+1) # Immigrants/Captured
        
        }
        
      } # end loop years
  } # end loop simul
  
  
  ratioMigrants[[paste0(iEXPE)]] <- apply(ratioMig,c(1,3),mean)
  ratioMigrantsVsExpl[[paste0(iEXPE)]] <- apply(ratioMigExpl,c(1,3),mean)
  


### Save results
save(ratioMigrants, ratioMigrantsVsExpl,file=paste0("results/MIGRANTS/MIGRANTS",iEXPE,".RData"))

q('no')
