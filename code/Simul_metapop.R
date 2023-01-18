#_________________ SCENARIOS _________________#

# ENVIRONMENTAL CONDITIONS
#CHOOSE ENV SCENARIO
Obs=FALSE #here you can choose to use observed environmental data (flow in mm/day and air temperature) or to simulate them
Diversity_env=FALSE #FALSE
scenarioEnvi=1 #scenario 1 pour absence de CC, scenario 2 pour CC
# if want to use the same for several scenarios
# 1. generate them from the generate_env_conditions.R code
# 2. if already generated, load them from the data folder: env_new_false.Rdata & env_new_true.Rdata
#load("data/environmental_conditions_Diversity_FALSE_30simul.RData")
#load("data/environmental_conditions_Diversity_TRUE_100simul.RData")
# if want to generate new environmental conditions for each scenario, use the river_climate_multi_Tair script in the simulation loop below


#h = c(1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7) # Philopatry (homing) rates #edit al - 22/03/21 5% and 15%
scenarioConnect = 4 #as.numeric(args[1]) #scenario 1 for h=1.00, scenario 2 for h=0.95, scenario 3 pour h=0.80

# 0: control /
# 1: no fishing on sink populations
# 2: no fishing on source populations
# 3: no fishing on sink populations & stronger fishing pressure on sources
# 4: no fishing on source populations & stronger fishing pressure on sinks
scenarioFishing = 1 #as.numeric(args[2])

#frates_vec =c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
scenarioFishingRate = 2#as.numeric(args[3])


#Sc_name <-""
#Sc_name <- readline(prompt="Enter name: ")
Sc_name <- paste0("FishStrat_",scenarioConnect,scenarioFishing,scenarioFishingRate)
if (Sc_name == "") {
  print('WARNING: no scenario name provided, temporary files and results will store into /tmp_ & /results_ folders')
  
}

# Create results folder
#ifelse(!dir.exists(paste0(dir_results,'results_',Sc_name)),dir.create(paste0(dir_results,'results_',Sc_name)), FALSE)

#_________________ DATA _________________#
#source("dataIbasam.R")
#_________________ PARAMETERS _________________#
source("parIbasam.R")
#source(paste0("parIbasam",args,".R"))
print(fishing.rates)


###############
nyear=30
npop=15

Alpha=rep(0.1,npop)
Rmax=rep(10,npop)

N<-phi<-A<-R <-ratio <-array(,dim=c(nyear,npop))
N[1,]<- Area*0.01
ratio[1,]<- 1

for (year in 2:nyear){
  for (pop in 1:npop){
    
    # Recruitment
    #phi[year,pop] <- Rmax[pop]/(N[year-1,pop]+Alpha[pop]*Rmax[pop]) # 	gSurv=Nmax_ / (Ntot + gBH_ * Nmax_);//relation of Beverton and Holt (1957) type for survival here.
    #phi[t] <- rmax/(s_theo[t]+rmax/alpha)
    #phi[year,pop] <- (Alpha[pop]*Rmax[pop])/(Alpha[pop]*N[year-1,pop]+Rmax[pop])
    phi[year,pop] <-1
    R[year,pop] <- N[year-1,pop] * phi[year,pop]
    
    # survival
    A[year,pop] <- R[year,pop] * 1
  }
  
    # dispersal
    tmp <- A[year,]
    disp <- tmp * pstray
    
    for (pop in 1:npop){
      
    ratio[year,pop]<- sum(disp[pop,])/sum(disp[,pop]) # Im / emi
      
    N[year,pop]<- sum(disp[,pop])#N[year-1,pop]* + img[year,pop] - emg[year,pop]
  }
}

plot(NULL,xlim=c(0,nyear),ylim=c(min(N),max(N)))
for (pop in 1:npop){
  points(1:nyear,N[,pop],col=pop,type='b')
}

plot(NULL,xlim=c(0,nyear),ylim=c(0.6,max(ratio)))
abline(h=1,lty=2)
for (pop in 1:npop){
  points(1:nyear,ratio[,pop],col=pop,type='b')
}
