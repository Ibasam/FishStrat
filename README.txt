This R project gathers all the data and R code needed to reproduce simulations and results analysis that are presented in the PhD chapter untitled "Management strategies of exploited Atlantic salmon metapopulation".

It is the third application study using MetaIBASAM, a demo-genetic agent-based model to simulate spatially structured salmon populations.

MetaIbasam is an extension of the existing IBASAM model (https://github.com/Ibasam/IBASAM/wiki) by incorporating a dispersal process to describe Atlantic salmon metapopulation and its eco-evolutionary dynamics. MetaIBASAM allows an investigation of the consequences of dispersal on local populations and network dynamics at the demographic, phenotypic, and genotypic levels. More generally, it allows to explore eco-evolutionary dynamics by taking into account complex interactions between ecological and evolutionary processes (plasticity, genetic adaptation and dispersal), feedbacks (e.g. genetic <-> demography) and trade-offs (e.g. growth vs survival). By doing so, one can investigate responses to changing environments and alternative management strategies of exploitation.

In this study, we ran different scenarios of dispersal rates and spatialiazed management strategies of exploitation that consider the metapopulation structure of populations, and looked at the evolutionary and demographic consequences of the strategies.

----

Contacts: 
amaia.lamarins@inrae.fr 
mathieu.buoro@inrae.fr

----

Below is described the workflow:


0. Install metaIbasam package v.0.0.6 through the .tar.gz file, define the directory in a terminal and create a "results" folder

> cd /folder/FishStrat
> mkdir results

1. Run the model for a defined scenario and number of simulations - launch simulations in a terminal:

> Rscript --vanilla metaIbasam.R 1 0 3 1 50 &

Arguments:
#1: scenarioConnect (e.g. 1 for 0 dispersal)
#2: scenarioFishing (e.g. 0 for fishing all populations)
#3: scenarioFishingrates (e.g. 3 for fishing at 7% 1SW and 15% MSW)
#4: first simulation (e.g. 1)
#5: last simulation (e.g. 50)

Each scenario (combination of arguments 1-8) has its own temporary folder and results folder. Several simulations can run in parallel, depending on the computer memory.



2. To launch simulations at a specific hour (using the linux package "at"):
# at: could be today, tomorrow,... friday,...
> at 04:00am tomorrow
> pkill -9 -u mbuoro R
> Rscript --vanilla metaIbasam.R 1 0 3 1 50 &
> Rscript --vanilla metaIbasam.R 2 0 3 1 50 &
> Rscript --vanilla metaIbasam.R 3 0 3 1 50 &

then crtl+D



3. Extract the results (e.g. DEMOGRAPHY.R, PHENOGENOTYPE.R, FITNESS.R) on scenarios defined beforehand in the file Rextract.sh :

> cd code/
> nohup bash Rextract.sh &


The Rextract.sh file launches the files "DEMOGRAPHY.R", "PHENOGENOTYPE.R", and "FITNESS.R" on each scenario defined in Rextract.sh and number of simulations defined in the files.


4. Plot the figures found in the paper.

Run the R code in the file "FIGURES.R" in the "code" folder.


5. If there is a need to move heavy files between folders.
 mv -v /folder/results/results_Dispersal_10311111/results_Dispersal_10311111/* /folder/results/results_Dispersal_10311111
 
 mv /folder/results/ /another_folder/results/
 mv /folder/results/PHENO* /folder/results/PHENOTYPE/
 
or between computers: source (147.100.121.93) -> destination (147.100.121.92)
 scp -r /folder/results/* mbuoro@147.100.121.92:/folder/results/
 
 
scn=105
mv /folder/results_FishStrat_$scn/* /folder/results/results_FishStrat_$scn/
rm -R /folder/results_FishStrat_$scn

 
 