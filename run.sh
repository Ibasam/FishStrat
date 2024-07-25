#!/bin/bash

#DIR="results/" # directory where to save results
#cd $DIR



# SIMULATION
echo "BEGINNING OF SIMULATIONS"
echo "PID du processus courant : $$"
STARTTIME=$(date +%s);

#for s in "${scn[@]}"
#do

Rscript --vanilla metaIbasam.R 7 3 1 10 1 50 &
#sleep 5;
#Rscript --vanilla metaIbasam.R 7 3 1 11 1 50 &
#sleep 5;
#Rscript --vanilla metaIbasam.R 4 3 1 11 1 50 &

#done # end loop 


ENDTIME=$(date +%s);
MINUTES=$(( ($ENDTIME - $STARTTIME) / 60 ));
echo "Simulations successful! Duration: $MINUTES minutes" 
echo "END OF SIMULATIONS"
