#!/bin/bash

#DIR="results/" # directory where to save results
#cd $DIR

#nSIMUL=10 # Nb simulations 100

# SIMULATION
echo "BEGINNING OF SIMULATIONS"
echo "PID du processus courant : $$"
STARTTIME=$(date +%s);

#scn=( 101 102 103 105 107 109 1011)

#scn=( 401 402 403 405 407 409 4011)
#scn=( 412 413 415 417 419 4111 4113)
#scn=( 422 423 425 427 429 4211 4213)

#scn=( 701 702 703 705 707 709 7011)
#scn=( 712 713 715 717 719 7111 7113)
#scn=( 722 723 725 727 729 7211 7213)

scn=(4113 4213 7113 7213)
for s in "${scn[@]}"
do

#STARTTIME=$(date +%s);
Rscript --vanilla DEMOGRAPHY.R $s &

echo "Extraction of scenario $s started!" 
sleep 2 # wait x seconds before starting the next extraction
done # end loop 


#ENDTIME=$(date +%s);
#MINUTES=$(( ($ENDTIME - $STARTTIME) / 60 ));
#echo "Simulations successful! Duration: $MINUTES minutes" 
#echo "END OF SIMULATIONS"
