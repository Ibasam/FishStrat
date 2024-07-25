#!/bin/bash

#DIR="results/" # directory where to save results
#cd $DIR


# EXTRACTION
echo "BEGINNING OF EXTRACTION"
#echo "PID du processus courant : $$"
#STARTTIME=$(date +%s);

#Example: 15% dispersal, Exploitation all populations, Selective on MSW
scn=(4011 4013 4015 4017 4019 40111 40113 40115 40117 40119 40121)


for s in "${scn[@]}"
do
Rscript --vanilla code/METAPOP.R $s &

done # end loop

# sleep 60m
# 
# scn=(4011 4012 4013 4014 4015 4016 4017 4018 4019 40110 40111) #exploit all
# 
# for s in "${scn[@]}"
# do
# Rscript --vanilla code/METAPOP.R $s &
# done # end loop

echo "END OF EXTRACTION"



