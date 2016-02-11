#!/bin/sh
i=0;
##############
# use OpenMP #
##############
export OMP_NUM_THREADS=2
################
#### CONFIG ####
################
end=20; # number of jobs
prefix=test
#################
#### /CONFIG ####
#################

echo "bsub -n 2 -W 24:00 -R 'rusage[mem=3000]' -J $prefix ./i2mart" >> submit
until [ $i -ge $end ]
do
echo "bsub -n 2 -W 24:00 -R 'rusage[mem=3000]' -J $prefix -w \"ended($prefix)\" ./i2mart" >> submit
i=`expr $i + 1`
done
chmod 744 submit
./submit
###############################################
# delete submit file after jobs are submitted #
###############################################
rm -r submit
