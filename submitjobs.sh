#!/bin/sh

# Run prefix
prefix=${PWD##*/}

# Use OpenMP
export OMP_NUM_THREADS=2

# Number of jobs, minimum 1
i=0;
end=1; 

# Generate initial conditions
echo "bsub -n 2 -W 0:30 -R 'rusage[mem=3000]' -J $prefix ./in2mart" >> submit

# Produce (chain-)jobs of desired length on 24h queue
until [ $i -ge $end ]
do
echo "bsub -n 2 -W 24:00 -R 'rusage[mem=3000]' -J $prefix -w \"ended($prefix)\" ./i2mart" >> submit
i=`expr $i + 1`
done
chmod 744 submit
./submit

# Delete submit file after jobs are submitted
rm -r submit
