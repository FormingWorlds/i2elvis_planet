# Run prefix
prefix=${PWD##*/}

echo "#!/bin/sh" >> submit.sh

# SLURM commands (man sbatch / http://www.arc.ox.ac.uk/content/slurm-job-scheduler)
echo "#SBATCH --job-name=$prefix" >> submit.sh
echo "#SBATCH --output=$prefix.out" >> submit.sh
echo "#SBATCH --error=$prefix.err" >> submit.sh
echo "#SBATCH --time=336:00:00" >> submit.sh
echo "#SBATCH --ntasks=1" >> submit.sh
echo "#SBATCH --mem=3000" >> submit.sh
echo "#SBATCH -p priority-rp" >> submit.sh
#echo "#SBATCH -p shared" >> submit.sh
echo "#SBATCH --cpus-per-task=1" >> submit.sh

# Print date and time before job runs
echo "date" >> submit.sh

# Load environment
echo "module load netcdf/4.4.1__intel2015" >> submit.sh
echo "module load intel-compilers/2015" >> submit.sh
echo "module load intel-mpi/2015" >> submit.sh
echo "module load intel-mkl/2015" >> submit.sh

# Generate initial conditions and run the code
#echo "./in2mart" >> submit.sh
echo "./i2mart" >> submit.sh

#Print date and time after job finished
echo "date" >> submit.sh

# Submit created script
sbatch ./submit.sh

# Delete submit file after jobs are submitted
rm -r submit.sh