#!/bin/bash
#
# The #SBATCH lines below are //not// commented out! These lines are read by the Slurm preprocessor. They need to start with a '#' character to work.
#
#SBATCH --job-name=tr_overlap_subset
#SBATCH --output=tr_overlap_subset.%j.%N.out
# Optional: Get an email notification on job completion:
#SBATCH --mail-type=end
#SBATCH --mail-user=delucmat@students.zhaw.ch
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --partition=single
#SBATCH --qos=single

module load slurm/17.11.8
module load gcc/7.3.0
module load bzip2/1.0.6
module load freetype/2.9.1 
module load curl/7.60.0
module load icu4c/60.1
module load libtiff/4.0.9
module load tcl/8.6.8
module load pcre/8.42
module load readline/7.0
module load glib/2.56.3-python-3.6.5-perl-5.26.2
module load cairo/1.16.0-a-python-3.6.5-perl-5.26.2
module load tk/8.6.8
module load jdk/11.0.2_9
module load ncurses/6.1
module load zlib/1.2.11
module load pango/1.41.0-a-python-3.6.5-perl-5.26.2
module load libx11/1.6.5
module load bison/3.0.5

# Not sure if the following is needed:
export R_HOME=$LSFM_CLUSTER_SCRATCH_USER_PATH/R-3.5.3

# Make R commands available in the path:
export PATH=$LSFM_CLUSTER_SCRATCH_USER_PATH/R-3.5.3/bin:$PATH

# Get the name of the compute node this job just got assigned:
NODE="$(hostlist -e $SLURM_NODELIST)"

# Copy all data from the local working dir to the local scratch dir on the node:
scp -rp * $NODE.c.hpc.zhaw.ch:/$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH

# Change current directory to the local scratch dir on the node:
cd $LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH

#srun Rscript ./compute_tr_all_overlap_DEBUG.R
#srun Rscript ./compute_tr_all_overlap.R
srun Rscript ./compute_tr_all_overlap_subset.R


