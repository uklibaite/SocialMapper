#!/bin/bash
#SBATCH -J runRE
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 64000
#SBATCH -t 0-12:00
#SBATCH --array=1-396

module load matlab/R2017a-fasrc02

export MATLABPATH=$DIRC'/code/'
cd $MATLABPATH

listOfFiles="/n/holylfs02/LABS/olveczky_lab/Ugne/bigRun/socList20230131.txt"
rowNumber=$SLURM_ARRAY_TASK_ID
fileName=$(sed "${rowNumber}q;d" $listOfFiles)

echo $MATLABPATH
echo $listOfFiles
echo $rowNumber
echo $fileName

matlab -nojvm -nodisplay -nosplash -nodesktop -r "cd('$MATLABPATH'); runReembedBigRunS_CCZ('$fileName','/n/holylfs02/LABS/olveczky_lab/Ugne/bigRun/trainingSocial_20230131.mat','/n/holylfs02/LABS/olveczky_lab/Ugne/bigRun/RE_S396/'); exit;"
