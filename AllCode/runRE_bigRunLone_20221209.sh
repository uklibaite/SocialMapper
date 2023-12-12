#!/bin/bash
#SBATCH -J runRE
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 64000
#SBATCH -t 0-12:00
#SBATCH --array=1-80

module load matlab/R2017a-fasrc02

export MATLABPATH=$DIRC'/code/'
cd $MATLABPATH

listOfFiles="/n/holylfs02/LABS/olveczky_lab/Ugne/bigRun/loneList_chd8_nrxn1.txt"
rowNumber=$SLURM_ARRAY_TASK_ID
fileName=$(sed "${rowNumber}q;d" $listOfFiles)

echo $MATLABPATH
echo $listOfFiles
echo $rowNumber
echo $fileName

matlab -nojvm -nodisplay -nosplash -nodesktop -r "cd('$MATLABPATH'); runReembedBigRun_HDZ('$fileName','/n/holylfs02/LABS/olveczky_lab/Ugne/bigRun/train1228.mat','/n/holylfs02/LABS/olveczky_lab/Ugne/bigRun/vecsVals_2022_12_22.mat','/n/holylfs02/LABS/olveczky_lab/Ugne/bigRun/RE_lone_chd8_nrxn1/'); exit;"
