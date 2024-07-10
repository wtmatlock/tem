#!/bin/bash
#SBATCH -J unicycler
#SBATCH -A bag.prj
#SBATCH -p long
#SBATCH --array 1-4:1

# This script is run in the Conda environment
# ...details...

echo "****************************************************"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Slurm Task ID: $SLURM_ARRAY_TASK_ID"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "****************************************************"

echo "****************************************************"
echo "Setting up paths and directories at: "`date`

cd /well/bag/fnd111/ecoli548
samples=/well/bag/fnd111/ecoli548/redo.txt
f=$(sed -n "$SLURM_ARRAY_TASK_ID"p $samples)
echo "Working on ${f}"
cd ./reads/$f

sr1=./short_1_qc.fastq.gz
sr2=./short_2_qc.fastq.gz
lr=./long_qc.fastq.gz

rasusa=/well/bag/fnd111/miniconda3/envs/unicycler/bin/rasusa
unicycler=/well/bag/fnd111/miniconda3/envs/unicycler/bin/unicycler
spades=/well/bag/fnd111/miniconda3/envs/unicycler/bin/spades.py
racon=/well/bag/fnd111/miniconda3/envs/unicycler/bin/racon

echo "Finished set up at: "`date`
echo "****************************************************"

echo "****************************************************"
echo "Running Unicycler assembly: "`date`

$rasusa --input $lr --coverage 200 --genome-size 5mb --output ./long_qc_sub.fastq.gz
lr=./long_qc_sub.fastq.gz 

$unicycler -1 $sr1 -2 $sr2 -l $lr -o ./unicycler_output --spades_path $spades --racon_path $racon -t 16 --min_component_size 500 --min_dead_end_size 500 --verbosity 1
$unicycler -1 $sr1 -2 $sr2 -l $lr -o ./unicycler_output_conservative --spades_path $spades --racon_path $racon -t 16 --min_component_size 500 --min_dead_end_size 500 --verbosity 1 --mode conservative
$unicycler -1 $sr1 -2 $sr2 -l $lr -o ./unicycler_output_bold --spades_path $spades --racon_path $racon -t 16 --min_component_size 500 --min_dead_end_size 500 --verbosity 1 --mode bold


echo "Finished working on ${f}"

echo "Finished at: "`date`
echo "*****************************************************"

echo "*****************************************************"
echo "Finished!"
echo "*****************************************************"
