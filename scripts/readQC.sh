#!/bin/bash
#SBATCH -J unicycler
#SBATCH -A bag.prj
#SBATCH -p long
#SBATCH --array 1-550:1

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
samples=/well/bag/fnd111/ecoli548/samples.txt
f=$(sed -n "$SLURM_ARRAY_TASK_ID"p $samples)
echo "Working on ${f}"
cd ./reads/$f

sr1=./short_1.fastq.gz
sr2=./short_2.fastq.gz
lr=./long.fastq.gz

fastp=/well/bag/fnd111/miniconda3/envs/readqc/bin/fastp
filtlong=/well/bag/fnd111/miniconda3/envs/readqc/bin/filtlong

echo "Finished set up at: "`date`
echo "****************************************************"

echo "****************************************************"
echo "Running Unicycler assembly: "`date`


$fastp --in1 $sr1 --in2 $sr2 --out1 short_1_qc.fastq.gz --out2 short_2_qc.fastq.gz --unpaired1 short_u.fastq.gz --unpaired2 short_u.fastq.gz
$filtlong --min_length 1000 --keep_percent 95 $lr | gzip > long_qc.fastq.gz

echo "Finished working on ${f}"

echo "Finished at: "`date`
echo "*****************************************************"

echo "*****************************************************"
echo "Finished!"
echo "*****************************************************"
