#!/bin/bash
#SBATCH -J flye
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
flye=/well/bag/fnd111/miniconda3/envs/flye/bin/flye
medaka=/well/bag/fnd111/miniconda3/envs/flye/bin/medaka_consensus
bwa=/well/bag/fnd111/miniconda3/envs/flye/bin/bwa
filter=/well/bag/fnd111/miniconda3/envs/flye/bin/polypolish_insert_filter.py
polypolish=/well/bag/fnd111/miniconda3/envs/flye/bin/polypolish

echo "Finished set up at: "`date`
echo "****************************************************"

echo "****************************************************"
echo "Running Flye assembly: "`date`

# $rasusa --input $lr --coverage 200 --genome-size 5mb --output ./long_qc_sub.fastq.gz
lr=./long_qc_sub.fastq.gz

$flye -o ./flye_output --plasmids --nano-raw $lr
$medaka -i $lr -d ./flye_output/assembly.fasta -o medaka_output -m r941_min_high_g360
$bwa index ./medaka_output/consensus.fasta
mkdir polypolish_output
$bwa mem -t 16 -a ./medaka_output/consensus.fasta $sr1 > ./polypolish_output/alignments_1.sam
$bwa mem -t 16 -a ./medaka_output/consensus.fasta $sr2 > ./polypolish_output/alignments_2.sam
$filter --in1 ./polypolish_output/alignments_1.sam --in2 ./polypolish_output/alignments_2.sam --out1 ./polypolish_output/filtered_1.sam --out2 ./polypolish_output/filtered_2.sam
$polypolish ./medaka_output/consensus.fasta ./polypolish_output/filtered_1.sam ./polypolish_output/filtered_2.sam > ./polypolish_output/polished.fasta


echo "Finished working on ${f}"

echo "Finished at: "`date`
echo "*****************************************************"

echo "*****************************************************"
echo "Finished!"
echo "*****************************************************"
