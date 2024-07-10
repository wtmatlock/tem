#!/bin/bash
#SBATCH -J copynumber
#SBATCH -A bag.prj
#SBATCH -p long
#SBATCH --array 1-435:1

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
samples=/well/bag/fnd111/ecoli548/assemblies.txt
f=$(sed -n "$SLURM_ARRAY_TASK_ID"p $samples)
echo "Working on ${f}"

bwa=/well/bag/fnd111/miniconda3/envs/copynumber/bin/bwa
samtools=/well/bag/fnd111/miniconda3/envs/copynumber/bin/samtools
sr1=/well/bag/fnd111/ecoli548/reads/"$f"/short_1_qc.fastq.gz
sr2=/well/bag/fnd111/ecoli548/reads/"$f"/short_2_qc.fastq.gz

echo "Finished set up at: "`date`
echo "****************************************************"
echo "****************************************************"
echo "Running annotation tools: "`date`

cd /well/bag/fnd111/ecoli548/assemblies/"$f"
mkdir ./copynumber

for contig in ./contigs_assembly_"$f"/*.fasta; do
	basename=$(basename $contig .fasta)
	
	$bwa index $contig
	$bwa mem -t 4 $contig $sr1 $sr2 | $samtools view -b > ./copynumber/"$basename"_mapped_reads.bam
	$samtools sort ./copynumber/"$basename"_mapped_reads.bam -o ./copynumber/"$basename"_mapped_reads.bam.sorted
	$samtools index ./copynumber/"$basename"_mapped_reads.bam.sorted	
	$samtools depth -a ./copynumber/"$basename"_mapped_reads.bam.sorted > ./copynumber/"$basename"_mapped_reads.bam.sorted.depth
	
	awk '{sum += $3} END {print "'"$f"'", "'"$basename"'", sum / NR}' ./copynumber/"$basename"_mapped_reads.bam.sorted.depth > ./copynumber/"$basename"_mean_depth.txt
 
done

echo "Finished working on ${f}"

echo "Finished at: "`date`
echo "*****************************************************"

echo "*****************************************************"
echo "Finished!"
echo "*****************************************************"
