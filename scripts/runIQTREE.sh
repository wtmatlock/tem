#!/bin/bash
#SBATCH -J iqtree
#SBATCH -A bag.prj
#SBATCH -p long

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

cd /well/bag/fnd111/ecoli548/phylos

iqtree=/well/bag/fnd111/miniconda3/envs/iqtree_env/bin/iqtree

echo "Finished set up at: "`date`
echo "****************************************************"
echo "****************************************************"
echo "Running Panaroo: "`date`

mkdir iqtree_all_chromos_98
cd iqtree_all_chromos_98

aln=/well/bag/fnd111/ecoli548/phylos/panaroo_all_chromos_98/core_gene_alignment_filtered.aln

$iqtree -s $aln -m GTR+F+I+R4 -keep-ident -B 1000 -mem 10GB --prefix GTR_F_I_R4
	
echo "Finished at: "`date`
echo "*****************************************************"

echo "*****************************************************"
echo "Finished!"
echo "*****************************************************"
