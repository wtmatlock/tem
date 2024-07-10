#!/bin/bash
#SBATCH -J panaroo
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

panaroo=/well/bag/fnd111/miniconda3/envs/panaroo_env/bin/panaroo

echo "Finished set up at: "`date`
echo "****************************************************"
echo "****************************************************"
echo "Running Panaroo: "`date`

mkdir panaroo_all_chromos_98

$panaroo --clean-mode sensitive -i ./chromos_no_blaEC/*/prokka_output/*.gff -o ./panaroo_all_chromos_98 --aligner mafft -a core --core_threshold 0.98 -t 16
	
echo "Finished at: "`date`
echo "*****************************************************"

echo "*****************************************************"
echo "Finished!"
echo "*****************************************************"
