#!/bin/bash
#SBATCH -J mapPLSDB
#SBATCH -A bag.prj
#SBATCH -p short
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
samples=/well/bag/fnd111/ecoli548/assemblies.txt
f=$(sed -n "$SLURM_ARRAY_TASK_ID"p $samples)
echo "Working on ${f}"

plsdb=/well/bag/fnd111/ecoli548/plsdb/contigs_plsdb
mash=/well/bag/fnd111/miniconda3/envs/mash_env/bin/mash

echo "Finished set up at: "`date`
echo "****************************************************"
echo "****************************************************"
echo "Running Mash screen: "`date`

cd /well/bag/fnd111/ecoli548/assemblies/"$f"
rm -r mash_output
mkdir mash_output

for contig_path in ./contigs_assembly_"$f"/*.fasta; do
	$mash sketch $contig_path
	contig=$(basename "$contig_path")
	find "$plsdb" -name '*.fasta' -print0 | while IFS= read -r -d '' plasmid; do 
		mash screen "$contig_path".msh "$plasmid" | awk -v plasmid_name="$plasmid" '{print $0"\t"plasmid_name}' >> ./mash_output/"$contig".screen.tab
	done
done

echo "Finished working on ${f}"

echo "Finished at: "`date`
echo "*****************************************************"

echo "*****************************************************"
echo "Finished!"
echo "*****************************************************"
