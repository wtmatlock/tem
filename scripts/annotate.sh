#!/bin/bash
#SBATCH -J annotate
#SBATCH -A bag.prj
#SBATCH -p long
#SBATCH --array 1-440:1

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

bakta=/well/bag/fnd111/miniconda3/envs/annotate/bin/bakta
bakta_path=/well/bag/fnd111/ecoli548/bakta_path/db
amrfinder=/well/bag/fnd111/miniconda3/envs/annotate/bin/amrfinder
rgi=/well/bag/fnd111/miniconda3/envs/rgi_env/bin/rgi
rgi_db=/well/bag/fnd111/ecoli548/card.json
makeblastdb=/well/bag/fnd111/miniconda3/envs/annotate/bin/makeblastdb
tblastn=/well/bag/fnd111/miniconda3/envs/annotate/bin/tblastn
tem1=/well/bag/fnd111/ecoli548/WP_000027057.1.fasta
blastn=/well/bag/fnd111/miniconda3/envs/annotate/bin/blastn
promoters=/well/bag/fnd111/ecoli548/TEM1-promoters.fasta
abricate=/well/bag/fnd111/miniconda3/envs/annotate/bin/abricate
mobtyper=/well/bag/fnd111/miniconda3/envs/mobtyper_env/bin/mob_typer
mlst=/well/bag/fnd111/miniconda3/envs/annotate/bin/mlst
ezclermont=/well/bag/fnd111/miniconda3/envs/annotate/bin/ezclermont

echo "Finished set up at: "`date`
echo "****************************************************"
echo "****************************************************"
echo "Running annotation tools: "`date`

cd /well/bag/fnd111/ecoli548/assemblies/"$f"
#mkdir ./{bakta_output,promoter_output,tem_output,plasmidfinder_output,mlst_output,ezclermont_output,amrfinder_output}
#mkdir ./mobtyper_output

#$rgi load --card_json $rgi_db --local

for contig in ./contigs_assembly_"$f"/*.fasta; do
	basename=$(basename $contig .fasta)
	
	#$bakta --db $bakta_path --force --output bakta_output --prefix "$basename" --genus Escherichia --species coli --keep-contig-headers "$contig"
	
	#$amrfinder -n $contig -o ./amrfinder_output/"$basename".tsv --plus --organism Escherichia
	
	#$rgi main -i $contig -o ./rgi_output/"$basename" --local --data wgs

	#$makeblastdb -in $contig -dbtype nucl -out ./tem_output/"$basename"/"$basename"
	#$tblastn -query $tem1 -db ./tem_output/"$basename"/"$basename" -out ./tem_output/"$basename"/"$basename".txt -outfmt 6

	#$makeblastdb -in $contig -dbtype nucl -out ./promoter_output/"$basename"/"$basename"
	#$blastn -query $promoters -db ./promoter_output/"$basename"/"$basename" -out ./promoter_output/"$basename"/"$basename".txt -outfmt 6

	#$abricate $contig --db plasmidfinder > ./plasmidfinder_output/"$basename".tsv

	#$mobtyper --infile $contig --out_file ./mobtyper_output/"$basename".txt

	#$mlst $contig > ./mlst_output/"$basename".tsv

	#$ezclermont $contig > ./ezclermont_output/"$basename".tsv
done

#awk -v f="$f" 'BEGIN{OFS="\t"} {$NF=f; print}' ./rgi_output/*.txt > ./rgi_output/report.tsv
#awk -v f="$f" 'BEGIN{OFS="\t"} {$NF=f; print}' ./amrfinder_output/*.tsv > ./amrfinder_output/report.tsv
#awk -v f="$f" 'BEGIN{OFS="\t"} {$NF=f; print}' ./tem_output/*/*.txt > ./tem_output/report.tsv 
#awk -v f="$f" 'BEGIN{OFS="\t"} {$NF=f; print}' ./promoter_output/*/*.txt > ./promoter_output/report.tsv
#awk '(NR == 1) || (FNR > 1)' ./plasmidfinder_output/*.tsv > ./plasmidfinder_output/report.tsv
#cat ./mlst_output/*.tsv > ./mlst_output/report.tsv
#cat ./ezclermont_output/*.tsv > ./ezclermont_output/report.tsv.tmp
#awk -v f="$f" 'BEGIN{OFS="\t"} {$3=f; print}' ./ezclermont_output/report.tsv.tmp > ./ezclermont_output/report.tsv
#rm ./ezclermont_output/report.tsv.tmp
awk -v f="$f" 'BEGIN{OFS="\t"} {$NF=f; print}' ./mobtyper_output/*.txt > ./mobtyper_output/"$f"_report.tsv

echo "Finished working on ${f}"

echo "Finished at: "`date`
echo "*****************************************************"

echo "*****************************************************"
echo "Finished!"
echo "*****************************************************"
