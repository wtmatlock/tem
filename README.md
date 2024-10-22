### *E. coli* phylogeny drives co-amoxiclav resistance through variable expression of *bla*<sub>TEM-1</sub>

Contained within the scripts directory is the code required to reproduce the analysis from the manuscript:

| Name       | Description | Notes |
|------------------|----------|----------|
| annotate.sh      | Runs a suite of genome annotation tools.         | Formatted to be run on a Slurm scheduler. |
| buildResults.R   | Reproduces all figure and sample statistics areported in the manuscript.         |          |
| copyNumber.sh    | Estimates the copy number of assembly contigs.        | Formatted to be run on a Slurm scheduler.         |
| extractGene.py   | Extracts a gene from a contig as a new FASTA file.         | Helper script.         | 
| modelCombined.R  |          |          |
| modelExpression.R|          |          |
| modelMIC.R       |          |          |
| readQC.sh        |          |          |
| removeGene.py    |          |          |
| runFlye.sh       |          | Formatted to be run on a Slurm scheduler.         |
| runIQTREE.sh     |          | Formatted to be run on a Slurm scheduler.         |
| runPanaroo.sh    |          | Formatted to be run on a Slurm scheduler.         |
| runUnicycler.sh  |          |  Formatted to be run on a Slurm scheduler.        |
| screenPlsdb.sh   |          |  Formatted to be run on a Slurm scheduler.        |
