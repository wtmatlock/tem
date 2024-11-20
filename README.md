## *E. coli* phylogeny drives co-amoxiclav resistance through variable expression of *bla*<sub>TEM-1</sub>

### Scripts

Contained within the scripts directory is the code required to reproduce the analysis from the manuscript:

| Name       | Description | Notes |
|------------------|----------|----------|
| `annotate.sh`      | Runs a suite of genome annotation tools.         | Formatted to be run on a Slurm scheduler. |
| `buildResults.R`   | Reproduces all figure and sample statistics reported in the manuscript.         |          |
| `copyNumber.sh`    | Estimates the copy number of assembly contigs.        | Formatted to be run on a Slurm scheduler.         |
| `extractGene.py`   | Extracts a gene from a contig as a new FASTA file.         | Helper script.         | 
| `modelCombined.R`  | Runs the combined expression and MIC model.          |          |
| `modelExpression.R`| Runs the expression model.         |          |
| `modelMIC.R`       | Runs the MIC model.         |          |
| `readQC.sh`        | Quality control for the short- and long-reads.         |          |
| `removeGene.py`    | Removes a gene from a contig as a new FASTA file.          | Helper script.          |
| `runFlye.sh`       | Runs Flye to assemble a genome from short- and long-reads.         | Formatted to be run on a Slurm scheduler.         |
| `runIQTREE.sh`     | Runs IQ-TREE on a set of assemblies.          | Formatted to be run on a Slurm scheduler.         |
| `runPanaroo.sh`    | Runs Panaroo on a set of assemblies.         | Formatted to be run on a Slurm scheduler.         |
| `runUnicycler.sh`  | Runs Unicycler to assemble a genome from short- and long-reads.           |  Formatted to be run on a Slurm scheduler.        |
| `screenPlsdb.sh`   | Finds the PLSDB plasmid sequence which best contains a query contig based on *k*-mers.        |  Formatted to be run on a Slurm scheduler.        |

### Running the models

The files provided within the data directory are used by the modelling scripts (`modelExpression.R`, `modelMIC.R`, and `modelCombined.R`).
