import os
from Bio.Seq import Seq

# open the TSV file and read in the data
with open('/Users/willmatlock/Desktop/chromosome-analysis/ec_positions.tsv', 'r') as f:
    data = [line.strip().split('\t') for line in f.readlines()]

# loop over each row in the TSV file
for row in data:
    # extract the relevant data
    assembly = row[0]
    contig = row[1]
    start = int(row[2])
    end = int(row[3])
    strand = row[4]

    # open the corresponding FASTA file and read in the sequence
    fasta_path = f'/Users/willmatlock/Desktop/chromosome-analysis/contigs_assembly_{assembly}/{contig}.fasta'
	
    os.mkdir(f'/Users/willmatlock/Desktop/chromosome-analysis/phylos/chromos_no_blaEC/{assembly}')
    output_path = f'/Users/willmatlock/Desktop/chromosome-analysis/phylos/chromos_no_blaEC/{assembly}/{assembly}.fasta'

    with open(fasta_path, 'r') as f:
        fasta_data = f.readlines()
    seq = ''.join(fasta_data[1:]).replace('\n', '')

    print(f'Before: {len(seq)}')

    # remove the gene sequence
    seq_without_gene = seq[:start-1] + seq[end:]

    print(f'After: {len(seq_without_gene)}')
    print(f'Gene length: {len(seq) - len(seq_without_gene)}')

    # write the modified sequence to a new FASTA file
    with open(output_path, 'w') as f:
        f.write(f'>{assembly}\n{seq_without_gene}\n')
