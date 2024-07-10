import os
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein

# open the TSV file and read in the data
with open('/well/bag/fnd111/ecoli548/ec_positions.tsv', 'r') as f:
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
    with open('/well/bag/fnd111/ecoli548/assemblies/{}/contigs_assembly_{}/{}.fasta'.format(assembly, assembly, contig), 'r') as f:
        fasta_data = f.readlines()
    seq = ''.join(fasta_data[1:]).replace('\n', '')

    # extract the desired subsequence
    subseq = seq[start-1:end]

    # reverse complement the subsequence if the strand is -1
    if strand == '-':
        subseq = str(Seq(subseq).reverse_complement())

    # write the subsequence to a new FASTA file
    with open('/well/bag/fnd111/ecoli548/blaEC/{}-{}-blaEC-{}-{}.fasta'.format(assembly, contig, start, end), 'w') as f:
        f.write('>{}_{}_EC_{}_{}\n{}\n'.format(assembly, contig, start, end, subseq))

    # translate the subsequence to protein and write it to a new FASTA file
    prot_seq = Seq(subseq, generic_dna).translate()
    with open('/well/bag/fnd111/ecoli548/blaEC/{}-{}-blaEC-{}-{}.fasta.prot'.format(assembly, contig, start, end), 'w') as f:
        f.write('>{}_{}_EC_{}_{}\n{}\n'.format(assembly, contig, start, end, prot_seq))
