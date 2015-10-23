#!/usr/local/bin/python3

import sys, os
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord

"""maskrRNAgbk.py: masking 16S regions in genbank files"""

__author__ = "Sarah Stevens"
__email__ = "sstevens2@wisc.edu"

def usage():
    print("Usage: maskrRNAgbk.py genbankfile.gbk")
    print("Finds each annotated rRNA sequence and replaces all nucleotides with N")
    print("Leaves annotation there so you know where it used to be")
    print("Makes fastafile for output")

if len(sys.argv) != 2:
	usage()
	sys.exit(2)

input=open(sys.argv[1], 'rU')
# Open the genbank file as SeqRecord,
records = list(SeqIO.parse(input, "genbank"))
for contig in records:
    # For every annotation
    for record in contig.features:
        # If the feature is rRNA
        if record.type == 'rRNA':
            #print(record.location) #testing
            #print('before:',contig.seq) #testing
            before_len=len(contig.seq) # for assert below
            # replace seqeunce with a new one with N's where the rRNA feature is
            contig.seq=Seq.Seq(str(contig.seq[0:record.location.start.position]).rstrip()+('N'*len(contig.seq[record.location.start.position:record.location.end.position])).rstrip()+str(contig.seq[record.location.end.position:]).rstrip())
            #print('after:',contig.seq) #testing
            assert before_len==len(contig.seq)  # Sanity check that the sequence is still the same size

# Write output to file
output_fna=open(os.path.splitext(sys.argv[1])[0]+'_norrna.fna','w')
SeqIO.write(records,output_fna,'fasta')
output_fna.close()
# Output genbank file too, currently disabled since the IMG locus tags are too long for proper genbank files
"""
output_gbk=open(os.path.splitext(sys.argv[1])[0]+'_norrna.gbk','w')
SeqIO.write(contig,output_gbk,'genbank')
output_gbk.close()
"""
