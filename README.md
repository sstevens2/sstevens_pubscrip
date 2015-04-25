# sstevens_pubscrip

###This repo contains scripts written by Sarah Stevens.


#### parse_phylosift_sts.py
	This script is meant to be used on the sequence_taxa_summary file
	output from the program phylosift run on a single genome or genome
	bin.  It goes though each linnean classification level and removes
	any marker gene hits below the set cutoffprob. value.  It then looks
	for a classification matching above the cutoffperc. value specificed.
#### filterbylen3.py
	This script takes a fasta file and filters out reads/contigs which are
	below your minlen cutoff.
