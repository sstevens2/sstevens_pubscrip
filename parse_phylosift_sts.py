import sys, argparse, pandas as pd

__author__ = "Sarah Stevens"
__email__ = "sstevens2@wisc.edu"

# Inputs
def parseArgs():
	parser = argparse.ArgumentParser(description='parse_phylosift_sts.py: takes a \
		sequence_taxa_summary.txt file from the program phyosift (run on only one \
		bin/genome) and uses a cutoff probability and percentage to determine the \
		taxonomy of that bin/genome')
	parser.add_argument('--stsfile','-sts' , action="store", dest='sts_in', type=str, \
		required=True, metavar='sequence_taxa_summary.txt', \
		help='This is the sequence_taxa_summary.txt file output by phylosift.')
	parser.add_argument('--cutoff_prob','-co_prob', action="store", dest='co_prob', \
		type=float, required=True, metavar='.90', help='The cutoffprob. removes any marker \
		genes with probabilies below that value, for each level of taxonomy. \
		Must be decimal form.')
	parser.add_argument('--cutoff_perc','-co_perc', action="store", dest='co_perc', \
		type=float, required=True, metavar='.70', help='The cutoffperc. only allows for \
		classification if the marker genes match at that level of taxonomy above this value.\
		Must be decimal form.')
	parser.add_argument('--out_name','-out', action="store", dest='outname', \
		type=str, required=True, metavar='sts_mygenome.txt', help='Allows you to specify the prefix \
		for the output filename')
	parser.add_argument('--bactArchMarkersOnly', '-bamo', action='store_true', \
		dest='concatonly', default=False, \
		help='Must use "True" or "False" exactly for using only the Bacterial/Archaeal Marker genes')
	args=parser.parse_args()
	return args.sts_in, args.co_prob, args.co_perc, args.outname, args.concatonly



#Bring in the inputs and set up output
sts_filename, co_prob, co_perc, outname, concatonly = parseArgs()
sts_file=pd.read_table(sts_filename,header=0)
# inputfile=open(sts_file, "rU")
# list1 =inputfile.readlines()
outputfile=open(outname+".prob+"+str(co_prob)+".perc"+str(co_perc)+".txt", "w")


# Remove all 'no rank' and hits less than co_prob
sts_file=sts_file[sts_file['Taxon_Rank']!='no rank']
sts_file=sts_file[sts_file['Cumulative_Probability_Mass']>=co_prob]
# If concatonly is true, remove all nonconcat
if concatonly:
    sts_file=sts_file[sts_file['Markers_Hit']=='concat']
taxon_ranks=['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
# For each taxon rank 
	for rank in taxon_ranks:
	#temp subset to only that level
	# group by Taxon names and count
	# get the biggest taxon name
		# if largest # hit / total > co_perc:
			# Write name to variable/file
		# else
			# close file
			# end program early
	# remove any marker genes that didn't match the biggest taxon name found in this taxon rank
