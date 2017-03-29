import sys, argparse

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
sts_file, co_prob, co_perc, outname, concatonly = parseArgs()
inputfile=open(sts_file, "rU")
list1 =inputfile.readlines()
outputfile=open(outname+".prob+"+str(co_prob)+".perc"+str(co_perc)+".txt", "w")


def prob_check(list2check): #Removes all lines below the co_prob or that have 'no rank', returns list with	all fines above or equal to the co_prob
	outlist=[]
	#checklist=[]
	for line in list2check:
		if not (float(line[5]) < co_prob):
			if not (line[3] == 'no rank'):
				#checklist.append(line[5])
				outlist.append(line)
	#print checklist
	return outlist

def match_check(list2check, rank2check): #Function for checking specific level of taxonomy for matching, returns if matches above the co_perc.  and if True, what the match is
	ranklist = []
	markerlist=[]
	#adds all the lines that match the rank to ranklist
	for line in list2check:
		if line[3] == rank2check:
			ranklist.append(line)
			markerlist.append(line[-1])
	rankcount=[]
	#checks all possible combos for matching ranks with that line
	for line in ranklist:
		count=(len(line[1].split("."))/2)
		#count=0
		for line2 in ranklist:
			if line != line2:
				if line[4] == line2[4]:
					#count+=1
					count=count+(len(line2[1].split("."))/2)
		rankcount.append(count)
	num_matches= max(rankcount)
	total=len(rankcount)
	bm_index=rankcount.index(max(rankcount))
	bm_name =ranklist[bm_index][4]
	index=0
	mmlist=[]
	for num in rankcount:
		if num != max(rankcount):
			mmlist.append(ranklist[index][1])
		index+=1
	#print num_matches, total,  (float(num_matches)/float(total)), co_perc
	if total <= 1:
		return False, mmlist
	if (float(num_matches)/float(total)) < co_perc:
		return False, mmlist
	else:
		outputfile.write(rank2check+": "+bm_name+"\n")
		print(bm_name)
		return True, mmlist

def trueMatch(matchTF): #terminates program if matchTF is false
	if matchTF == False:
		outputfile.close()
		sys.exit()
		
def mremover(list2check, mlist): #removes markers ruled out from previous levels
	outputlist=[]
	for line in list2check:
			if (line[1] not in mlist) and (line[1] not in outputlist):
				if line not in outputlist:
					outputlist.append(line)
	return outputlist

def checkrank(list2check, rank): #checks that all the makers have that rank, returns list of markers to remove
	missmarklist=[]
	rankmatchlist=[]
	for line in list2check:
		if (line[3] == rank):
			rankmatchlist.append(line[1])
	for line in list2check:
		if (line[1] not in rankmatchlist) and (line[1] not in missmarklist):
			missmarklist.append(line[1])
	return missmarklist

def bacterialmonly(list2check):
	outputlist=[]
	for line in list2check:
		if line[-1].split("\n")[0] == "concat" or line[-1].split("\n")[0] == "16s_reps_bac":
			outputlist.append(line)
	return outputlist

#Sort lines of inputfile into lists to be used by prob_check
list2=[]
for line in list1:
	list2.append(line.split("\t"))
inputfile.close()
#print(list2[0])
list2.pop(0)


#run prob_check on the list
prob_list=prob_check(list2)
#print prob_list
#run match-check on each level
taxon_ranks=['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
for rank in taxon_ranks:
	if len(prob_list) != 0:
		if concatonly == True:
			prob_list=bacterialmonly(prob_list)
		norankmlist=checkrank(prob_list, rank)
		prob_list=mremover(prob_list,norankmlist)
		matchTF, mismlist =match_check(prob_list, rank)
	if len(mismlist) != 0:
		prob_list=mremover(prob_list, mismlist)
	trueMatch(matchTF)
