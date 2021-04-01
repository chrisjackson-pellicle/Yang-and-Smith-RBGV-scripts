# Author: Alexander Schmidt-Lebuhn


from Bio import AlignIO
from Bio import Phylo
import Bio.Align
import sys
import glob

# first parameter is folder with subtrees; second is extension of subtree files for glob; third is alignment folder; all three are required
treefilefolder = sys.argv[1]
if (treefilefolder[len(treefilefolder)-1] != '/'):
	treefilefolder = treefilefolder + '/'

treefileextension = sys.argv[2]

alignmentfoldername = sys.argv[3]
if (alignmentfoldername[len(alignmentfoldername)-1] != '/'):
	alignmentfoldername = alignmentfoldername + '/'

# optional fourth parameter is starting characters that uniquely identify the outgroup sequence(s) that may have been cut off by MO and RT scripts
# may want to enter a nonsense string to avoid grabbing the outgroup twice when having used MI, which often keeps it
# works only with one outgroup taxon, unfortunately
if len(sys.argv)>4:
	targetname = sys.argv[4]
	targetname = targetname.split(",")
else:
	targetname = ["sunf"]			# if not provided, it is assumed that the target are the sunflower sequences from the Mandel et al 2014 Compositae bait set

# optional fifth parameter is the folder where to output the subselected fasta files. Default is folder with alignment files
if len(sys.argv)>5:
	outfoldername = sys.argv[5]
	if (outfoldername[len(outfoldername)-1] != '/'):
		outfoldername = outfoldername + '/'
else:
	outfoldername = alignmentfoldername

treefilename = glob.glob(treefilefolder + '*' + treefileextension)

for i in range(0,len(treefilename)):
	thistree = Phylo.read(treefilename[i], "newick")
	terminals = thistree.get_terminals()
	# get last part of treefile name, i.e. the actual file name without path
	nameholder = treefilename[i].split("/")
	outfilename = nameholder[len(nameholder)-1]
	# split first on . and take first item, then on underscores and take first item
	# THIS MEANS YOU CANNOT HAVE UNDERSCORES IN GENE NAMES
	# along the way we store the
	nameholder = outfilename.split(".")[0]
	nameholder = nameholder.split("_")[0]
	# now glob for the alignment file name
	# THIS REQUIRES NON-NESTED GENE NAMES! Avoid e.g. LEAFY and LEAFYLIKE
	# You also MUST NOT have multiple fasta files for the same gene in the alignment file folder
	# all of this may be annoying, but it makes the script maximally flexible in dealing with different
	# kinds of tree file names and alignment file names
	alignmentfilename = glob.glob(alignmentfoldername+nameholder+"*.fasta")
	alignmentfilename = alignmentfilename[0]
	outfilename = outfilename[:(len(outfilename)-len(treefileextension))] #slice off extension from output file name
	print(alignmentfilename)
	alignment = AlignIO.read(alignmentfilename, "fasta")
	subalignment = Bio.Align.MultipleSeqAlignment([])
	for k in range(0,len(alignment)):
		for j in range(0,len(terminals)):
			if (terminals[j].name == alignment[k].id):
				subalignment.append(alignment[k])
		for l in range(0,len(targetname)):
			if alignment[k].id[:len(targetname[l])] == targetname[l]:
				subalignment.append(alignment[k])
	AlignIO.write(subalignment, outfoldername+outfilename+".selected.fa", "fasta")
	#AlignIO.write(subalignment, alignmentfoldername+nameholder+".selected.fa", "fasta")

# how to derive alignmentfilename: At1g01470_1.subtree; I assume it is _2 for the second iteration, etc.
# so just slice off "_1.subtree" to derive alignment file name (At1g01470.fa)
# outputfilename: remove last three letters (".fa"), attach "_rr.fa"
