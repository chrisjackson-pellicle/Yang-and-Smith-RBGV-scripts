# Author: Alexander Schmidt-Lebuhn (modified by Chris Jackson Julu 2021)

from Bio import AlignIO
from Bio import Phylo
import Bio.Align
import sys
import glob
import os  # CJJ
import argparse  # CJJ
import logging  # CJJ
import socket  # CJJ
import fnmatch  # CJJ

# # first parameter is folder with subtrees; second is extension of subtree files for glob; third is alignment folder; all three are required
# treefilefolder = sys.argv[1]
# if (treefilefolder[len(treefilefolder)-1] != '/'):
# 	treefilefolder = treefilefolder + '/'
#
# treefileextension = sys.argv[2]
#
# alignmentfoldername = sys.argv[3]
# if (alignmentfoldername[len(alignmentfoldername)-1] != '/'):
# 	alignmentfoldername = alignmentfoldername + '/'
#
# # optional fourth parameter is starting characters that uniquely identify the outgroup sequence(s) that may have been cut off by MO and RT scripts
# # may want to enter a nonsense string to avoid grabbing the outgroup twice when having used MI, which often keeps it
# # works only with one outgroup taxon, unfortunately
# if len(sys.argv)>4:
# 	targetname = sys.argv[4]
# 	targetname = targetname.split(",")
# else:
# 	targetname = ["sunf"]			# if not provided, it is assumed that the target are the sunflower sequences from the Mandel et al 2014 Compositae bait set
#
# # optional fifth parameter is the folder where to output the subselected fasta files. Default is folder with alignment files
# if len(sys.argv)>5:
# 	outfoldername = sys.argv[5]
# 	if (outfoldername[len(outfoldername)-1] != '/'):
# 		outfoldername = outfoldername + '/'
# else:
# 	outfoldername = alignmentfoldername

########################################################################################################################
########################################################################################################################
# Get current working directory and host name

cwd = os.getcwd()
host = socket.gethostname()

# Configure logger:

# Create a custom logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Create handlers
c_handler = logging.StreamHandler(sys.stdout)
existing_log_file_numbers = [int(file.split('_')[-1]) for file in os.listdir('.') if fnmatch.fnmatch(file, '*.mylog*')]
if not existing_log_file_numbers:
    new_log_number = 1
else:
    new_log_number = sorted(existing_log_file_numbers)[-1] + 1
f_handler = logging.FileHandler(f'logging_file.mylog_{new_log_number}', mode='w')

# Create formatters and add it to handlers
c_format = logging.Formatter('%(message)s')
f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
c_handler.setFormatter(c_format)
f_handler.setFormatter(f_format)

# Add handlers to the logger
logger.addHandler(c_handler)
logger.addHandler(f_handler)

########################################################################################################################
########################################################################################################################


def createfolder(directory):
    """
    Attempts to create a directory named after the name provided, and provides an error message on failure
    """
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        logger.info(f'Error: Creating directory: {directory}')


def subsample_alignments(tree_folder, tree_suffix, alignment_folder, output_folder, from_cut_internal_branches=False):
	"""
	Takes a pruned/QC'd tree file, finds the original matching alignment, and sub-samples that alignment to recover
	only sequences corresponding to tree termini.
	"""
	createfolder(output_folder)

	for tree in glob.glob(f'{tree_folder}/*{tree_suffix}'):
		read_tree = Phylo.read(tree, "newick")
		tree_terminals = read_tree.get_terminals()
		tree_basename = os.path.basename(tree)

		# Derive the matching alignment file name depending on input tree file name:
		if from_cut_internal_branches:
			alignment_prefix = '_'.join(tree_basename.split('_')[0:-1])
			output_alignment_prefix = tree_basename.split('.')[0]
			print(f'alignment_prefix is: {alignment_prefix}')
			matching_alignment = f'{alignment_folder}/{alignment_prefix}.paralogs.aln.hmm.trimmed.fasta'
			print(matching_alignment)
		else:
			alignment_prefix = tree_basename.split('.')[0]
			output_alignment_prefix = alignment_prefix
			print(f'alignment_prefix is: {alignment_prefix}')
			matching_alignment = f'{alignment_folder}/{alignment_prefix}.outgroup_added.aln.trimmed.fasta'
			print(f'matching_alignment is: {matching_alignment}')

		# Read in original alignments and select seqs matching tree termini:
		alignment = AlignIO.read(matching_alignment, "fasta")
		subalignment = Bio.Align.MultipleSeqAlignment([])
		for k in range(0, len(alignment)):
			for j in range(0, len(tree_terminals)):
				if tree_terminals[j].name == alignment[k].id:
					subalignment.append(alignment[k])

		# Write an alignment of the sub-selected sequences:
		AlignIO.write(subalignment, f'{output_folder}/{output_alignment_prefix}.selected.fa', "fasta")


def parse_arguments():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('tree_file_folder', type=str, help='directory containing tree files')
	parser.add_argument('tree_file_extension', type=str, help='extension of tree files')
	parser.add_argument('alignment_file_folder', type=str, help='directory containing alignments to subsample from')
	parser.add_argument('output_alignment_folder', type=str, help='directory to contain output alignments')
	parser.add_argument('-from_cut_internal_branches', action='store_true', default=False,
						help='If set, trees are from the final QC stage during Yang and Smith pipeline i.e., the '
							 '<folder 06_cut_internal_branches>')
	results = parser.parse_args()
	return results


########################################################################################################################
########################################################################################################################
# Run script:

def main():
	results = parse_arguments()
	print(results)
	subsample_alignments(results.tree_file_folder, results.tree_file_extension, results.alignment_file_folder,
						 results.output_alignment_folder, from_cut_internal_branches=results.from_cut_internal_branches)


########################################################################################################################
########################################################################################################################

if __name__ == '__main__':
	if not len(sys.argv) >= 1:
		print(__doc__)
		sys.exit()
	sys.exit(main())

########################################################################################################################
########################################################################################################################