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
import shutil


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


def batch_input_files(gene_fasta_directory, output_directory, batch_size=20):
    """
    Takes a folder of fasta files, and splits them in to batch folders according to the number provided by
    parameter batch_size.

    Parameters
    ----------
    gene_fasta_directory : path to the folder containing gene/paralog fasta files
    batch_size : number of gene/paralog fasta files to output in each batch folder

    Returns
    -------

    """
    createfolder(output_directory)

    fasta_file_list = glob.glob(f'{gene_fasta_directory}/*.fa')

    def chunks(lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    batches = list(chunks(fasta_file_list, batch_size))
    batch_num = 1
    for batch in batches:
        createfolder(f'{output_directory}/batch_{batch_num}')
        for fasta_file in batch:
            shutil.copy(fasta_file, f'{output_directory}/batch_{batch_num}')
        batch_num += 1


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('tree_file_folder', type=str, help='directory containing tree files')
    parser.add_argument('tree_file_extension', type=str, help='extension of tree files')
    parser.add_argument('alignment_file_folder', type=str, help='directory containing alignments to subsample from')
    parser.add_argument('output_alignment_folder', type=str, help='directory to contain output alignments')
    parser.add_argument('-from_cut_internal_branches', action='store_true', default=False,
                        help='If set, trees are from the final QC stage during Yang and Smith pipeline i.e., the '
                             '<folder 06_cut_internal_branches>')
    parser.add_argument('-batch_size', type=int, default=20,
                        help='Number of fasta files in each batch, from input paralog fasta files')
    results = parser.parse_args()
    return results


########################################################################################################################
########################################################################################################################
# Run script:

def main():
    folder_01 = f'{cwd}/08_selected_alignments_batch_folders'

    args = parse_arguments()
    print(args)
    subsample_alignments(args.tree_file_folder, args.tree_file_extension, args.alignment_file_folder,
                         args.output_alignment_folder, from_cut_internal_branches=args.from_cut_internal_branches)

    # Batch fasta files for alignment and tree-building steps;
    batch_input_files(args.output_alignment_folder, folder_01, batch_size=args.batch_size)


########################################################################################################################
########################################################################################################################

if __name__ == '__main__':
    if not len(sys.argv) >= 1:
        print(__doc__)
        sys.exit()
    sys.exit(main())

########################################################################################################################
########################################################################################################################