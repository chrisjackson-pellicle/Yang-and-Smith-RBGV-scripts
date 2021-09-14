#!/usr/bin/env python

# Author: Chris Jackson

"""
- Checks gene names in paralog files and the external outgroup file (if provided) for dots, and converts them to
  underscores.
- Checks if there an outgroup (internal or external) for each gene paralog file.
- Takes a folder of fasta files, and splits them in to batch folders according to the number provided by parameter
  batch_size.
"""

import logging
import sys
import argparse
import os
import socket
import gzip
import re
import fnmatch
import glob
import shutil
from collections import defaultdict
from Bio import SeqIO


if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")


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
# Define functions:


def createfolder(directory):
    """
    Attempts to create a directory named after the name provided, and provides an error message on failure
    """
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        logger.info(f'Error: Creating directory: {directory}')


def file_exists_and_not_empty(file_name):
    """
    Check if file exists and is not empty by confirming that its size is not 0 bytes
    """
    return os.path.isfile(file_name) and not os.path.getsize(file_name) == 0


def concatenate_small(outfile, *args):
    """
    e.g. concatenate_small('test.fastq', 'IDX01_S1_L001_R1_001.fastq', 'IDX01_S1_L001_R2_001.fastq')
    """
    with open(outfile, 'a+') as outfile:
        for filename in args:
            with open(filename, 'r') as infile:
                outfile.write(infile.read())


def gunzip(file):
    """
    Unzips a .gz file unless unzipped file already exists
    """
    expected_unzipped_file = re.sub('.gz', '', file)
    if not file_exists_and_not_empty(expected_unzipped_file):
        with open(expected_unzipped_file, 'w') as outfile:
            with gzip.open(file, 'rt') as infile:
                outfile.write(infile.read())
        os.remove(file)


def check_outgroup_coverage(folder_of_paralog_files, list_of_internal_outgroups, file_of_external_outgroups,
                            list_of_external_outgroups=None):
    """
    Check the number of genes that have an outgroup sequence in either the list_of_internal_outgroups (i.e.
    corresponding to samples within the existing paralog fasta file), or within a file of external outgroup sequences
    (i.e. new taxa to add as outgroups). Writes a report of gene coverage and corresponding outgroup(s).
    """
    # print(list_of_internal_outgroups)
    # print(list_of_external_outgroups)

    # Read in paralog fasta files, and create a dictionary of gene_id:list_of_seq_names:
    paralog_dict = defaultdict(list)
    for fasta in glob.glob(f'{folder_of_paralog_files}/*'):
        gene_id = os.path.basename(fasta).split('.')[0]  # CJJ get prefix e.g. '4471'
        seqs = SeqIO.parse(fasta, 'fasta')
        for seq in seqs:
            paralog_dict[gene_id].append(seq.name)

    # Read in external outgroups file, and create a dictionary of gene_id:list_of_seq_names:
    if file_of_external_outgroups:  # CJJ dict not created if no outgroups file provided
        external_outgroup_dict = defaultdict(list)
        seqs = SeqIO.parse(file_of_external_outgroups, 'fasta')
        for seq in seqs:
            gene_id = seq.name.split('-')[-1]  # CJJ get gene id e.g. '4471'
            external_outgroup_dict[gene_id].append(seq.name)

    # Check whether there is a seq for each taxon in list_of_internal_outgroups for each gene, and create a dictionary
    internal_outgroup_coverage_dict = defaultdict(list)
    if list_of_internal_outgroups:
        for gene, sequence_names in paralog_dict.items():
            internal_outgroup = 0
            for taxon in list_of_internal_outgroups:
                if taxon in sequence_names or f'{taxon}.main' in sequence_names:  # CJJ i.e. if there are paralogs
                    internal_outgroup += 1
                    internal_outgroup_coverage_dict[gene].append(taxon)
                else:
                    # print(f'No seq for gene {gene} for taxon {taxon}')
                    pass
            if internal_outgroup == 0:
                internal_outgroup_coverage_dict[gene].append('No internal outgroup')

    # Check whether there is a seq for each taxon in list_of_external_outgroups for each gene, and create a dictionary
    external_outgroup_coverage_dict = defaultdict(list)
    if file_of_external_outgroups and list_of_external_outgroups:  # CJJ i.e. if filtering the external outgroups for
        # CJJ specified taxa
        for gene, sequences in paralog_dict.items():
            gene_lookup = external_outgroup_dict[gene]
            if len(gene_lookup) == 0:
                external_outgroup_coverage_dict[gene].append('No external outgroup')
            for taxon in list_of_external_outgroups:
                if taxon in ['-'.join(name.split('-')[:-1]) for name in gene_lookup]:  # CJJ i.e the names without
                    # CJJ the gene id suffix)
                    external_outgroup_coverage_dict[gene].append(taxon)
                else:
                    # print(f'Taxon {taxon} not found for gene {gene}!')
                    pass

    if file_of_external_outgroups and not list_of_external_outgroups:  # CJJ i.e. if NOT filtering the external
        # outgroups for specified taxa
        for gene, sequences in paralog_dict.items():
            gene_lookup = external_outgroup_dict[gene]
            if len(gene_lookup) == 0:
                external_outgroup_coverage_dict[gene].append('No external outgroup')
            for seq in gene_lookup:
                external_outgroup_coverage_dict[gene].append('-'.join(seq.split('-')[:-1]))

    # Iterate over all genes from paralogs dict, and check for internal and/or external outgroups. Write a tsv file
    # of results:
    with open(f'outgroup_coverage_report.tsv', 'w') as tsv_report:
        number_of_paralog_files = len(paralog_dict)
        num_paralog_files_with_internal_outgroup = 0
        num_paralog_files_with_external_outgroup = 0
        num_paralog_files_with_both_internal_and_external_outgroup = 0
        tsv_report.write(f'Gene_name\tInternal_outgroup_taxa\tExternal_outgroup_taxa\n')
        for gene, sequences in paralog_dict.items():
            internal_outgroups = internal_outgroup_coverage_dict[gene]
            external_outgroups = external_outgroup_coverage_dict[gene]
            if len(internal_outgroups) != 0 and 'No internal outgroup' not in internal_outgroups:
                # print("LUCY")
                num_paralog_files_with_internal_outgroup += 1
            if len(external_outgroups) != 0 and 'No external outgroup' not in external_outgroups:
                num_paralog_files_with_external_outgroup += 1
            if len(internal_outgroups) != 0 and 'No internal outgroup' not in internal_outgroups and len(
                    external_outgroups) != 0 and 'No external outgroup' not in external_outgroups:
                num_paralog_files_with_both_internal_and_external_outgroup += 1
            tsv_report.write(f'{gene}\t{";".join(internal_outgroups)}\t{";".join(external_outgroups)}\n')

        # Print to stdout, to be captured by Nextflow process to e.g. print a warning
        print(f'*** OUTGROUP COVERAGE STATISTICS ***\n')
        print(f'Number of paralog files with at least one internal outgroup sequence: '
              f'{num_paralog_files_with_internal_outgroup} of {number_of_paralog_files}')
        print(f'Number of paralog files with at least one external outgroup sequence: '
              f'{num_paralog_files_with_external_outgroup} of {number_of_paralog_files}')
        print(f'Number of paralog files with at least one internal AND external outgroup sequence: '
              f'{num_paralog_files_with_both_internal_and_external_outgroup} of {number_of_paralog_files}\n')


def sanitise_gene_names(paralogs_folder, file_of_external_outgroups, sanitised_paralog_output_folder):
    """
    Checks gene names in paralog files and the external outgroup file (if latter provided) for dots, and converts
    them to underscores.
    """

    # Sanitise filenames in input paralogs folder:
    createfolder(sanitised_paralog_output_folder)
    for file in glob.glob(f'{paralogs_folder}/*'):
        basename = os.path.basename(file)
        if not re.search('.paralogs.fasta', basename):
            logger.info(f'\nFile "{basename}" appears not to follow the expected naming convention '
                        f'"geneName.paralogs.fasta". Please check your input files!\n')
            sys.exit(1)
        gene_name = basename.split('.paralogs.fasta')[0]
        gene_name_sanitised = re.sub('[.]', '_', gene_name)
        paralog_filename_sanitised = f'{gene_name_sanitised}.paralogs.fasta'
        shutil.copy(file, f'{sanitised_paralog_output_folder}/{paralog_filename_sanitised}')

    # Sanitise gene names in the external outgroup fasta file, if provided:
    if file_of_external_outgroups:
        basename = os.path.basename(file_of_external_outgroups)
        filename, extension = os.path.splitext(basename)
        sanitised_seqs_to_write = []
        seqs = SeqIO.parse(file_of_external_outgroups, 'fasta')
        for seq in seqs:
            gene_name = seq.id.split('-')[-1]
            sample_name = seq.id.split('-')[0]
            gene_name_sanitised = re.sub('[.]', '_', gene_name)
            seq_name_sanitised = f'{sample_name}-{gene_name_sanitised}'

            # Re-name sequence:
            seq.id = seq_name_sanitised
            seq.name = seq_name_sanitised
            seq.description = seq_name_sanitised

            sanitised_seqs_to_write.append(seq)
        with open(f'{filename}_sanitised{extension}', 'w') as sanitised_outgroups_files:
            SeqIO.write(sanitised_seqs_to_write, sanitised_outgroups_files, 'fasta')

        return sanitised_paralog_output_folder, f'{filename}_sanitised{extension}'

    return sanitised_paralog_output_folder, None


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

    fasta_file_list = glob.glob(f'{gene_fasta_directory}/*.fasta')

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
    parser.add_argument('gene_fasta_directory', type=str, help='directory contains fasta files including paralogs')
    parser.add_argument('-external_outgroups_file', type=str,
                        help='file in fasta format with additional outgroup sequences to add to each gene')
    parser.add_argument('-external_outgroup', action='append', type=str, dest='external_outgroups',
                        help='<Required> Set flag. If one or more taxon names are provided, only use these sequences '
                             'from the user-provided external_outgroups_file', required=False)
    parser.add_argument('-internal_outgroup', action='append', type=str, dest='internal_outgroups',
                        help='<Required> Set flag', required=False, default=None)
    parser.add_argument('-batch_size', type=int, default=20,
                        help='Number of fasta files in each batch, from input paralog fasta files')

    results = parser.parse_args()
    return results


########################################################################################################################
########################################################################################################################
# Run script:

def main():
    args = parse_arguments()
    folder_00 = f'{cwd}/00_paralogs_gene_names_sanitised'
    folder_01 = f'{cwd}/01_batch_folders'


    # Check gene names in paralog files and the external outgroup file (if provided) for dots, and convert to
    # underscores:
    paralogs_folder_sanitised, external_outgroups_file_sanitised = sanitise_gene_names(args.gene_fasta_directory,
                                                                                       args.external_outgroups_file,
                                                                                       folder_00)

    # Check coverage of outgroup sequences:
    check_outgroup_coverage(paralogs_folder_sanitised,
                            args.internal_outgroups,
                            external_outgroups_file_sanitised,
                            list_of_external_outgroups=args.external_outgroups)

    # Batch files into separate folders:
    batch_input_files(paralogs_folder_sanitised,
                      folder_01,
                      batch_size=args.batch_size)


########################################################################################################################
########################################################################################################################

if __name__ == '__main__':
    if not len(sys.argv) >= 1:
        print(__doc__)
        sys.exit()
    sys.exit(main())

########################################################################################################################
########################################################################################################################
