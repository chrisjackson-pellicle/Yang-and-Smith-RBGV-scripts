#!/usr/bin/env python

# Author: Chris Jackson

"""
- Checks if there an outgroup (internal or external) for each gene paralog file.
- Aligns the paralog fasta file using mafft, and if the option -no_supercontigs is provided,
  realigns using Clustal Omega (which can do a better job when alignment contains contigs from different regions of
  the full-length reference e.g. split between 5' and 3' halves)
- Trims alignments with Trimal
- Runs HmmCleaner.pl on the alignments.
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
import subprocess
import shutil
import copy
from collections import defaultdict
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline, ClustalOmegaCommandline
from concurrent.futures.process import ProcessPoolExecutor
from multiprocessing import Manager
from concurrent.futures import wait


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
    # Check if file exist and is not empty
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


def done_callback(future_returned):
    """
    Callback function for ProcessPoolExecutor futures; gets called when a future is cancelled or 'done'.
    """
    if future_returned.cancelled():
        logger.info(f'{future_returned}: cancelled')
    elif future_returned.done():
        error = future_returned.exception()
        if error:
            logger.info(f'{future_returned}: error returned: {error}')
        else:
            result = future_returned.result()
    return result


# def get_outgroup_seqs(target_fasta_file, outgroup_list):
#     """
#     Searches a list of outgroup taxon names against a fasta file of target sequences
#     (e.g. the 353 angiosperm target file), and checks if there's an outgroup taxon for each gene. Exits if not.
#     Returns a defaultdict of gene_id:outgroup_Seq_object_list
#     """
#     unique_gene_names = set()
#     outgroup_lists = defaultdict(list)
#     with open(target_fasta_file, 'r') as target_file:
#         seqs = SeqIO.parse(target_file, 'fasta')
#         for seq in seqs:
#             taxon_name = re.split('-', seq.name)[0]
#             gene_id = re.split('-', seq.name)[-1]
#             unique_gene_names.add(gene_id)
#             if taxon_name in outgroup_list:
#                 outgroup_lists[gene_id].append(seq)
#     if len(outgroup_lists) == len(unique_gene_names):
#         logger.info(f'Found an outgroup for each of {len(unique_gene_names)} genes in the targetfile provided...')
#     else:
#         # sys.exit(f'Only found an outgroup for {len(outgroup_lists)} of {len(unique_gene_names)} genes in the '
#         #          f'targetfile provided! Perhaps add more species name for outgroup selection?')
#         logger.info(f'Only found an outgroup for {len(outgroup_lists)} of {len(unique_gene_names)} genes in the '
#                     f'targetfile provided! Perhaps add more species name for outgroup selection?')
#     return outgroup_lists


# def no_paralogs(paralog_directory, output_folder):
#     """
#     Recovers main or single sequences only - no paralogs (i.e seq.0, seq.1)
#     """
#     createfolder(output_folder)
#     for fasta_file in glob.glob(f'{paralog_directory}/*.fasta'):
#         gene_id = os.path.basename(fasta_file).split('.')[0]
#         seqs_to_recover = []
#         seq_dict = defaultdict(list)
#         seqs = SeqIO.parse(fasta_file, 'fasta')
#         for seq in seqs:
#             seq_dict[seq.name.split('.')[0]].append(seq)
#         for key, value in seq_dict.items():
#             for sequence in value:
#                 try:
#                     paralog_id = sequence.name.split('.')[1]
#                     if paralog_id == 'main':
#                         seqs_to_recover.append(sequence)
#                 except IndexError:
#                     seqs_to_recover.append(sequence)
#         with open(f'{output_folder}/{gene_id}.no_paralogs.fasta', 'w') as no_paralog_file:
#             SeqIO.write(seqs_to_recover, no_paralog_file, 'fasta')


def run_hmm_cleaner(input_folder, output_folder):
    """
    Runs HmmCleaner.pl on each alignments within a provided folder.
    """
    createfolder(output_folder)
    for alignment in glob.glob(f'{input_folder}/*.aln.trimmed.fasta'):
        # command = f'HmmCleaner.pl {alignment} -costs -0.99 -0.01 0.01 0.1 -symfrac 0.0'
        # command = f'HmmCleaner.pl {alignment} --large'

        command = f'/usr/bin/perl /usr/local/bin/HmmCleaner.pl {alignment}'

        if host == 'RBGs-MacBook-Air.local':
            command = f'/Users/chrisjackson/perl5/perlbrew/perls/perl-5.26.2/bin/perl ' \
                      f'/Users/chrisjackson/perl5/perlbrew/perls/perl-5.26.2/bin/HmmCleaner.pl {alignment}'

        # print(command)
        # run = subprocess.run(command, shell=True, check=True, capture_output=True)
        # print(run)
        try:
            run = subprocess.run(command, shell=True, check=True, capture_output=True)

            # Filter out empty sequences comprised only of dashes
            try:
                logger.debug(f'trying command {command}')
                seqs_to_retain = []
                hmm_file = re.sub('aln.trimmed.fasta', 'aln.trimmed_hmm.fasta', str(alignment))
                hmm_file_output = re.sub('aln.trimmed.fasta', 'aln.hmm.trimmed.fasta', str(alignment))
                # print(f'hmm_file is {hmm_file}')
                with open(hmm_file, 'r') as hmm_fasta:
                    seqs = SeqIO.parse(hmm_fasta, 'fasta')
                    for seq in seqs:
                        characters = set(character for character in seq.seq)
                        if len(characters) == 1 and '-' in characters:
                            pass
                        else:
                            seqs_to_retain.append(seq)
                with open(hmm_file_output, 'w') as filtered_hmm_fasta:
                    SeqIO.write(seqs_to_retain, filtered_hmm_fasta, 'fasta')
            except:
                pass
        except subprocess.CalledProcessError:
            hmm_file_output = re.sub('aln.trimmed.fasta', 'aln.hmm.trimmed.fasta', str(alignment))
            logger.info(f"Couldn't run HmmCleaner for alignment {alignment} using command {command}")
            logger.info(f'Copying alignment {alignment} to {hmm_file_output} anyway...')
            shutil.copy(alignment, hmm_file_output)

    # for file in glob.glob(f"{input_folder}/*hmm*"):
    for file in glob.glob(f"{input_folder}/*aln.hmm.trimmed*"):
        try:
            shutil.move(file, output_folder)
        except shutil.Error:
            # raise
            pass


# def append_outgroup_seqs(paralog_directory, outgroups_dict, output_folder):
#     """
#     Adds an outgroup sequence recovered from the target file to the fasta gene list recovered for each gene via
#     paralog_retriever.
#     """
#     createfolder(output_folder)
#     # print(outgroups_dict)
#     for fasta_file in glob.glob(f'{paralog_directory}/*.fasta'):
#         basename = os.path.basename(fasta_file)
#         gene_id = basename.split('.')[0]
#         # seqs = list(SeqIO.parse(fasta_file, 'fasta'))  # Turn this into a dictionary instead, so I can check search
#         seqs_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
#         # if seqs:  # Don't process HybPiper paralog fasta files without sequences
#         if seqs_dict:  # Don't process HybPiper paralog fasta files without sequences
#             for seq in outgroups_dict[gene_id]:
#                 seq_taxon_name = seq.name.split('-')[0]
#                 try:
#                     seqs_dict[seq_taxon_name].name = f'{seq.name.split("-")[0]}.outgroup'
#                     seqs_dict[seq_taxon_name].id = f'{seq.id.split("-")[0]}.outgroup'
#                     seqs_dict[seq_taxon_name].description = ''
#                 except KeyError:  # i.e. the outgroup sequence isn't already in the fast file
#                     print(f'Sequence {seq_taxon_name} is not already in the fasta file')
#                     seq.name = f'{seq.name.split("-")[0]}.outgroup'
#                     seq.id = f'{seq.id.split("-")[0]}.outgroup'
#                     seq.description = ''
#                     seqs_dict[seq.name] = seq
#             with open(f'{output_folder}/{gene_id}.outgroup_added.fasta', 'w') as outgroup_added:
#                 SeqIO.write(seqs_dict.values(), outgroup_added, 'fasta')


def mafft_align(fasta_file, algorithm, output_folder, counter, lock, num_files_to_process, threads=2,
                no_supercontigs=False):
    """
    Uses mafft to align a fasta file of sequences, using the algorithm and number of threads provided. Trims
    alignment with Trimal if no_supercontigs=True. Returns filename of the alignment produced.
    """
    createfolder(output_folder)
    fasta_file_basename = os.path.basename(fasta_file)
    expected_alignment_file = f'{output_folder}/{re.sub(".fasta", ".aln.fasta", fasta_file_basename)}'

    try:
        assert file_exists_and_not_empty(expected_alignment_file)
        logger.debug(f'Alignment exists for {fasta_file_basename}, skipping...')
        with lock:
            counter.value += 1
        return os.path.basename(expected_alignment_file)
    except AssertionError:
        mafft_cline = (MafftCommandline(algorithm, adjustdirection='true', thread=threads, input=fasta_file))
        stdout, stderr = mafft_cline()
        with open(expected_alignment_file, 'w') as alignment_file:
            alignment_file.write(stdout)

        if not no_supercontigs:
            trimmed_alignment = re.sub('.aln.fasta', '.aln.trimmed.fasta', expected_alignment_file)
            run_trim = subprocess.run(['trimal', '-in', expected_alignment_file, '-out', trimmed_alignment,
                                       '-gapthreshold', '0.12', '-terminalonly', '-gw', '1'], check=True)
        with lock:
            counter.value += 1
        logger.debug(f'Aligned file {fasta_file_basename}')
        return os.path.basename(expected_alignment_file)
    finally:
        logger.debug(f'\rFinished generating alignment for file {fasta_file_basename}, '
              f'{counter.value}/{num_files_to_process}', end='')


def mafft_align_multiprocessing(fasta_to_align_folder, alignments_output_folder, algorithm='linsi', pool_threads=1,
                                mafft_threads=2, no_supercontigs=False):
    """
    Generate alignments via function <align_targets> using multiprocessing.
    """
    createfolder(alignments_output_folder)
    logger.debug('Generating alignments for fasta files using mafft...\n')
    target_genes = [file for file in sorted(glob.glob(f'{fasta_to_align_folder}/*.fasta'))]

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(mafft_align, fasta_file, algorithm, alignments_output_folder, counter, lock,
                                      num_files_to_process=len(target_genes), threads=mafft_threads,
                                      no_supercontigs=no_supercontigs)
                          for fasta_file in target_genes]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")
    alignment_list = [alignment for alignment in glob.glob(f'{alignments_output_folder}/*.aln.fasta') if
                      file_exists_and_not_empty(alignment)]
    logger.debug(f'\n{len(alignment_list)} alignments generated from {len(future_results)} fasta files...\n')


def clustalo_align(fasta_file, output_folder, counter, lock, num_files_to_process, threads=2):
    """
    Uses clustal omega to align a fasta file of sequences, using the algorithm and number of threads provided.
    Trims alignment with Trimal. Returns filename of the alignment produced.
    """
    createfolder(output_folder)
    fasta_file_basename = os.path.basename(fasta_file)
    expected_alignment_file = f'{output_folder}/{fasta_file_basename}'

    try:
        assert file_exists_and_not_empty(expected_alignment_file)
        logger.debug(f'Alignment exists for {fasta_file_basename}, skipping...')
    except AssertionError:
        # clustalomega_cline = ClustalOmegaCommandline(infile=fasta_file, outfile=expected_alignment_file,
        #                                              verbose=True, auto=True, threads=threads, pileup=True)
        clustalomega_cline = ClustalOmegaCommandline(infile=fasta_file, outfile=expected_alignment_file,
                                                     verbose=True, auto=True, threads=threads)
        clustalomega_cline()
        trimmed_alignment = re.sub('.aln.fasta', '.aln.trimmed.fasta', expected_alignment_file)
        run_trim = subprocess.run(['trimal', '-in', expected_alignment_file, '-out', trimmed_alignment, '-gapthreshold',
                                   '0.12', '-terminalonly', '-gw', '1'], check=True)
    finally:
        with lock:
            counter.value += 1
            logger.debug(f'\rFinished generating alignment for file {fasta_file_basename}, '
                  f'{counter.value}/{num_files_to_process}', end='')
        return os.path.basename(expected_alignment_file)


def clustalo_align_multiprocessing(fasta_to_align_folder, alignments_output_folder, pool_threads=1,
                                   clustalo_threads=2):
    """ 
    Generate alignments via function <clustalo_align> using multiprocessing.
    """
    createfolder(alignments_output_folder)
    logger.debug('Generating alignments for fasta files using clustal omega...\n')
    target_genes = [file for file in sorted(glob.glob(f'{fasta_to_align_folder}/*.fasta'))]

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(clustalo_align, fasta_file, alignments_output_folder, counter, lock,
                                      num_files_to_process=len(target_genes), threads=clustalo_threads)
                          for fasta_file in target_genes]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")
    alignment_list = [alignment for alignment in glob.glob(f'{alignments_output_folder}/*.aln.fasta') if
                      file_exists_and_not_empty(alignment)]
    logger.debug(f'\n{len(alignment_list)} alignments generated from {len(future_results)} fasta files...\n')


# def macse_align(fasta_file, counter, lock, num_files_to_process, output_folder=None):
#     """uses macse to align a fasta file of nucleotide sequences via protein translation, using the number of threads
#     provided. If file exists, returns without aligning"""
#     createfolder(output_folder)
#     fasta_file_basename = os.path.basename(fasta_file)
#     expected_alignment_file = f'{output_folder}/{re.sub(".fasta", ".aln.fasta", fasta_file_basename)}'
#
#     try:
#         assert file_exists_and_not_empty(expected_alignment_file)
#         logger.debug(f'Alignment exists for {fasta_file_basename}, skipping...')
#     except AssertionError:
#         run = subprocess.run(['java', '-jar', '/home/chris/01_installs/macse_v2.04.jar', '-prog', 'alignSequences',
#                               '-seq', fasta_file, '-out_NT', expected_alignment_file], check=True)
#     finally:
#         with lock:
#             counter.value += 1
#             print(f'\rFinished generating alignment for file {fasta_file_basename}, '
#                   f'{counter.value}/{num_files_to_process}', end='')
#             return os.path.basename(expected_alignment_file)


# def macse_align_multiprocessing(fasta_to_align_folder, alignments_output_folder, pool_threads=1):
#     """
#     Generate alignments via function <macse_align> using multiprocessing.
#     """
#     createfolder(alignments_output_folder)
#     logger.info('Generating alignments for fasta file without paralogs...\n')
#     target_genes = [file for file in sorted(glob.glob(f'{fasta_to_align_folder}/*.fasta'))]
#
#     with ProcessPoolExecutor(max_workers=pool_threads) as pool:
#         manager = Manager()
#         lock = manager.Lock()
#         counter = manager.Value('i', 0)
#         future_results = [pool.submit(macse_align, fasta_file, counter, lock, num_files_to_process=len(target_genes),
#                                       output_folder=alignments_output_folder) for fasta_file in target_genes]
#         for future in future_results:
#             future.add_done_callback(done_callback)
#         wait(future_results, return_when="ALL_COMPLETED")
#     alignment_list = [alignment for alignment in glob.glob(f'{alignments_output_folder}/*.aln.fasta') if
#                       file_exists_and_not_empty(alignment)]
#     logger.info(f'\n{len(alignment_list)} alignments generated from {len(future_results)} fasta files...\n')


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
        # specified taxa
        for gene, sequences in paralog_dict.items():
            gene_lookup = external_outgroup_dict[gene]
            if len(gene_lookup) == 0:
                external_outgroup_coverage_dict[gene].append('No external outgroup')
            for taxon in list_of_external_outgroups:
                if taxon in ['-'.join(name.split('-')[:-1]) for name in gene_lookup]:  # CJJ i.e the names without
                    # the gene id suffix)
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
            if 'No internal outgroup' not in internal_outgroups:
                num_paralog_files_with_internal_outgroup += 1
            if 'No external outgroup' not in external_outgroups:
                num_paralog_files_with_external_outgroup += 1
            if 'No internal outgroup' not in internal_outgroups and 'No external outgroup' not in external_outgroups:
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
    parser.add_argument('-pool', type=int, help='Number of threads to use for the multiprocessing pool', default=1)
    parser.add_argument('-threads', type=int, help='Number of threads to use for the multiprocessing pool function',
                        default=1)
    parser.add_argument('-skip_hmmcleaner', action='store_true', help='skip the hmmcleaner step')
    parser.add_argument('-no_supercontigs', action='store_true', default=False,
                        help='If specified, realign mafft alignments with clustal omega')

    results = parser.parse_args()
    return results


########################################################################################################################
########################################################################################################################
# Run script:

def main():
    results = parse_arguments()
    # folder_01 = f'{cwd}/01_outgroup_added'
    folder_01a = f'{cwd}/01a_mafft_alignments'
    folder_01b = f'{cwd}/01_alignments'
    folder_02 = f'{cwd}/02_alignments_hmmcleaned'

    # Check coverage of outgroup sequences
    check_outgroup_coverage(results.gene_fasta_directory, results.internal_outgroups, results.external_outgroups_file,
                            list_of_external_outgroups=results.external_outgroups)

    if not results.no_supercontigs:  # i.e. if it's a standard run with supercontigs produced.
        logger.debug(f'Running without no_supercontigs option - aligning with mafft only')
        mafft_align_multiprocessing(results.gene_fasta_directory, folder_01b, algorithm='linsi',
                                    pool_threads=results.pool,
                                    mafft_threads=results.threads, no_supercontigs=False)
        run_hmm_cleaner(folder_01b, folder_02)
    elif results.no_supercontigs:  # Re-align with Clustal Omega.
        logger.debug(f'Running with no_supercontigs option - realigning with clustal omega')
        mafft_align_multiprocessing(results.gene_fasta_directory, folder_01a, algorithm='linsi',
                                    pool_threads=results.pool,
                                    mafft_threads=results.threads, no_supercontigs=True)
        clustalo_align_multiprocessing(folder_01a, folder_01b, pool_threads=results.pool,
                                       clustalo_threads=results.threads)
        run_hmm_cleaner(folder_01b, folder_02)


########################################################################################################################
########################################################################################################################

if __name__ == '__main__':
    if not len(sys.argv) >= 1:
        print(__doc__)
        sys.exit()
    sys.exit(main())

########################################################################################################################
########################################################################################################################
