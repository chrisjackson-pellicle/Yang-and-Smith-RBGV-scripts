#!/usr/bin/env python

# Author: Chris Jackson

"""
Aligns the fasta file using mafft.

########################################################################################################################
Additional information:

Text Here
########################################################################################################################

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
from collections import defaultdict
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline, ClustalOmegaCommandline, MuscleCommandline
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


def multiprocessing(input_folder, output_folder, pool_threads=1):
    """
    General function for multiprocessing. Adds a callback function to each future.
    """

    createfolder(output_folder)
    to_process_list = [file for file in sorted(glob.glob(f'{input_folder}'))]

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(function_x, function_arg_1, counter, lock,
                                      num_files_to_process=len(to_process_list))
                          for function_arg_1 in to_process_list]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")


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


def mafft_align(fasta_file, algorithm, output_folder, counter, lock, num_files_to_process, threads=2,
                no_supercontigs=False, use_muscle=False):
    """
    Uses mafft to align a fasta file of sequences, using the algorithm and number of threads provided. Returns filename
    of the alignment produced.
    """

    fasta_file_basename = os.path.basename(fasta_file)
    expected_alignment_file = f'{output_folder}/{re.sub(".fasta", ".aln.fasta", fasta_file_basename)}'

    try:
        assert file_exists_and_not_empty(expected_alignment_file)
        logger.debug(f'Alignment exists for {fasta_file_basename}, skipping...')
        with lock:
            counter.value += 1
        return os.path.basename(expected_alignment_file)
    except AssertionError:
        if use_muscle:
            # print(fasta_file_basename)
            logger.info('Alignment will be performed using MUSCLE rather than MAFFT!')
            muscle_cline = MuscleCommandline(input=fasta_file, out=expected_alignment_file)
            stdout, stderr = muscle_cline()
        else:
            mafft_cline = (MafftCommandline(algorithm, adjustdirection='true', thread=threads, input=fasta_file))
            stdout, stderr = mafft_cline()
            with open(expected_alignment_file, 'w') as alignment_file:
                alignment_file.write(stdout)
        with lock:
            counter.value += 1
        logger.debug(f'Aligned file {fasta_file_basename}')
        return os.path.basename(expected_alignment_file)
    finally:
        print(f'\rFinished generating alignment for file {fasta_file_basename}, '
              f'{counter.value}/{num_files_to_process}', end='')


def mafft_align_multiprocessing(fasta_to_align_folder, alignments_output_folder, algorithm='linsi', pool_threads=1,
                                mafft_threads=2, no_supercontigs=False, use_muscle=False):
    """
    Generate alignments via function <align_targets> using multiprocessing.
    """

    input_folder_basename = os.path.basename(fasta_to_align_folder)
    output_folder = f'{input_folder_basename}_alignments'
    # print(f'output_folder from mafft_align_multiprocessing: {output_folder}')
    createfolder(output_folder)

    logger.info(f'Generating alignments for fasta files in folder {fasta_to_align_folder}...\n')
    target_genes = [file for file in sorted(glob.glob(f'{fasta_to_align_folder}/*.fasta'))]

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(mafft_align,
                                      fasta_file,
                                      algorithm,
                                      alignments_output_folder,
                                      counter,
                                      lock,
                                      num_files_to_process=len(target_genes),
                                      threads=mafft_threads,
                                      no_supercontigs=no_supercontigs,
                                      use_muscle=use_muscle)
                          for fasta_file in target_genes]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")
    alignment_list = [alignment for alignment in glob.glob(f'{alignments_output_folder}/*.aln.fasta') if
                      file_exists_and_not_empty(alignment)]
    logger.info(f'\n{len(alignment_list)} alignments generated from {len(future_results)} fasta files...\n')

    return output_folder


def clustalo_align(fasta_file, output_folder, counter, lock, num_files_to_process, threads=2):
    """
    Uses clustal omega to align a fasta file of sequences, using the algorithm and number of threads provided. Returns
    filename
    of the alignment produced.
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
    finally:
        with lock:
            counter.value += 1
            print(f'\rFinished generating alignment for file {fasta_file_basename}, '
                  f'{counter.value}/{num_files_to_process}', end='')
        return os.path.basename(expected_alignment_file)


def clustalo_align_multiprocessing(fasta_to_align_folder, alignments_output_folder, pool_threads=1,
                                clustalo_threads=2):
    """
    Generate alignments via function <clustalo_align> using multiprocessing.
    """

    input_folder_basename = os.path.basename(fasta_to_align_folder)
    output_folder = f'{input_folder_basename}_clustal'
    # print(f'output_folder from clustal_align_multiprocessing: {output_folder}')
    createfolder(output_folder)

    logger.info('Generating alignments for fasta files using clustal omega...\n')
    target_genes = [file for file in sorted(glob.glob(f'{fasta_to_align_folder}/*.fasta'))]

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(clustalo_align,
                                      fasta_file,
                                      alignments_output_folder,
                                      counter,
                                      lock,
                                      num_files_to_process=len(target_genes),
                                      threads=clustalo_threads)
                          for fasta_file in target_genes]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")
    alignment_list = [alignment for alignment in glob.glob(f'{alignments_output_folder}/*.aln.fasta') if
                      file_exists_and_not_empty(alignment)]
    logger.info(f'\n{len(alignment_list)} alignments generated from {len(future_results)} fasta files...\n')

    return output_folder


def strip_names_for_concat(input_folder):
    """
    Strips everything after a dot ('.') from the name of each sequence in an alignment file.
    Returns the name of the output folder.
    """

    input_folder_basename = os.path.basename(input_folder)
    output_folder = f'{input_folder_basename}_stripped_names'
    # print(f'output_folder from mafft_align_multiprocessing: {output_folder}')
    createfolder(output_folder)

    for alignment in glob.glob(f'{input_folder}/*.fa'):
        alignment_basename = os.path.basename(alignment)
        seqs = AlignIO.read(alignment, 'fasta')
        for seq in seqs:
            seq.name = f'{seq.name.split(".")[0]}'
            seq.id = f'{seq.id.split(".")[0]}'
            seq.description = ''
        with open(f'{output_folder}/{re.sub(".fa", "_stripped.fasta", alignment_basename)}', 'w') as stripped:
            AlignIO.write(seqs, stripped, 'fasta')
    return output_folder


def parse_arguments():
    parser = argparse.ArgumentParser(description=print(__doc__),
                                     # formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog='text here')
    parser.add_argument('gene_fasta_directory', type=str, help='directory contains fasta files including paralogs')
    parser.add_argument('-threads_pool', type=int, help='Number of threads to use for the multiprocessing pool',
                        required=True)
    parser.add_argument('-threads_mafft', type=int, help='Number of threads to use for mafft',
                        required=True)
    parser.add_argument('-no_supercontigs', action='store_true', default=False, help='If specified, realign mafft '
                                                                                     'alignments with clustal omega')
    parser.add_argument('-use_muscle', action='store_true', default=False,
                        help='If specified, use muscle rather than mafft')
    parser.add_argument('-mafft_algorithm', default='auto', help='Algorithm to use for mafft alignments')

    results = parser.parse_args()
    return results


########################################################################################################################
########################################################################################################################
# Run script:

def main():
    results = parse_arguments()

    output_folder = strip_names_for_concat(results.gene_fasta_directory)

    if not results.no_supercontigs:  # i.e. if it's a standard run.
        alignments_output_folder = mafft_align_multiprocessing(output_folder,
                                                               algorithm=results.mafft_algorithm,
                                                               pool_threads=results.threads_pool,
                                                               mafft_threads=results.threads_mafft,
                                                               no_supercontigs=results.no_supercontigs,
                                                               use_muscle=results.use_muscle)

    elif results.no_supercontigs:  # re-align with Clustal Omega.
         alignments_output_folder = mafft_align_multiprocessing(output_folder,
                                                                algorithm=results.mafft_algorithm,
                                                                pool_threads=results.threads_pool,
                                                                mafft_threads=results.threads_mafft,
                                                                no_supercontigs=results.no_supercontigs,
                                                                use_muscle=results.use_muscle)

        clustal_alignment_output_folder = clustalo_align_multiprocessing(alignments_output_folder,
                                                                         pool_threads=results.threads_pool,
                                                                         clustalo_threads=results.threads_mafft)


########################################################################################################################
########################################################################################################################

if __name__ == '__main__':
    if not len(sys.argv) >= 1:
        print(__doc__)
        sys.exit()
    sys.exit(main())

########################################################################################################################
########################################################################################################################
