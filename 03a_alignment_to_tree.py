#!/usr/bin/env python

# Author: Chris Jackson

"""
help __doc__
Text Here
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
from collections import defaultdict
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline
from concurrent.futures.process import ProcessPoolExecutor
from multiprocessing import Manager
from concurrent.futures import wait
import subprocess
import shutil

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
logger.setLevel(logging.DEBUG)

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


def iqtree(alignment_file, output_folder, iqtree_threads, counter, lock, num_files_to_process, bootstraps=False):
    """
    Generate trees from alignments using iqtree
    """
    alignment_file_basename = os.path.basename(alignment_file)
    expected_output_file = f'{output_folder}/{alignment_file_basename}.treefile'

    try:
        assert file_exists_and_not_empty(expected_output_file)
        logger.debug(f'Output exists for {expected_output_file}, skipping...')
        with lock:
            counter.value += 1
        return os.path.basename(expected_output_file)
    except AssertionError:
        try:
            if bootstraps:
                subprocess.run(['iqtree', '-redo', '-pre', f'{output_folder}/{alignment_file_basename}', '-s',
                                alignment_file, '-m', 'GTR+G', '-bb', '1000', '-bnni', '-nt',
                                str(iqtree_threads), '-quiet'], check=True)
            else:
                subprocess.run(
                    ['iqtree', '-redo', '-pre', f'{output_folder}/{alignment_file_basename}', '-s', alignment_file,
                     '-m', 'GTR+G', '-nt',
                     str(iqtree_threads), '-quiet'], check=True)

        except:
            logger.info(f'\nNo tree produced for {alignment_file}- fewer than 3 sequences in alignment?\n')
        with lock:
            counter.value += 1
        return os.path.basename(expected_output_file)
    finally:
        print(f'\rFinished generating output {os.path.basename(expected_output_file)}, {counter.value}/{num_files_to_process}', end='')


def iqtree_multiprocessing(alignments_folder, pool_threads=1, iqtree_threads=2, bootstraps=False):
    """
    Generate iqtree trees using multiprocessing.
    """

    input_folder_basename = os.path.basename(alignments_folder)
    output_folder = f'{input_folder_basename}_tree_files'
    createfolder(output_folder)

    logger.info('Generating IQTREE trees from alignments...\n')
    alignments = [file for file in sorted(glob.glob(f'{alignments_folder}/*trimmed.fasta'))]
    alignments.extend([file for file in sorted(glob.glob(f'{alignments_folder}/*trimmed_hmm.fasta'))])
    # print(alignments)

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(iqtree, alignment, output_folder, iqtree_threads, counter, lock,
                                      num_files_to_process=len(alignments), bootstraps=bootstraps)
                          for alignment in alignments]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")
    tree_list = [tree for tree in glob.glob(f'{output_folder}/*.treefile') if file_exists_and_not_empty(tree)]
    logger.info(f'\n{len(tree_list)} trees generated from {len(future_results)} fasta files...\n')


def fasttree_multiprocessing(alignments_folder, pool_threads=1, bootstraps=False):
    """
    Generate FastTree trees using multiprocessing.
    """

    input_folder_basename = os.path.basename(alignments_folder)
    output_folder = f'{input_folder_basename}_tree_files'
    createfolder(output_folder)

    logger.info('Generating FastTree trees from alignments...\n')
    alignments = [file for file in sorted(glob.glob(f'{alignments_folder}/*trimmed.fasta'))]
    alignments.extend([file for file in sorted(glob.glob(f'{alignments_folder}/*trimmed_hmm.fasta'))])
    # print(alignments)

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(fasttree, alignment, output_folder, counter, lock,
                                      num_files_to_process=len(alignments), bootstraps=bootstraps)
                          for alignment in alignments]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")
    tree_list = [tree for tree in glob.glob(f'{output_folder}/*.treefile') if file_exists_and_not_empty(tree)]
    logger.info(f'\n{len(tree_list)} trees generated from {len(future_results)} fasta files...\n')


def fasttree(alignment_file, output_folder, counter, lock, num_files_to_process, bootstraps=False):
    """
    Generate trees from alignments using FastTree
    """
    alignment_file_basename = os.path.basename(alignment_file)
    expected_output_file = f'{output_folder}/{alignment_file_basename}.treefile'

    print('HERE I AM')

    try:
        assert file_exists_and_not_empty(expected_output_file)
        logger.debug(f'Output exists for {expected_output_file}, skipping...')
        with lock:
            counter.value += 1
        return os.path.basename(expected_output_file)
    except AssertionError:
        try:
            if bootstraps:
                fasttree_command = f'FastTreeMP -quote -gtr -nt < {alignment_file} > {expected_output_file}'
                result = subprocess.run(fasttree_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True)
                logger.debug(f'FastTreeMP check_returncode() is: {result.check_returncode()}')
                logger.debug(f'FastTreeMP stdout is: {result.stdout}')
                logger.debug(f'FastTreeMP stderr is: {result.stderr}')

            else:
                fasttree_command = f'FastTreeMP -quote -gtr -nt < {alignment_file} > {expected_output_file}'
                result = subprocess.run(fasttree_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True)
                logger.debug(f'FastTreeMP check_returncode() is: {result.check_returncode()}')
                logger.debug(f'FastTreeMP stdout is: {result.stdout}')
                logger.debug(f'FastTreeMP stderr is: {result.stderr}')

        except subprocess.CalledProcessError as exc:
            logger.error(f'FastTreeMP FAILED. Output is: {exc}')
            logger.error(f'FastTreeMP stdout is: {exc.stdout}')
            logger.error(f'FastTreeMP stderr is: {exc.stderr}')

        # except:
        #     logger.info(f'\nNo tree produced for {alignment_file}- fewer than 3 sequences in alignment?\n')
        with lock:
            counter.value += 1
        return os.path.basename(expected_output_file)
    finally:
        print(f'\rFinished generating output {os.path.basename(expected_output_file)}, {counter.value}'
              f'/{num_files_to_process}', end='')


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('alignment_directory', type=str, help='directory contains fasta files including paralogs')
    parser.add_argument('-threads_pool', type=int, help='<Required> Set flag',
                        required=True)
    parser.add_argument('-threads_iqtree', type=str, help='<Required> Set flag',
                        required=True)
    parser.add_argument('-generate_bootstraps', action='store_true', default=False,
                        help='Create bootstraps for trees using UFBoot')
    parser.add_argument('-use_fasttree', action='store_true', default=False,
                        help='Use FastTree instead of IQTREE')

    results = parser.parse_args()
    return results


########################################################################################################################
########################################################################################################################
# Run script:

def main():
    args = parse_arguments()
    print(args)

    if args.use_fasttree:
        print(f'Running using FastTree!')
        fasttree_multiprocessing(args.alignment_directory,
                                 pool_threads=args.threads_pool,
                                 bootstraps=args.generate_bootstraps)  # Uses OpenMP with max threads default
    else:
        iqtree_multiprocessing(args.alignment_directory,
                               pool_threads=args.threads_pool,
                               iqtree_threads=args.threads_iqtree,
                               bootstraps=args.generate_bootstraps)


########################################################################################################################
########################################################################################################################

if __name__ == '__main__':
    if not len(sys.argv) >= 1:
        print(__doc__)
        sys.exit()
    sys.exit(main())

########################################################################################################################
########################################################################################################################
