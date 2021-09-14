#!/usr/bin/env python

# Author: Chris Jackson

"""
- Checks gene names in paralog files and the external outgroup file (if provided) for dots, and converts them to
  underscores.
- Checks if there an outgroup (internal or external) for each gene paralog file.
- Aligns the paralog fasta file using mafft, and if the option -no_supercontigs is provided,
  realigns using Clustal Omega (which can do a better job when alignment contains contigs from different regions of
  the full-length reference e.g. split between 5' and 3' halves).
- Trims alignments with Trimal.
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
        return
    elif future_returned.done():
        error = future_returned.exception()
        if error:
            logger.info(f'{future_returned}: error returned: {error}')
        else:
            result = future_returned.result()
            return result
    # return result


def run_hmm_cleaner(input_folder, output_folder):
    """
    Runs HmmCleaner.pl on each alignments within a provided folder.
    """
    createfolder(output_folder)
    for alignment in glob.glob(f'{input_folder}/*.aln.trimmed.fasta'):
        command = f'/usr/bin/perl /usr/local/bin/HmmCleaner.pl {alignment}'

        # print(f'hostname is: {host}')
        # if host == '192-168-1-111.tpgi.com.au':
        if host == 'RBGs-MacBook-Air.local':
            command = f'/Users/chrisjackson/perl5/perlbrew/perls/perl-5.26.2/bin/perl ' \
                      f'/Users/chrisjackson/perl5/perlbrew/perls/perl-5.26.2/bin/HmmCleaner.pl {alignment}'
        try:
            run = subprocess.run(command, shell=True, check=True, capture_output=True)

            # Filter out empty sequences comprised only of dashes
            try:
                logger.debug(f'trying command {command}')
                seqs_to_retain = []
                hmm_file = re.sub('aln.trimmed.fasta', 'aln.trimmed_hmm.fasta', str(alignment))
                hmm_file_output = re.sub('aln.trimmed.fasta', 'aln.hmm.trimmed.fasta', str(alignment))
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

    for file in glob.glob(f"{input_folder}/*aln.hmm.trimmed*"):
        try:
            shutil.move(file, output_folder)
        except shutil.Error:
            # raise
            pass


def remove_r_prefix(alignment):
    """
    Takes a fasta alignment, removes any '_R_' prefix in fasta headers (inserted by mafft if a sequences was
    reversed) and writes a new alignment to the same filename.
    """
    with open(alignment) as alignment_handle:
        alignment_obj = AlignIO.read(alignment_handle, 'fasta')
        for seq in alignment_obj:
            if seq.name.startswith('_R_'):
                seq.name = seq.name.lstrip('_R_')
                seq.id = seq.id.lstrip('_R_')
        with open(alignment, 'w') as new_alignment_handle:
            AlignIO.write(alignment_obj, new_alignment_handle, 'fasta')


def mafft_align(fasta_file, algorithm, output_folder, counter, lock, num_files_to_process, threads=2,
                no_supercontigs=False, use_muscle=False):
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
        if use_muscle:
            # print(fasta_file_basename)
            logger.info('Alignment will be performed using MUSCLE rather than MAFFT!')
            muscle_cline = MuscleCommandline(input=fasta_file, out=expected_alignment_file)
            stdout, stderr = muscle_cline()
        else:
            if algorithm == 'auto':
                mafft_cline = (MafftCommandline(auto='true', adjustdirection='true', thread=threads, input=fasta_file))
            else:
                mafft_cline = (MafftCommandline(algorithm, adjustdirection='true', thread=threads, input=fasta_file))
            logger.info(f'Performing MAFFT alignment with command" {mafft_cline}')
            stdout, stderr = mafft_cline()
            with open(expected_alignment_file, 'w') as alignment_file:
                alignment_file.write(stdout)
            remove_r_prefix(expected_alignment_file)

        if not no_supercontigs:
            trimmed_alignment = re.sub('.aln.fasta', '.aln.trimmed.fasta', expected_alignment_file)
            run_trim = subprocess.run(['trimal', '-in', expected_alignment_file, '-out', trimmed_alignment,
                                       '-gapthreshold', '0.12', '-terminalonly', '-gw', '1'], check=True)
        with lock:
            counter.value += 1
        logger.debug(f'Aligned file {fasta_file_basename}')
        return os.path.basename(expected_alignment_file)
    finally:
        logger.debug(f'\rFinished generating alignment for file {fasta_file_basename}, {counter.value}'
                     f'/{num_files_to_process}', end='')


def mafft_align_multiprocessing(fasta_to_align_folder, alignments_output_folder, algorithm='linsi', pool_threads=1,
                                mafft_threads=2, no_supercontigs=False, use_muscle=False, mafft_algorithm_auto=False):
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
        future_results = [pool.submit(mafft_align,
                                      fasta_file,
                                      algorithm,
                                      alignments_output_folder,
                                      counter, lock,
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
        remove_r_prefix(expected_alignment_file)
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


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gene_fasta_directory', type=str, help='directory containing fasta files (with '
                                                               'sanitised gene names) including paralogs')
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
    parser.add_argument('-use_muscle', action='store_true', default=False,
                        help='If specified, use muscle rather than mafft')
    # group_1 = parser.add_mutually_exclusive_group()
    parser.add_argument('-mafft_algorithm', default='auto', help='Algorithm to use for mafft alignments')
    # group_1.add_argument('-mafft_algorithm_auto', action='store_true', default=False, help='use --auto for mafft '
    #                                                                                        'algorithm')

    results = parser.parse_args()
    return results


########################################################################################################################
########################################################################################################################
# Run script:

def main():
    args = parse_arguments()
    folder_01a = f'{cwd}/01a_mafft_alignments'
    folder_01b = f'{cwd}/01_alignments'
    folder_02 = f'{cwd}/02_alignments_hmmcleaned'

    if not args.no_supercontigs:  # i.e. if it's a standard run with supercontigs produced.
        logger.debug(f'Running without no_supercontigs option - aligning with mafft or muscle only')
        mafft_align_multiprocessing(args.gene_fasta_directory,
                                    folder_01b,
                                    algorithm=args.mafft_algorithm,
                                    pool_threads=args.pool,
                                    mafft_threads=args.threads,
                                    no_supercontigs=False,
                                    use_muscle=args.use_muscle)

        run_hmm_cleaner(folder_01b,
                        folder_02)

    elif args.no_supercontigs:  # Re-align with Clustal Omega.
        logger.debug(f'Running with no_supercontigs option - realigning with clustal omega')
        mafft_align_multiprocessing(args.gene_fasta_directory,
                                    folder_01a,
                                    algorithm=args.mafft_algorithm,
                                    pool_threads=args.pool,
                                    mafft_threads=args.threads,
                                    no_supercontigs=True,
                                    use_muscle=args.use_muscle)

        clustalo_align_multiprocessing(folder_01a,
                                       folder_01b,
                                       pool_threads=args.pool,
                                       clustalo_threads=args.threads)

        run_hmm_cleaner(folder_01b,
                        folder_02)


########################################################################################################################
########################################################################################################################

if __name__ == '__main__':
    if not len(sys.argv) >= 1:
        print(__doc__)
        sys.exit()
    sys.exit(main())

########################################################################################################################
########################################################################################################################
