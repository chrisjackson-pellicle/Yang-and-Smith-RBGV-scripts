#!/usr/bin/env python

# Author: Chris Jackson

"""
Takes a folder of tree files as input. Reads each tree newick file, searches for any polytomies, resolves them
arbitrarily using ete3, unroots the tree, and writes resolved trees to a new folder.
"""

import logging
import sys
import argparse
import os
import socket
from ete3 import Tree
import datetime
import glob

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")


########################################################################################################################
########################################################################################################################
# Get current working directory and host name

cwd = os.getcwd()
host = socket.gethostname()

# Configure logger:
def setup_logger(name, log_file, console_level=logging.INFO, file_level=logging.DEBUG,
                 logger_object_level=logging.DEBUG):
    """
    Function to create a logger instance.

    By default, logs level DEBUG and above to file.
    By default, logs level INFO and above to stderr and file.

    :param string name: name for the logger instance
    :param string log_file: filename for log file
    :param string console_level: level for logging to console
    :param string file_level: level for logging to file
    :param string logger_object_level: level for logger object
    :return: a logger object
    """

    # Get date and time string for log filename:
    date_and_time = datetime.datetime.now().strftime("%Y-%m-%d-%H_%M_%S")

    # Log to file:
    file_handler = logging.FileHandler(f'{log_file}_{date_and_time}.log', mode='w')
    file_handler.setLevel(file_level)
    file_format = logging.Formatter('%(asctime)s - %(filename)s - %(name)s - %(funcName)s - %(levelname)s - %('
                                    'message)s')
    file_handler.setFormatter(file_format)

    # Log to Terminal (stdout):
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(console_level)
    console_format = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_format)

    # Setup logger:
    logger_object = logging.getLogger(name)
    logger_object.setLevel(logger_object_level)  # Default level is 'WARNING'

    # Add handlers to the logger
    logger_object.addHandler(console_handler)
    logger_object.addHandler(file_handler)

    return logger_object


# Create logger(s):
logger = setup_logger(__name__, 'resolve_polytomies')

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


def resolve_polytomies(treefile_directory, resolved_treefile_output_directory):
    """
    Iterates over tree newick files in a directory. For each tree, checks to polyotomies and arbitrarily resolves them.

    """

    createfolder(resolved_treefile_output_directory)

    for treefile in glob.glob(f'{treefile_directory}/*.treefile'):
        tree_basename = os.path.basename(treefile)
        file_name, ext = os.path.splitext(tree_basename)

        # Randomly resolve any polytomies (as FastTree can generate polytomies):
        ete3_tree = Tree(treefile)
        polytomy = False

        for node in ete3_tree.iter_descendants():
            if len(node.children) > 2:
                logger.info(f'Tree {tree_basename} has a polytomy at node {node}')
                polytomy = True

        if polytomy:
            logger.info(f'Polytomies in tree {tree_basename} will be randomly resolved!')

            ete3_tree.resolve_polytomy(recursive=True, default_dist='0.01')
            for node in ete3_tree.iter_descendants():
                if len(node.children) > 2:
                    raise ValueError(f'Tree {tree_basename} still has a polytomy at node {node}!')
            ete3_tree.unroot()
            ete3_tree.write(outfile=f'{resolved_treefile_output_directory}/{file_name}_polytomies_resolved'
                                    f'{ext}', format=0)

        else:
            logger.info(f'Tree {tree_basename} has no polytomies - writing unchanged to file'
                        f' {file_name}_polytomies_resolved{ext}')
            ete3_tree.write(outfile=f'{resolved_treefile_output_directory}/{file_name}_polytomies_resolved'
                                    f'{ext}', format=0)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('treefile_directory', type=str, help='directory contains newick tree files')
    parser.add_argument('--resolved_treefile_output_directory',
                        default='treefiles_polytomies_resolved',
                        type=str, help='Output directory for resolved newick =tree files')

    results = parser.parse_args()
    return results


########################################################################################################################
########################################################################################################################
# Run script:

def main():
    args = parse_arguments()
    print(args)

    resolve_polytomies(args.treefile_directory, args.resolved_treefile_output_directory)


########################################################################################################################
########################################################################################################################

if __name__ == '__main__':
    if not len(sys.argv) >= 1:
        print(__doc__)
        sys.exit()
    sys.exit(main())

########################################################################################################################
########################################################################################################################
