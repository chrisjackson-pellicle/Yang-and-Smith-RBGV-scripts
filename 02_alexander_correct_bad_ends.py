# Author: Alexander Schmidt-Lebuhnl modified by Chris Jackson

"""
__doc__ text here
"""

from Bio import AlignIO, SeqIO
import sys
import string
import glob
import os
import argparse


def createfolder(directory):
    """
    Attempts to create a directory named after the name provided, and provides an error message on failure
    """
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print(f'Error: Creating directory: {directory}')


def parse_arguments():
    parser = argparse.ArgumentParser(description=print(__doc__), formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-gene_folder', type=str, help='path to folder with alignments', default='./')
    parser.add_argument('-matchvalue', type=int, help='how many bases need to match to stop cutting from ends',
                        default=5)
    parser.add_argument('-fastaext', type=str, help='extension of fasta files', default='.fasta')
    parser.add_argument('-targetname', action='append', type=str, dest='targetname',
                        help='beginning of bait / target / outgroup name used for comparison')
    parser.add_argument('-mingaplength', type=int,
                        help='minimum size of gap to examine match to target of sequence areas around it', default=15)
    parser.add_argument('-outputfolder', type=str, help='output folder for cut alignments', default='output_folder')

    results = parser.parse_args()
    return results


results = parse_arguments()
inputfilenames = glob.glob(f'{results.gene_folder}/*{results.fastaext}')
createfolder(results.outputfolder)

for i in range(0, len(inputfilenames)):
    # read alignment
    alignmentfilename = inputfilenames[i]
    alignment = AlignIO.read(alignmentfilename, "fasta")
    # only proceed if alignment has more than one sequence
    if len(alignment) > 1:
        target = -99  # CJJ I don't understand the significance of this, ask Alexander
        for seq in alignment:  # CJJ remove the _R_ appended to sequences by mafft
            if seq.id[:3] == "_R_":
                seq.id = seq.id[3:]
            seq.name = ''
            seq.description = ''
        alignment_sequence_dict = {}
        for seq in alignment:
            alignment_sequence_dict[seq.id.split('.')[0]] = seq
        for target in results.targetname:
            try:
                trim_to_target_seq = alignment_sequence_dict[target]
                break
            except KeyError:
                print(f'Target {target} not found in alignment {alignmentfilename}')

        if target != -99:  # CJJ I don't understand the significance of this, ask Alexander
            didicut = False
            # for i in range(0, len(alignment)):
            for sequence in alignment:
                if sequence.id != trim_to_target_seq.id:
                    # print(sequence.id)

                    # if i != target:
                    # walk through alignment from the left
                    # ensure that cuts will be made directly from beginning of alignment and after longer gaps,
                    # but not after small gaps
                    x = 0
                    longgap = True
                    # print(len(sequence.seq))
                    while x < (len(sequence.seq) - 1):  # CJJ increment along each base of the sequence
                        gaplengthcount = 0
                        while (sequence.seq[x] == "-") & (x < (len(sequence.seq) - 1)):  # CJJ record gaps in sequence
                            x = x + 1
                            # print(x)
                            gaplengthcount = gaplengthcount + 1
                        if gaplengthcount > results.mingaplength:  # change this to adjust minimum gap length for cutting at the edges
                            # print(f'CJJ gaplengthcount: {gaplengthcount}')
                            longgap = True
                        if longgap & (x < (len(sequence.seq) - 1 - results.matchvalue)):
                            # print(f'CJJ longap is {longgap}')
                            # keep turning bases into gaps until at least "matchvalue" consecutive ones were identical to target
                            shouldicut = False
                            startingmatch = x
                            # print(f'startingmatch: {startingmatch}')
                            x = x + results.matchvalue  # CJJ add the number of matchingvalue bases for the purpose of taking a slice
                            # print(x)
                            while (x < len(sequence.seq) - 1) & (
                                    sequence.seq[(x - results.matchvalue):x].upper() != trim_to_target_seq.seq[(
                                    x - results.matchvalue):x].upper()):  # CJJ take a slice of the alignment after a gap, and compare it to the target sequence slice
                                # print(f'CJJ not a match yet!')
                                x = x + 1
                                shouldicut = True
                            if shouldicut:
                                # print(f'CJJ about to parse trimmed sequence with x value {x}...')
                                sequence.seq = sequence.seq[:startingmatch] + '-' * (
                                            x - startingmatch - results.matchvalue) + sequence.seq[(
                                                                                                               x - results.matchvalue):]  # CJJ replace non-matching bases with gaps
                                dididcut = True
                            # print(f'CJJ: didicut from left end: {dididcut}')
                        # move to end of remaining sequence section
                        x = x + 1
                        if x < (len(sequence.seq) - 1):
                            while (sequence.seq[x] != "-") & (x < (len(sequence.seq) - 2)):  # CJJ why -2 here?
                                x = x + 1
                        longgap = False
                    # now the same from the right end
                    # print(f'CJJ starting from right hand side now...')
                    x = alignment.get_alignment_length() - 1
                    longgap = True
                    while x > 0:
                        gaplengthcount = 0
                        while (sequence.seq[x] == "-") & (x > 1):
                            x = x - 1
                            gaplengthcount = gaplengthcount + 1
                        if gaplengthcount > results.mingaplength:  # change this to adjust minimum gap length for cutting at the edges
                            longgap = True
                        if longgap & (x > results.matchvalue):
                            # keep turning bases into gaps until at least "matchvalue" consecutive ones were identical to target
                            shouldicut = False
                            startingmatch = x
                            x = x - results.matchvalue
                            while sequence.seq[x:(x + results.matchvalue)].upper() != trim_to_target_seq.seq[x:(
                                    x + results.matchvalue)].upper():
                                x = x - 1
                                shouldicut = True
                            if shouldicut:
                                # print(f'CJJ about to parse sequence from right hand side!')
                                sequence.seq = sequence.seq[:(x + results.matchvalue)] + '-' * (
                                            startingmatch - x - results.matchvalue + 1) + sequence.seq[
                                                                                          (startingmatch + 1):]
                                dididcut = True
                            # print(f'CJJ dididcut from right hand side {dididcut} ')
                        # move to end of remaining sequence section
                        if x > 0:
                            x = x - 1
                            while (sequence.seq[x] != "-") & (x > 1):
                                x = x - 1
                        longgap = False
                # print(f'CJJ dididcut {dididcut}')
            if dididcut:
                print("Removed bad ends in " + alignmentfilename)
        else:
            print("Target not found in " + alignmentfilename)
        # write alignment
        alignmentfilename = os.path.basename(alignmentfilename)
        AlignIO.write(alignment,
                      f'{results.outputfolder}/{alignmentfilename[:len(alignmentfilename) - len(results.fastaext)]}.internalcut{results.fastaext}',
                      "fasta")

        # Filter out empty sequences comprised only of dashes
        try:
            seqs_to_retain = []
            trimmed_file = f'{results.outputfolder}/{alignmentfilename[:len(alignmentfilename) - len(results.fastaext)]}.internalcut{results.fastaext}'
            with open(trimmed_file, 'r') as trimmed_fasta:
                seqs = SeqIO.parse(trimmed_fasta, 'fasta')
                for seq in seqs:
                    characters = set(character for character in seq.seq)
                    if len(characters) == 1 and '-' in characters:
                        pass
                    else:
                        seqs_to_retain.append(seq)
            with open(trimmed_file, 'w') as filtered_hmm_fasta:
                SeqIO.write(seqs_to_retain, filtered_hmm_fasta, 'fasta')
        except:
            pass
