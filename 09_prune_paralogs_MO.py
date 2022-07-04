# Author: Yang and Smith, modified by Alexander Schmidt-Lebuhn

"""
Taxon duplication? --No--> output one-to-one orthologs
		|
	   Yes
	    |
Outgroup present? --No--> ignore this homolog
		|
	   Yes
	    |
Outgroup taxon duplication? --Yes--> ignore this homolog
		|
	    No
	    |
Outgroup monophyletic? --No--> ignore this homolog
		|
	   Yes
	    |
Infer orthologs by using monophyletic, non-repeating outgroups

If not to output 1-to-1 orthologs, for example, already analysed these
set OUTPUT_1to1_ORTHOLOGS to False
"""

import phylo3
import newick3
import os
import sys
import argparse
import glob
import shutil


# if pattern changes, change it here
# given tip label, return taxon name identifier
def get_name(label):
    return label.split(".")[0]


def get_clusterID(filename):
    return filename.split(".")[0]


def get_front_labels(node):
    leaves = node.leaves()
    return [i.label for i in leaves]


def get_back_labels(node, root):
    all_labels = get_front_labels(root)
    front_labels = get_front_labels(node)
    return set(all_labels) - set(front_labels)  # labels do not repeat


def get_front_names(node):  # may include duplicates
    labels = get_front_labels(node)
    return [get_name(i) for i in labels]


def get_front_outgroup_names(node):
    names = get_front_names(node)
    return [i for i in names if i in OUTGROUPS]


def get_back_names(node, root):  # may include duplicates
    back_labels = get_back_labels(node, root)
    return [get_name(i) for i in back_labels]


def remove_kink(node, curroot):
    if node == curroot and curroot.nchildren == 2:
        # move the root away to an adjacent none-tip
        if curroot.children[0].istip:  # the other child is not tip
            curroot = phylo3.reroot(curroot, curroot.children[1])
        else:
            curroot = phylo3.reroot(curroot, curroot.children[0])
    # ---node---< all nodes should have one child only now
    length = node.length + (node.children[0]).length
    par = node.parent
    kink = node
    node = node.children[0]
    # parent--kink---node<
    par.remove_child(kink)
    par.add_child(node)
    node.length = length
    return node, curroot


# check if outgroups are monophyletic and non-repeating and reroot
# otherwise return None
def reroot_with_monophyletic_outgroups(root):
    lvs = root.leaves()
    outgroup_matches = {}  # key is label, value is the tip node object
    # Since no taxon repeat in outgroups name and leaf is one-to-one
    outgroup_labels = []
    for leaf in lvs:
        label = leaf.label  # CJJ e.g. 376678.main or 376728.0, etc
        name = get_name(label)  # CJJ e.g. 376678 or 376728, etc
        if name in OUTGROUPS:
            outgroup_matches[label] = leaf
            outgroup_labels.append(leaf.label)
    if len(outgroup_labels) == 1:  # one single outgroup
        # cannot reroot on a tip so have to go one more node into the ingroup
        new_root = outgroup_matches[outgroup_labels[0]].parent
        return phylo3.reroot(root, new_root)
    else:  # has multiple outgroups. Check monophyly and reroot
        newroot = None
        for node in root.iternodes():  # CJJ iterate over nodes and try to find one with monophyletic outgroup
            if node == root:
                continue  # skip the root
            front_names = get_front_names(node)
            back_names = get_back_names(node, root)
            front_in_names, front_out_names, back_in_names, back_out_names = 0, 0, 0, 0
            for i in front_names:
                if i in OUTGROUPS:
                    front_out_names += 1
                else:
                    front_in_names += 1
            for j in back_names:
                if j in OUTGROUPS:
                    back_out_names += 1
                else:
                    back_in_names += 1
            if front_in_names == 0 and front_out_names > 0 and back_in_names > 0 and back_out_names == 0:
                newroot = node  # ingroup at back, outgroup in front
                break
            if front_in_names > 0 and front_out_names == 0 and back_in_names == 0 and back_out_names > 0:
                newroot = node.parent  # ingroup in front, outgroup at back
                break
        if newroot != None:
            return phylo3.reroot(root, newroot)
        else:
            return None


def prune_paralogs_from_rerooted_homotree(root):
    if len(get_front_names(root)) == len(set(get_front_names(root))):
        return root  # no pruning needed
    # check for duplications at the root first
    # one or two of the trifurcating root clades are ingroup clades
    node0, node1, node2 = root.children[0], root.children[1], root.children[2]
    out0, out1, out2 = len(get_front_outgroup_names(node0)), len(get_front_outgroup_names(node1)), len(
        get_front_outgroup_names(node2))
    if out0 == 0 and out1 == 0:  # 0 and 1 are the ingroup clades
        name_set0 = set(get_front_names(node0))
        name_set1 = set(get_front_names(node1))
        if len(name_set0.intersection(name_set1)) > 0:
            if len(name_set0) > len(name_set1):  # cut the side with less taxa
                root.remove_child(node1)
                node1.prune()
            else:
                root.remove_child(node0)
                node0.prune()
    elif out1 == 0 and out2 == 0:  # 1 and 2 are the ingroup clades
        name_set1 = set(get_front_names(node1))
        name_set2 = set(get_front_names(node2))
        if len(name_set1.intersection(name_set2)) > 0:
            if len(name_set1) > len(name_set2):  # cut the side with less taxa
                root.remove_child(node2)
                node2.prune()
            else:
                root.remove_child(node1)
                node1.prune()
    elif out0 == 0 and out2 == 0:  # 0 and 2 are the ingroup clades
        name_set0 = set(get_front_names(node0))
        name_set2 = set(get_front_names(node2))
        if len(name_set0.intersection(name_set2)) > 0:
            if len(name_set0) > len(name_set2):  # cut the side with less taxa
                root.remove_child(node2)
                node2.prune()
            else:
                root.remove_child(node0)
                node0.prune()

    while len(get_front_names(root)) > len(set(get_front_names(root))):
        for node in root.iternodes(order=0):  # PREORDER, root to tip  CJJ: this tree includes outgroup taxa

            if node.istip:
                continue
            elif node == root:
                continue

            child0, child1 = node.children[0], node.children[1]
            name_set0 = set(get_front_names(child0))
            name_set1 = set(get_front_names(child1))
            if len(name_set0.intersection(name_set1)) > 0:
                if len(name_set0) > len(name_set1):  # cut the side with less taxa
                    node.remove_child(child1)
                    child1.prune()
                else:
                    node.remove_child(child0)
                    child0.prune()
                node, root = remove_kink(node, root)  # no re-rooting here
                break
    return root


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('treefile_folder', type=str, help='A folder containing tree files')
    parser.add_argument('treefile_suffix', type=str, default='.treefile',
                        help='suffix for the tree files, default is ".treefile".')
    parser.add_argument('minimum_taxa', type=int, help='Minimum number of taxa required')
    parser.add_argument('output_directory', type=str, help='Name for output directory')
    parser.add_argument('in_and_outgroup_list', type=str, help='Text file listing in and out-group taxa')

    args = parser.parse_args()

    OUTPUT_1to1_ORTHOLOGS = True
    INGROUPS = []
    OUTGROUPS = []

    # taxon_code_file_file = sys.argv[5]
    with open(args.in_and_outgroup_list, "r") as infile:
        for line in infile:
            if len(line) < 3:
                continue
            sample = line.strip().split("\t")
            if sample[0] == "IN":
                INGROUPS.append(sample[1])
            elif sample[0] == "OUT":
                OUTGROUPS.append(sample[1])
            else:
                print("Check taxon_code_file file format")
                sys.exit()

    if len(set(INGROUPS) & set(OUTGROUPS)) > 0:
        print("Taxon ID", set(INGROUPS) & set(OUTGROUPS), "in both ingroup and outgroup")
        sys.exit(0)

    print(len(INGROUPS), "ingroup taxa and", len(OUTGROUPS), "outgroup taxa read")

    for treefile in glob.glob(f'{args.treefile_folder}/*{args.treefile_suffix}'):
        treefile_basename = os.path.basename(treefile)

        # read in the tree and check number of taxa
        outID = f'{args.output_directory}/{get_clusterID(treefile_basename)}'

        # with open(inDIR + treefile, "r") as infile:
        with open(treefile, "r") as infile:
            intree = newick3.parse(infile.readline())
            curroot = intree
            names = get_front_names(curroot)
            num_tips, num_taxa = len(names), len(set(names))
            if num_taxa < args.minimum_taxa:
                continue  # not enough taxa

        # If the homolog has no taxon duplication, no cutting is needed
        if num_tips == num_taxa:
            if OUTPUT_1to1_ORTHOLOGS:
                shutil.copy(treefile, f'{outID}.1to1ortho.tre')
        else:
            # now need to deal with taxon duplications
            # check to make sure that the ingroup and outgroup names were set correctly
            unrecognised_names = []
            for name in names:
                if name not in INGROUPS and name not in OUTGROUPS:
                    print("check name", name)
                    unrecognised_names.append(name)
                    with open('prune_MO_trees_skipped_incorrect_names.txt', 'a+') as mafft_reversed_log:
                        mafft_reversed_log.write(f'Check outgroup sequence for tree {treefile}, sequence {name}\n')
            if len(unrecognised_names) != 0:
                print(f'unrecognised_names: {unrecognised_names}')
                continue  # CJJ should just skip these trees, not stop completely

            outgroup_names = get_front_outgroup_names(curroot)

            # if no outgroup at all, do not attempt to resolve gene duplication
            if len(outgroup_names) == 0:
                print(f"duplicated taxa in unrooted tree {treefile}")
                # CJJ I don't understand why the message above is printed - shouldn't it say 'no outgroup taxa for tree?

            # skip the homolog if there are duplicated outgroup taxa
            elif len(outgroup_names) > len(set(outgroup_names)):
                print(f"Outgroup contains taxon repeats, skipping tree {treefile_basename}")

            else:  # at least one outgroup present and there's no outgroup duplication
                if curroot.nchildren == 2:  # need to reroot
                    temp, curroot = remove_kink(curroot, curroot)
                curroot = reroot_with_monophyletic_outgroups(curroot)

                # only return one tree after pruning
                if curroot != None:  # CJJ i.e. the outgroup was monophyletic
                    with open(f'{outID}.reroot', "w") as outfile:
                        outfile.write(newick3.tostring(curroot) + ";\n")  # CJJ write re-rooted trees to file
                    ortho = prune_paralogs_from_rerooted_homotree(curroot)
                    if len(set(get_front_names(curroot))) >= args.minimum_taxa:
                        with open(f'{outID}.ortho.tre', "w") as outfile:
                            outfile.write(newick3.tostring(ortho) + ";\n")
                    else:
                        print(f"Not enough taxa after pruning for tree {treefile_basename}")
                else:
                    print(f"Outgroup non-monophyletic for tree {treefile_basename}")

