#! /usr/bin/env python
'''
############################################
# Push entries of ntuple into an hdf5 file #
############################################
'''
import argparse
import os
import sys
from array import array
import numpy
import h5py
from datetime import datetime
import getpass
try:
    import ROOT
    ROOT.gROOT.SetBatch(True) # ROOT in batch mode
except ImportError:
    print "Cant find root lib is it in $PYTHONPATH?"
    sys.exit()

def process(file_name, tree_name, branches, out_file, n_events):
    '''
    read ntp make hdf5
    '''
    # Get the event tree
    tree = ROOT.TChain(tree_name)
    tree.Add(file_name)
    #tree.Print()
    if not tree:
        print "Error: No tree named %s in %s" %(tree_name, file_name)
        sys.exit()

    # Check branches exist
    branches_avail = [x.GetName() for x in tree.GetListOfBranches()]
    print "branches_available:", branches_avail
    branches = [branch.strip() for branch in branches]
    print "branches specified:", branches
    for branch in branches:
        if not branch in branches_avail:
            print "Error branch '%s' not a branch in input tree" %(branch)
            print "Branches available are: \n"
            print "\t".join(branches_avail)
            sys.exit()

    # output
    if n_events < 0:
        n_events = tree.GetEntries()

    # loop over events and fill the branches of new ntuple
    vals = [] #empty array for the case the tree has zero entries
    for index, entry in enumerate(tree):
        vals.append([entry.__getattr__(branch) for branch in branches])
        if index % 10000 == 0:
            print index, "/", n_events
    vals = numpy.array(vals)

    # vals = numpy.array([]) #empty array for the case the tree has zero entries
    # for index, entry in enumerate(tree):
        # if index > n_events:
            # break
        # branch_vals = array('f', [entry.__getattr__(branch) for branch in branches])
        # if index == 0:
            # vals = numpy.array([branch_vals])
        # else:
            # vals = numpy.concatenate((vals, [branch_vals]), axis=0)

    # Save
    h5_file_out = h5py.File(out_file, 'w')
    # meta information
    h5_file_out.attrs[u'file_name'] = os.path.basename(out_file)
    h5_file_out.attrs[u'file_time'] = str(datetime.now())
    h5_file_out.attrs[u'author'] = getpass.getuser()
    h5_file_out.attrs[u'HDF5_Version'] = str(h5py.version.hdf5_version)
    h5_file_out.attrs[u'h5py_version'] = str(h5py.version.version)

    # within the hdf5 file, save data into a 'dataset' group
    dataset = h5_file_out.create_group(u'dataset')
    dataset.attrs[u'branches'] = branches # meta info
    # save each branch as dataset (using branch names)
    for i, branch in enumerate(branches):
      dataset.create_dataset(str(branch), data=vals[:, i])

    h5_file_out.close()

    print "Written %i entries of branch(es) '%s' \nin file %s"\
        %(n_events, ":".join(branches), out_file)

def main(args):
    '''
    load args
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str, )
    parser.add_argument('-treename', metavar='-t', type=str, default="output")
    parser.add_argument('branches', type=str)
    parser.add_argument('-outfile', metavar='-o', type=str, default="")
    parser.add_argument('-nevents', metavar='-nev', type=int, default=-1)
    args = parser.parse_args(args)

    if args.outfile == "":
        outfile = os.path.splitext(args.filename)[0]
        outfile += ".h5"
        args.outfile = outfile
    process(args.filename, args.treename, args.branches, args.outfile, args.nevents)

if __name__ == "__main__":
    main(sys.argv[1:])
