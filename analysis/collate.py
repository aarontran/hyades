#!/usr/bin/env python
"""
collate.py

Script to collate multiple HDF5 files, with datasets d_0, d_1, d_2, etc
of ranks N_0, N_1, N_2, etc into a new file with same datasets, and ranks +=1
for all. The new dimension is the front axis of the collated datasets.

Requires HDF5 and h5py.
"""

from __future__ import division, print_function

import argparse
from datetime import datetime
from glob import glob
import numpy as np
import os.path as path

import h5py


def dump_hdf5(fname, **kwds):
    """
    Dump numpy array to HDF5 file.
    Interface mimics numpy.savez(...)
    Input
    Inputs:
        fname = filename path to write
        kwds = key-value dict, key = dataset name, val = array or data to save
    """
    with h5py.File(fname, "w") as f:
        for k, v in kwds.items():
            f.create_dataset(k, data=v)
    return


def main():
    """User-facing interface for command-line usage"""

    started = datetime.now()

    description = "Collate multiple HDF5 files into one"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('dir', type=str,
                        help='path to HDF5 file directory')
    parser.add_argument('stem', type=str,
                        help='file stem, expect form stem.xxx.hdf5')
    parser.add_argument('--dry-run', action='store_true',
                        help='print file paths to write but do no computation')
    args = parser.parse_args()

    # List all the file paths
    dumps = glob(path.join(args.dir, f"{args.stem:s}.*.hdf5"))

    # Which time steps do we need to process?
    steps = []
    for d in dumps:
        numstr = d.split('.')[-2]
        steps.append(int(numstr))
    steps = np.array(sorted(steps))

    # hardcoded dump into cwd
    fname_out = f"{args.stem:s}.all.hdf5"

    collated = dict()

    for ii, step in enumerate(steps):

        fname = path.join(args.dir, f"{args.stem:s}.{step:d}.hdf5")

        if args.dry_run:
            print(f"Would process {fname:s}")
            continue

        with h5py.File(fname, "r") as f:
            for k, v in f.items():
                if k not in collated:
                    assert ii == 0
                    collated[k] = np.empty((len(steps),*v.shape), dtype=v.dtype)
                collated[k][ii,:] = v[()]

        print(f"Done step {step:d} elapsed", datetime.now()-started)

    if args.dry_run:
        print(f"Would collate into {fname_out:s}")
    else:
        dump_hdf5(fname_out, **collated)

    print("collate.py all done")


if __name__ == '__main__':
    main()
