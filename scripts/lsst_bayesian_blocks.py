#!/usr/bin/env python
import argparse

from lsst_transients.run_bayesian_blocks import BayesianBlocks
from lsst_transients.utils.file_io import sanitize_filename

if __name__ == "__main__":

    # Collect input from the user
    parser = argparse.ArgumentParser(description="Collect information from each visit")

    parser.add_argument('--db', type=str, help='Name of the _database', required=True)
    parser.add_argument('--outfile', type=str, help='Name of the file to store the blocks in', required=True)
    parser.add_argument('--min_edges', type=int, help='Minimum number of block edges for a region to be saved',
                        required=False, default=4)
    parser.add_argument('--p0', type=float, help='Null-hyp probability',
                        required=False, default=1e-5)

    args = parser.parse_args()

    bb = BayesianBlocks(sanitize_filename(args.db), sanitize_filename(args.outfile), args.min_edges)
    bb.run_algorithm(args.p0)


