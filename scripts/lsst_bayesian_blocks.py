#!/usr/bin/env python
import argparse

from lsst_transients.run_bayesian_blocks import BayesianBlocks

if __name__ == "__main__":

    # Collect input from the user
    parser = argparse.ArgumentParser(description="Collect information from each visit")
    parser.add_argument('-s', '--sort', type=bool, help='Sort the visits in the _database (True or False)',
                        required=False, default=False)
    parser.add_argument('-n', '--name', type=str, help='Name of the _database (without .db)', required=True)
    parser.add_argument('-f', '--file', type=str, help='Name of the file to store the blocks in', required=True)
    parser.add_argument('-b', '--blocks', type=int, help='Minimum number of blocks for a region to be saved',
                        required=False, default=4)

    args = parser.parse_args()

    bb = BayesianBlocks(args.name, args.file, args.blocks)
    bb.run_algorithm(args.sort)


