#!/usr/bin/env python
import argparse

from lsst_transients.run_bayesian_blocks import BayesianBlocks

if __name__ == "__main__":

    # Collect input from the user
    parser = argparse.ArgumentParser(description="Collect information from each visit")
    parser.add_argument('-s', '--sort', type=bool, help='Sort the visits in the database (True or False)',
                        required=False, default=False)
    parser.add_argument('-n', '--name', type=str, help='Name of the database (without .db)', required=True)

    args = parser.parse_args()

    bb = BayesianBlocks(args.name)
    bb.run_algorithm(args.sort)


