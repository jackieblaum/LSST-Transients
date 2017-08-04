#!/usr/bin/env python
import argparse

from lsst_transients.examine_transient_candidates import examine_transient_candidates

if __name__ == "__main__":

    # Collect input from the user
    parser = argparse.ArgumentParser(description="Collect information from each visit")
    parser.add_argument('-n', '--name', type=str, help='Name of the _database (without .db)', required=True)
    parser.add_argument('-f', '--file', type=str, help='Name of the yaml file in which the blocks are stored '
                                                        '(without .yml)', required=True)

    args = parser.parse_args()

    examine_transient_candidates(args.name, '%s.yml' % args.file)