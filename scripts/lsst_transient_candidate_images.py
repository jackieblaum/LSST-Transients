#!/usr/bin/env python
import argparse

from lsst_transients.examine_transient_candidates import examine_transient_candidates

if __name__ == "__main__":

    # Collect input from the user
    parser = argparse.ArgumentParser(description="Collect information from each visit")
    parser.add_argument('-n', '--name', type=str, help='Name of the database (without .db)', required=True)
    parser.add_argument('-f', '--file', type=str, help='Name of the yaml file in which the blocks are stored '
                                                        '(without .yml)', required=True)
    parser.add_argument('-m', '--multiply', type=int, help='Size of the images will be this many times the '
                                                           'size of the regions', required=False, default=10)
    parser.add_argument('-v', '--visits', type=str, help='Name of the directory storing the visit files', required=True)
    parser.add_argument('-b', '--filter', type=str, help='Name of the filter to analyze (r, u, v...)', required=True)

    args = parser.parse_args()

    examine_transient_candidates(args.name, '%s.yml' % args.file, args.multiply, args.visits, args.filter)