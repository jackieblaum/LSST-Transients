#!/usr/bin/env python
import argparse

from lsst_transients.examine_transient_candidates import examine_transient_candidates
from lsst_transients.utils.file_io import sanitize_filename
from lsst_transients.grid import Grid

if __name__ == "__main__":

    # Collect input from the user
    parser = argparse.ArgumentParser(description="Collect information from each visit")
    parser.add_argument('--db', type=str, help='Name of the _database (without .db)', required=True)
    parser.add_argument('--grid', type=str, help='Database with the grid', required=True)
    parser.add_argument('-f', '--results', type=str, help='Name of the json file in which the blocks are stored',
                        required=True)
    parser.add_argument('-m', '--multiply', type=float, help='Size of the images will be this many times the '
                                                           'size of the regions', required=False, default=2)
    parser.add_argument('-v', '--visits', type=str, help='Name of the directory storing the visit files', required=True)
    parser.add_argument('-b', '--filter', type=str, help='Name of the filter to analyze (r, u, v...)', required=True)
    parser.add_argument('-r', '--regid', type=str, help='Region ID of the region to be examined (ie. reg1)', required=True)

    args = parser.parse_args()

    db = sanitize_filename(args.db)
    results = sanitize_filename(args.results)
    visits = sanitize_filename(args.visits)

    grid = Grid.read_from(sanitize_filename(args.grid))

    examine_transient_candidates(db, grid, args.regid, results, args.multiply, visits, args.filter)