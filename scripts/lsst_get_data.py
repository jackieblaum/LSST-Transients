#!/usr/bin/env python
import argparse

from lsst_transients.data_database import DataDatabase
from lsst_transients.grid import Grid
from lsst_transients.utils.file_io import sanitize_filename

if __name__ == "__main__":

    # Collect input from the user
    parser = argparse.ArgumentParser(description="Collect information from each visit")
    parser.add_argument('--visits_dir', type=str, help='Input path to folder with all visits', required=True)
    parser.add_argument('--grid', type=str, help='File containing the grid', required=True)
    parser.add_argument('--filter', type=str, help='Name of the filter to analyze (r, u, v...)', required=True)
    parser.add_argument("--db", type=str, help='File name for the output database', required=True)

    parser.add_argument("--ncpus", type=int, help='Number of cpus to use', required=False, default=1)

    args = parser.parse_args()

    # First get the grid
    grid = Grid.read_from(sanitize_filename(args.grid))

    with DataDatabase(sanitize_filename(args.db)) as db:

        db.fill_visits(sanitize_filename(args.visits_dir), args.filter, grid, args.ncpus)
