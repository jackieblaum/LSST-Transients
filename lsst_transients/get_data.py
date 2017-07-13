from grid import Grid
from create_db import Database

import os
import argparse
import astropy.io.fits as pyfits

if __name__ == "__main__":

    # Collect input from the user
    parser = argparse.ArgumentParser(description="Collect information from each visit")
    parser.add_argument('-p', '--path', type=str, help='Input path to folder with all visit files', required=True)
    parser.add_argument('-r', '--regfile', type=str, help='Name of the region file, include .reg', required=True)
    
    args = parser.parse_args()
        
    db = Database("lsst_transients")
    db.fill_visits(args.regfile, args.path)
    db.close()
