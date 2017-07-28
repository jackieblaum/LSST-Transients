#!/usr/bin/env python
import argparse

from lsst_transients.create_db import Data_Database

if __name__ == "__main__":

    # Collect input from the user
    parser = argparse.ArgumentParser(description="Collect information from each visit")
    parser.add_argument('-p', '--path', type=str, help='Input path to folder with all visit files', required=True)
    parser.add_argument('-n', '--name', type=str, help='Name for the database', required=True)
    parser.add_argument('-f', '--flux', type=bool, help='Get the fluxes from the files (True or False)',
                        required=False, default=True)
    parser.add_argument('-c', '--conditions', type=bool, help='Get the conditions from the files (True or False)',
                        required=False, default=True)
    
    args = parser.parse_args()
        
    db = Data_Database("%s.db" % args.name, first=False)
    db.fill_visits(args.path, args.flux, args.conditions)
    db.close()