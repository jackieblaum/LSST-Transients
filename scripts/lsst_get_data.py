#!/usr/bin/env python
import argparse

from lsst_transients.create_db_notpyregion import Data_Database

if __name__ == "__main__":

    # Collect input from the user
    parser = argparse.ArgumentParser(description="Collect information from each visit")
    parser.add_argument('-p', '--path', type=str, help='Input path to folder with all visit files', required=True)
    
    args = parser.parse_args()
        
    db = Data_Database("lsst_transients.db")
    db.fill_visits(args.regfile, args.path)
    db.close()
