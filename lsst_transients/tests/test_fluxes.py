import argparse
import numpy as np

from create_db import Database
from lsst_transients.data_database import DataDatabase

if __name__ == "__main__":
    # Collect input from the user
    parser = argparse.ArgumentParser(description="Collect information from each visit")
    parser.add_argument('-p', '--path', type=str, help='Input path to folder with all visit files', required=True)
    parser.add_argument('-b', '--filter', type=str, help='Name of the filter to analyze (r, u, v...)', required=True)
    parser.add_argument('-f', '--flux', type=bool, help='Get the fluxes from the files (True or False)',
                        required=False, default=True)
    parser.add_argument('-c', '--conditions', type=bool, help='Get the conditions from the files (True or False)',
                        required=False, default=False)

    args = parser.parse_args()

    db = Database("test1.db")
    db.fill_visits(args.path, args.filter, args.flux, args.conditions)

    db2 = DataDatabase("test2.db", first=False)
    db2.fill_visits(args.path, args.filter, args.flux, args.conditions)
    fluxes1 = db.get_flux_list()
    fluxes2 = db2.get_flux_list()
    assert np.allclose(fluxes1, fluxes2, rtol=0.15), "New fluxes do not match pyregion fluxes:" \
                                                                "\nPyregion: %s\n " \
                                                                "New: %s" \
                                                                % (str(fluxes1), str(fluxes2))
    db.close()