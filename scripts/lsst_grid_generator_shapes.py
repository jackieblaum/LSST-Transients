#!/usr/bin/env python
import argparse
import os

import astropy.io.fits as pyfits

from lsst_transients.grid import Grid
from lsst_transients.data_database import DataDatabase
from lsst_transients.utils.logging_system import get_logger

log = get_logger(os.path.basename(__file__))


if __name__ == "__main__":

    # Collect input from the user
    parser = argparse.ArgumentParser(description="Grid CircularRegion Generator")
    parser.add_argument('-i', '--input_example', type=str, help='FITS file for one of the visits', required=True)
    parser.add_argument('-d', '--diameter', type=float, help='Length of the sides or diameter of the regions in pixels',
                        required=True)
    parser.add_argument('-f', '--fraction', type=float, help='Fraction of the region for which they _overlap',
                        required=True)
    parser.add_argument('-s', '--shape', type=str, help='Either circle or square regions', required=True)
    parser.add_argument('-n', '--name', type=str, help='Name of the database', required=True)

    args = parser.parse_args()

    # Make sure the database file does not exist
    if os.path.exists("%s.db" % args.name):

        raise IOError("Database %s already exists" % args.name)

    # Open the headers from the fits file inputted by the user
    with pyfits.open(args.input_example) as f:

        # Primary header
        header = f[0].header

        # Header of the background-subtracted image
        header2 = f[1].header

    # Get the pixels numbers and the center pixel
    min_pix_x = 1.0  # Fortran convention as in ds9
    max_pix_x = header2['NAXIS1']

    min_pix_y = 1.0  # Fortran convention as in ds9
    max_pix_y = header2['NAXIS2']

    # Some images might be rotated
    if header['ROTANG']:

        rotation_angle = header['ROTANG']

    else:

        rotation_angle = 0

    # Generate the grid of regions
    grid = Grid(min_pix_x, max_pix_x, min_pix_y, max_pix_y, args.fraction, args.diameter, args.shape)
    region_centers_wcs, angular_distance_arcsec = grid.get_grid(args.input_example)

    log.info("Generated %i regions" % len(region_centers_wcs))

    # Write the regions to a database
    db = DataDatabase("%s.db" % args.name)
    num_regs = db.fill_reg(region_centers_wcs, args.shape, angular_distance_arcsec, rotation_angle)
    db.close()

    log.info("Input file: %s" % args.input_example)
    log.info("Database file: %s.db" % args.name)
    log.info("Side length of regions in pixels/arcsec: %.2f/ %.2f" % (args.diameter, angular_distance_arcsec))
    log.info("Overlap (fraction of region): %f" % args.fraction)
