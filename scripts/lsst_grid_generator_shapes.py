#!/usr/bin/env python
import argparse
import os

import astropy.io.fits as pyfits

from lsst_transients.grid import GridFactory
from lsst_transients.utils.file_io import sanitize_filename
from lsst_transients.utils.logging_system import get_logger

log = get_logger(os.path.basename(__file__))


if __name__ == "__main__":

    # Collect input from the user
    parser = argparse.ArgumentParser(description="Grid CircularRegion Generator")
    parser.add_argument('--input_example', type=str, help='FITS file for one of the visits', required=True)
    parser.add_argument('--diameter', type=float, help='Length of the sides or diameter of the regions in pixels',
                        required=True)
    parser.add_argument('--fraction', type=float, help='Fraction of the region for which they _overlap',
                        required=True)
    parser.add_argument('--outfile', type=str, help='Name of the output grid file', required=True)

    args = parser.parse_args()

    grid_file = sanitize_filename(args.outfile)

    # Make sure the _database file does not exist
    if os.path.exists(grid_file):

        raise IOError("Database %s already exists" % grid_file)

    # Open the headers from the fits file inputted by the user
    input_example = sanitize_filename(args.input_example)
    with pyfits.open(input_example) as f:

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
    # if header['ROTANG']:
    #
    #     rotation_angle = header['ROTANG']
    #
    # else:
    #
    #     rotation_angle = 0

    # Generate the grid of regions
    grid_factory = GridFactory(min_pix_x, max_pix_x, min_pix_y, max_pix_y, args.fraction, args.diameter)
    grid = grid_factory.get_grid(input_example)

    log.info("Generated %i regions" % grid.n_regions)
    log.info("Side length of regions in pixels/arcsec: %.2f/ %.2f" % (args.diameter, grid["reg0"]['radius']))
    log.info("Overlap (fraction of region): %.2f" % args.fraction)

    grid.write_to(grid_file)

    log.info("Grid written into %s" % grid_file)