from grid import Grid

import astropy.io.fits as pyfits
import argparse

if __name__ == "__main__":

    # Collect input from the user
    parser = argparse.ArgumentParser(description="Grid Region Generator")
    parser.add_argument('-i', '--input', type=str, help='Input filename', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output filename', required=True)
    parser.add_argument('-d', '--diameter', type=float, help='Length of the sides or diameter of the regions in pixels', required=True)
    parser.add_argument('-f', '--fraction', type=float, help='Fraction of the region for which they overlap',
                        required=True)
    parser.add_argument('-s', '--shape', type=str, help='Either circle or square regions', required=True)

    args = parser.parse_args()

    # Open the headers from the fits file inputted by the user
    header = pyfits.getheader(args.infile, 0)
    header2 = pyfits.getheader(args.infile, 1)

    # Set variables by reading from the header
    center_pix_x = header2['CRPIX1']
    center_pix_y = header2['CRPIX2']
    max_pix_x = header2['NAXIS1']
    max_pix_y = header2['NAXIS2']
	
    if header['ROTANG']:
        rotation_angle = header['ROTANG']

    # Generate the grid
    grid = Grid(center_pix_x, center_pix_y, 0, max_pix_x, 0, max_pix_y, args.fraction, args.diameter, args.shape)
    angular_distance_arcsec = grid.write_grid(args.infile, args.output, args.shape)

    print("Input file: %s" % args.input)
    print("Output file: %s" % args.output)
    print("Side length of regions in pixels/arcsec: %f/ %f" % (args.diameter, angular_distance_arcsec))
    print("Overlap (fraction of region): %f" % args.fraction)
