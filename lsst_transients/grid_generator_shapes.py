from region import Region
from utils.cartesian_product import cartesian_product
from utils.pix2world import pix2world

import numpy as np
import astropy.io.fits as pyfits
import argparse
from astropy.coordinates.angle_utilities import angular_separation
import astropy.units as u

class Grid(object):
    '''
    The constructor takes the _x and _y coordinates of the center of the image (xcenter, ycenter)
    and the maximum and minimum _x and _y coordinates of the image (xcenter, ycenter).

    It also sets the overlap factor (overlapfactor) and diameter/side length (r) of the regions in the grid, and then calculates the overlap by multiplying these two values.

    An empty array of Region objects is also instantiated.
    '''

    def __init__(self, xcenter, ycenter, xmin, xmax, ymin, ymax, overlapfactor, d, shape):
        '''
        Constructs a grid object. A region list is instantiated, and the overlap length is set to be the product of the overlap factor and diameter/r length.

        :param xcenter: The _x-coordinate of the center of the image
        :param ycenter: The _y-coordinate of the center of the image
        :param xymax: The maximum _x- and _y-coordinates of the image
        :param xymin: The minimum _x- and _y-coordinates of the image
        :param overlapfactor: The fraction of the length of the region by which the regions overlap one another on the grid
        :param r: The length of the sides or diameter of the regions in the grid
	:param shape: The shape of the regions; either circle or square
        '''

        self.xcenter = xcenter
        self.ycenter = ycenter

        self._xmin = xmin
        self._xmax = xmax
        self._ymin = ymin
        self._ymax = ymax

        self.overlapfactor = overlapfactor
        self.d = d
        self.overlap = overlapfactor * d

        self.shape = shape

        self.debug = False

    def _in_range(self, a_region):
        '''
        A helper method used by do_steps() in order to determine if the Region is in range of the region of
        interest and should be added to the Region array.

        Returns true if out of range, returns false if in range.

        *Not yet fully implemented to account for rotation angle.

        :param a_region: the region to be checked
        :return: False if the region is out of the region of interest, True otherwise
        '''

        x = a_region.x
        y = a_region.y

        if self._xmin < x < self._xmax and self._ymin < y < self._ymax:

            return True

        else:

            return False

    def _do_round(self, delta):

        # Let's figure out the iteration we are in
        iteration = delta / (self.d - self.overlap)

        # Let's compute the start point of the round
        start_x = self.xcenter + delta
        start_y = self.ycenter + delta

        # How many steps do we need to complete one side of the round?
        steps = iteration * 2

        xs = np.linspace(start_x - steps * (self.d - self.overlap), start_x, steps + 1)
        ys = np.linspace(start_y - steps * (self.d - self.overlap), start_y, steps + 1)

        this_regions = map(lambda (x, y): Region(x, y, self.d, self.shape), cartesian_product([xs, ys]))

        return this_regions

    @staticmethod
    def distance(x1, y1, x2, y2):

        return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

    def get_grid(self):
        '''
        Calls gen_grid() using the location of the center of the region of interest and
        the starting step value 0.

        Returns the array of Region objects.

        :return:
        '''

        # Figure out maximum diagonal
        d1 = self.distance(self._xmin, self._ymin, self.xcenter, self.ycenter)
        d2 = self.distance(self._xmin, self._ymax, self.xcenter, self.ycenter)

        longest_diagonal = max(d1, d2)

        all_regions = self._do_round(longest_diagonal)

        survived_regions = filter(self._in_range, all_regions)

        return survived_regions

    def pix_to_wcs(self, filename, regions):
        '''
        Converts the center coordinates for each region in the array from pixels to WCS and
        stores them in a new array.

        :param filename: The name of the file from which the HDUlist will be loaded
        :param regions: The array of regions with locations in pixel coordinates
        :return: The array of regions with locations in WCS
        '''

        arr = []
        for region in regions:
            arr.append([region.x, region.y])

        regions_wcs = pix2world(filename, arr)

        return regions_wcs


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

    # Open the header from the fits file inputted by the user
    header = pyfits.getheader(args.input, 0)

    # Set variables by reading from the header
    center_pix_x = header['CRPIX1']
    center_pix_y = header['CRPIX2']
    max_pix_x = header['NAXIS1']
    max_pix_y = header['NAXIS2']
    rotation_angle = header['ROTANG']

    # Generate the grid
    grid = Grid(center_pix_x, center_pix_y, 0, max_pix_x, 0, max_pix_y, args.fraction, args.diameter, args.shape)
    regions = grid.get_grid()

    # Convert the locations of the regions to WCS
    regions_wcs = grid.pix_to_wcs(args.input, regions)

    # Convert the side length/diameter of the regions to WCS
    # NOTE: we assume that on this angular scale there is no distortion due to the projection
    # TODO: fix the projection for the entire LSST field of view (3.5 deg)

    # Get the angular distance between the upper left corner and the upper right corner of the
    # center region
    center_region = regions[0]
    upper_left_x = center_region.x - center_region.d / 2.0
    upper_y = center_region.y - center_region.d / 2.0
    upper_right_x = center_region.x + center_region.d / 2.0

    upper_left_sky_coords, upper_right_sky_coords = pix2world(args.input,
                                                              [[upper_left_x, upper_y], [upper_right_x, upper_y]])

    angular_distance = angular_separation(upper_left_sky_coords[0] * u.deg, upper_left_sky_coords[1] * u.deg,
                                          upper_right_sky_coords[0] * u.deg, upper_right_sky_coords[1] * u.deg)

    angular_distance_arcsec = angular_distance.to("arcsec").value

    # Write the grid to a DS9 region file
    with open(args.output, "w+") as f:
        f.write("icrs\n")

        for (ra, dec) in regions_wcs:
            if args.shape == "square":

                f.write('box(%f, %f, %f", %f", %f)\n' % (ra, dec, angular_distance_arcsec, angular_distance_arcsec,
                                                     360 - rotation_angle))
            elif args.shape == "circle":

                f.write('circle(%f, %f, %f)\n' % (ra, dec, angular_distance_arcsec/2.0))

    print("Input file: %s" % args.input)
    print("Output file: %s" % args.output)
    print("Side length of regions in pixels/arcsec: %f/ %f" % (args.diameter, angular_distance_arcsec))
    print("Overlap (fraction of region): %f" % args.fraction)
