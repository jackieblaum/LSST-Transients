import numpy as np
import astropy.units as u

from astropy.coordinates.angle_utilities import angular_separation
from region import Region
from create_db_notpyregion import Data_Database
from utils.cartesian_product import cartesian_product
from utils.pix2world import pix2world

class Grid(object):
    '''
    A Grid object contains overlapping circular or square regions that can overlay an image. The regions can be stored
    in a .reg file and in a database.
    '''

    def __init__(self, xcenter, ycenter, xmin, xmax, ymin, ymax, overlapfactor, d, shape):
        '''
        Constructs a grid object. A region list is instantiated, and the overlap length is set to be the product of the
        overlap factor and diameter/side length.

        :param xcenter: The _x-coordinate of the center of the image
        :param ycenter: The _y-coordinate of the center of the image
        :param xmin: The minimum x-coordinate of the image
        :param xmax: The maximum x-coordinate of the image
        :param ymin: The minimum y-coordinate of the image
        :param ymax: The maximum y-coordinate of the image
        :param overlapfactor: The fraction of the length of the region by which the regions overlap one another on the grid
        :param d: The length of the sides or diameter of the regions in the grid
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


    def _in_range(self, a_region):
        '''
        A helper method used by do_steps() in order to determine if the Region is in range of the region of
        interest and should be added to the Region array.

        Returns true if out of range, returns false if in range.

        :param a_region: the region to be checked
        :return: False if the region is out of the region of interest, True otherwise
        '''

        # TODO: account for rotation angle

        x = a_region.x
        y = a_region.y

        # Set the limits such that all parts of the image will be covered by at least one region
        if self._xmin-1 - (self.d/2) < x < self._xmax+1 + (self.d/2) and self._ymin-1 - (self.d/2) < y < self._ymax+1 + (self.d/2):

            return True

        else:

            return False


    def _get_all_coords(self, delta):
        '''
        Compute the x- and y-coordinate for the centers of all regions in the grid by taking the cartesian product, then
        create all the corresponding Region objects.

        :param delta: The longest diagonal from the center of the image to the edge
        :return: The Region objects for all regions in the grid
        '''

        # Let's figure out the iteration we are in
        iteration = np.floor(delta / (self.d - self.overlap))
        
        # How many steps do we need to complete one side of the round?
        steps = iteration * 2

        # Let's compute the start point of the round
        start_x = self.xcenter + (self.d - self.overlap)*iteration
        start_y = self.ycenter + (self.d - self.overlap)*iteration

        # Get all the x- and y- coordinates for the centers of the regions
        xs = np.linspace(start_x - steps * (self.d - self.overlap), start_x, steps + 1)
        ys = np.linspace(start_y - steps * (self.d - self.overlap), start_y, steps + 1)

        # Get all the coordinate pairs and make the Region objects
        this_regions = map(lambda (x, y): Region(x, y, self.d, self.shape), cartesian_product([xs, ys]))

        return this_regions


    @staticmethod
    def distance(x1, y1, x2, y2):
        '''
        Compute the distance between two points.

        :param x1: The x-coordinate of the first point
        :param y1: The y-coordinate of the first point
        :param x2: The x-coordinate of the second point
        :param y2: The y-coordinate of the second point
        :return: The distance between the two points
        '''

        return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)


    def get_grid(self):
        '''
        After finding the longest diagonal between the center and edge of the image, call _get_all_coords to get all the
        Region objects, then make sure they are all within the desired bounds by calling _in_range.

        :return: The array of Region objects that are within the bounds of the image.
        '''

        # Figure out maximum diagonal
        d1 = self.distance(self._xmin, self._ymin, self.xcenter, self.ycenter)
        d2 = self.distance(self._xmin, self._ymax, self.xcenter, self.ycenter)

        longest_diagonal = max(d1, d2)

        all_regions = self._get_all_coords(longest_diagonal)

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


    def write_grid(self, infile, outfile, shape, rotation_angle, name):
        '''
        Converts the grid coordinates from pixels to WCS and writes the grid to a file.

        :param infile: The name of the fits file inputted by the user
        :param outfile: The name of the region file to be outputted, specified by the user
        :param shape: Either circle or square, as specified by the user
        :param rotation_angle: The angle of rotation of the image in the fits file
        :return: The angular distance of the diameter/length of the side of the region
        '''

        # Get all the region objects
        regions = self.get_grid()

        # Convert the locations of the regions to WCS
        regions_wcs = self.pix_to_wcs(infile, regions)

        # Convert the side length/diameter of the regions to WCS
        # NOTE: we assume that on this angular scale there is no distortion due to the projection
        # TODO: fix the projection for the entire LSST field of view (3.5 deg)

        # Get the angular distance between the upper left corner and the upper right corner of the
        # center region (doesn't make a difference if the region is a circle or square)
        center_region = regions[0]
        upper_left_x = center_region.x - center_region.d / 2.0
        upper_y = center_region.y - center_region.d / 2.0
        upper_right_x = center_region.x + center_region.d / 2.0

        upper_left_sky_coords, upper_right_sky_coords = pix2world(infile,
                                                              [[upper_left_x, upper_y], [upper_right_x, upper_y]])

        angular_distance = angular_separation(upper_left_sky_coords[0] * u.deg, upper_left_sky_coords[1] * u.deg,
                                          upper_right_sky_coords[0] * u.deg, upper_right_sky_coords[1] * u.deg)

        angular_distance_arcsec = angular_distance.to("arcsec").value

        # Write the grid to a DS9 region file
        with open(outfile, "w+") as f:
            f.write("icrs\n")

            for (ra, dec) in regions_wcs:
                if shape == "square":

                    f.write('box(%f, %f, %f", %f", %f)\n' % (ra, dec, angular_distance_arcsec, angular_distance_arcsec,
                                                     360 - rotation_angle))
                elif shape == "circle":

                    f.write('circle(%f, %f, %f")\n' % (ra, dec, angular_distance_arcsec/2.0))
        
        # Write the regions to a database
        db = Data_Database("%s.db" % name)
        num_regs = db.fill_reg(outfile)

        # Initialize empty flux tables and condition table
        db.init_flux_tables(num_regs)
        db.init_cond_table()
        db.close()

        return angular_distance_arcsec

