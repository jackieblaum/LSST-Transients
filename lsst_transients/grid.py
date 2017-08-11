import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates.angle_utilities import angular_separation

from circular_region import CircularRegion
from utils.cartesian_product import cartesian_product
from utils.pix2world import pix2world

GRID_FRAME_NAME = 'regions'


class Grid(object):

    def __init__(self, region_ids, centers_ra, centers_dec, radiuses, df=None):

        if df is None:

            # Normal construction

            self._df = pd.DataFrame(columns=['ra', 'dec', 'radius'], index=region_ids)

            self._df['ra'] = centers_ra
            self._df['dec'] = centers_dec
            self._df['radius'] = radiuses

        else:

            # Used by from_hdf

            self._df = df

    @classmethod
    def read_from(cls, filename):

        df = pd.read_hdf(filename, GRID_FRAME_NAME)

        return cls(None, None, None, None, df=df)

    def write_to(self, filename):

        self._df.to_hdf(filename, GRID_FRAME_NAME, format='table')

    @property
    def n_regions(self):

        return self._df.shape[0]

    @property
    def centers_ra(self):
        return self._df['ra']

    @property
    def centers_dec(self):
        return self._df['dec']

    @property
    def radiuses(self):
        return self._df['radius']

    @property
    def region_ids(self):
        return self._df.index

    def __getitem__(self, region_id):

        return self._df.loc[region_id]

    def get_ds9_region(self, region_id):

        region = self[region_id]

        return 'circle(%s,%s,%s")' % (region['ra'], region['dec'], region['radius'])

    def write_regions(self):

        with open('regfile.reg', 'w+') as f:
            f.write('icrs\n')
            for region_id in self._df.index:
                region = self[region_id]
                f.write('circle(%s,%s,%s")\n' % (region['ra'], region['dec'], region['radius']))



class GridFactory(object):
    '''
    A Grid object contains overlapping circular or square regions that can overlay an image. The regions can be stored
    in a .reg file and in a _database.
    '''

    def __init__(self, xmin, xmax, ymin, ymax, overlapfactor, r):
        '''
        Constructs a grid object. A region list is instantiated, and the _overlap length is set to be the product of the
        _overlap factor and radius/side length.

        :param _xcenter: The _x-coordinate of the center of the image
        :param _ycenter: The _y-coordinate of the center of the image
        :param xmin: The minimum x-coordinate of the image
        :param xmax: The maximum x-coordinate of the image
        :param ymin: The minimum y-coordinate of the image
        :param ymax: The maximum y-coordinate of the image
        :param overlapfactor: The fraction of the length of the region by which the regions _overlap one another on the grid
        :param r: The radius of the regions in the grid
        '''

        self._xmin = xmin
        self._xmax = xmax
        self._ymin = ymin
        self._ymax = ymax

        # The center of the image in pixel coordinates

        self._xcenter = self._xmin + (self._xmax - self._xmin + 1) / 2.0
        self._ycenter = self._ymin + (self._ymax - self._ymin + 1) / 2.0

        assert float(overlapfactor) < 1, "The _overlap factor must be < 1"

        self._overlapfactor = float(overlapfactor)
        self._r = r
        self._overlap = overlapfactor * r * 2

        self._regions_wcs = None
        self._angular_distance_arcsec = None

    def _in_range(self, a_region):
        '''
        A helper method used by do_steps() in order to determine if the CircularRegion is in range of the region of
        interest and should be added to the CircularRegion array.

        Returns true if out of range, returns false if in range.

        :param a_region: the region to be checked
        :return: False if the region is out of the region of interest, True otherwise
        '''

        # TODO: account for rotation angle

        x = a_region.x
        y = a_region.y

        # Set the limits such that all parts of the image will be covered by at least one region
        if self._xmin - 1 - self._r < x < self._xmax + 1 + self._r and self._ymin - 1 - self._r < y < self._ymax + 1 + self._r:

            return True

        else:

            return False

    def _get_all_coords(self, delta):
        '''
        Compute the x- and y-coordinate for the centers of all regions in the grid by taking the cartesian product, then
        create all the corresponding CircularRegion objects.

        :param delta: The longest diagonal from the center of the image to the edge
        :return: The CircularRegion objects for all regions in the grid
        '''

        # Let's figure out the iteration we are in
        iteration = np.floor(delta / ((self._r * 2) - self._overlap))

        # How many steps do we need to complete one side of the round?
        steps = iteration * 2

        # Let's compute the start point of the round
        start_x = self._xcenter + ((self._r * 2) - self._overlap) * iteration
        start_y = self._ycenter + ((self._r * 2) - self._overlap) * iteration

        # Get all the x- and y- coordinates for the centers of the regions
        xs = np.linspace(start_x - steps * ((self._r * 2) - self._overlap), start_x, steps + 1)
        ys = np.linspace(start_y - steps * ((self._r * 2) - self._overlap), start_y, steps + 1)

        # Get all the coordinate pairs and make the CircularRegion objects
        this_regions = map(lambda (x, y): CircularRegion(x, y, self._r), cartesian_product([xs, ys]))

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

    def _compute_grid(self):
        '''
        After finding the longest diagonal between the center and edge of the image, call _get_all_coords to get all the
        CircularRegion objects, then make sure they are all within the desired bounds by calling _in_range.

        :return: The array of CircularRegion objects that are within the bounds of the image.
        '''

        # Figure out maximum diagonal
        d1 = self.distance(self._xmin, self._ymin, self._xcenter, self._ycenter)
        d2 = self.distance(self._xmin, self._ymax, self._xcenter, self._ycenter)

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


    def get_grid(self, infile):
        '''
        Converts the grid coordinates from pixels to WCS and writes the grid to a file.

        :param infile: The name of the fits file inputted by the user
        :return: List of regions in sky coordinates and the angular distance in arcsec
        '''

        # Get all the region objects
        regions = self._compute_grid()

        # Convert the locations of the regions to WCS
        regions_wcs = self.pix_to_wcs(infile, regions)

        # Convert the radius of the regions to WCS
        # NOTE: we assume that on this angular scale there is no distortion due to the projection
        # TODO: fix the projection for the entire LSST field of view (3.5 deg)

        # Get the angular distance between the left edge and the center of the
        # center region
        center_region = regions[0]
        left_x = center_region.x - center_region.r
        left_y = center_region.y

        left_sky_coords, center_sky_coords = pix2world(infile,[[left_x, left_y], [center_region.x, center_region.y]])

        angular_distance = angular_separation(left_sky_coords[0] * u.deg, left_sky_coords[1] * u.deg,
                                              center_sky_coords[0] * u.deg, center_sky_coords[1] * u.deg)

        angular_distance_arcsec = angular_distance.to("arcsec").value

        ras = np.zeros(len(regions))
        decs = np.zeros(len(regions))
        ids = []

        for i, (ra, dec) in enumerate(regions_wcs):

            ras[i] = ra
            decs[i] = dec
            ids.append("reg%i" % i)

        return Grid(ids, ras, decs, [angular_distance_arcsec] * len(ras))
