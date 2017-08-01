import numpy as np

from lsst_transients.utils.cartesian_product import cartesian_product

class CircularRegion(object):
    '''
    A CircularRegion object is a square or circular region, many of which form a grid.
    '''

    def __init__(self, x, y, d):
        '''
        Constructs a CircularRegion object.

        :param x: The horizontal location of the center of the region
        :param y: The vertical location of the center of the region
        :param d: The length of the sides of the region if it is a square, otherwise the diameter
		  of the region if it is a circle (in pixels)
        '''

        self._x = x
        self._y = y
        self._diameter_pixels = d

        self._mask = None

    @property
    def x(self):

        return self._x

    @property
    def y(self):

        return self._y

    @property
    def d(self):

        return self._diameter_pixels

    def _get_mask(self):

        return self._mask

    def _set_mask(self, new_mask):

        self._mask = new_mask  # type: np.ndarray

    mask = property(_get_mask, _set_mask, doc="Sets/gets the mask for this region")

    def __str__(self):
        '''
        Overrides the string function for a Box object.

        :return: a string in the form "box(_x, _y, width, height, angle)" to create a box region in DS9
        '''
        return "circle(" + str(self.x) + ", " + str(self._y) + ", " + str(self._diameter_pixels) + ")"

    def _fix_x(self, x, data_shape):

        return min(max(0, x), data_shape[1])

    def _fix_y(self, y, data_shape):

        return min(max(0, y), data_shape[0])

    def _fix_point(self, point, data_shape):

        return [self._fix_x(point[0], data_shape), self._fix_y(point[1], data_shape)]

    def get_boundingbox(self, data_shape, padding=5):
        '''
        :param padding: padding for the bounding box (default: 5)
        :return corners: An array of sky coordinate pairs for the four corners of the bounding box with the provided
        padding
        '''

        corner1 = [int(np.floor(float(self.x) - (float(self.d) / 2.0) - padding)),
                   int(np.floor(float(self.y) - (float(self.d) / 2.0) - padding))]
        corner2 = [int(np.floor(float(self.x) - (float(self.d) / 2.0) - padding)),
                   int(np.ceil(float(self.y) + (float(self.d) / 2.0) + padding))]
        corner3 = [int(np.ceil(float(self.x) + (float(self.d) / 2.0) + padding)),
                   int(np.floor(float(self.y) - (float(self.d) / 2.0) - padding))]
        corner4 = [int(np.ceil(float(self.x) + (float(self.d) / 2.0) + padding)),
                   int(np.ceil(float(self.y) + (float(self.d) / 2.0) + padding))]

        return [self._fix_point(corner1, data_shape),
                self._fix_point(corner2, data_shape),
                self._fix_point(corner3, data_shape),
                self._fix_point(corner4, data_shape)]

    def _get_distances(self, x, y, x2, y2):
        '''
        Get the distance between two points.

        :param x: x-coordinate(s) of the first point(s)
        :param y: y-coordinate(s) of the first point(s)
        :param x2: x-coordinate(s) of the second point(s)
        :param y2: y-coordinate(s) of the second point(s)

        :return distances: An array of the distances between the pairs of points given
        '''

        distances = np.sqrt((x - x2) ** 2 + (y - y2) ** 2)

        return distances

    def compute_mask(self, data_array):

        nx, ny = data_array.shape
        y, x = np.ogrid[-self.x:nx - self.x, -self.y:ny - self.y]
        mask = (x * x + y * y <= (self._diameter_pixels / 2.0)**2)  # type: np.ndarray

        self._mask = mask

        return mask

    def sum_within_region(self, array):

        if self._mask is None:

            self.compute_mask(array)

        if np.sum(self._mask) == 0:

            return np.nan

        else:

            return np.sum(array[self._mask])

