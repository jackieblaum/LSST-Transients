import numpy as np
import re

import pyregion
from astropy.wcs.utils import proj_plane_pixel_scales
import astropy.units as u

class CircularRegion(object):
    '''
    A CircularRegion object is a square or circular region, many of which form a grid.
    '''

    def __init__(self, x, y, r, ds9_string=None):
        '''
        Constructs a CircularRegion object.

        :param x: The horizontal location of the center of the region
        :param y: The vertical location of the center of the region
        :param r: Radius
		  of the region if it is a circle (in pixels)
        '''

        self._x = x
        self._y = y
        self._radius_pixels = r

        self._mask = None

        self._ds9_string = ds9_string

        if self._ds9_string is not None:

            self._pyregion = pyregion.parse(self._ds9_string)

    @classmethod
    def from_ds9_region(cls, w, ds9_string):
        """
        Returns a mask which selects all the pixels within the provided region

        :param w: The WCS from the header
        :param ds9_string: The information for a region in a string format

        :return mask array, corners of bounding box

        """

        # Get info about the region from the DS9 string
        split = re.split("[, (\")]+", ds9_string)
        shape = split[0]

        assert shape == "circle", "Only circular regions are supported"

        ra = float(split[1])
        dec = float(split[2])
        radius = float(split[3])

        # Get the radius in pixels assuming that the radius is the same in WCS for all regions

        pixel_scale = proj_plane_pixel_scales(w)

        assert np.isclose(pixel_scale[0], pixel_scale[1],
                          rtol=1e-2), "Pixel scale is different between X and Y direction"

        # We take the geometric average of the pixel scales in the X and Y direction, as done
        # in pyregion

        pixel_scale_with_units = (pixel_scale[0] * pixel_scale[1])**0.5 * u.Unit(w.wcs.cunit[0])

        pixel_scale_arcsec = pixel_scale_with_units.to("arcsec").value

        radius_pix = radius / pixel_scale_arcsec

        # Find the smallest square around the region
        x, y = w.all_world2pix([[ra, dec]], 0)[0]

        # Make a region object given the information found previously
        reg = cls(x, y, radius_pix, ds9_string=ds9_string)

        return reg

    @property
    def x(self):

        return self._x

    @property
    def y(self):

        return self._y

    @property
    def r(self):

        return self._radius_pixels

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

        r = 'image;circle(%s,%s,%s)' % (self.x, self.y, self._radius_pixels)

        return r

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

        corner1 = [int(np.floor(float(self.x) - float(self.r) - padding)),
                   int(np.floor(float(self.y) - float(self.r) - padding))]
        corner2 = [int(np.floor(float(self.x) - float(self.r) - padding)),
                   int(np.ceil(float(self.y) + float(self.r) + padding))]
        corner3 = [int(np.ceil(float(self.x) + float(self.r) + padding)),
                   int(np.floor(float(self.y) - float(self.r) - padding))]
        corner4 = [int(np.ceil(float(self.x) + float(self.r) + padding)),
                   int(np.ceil(float(self.y) + float(self.r) + padding))]

        return [self._fix_point(corner1, data_shape),
                self._fix_point(corner2, data_shape),
                self._fix_point(corner3, data_shape),
                self._fix_point(corner4, data_shape)]

    def compute_mask(self, data_array):

        nx, ny = data_array.shape

        # We need to start at 1 because ds9 regions follow the FORTRAN convention (arrays start at 1)

        y, x = np.ogrid[1:nx+1, 1:ny+1]

        mask = ((x-self.x)**2 + (y-self.y)**2 <= self._radius_pixels**2)  # type: np.ndarray

        self._mask = mask

        return mask

    def _compute_mask_pyregion(self, data, header):

        #reg_definition = 'icrs;%s' % ds9_string.replace(" ", "")

        filt = self._pyregion.as_imagecoord(header).get_filter()

        mask = filt.mask(data)

        return mask

    def sum_within_region(self, array):

        if self._mask is None:

            self.compute_mask(array)

        if np.sum(self._mask) == 0:

            return np.nan

        else:

            return np.sum(array[self._mask])

