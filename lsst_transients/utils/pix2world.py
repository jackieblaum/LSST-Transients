import numpy
from astropy import wcs
from astropy.io import fits as pyfits

def pix2world(filename, array):
    '''
    Converts the pixels of the array into the WCS, which is read from the given file.

    :param filename: The file to be read
    :param array: The array to be converted from pixels to WCS
    :return: The converted array (in WCS)
    '''

    # Load the FITS hdulist
    hdulist = pyfits.open(filename)

    # Parse the WCS keywords in the primary HDU
    w = wcs.WCS(hdulist[1].header)

    # Pixel coordinates
    pixels = numpy.array(array, numpy.float)

    # Convert the pixels to WCS
    world = w.wcs_pix2world(pixels, 0)

    return world
