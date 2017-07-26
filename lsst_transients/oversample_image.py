from scipy import ndimage
import scipy.misc
import numpy


def oversample(imageData, imageWCS, scaleFactor, threshold=None):
    '''
    
    :param imageData: Input array
    :param imageWCS: Header?
    :param scaleFactor: Zoom factor
    :param threshold:
    
    :return scaledData:
    :return scaledWCS:
    '''

    scaledData = scipy.misc.imresize(imageData, float(scaleFactor), interp='bilinear', mode='F')

    if threshold != None:
        # Resample with constant interpolation
        mask = ndimage.zoom(imageData, scaleFactor, order=0, mode='nearest')

        # Make a mask
        idx = mask <= threshold
        mask[idx] = 0
        mask[~idx] = 1

        # Restore zeros
        scaledData *= mask

        del mask

    # Take care of offset due to rounding in scaling image to integer pixel dimensions
    properDimensions = numpy.array(imageData.shape) * scaleFactor
    offset = properDimensions - numpy.array(scaledData.shape)

    # Rescale WCS
    try:
        oldCRPIX1 = imageWCS['CRPIX1']
        oldCRPIX2 = imageWCS['CRPIX2']
        CD11 = imageWCS['CD1_1']
        CD21 = imageWCS['CD2_1']
        CD12 = imageWCS['CD1_2']
        CD22 = imageWCS['CD2_2']
    except KeyError:
        # Try the older FITS header format
        try:
            oldCRPIX1 = imageWCS['CRPIX1']
            oldCRPIX2 = imageWCS['CRPIX2']
            CD11 = imageWCS['CDELT1']
            CD21 = 0
            CD12 = 0
            CD22 = imageWCS['CDELT2']
        except KeyError:
            scaledWCS = imageWCS.copy()
            return {'data': scaledData, 'wcs': scaledWCS}

    CDMatrix = numpy.array([[CD11, CD12], [CD21, CD22]], dtype=numpy.float64)
    scaleFactorMatrix = numpy.array([[1.0 / scaleFactor, 0], [0, 1.0 / scaleFactor]])
    scaledCDMatrix = numpy.dot(scaleFactorMatrix, CDMatrix)

    scaledWCS = imageWCS.copy()
    scaledWCS['NAXIS1'] = scaledData.shape[1]
    scaledWCS['NAXIS2'] = scaledData.shape[0]
    scaledWCS['CRPIX1'] = oldCRPIX1 * scaleFactor + offset[1]
    scaledWCS['CRPIX2'] = oldCRPIX2 * scaleFactor + offset[0]
    scaledWCS['CD1_1'] = scaledCDMatrix[0][0]
    scaledWCS['CD2_1'] = scaledCDMatrix[1][0]
    scaledWCS['CD1_2'] = scaledCDMatrix[0][1]
    scaledWCS['CD2_2'] = scaledCDMatrix[1][1]

    return scaledData, scaledWCS