from scipy import ndimage
import scipy.misc
import numpy


def oversample(image_data, image_header, scale_factor, threshold=None, method='bilinear'):
    '''
    
    :param image_data: Input array
    :param image_header: Header?
    :param scale_factor: Zoom factor
    :param threshold:
    
    :return scaledData:
    :return scaledWCS:
    '''

    if scale_factor == 1:

        return image_data, image_header

    scaledData = scipy.misc.imresize(image_data, float(scale_factor), interp=method, mode='F')

    if threshold != None:
        # Resample with constant interpolation
        mask = ndimage.zoom(image_data, scale_factor, order=0, mode='nearest')

        # Make a mask
        idx = mask <= threshold
        mask[idx] = 0
        mask[~idx] = 1

        # Restore zeros
        scaledData *= mask

        del mask

    # Take care of offset due to rounding in scaling image to integer pixel dimensions
    properDimensions = numpy.array(image_data.shape) * scale_factor
    offset = properDimensions - numpy.array(scaledData.shape)

    # Rescale WCS
    try:
        oldCRPIX1 = image_header['CRPIX1']
        oldCRPIX2 = image_header['CRPIX2']
        CD11 = image_header['CD1_1']
        CD21 = image_header['CD2_1']
        CD12 = image_header['CD1_2']
        CD22 = image_header['CD2_2']
    except KeyError:
        # Try the older FITS header format
        try:
            oldCRPIX1 = image_header['CRPIX1']
            oldCRPIX2 = image_header['CRPIX2']
            CD11 = image_header['CDELT1']
            CD21 = 0
            CD12 = 0
            CD22 = image_header['CDELT2']
        except KeyError:
            scaledWCS = image_header.copy()
            return {'data': scaledData, 'wcs': scaledWCS}

    CDMatrix = numpy.array([[CD11, CD12], [CD21, CD22]], dtype=numpy.float64)
    scaleFactorMatrix = numpy.array([[1.0 / scale_factor, 0], [0, 1.0 / scale_factor]])
    scaledCDMatrix = numpy.dot(scaleFactorMatrix, CDMatrix)

    scaledWCS = image_header.copy()
    scaledWCS['NAXIS1'] = scaledData.shape[1]
    scaledWCS['NAXIS2'] = scaledData.shape[0]
    scaledWCS['CRPIX1'] = oldCRPIX1 * scale_factor + offset[1]
    scaledWCS['CRPIX2'] = oldCRPIX2 * scale_factor + offset[0]
    scaledWCS['CD1_1'] = scaledCDMatrix[0][0]
    scaledWCS['CD2_1'] = scaledCDMatrix[1][0]
    scaledWCS['CD1_2'] = scaledCDMatrix[0][1]
    scaledWCS['CD2_2'] = scaledCDMatrix[1][1]

    return scaledData, scaledWCS