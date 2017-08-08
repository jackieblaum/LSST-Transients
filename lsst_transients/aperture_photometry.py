import pandas as pd
import numpy as np
import os
import multiprocessing

import astropy.io.fits as pyfits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
import photutils

from lsst_transients.utils.loop_with_progress import loop_with_progress
from lsst_transients.utils.logging_system import get_logger

log = get_logger(__name__)

PHOTOMETRY_NULL_VALUE = -1e9


def _get_background_uncertainty(oversampled_bkgsub_image, oversampled_counts_image, mask_data):
    '''
    Computes the error of the flux for each region of the image.

    :param oversampled_bkgsub_image: A list of fluxes for each region for the background subtracted image
    :param oversampled_counts_image: A list of counts for each region for the original image (before subtraction)
    :param mask_data: A list of counts for each region for the mask image

    :return The error on the background for this visit
    '''

    # Transform the original mask (which contains not only 0 and 1, but also other numbers)
    # to a mask of booleans which selects the pixels used for the background estimation by the LSST software
    mask = (np.array(mask_data) == 0)

    # The flat background level is just the median of all the values of the background pixels

    bkgd_level = np.median(oversampled_counts_image[mask])  # type: float

    bkgd_level_check = np.median(oversampled_counts_image[mask] - oversampled_bkgsub_image[mask])  # type: float

    assert np.isclose(bkgd_level, bkgd_level_check, rtol=0.05), "Background level is too different:" \
                                                                "\nOriginal with mask: %f\n " \
                                                                "Original with mask minus background-subtracted: %f" \
                                                                % (bkgd_level, bkgd_level_check)

    bkgd_error = np.std(oversampled_counts_image[mask])

    return bkgd_error


def _get_region_fluxes(header_nobkgd, data_nobkgd, data_orig, data_masks, apertures):
    """
    Get the flux and its error for each region in this visit

    :param header_nobkgd:
    :param data_nobkgd:
    :param data_orig:
    :param data_masks:
    :param apertures:
    :return:
    """

    log.debug("Getting background uncertainty")

    # Get the background error (which is the same for all regions)
    background_uncertainty = _get_background_uncertainty(data_nobkgd, data_orig, data_masks)

    # Add the fluxes within each region by calling _sum_flux
    log.debug("Measuring flux for each region...")

    # We need the WCS

    w = WCS(header_nobkgd)

    # Transforms all the regions in pixels coordinates

    fluxes_table = photutils.aperture_photometry(data_nobkgd, apertures, error=np.sqrt(data_orig),
                                                 method='exact', wcs=w)

    log.debug("done")

    fluxes = fluxes_table['aperture_sum'].data

    log.debug("Computing errors...")

    # Adding the background uncertainty to the photometric uncertainty

    fluxes_errors = np.sqrt(fluxes_table['aperture_sum_err'].data ** 2 + background_uncertainty ** 2)

    log.debug("done")

    log.debug("Median value for the regions: %.3f" % np.nanmedian(fluxes))
    log.debug("Median error for the regions: %.3f" % np.nanmedian(fluxes_errors))

    return fluxes, fluxes_errors

def _get_apertures(grid):
    # NOTE: we assume all apertures have the same radius
    r = np.unique(grid.radiuses)
    assert r.shape[0] == 1, "Regions have different radiuses!"

    apertures = photutils.SkyCircularAperture(SkyCoord(ra=grid.centers_ra * u.deg,
                                                       dec=grid.centers_dec * u.deg,
                                                       frame='icrs'),
                                              r[0] * u.arcsec)

    return apertures


def find_visits_files(path, selected_filter):
    """
    Finds all visits files and return them in a list ordered by observation time

    :param path:
    :return:
    """

    # Collect all the visit files from the directory
    times_with_visits = []

    for root, dirs, files in os.walk(path):

        for name in files:

            if 'bkgd' not in name and '.fits' in name:

                filename = os.path.join(root, name)

                # Get filter name

                try:

                    with pyfits.open(filename) as f:

                        this_filter = f[0].header["FILTER"]
                        observation_time = f[0].header['MJD-OBS']
                        obsid = f[0].header['OBSID']

                except:

                    raise

                else:

                    # Add if it matches the filter we are looking for

                    if this_filter.replace(" ", "") == selected_filter:
                        times_with_visits.append((observation_time, filename, obsid))

    # Sort in tuples works by using the first element in each tuple
    times_with_visits.sort()

    visits_sorted = map(lambda x: x[1], times_with_visits)
    obsid_sorted = map(lambda x: x[2], times_with_visits)

    assert np.all(np.array(obsid_sorted) == np.unique(obsid_sorted)), "Observation IDs are not unique"

    return visits_sorted, obsid_sorted


def bulk_aperture_photometry(path, selected_filter, grid):
    '''
    Fills the dataframes that are indexed by visits. It first fills the flux tables and then fills the conditions table.

    :param path: The path to folder with all visit files
    :param selected_filter: The filter to be used
    :param grid: a Grid object

    :return None
    '''

    log.info("Selecting visits...")

    visits_sorted, obsid_sorted = find_visits_files(path, selected_filter)

    log.info("Found %i visits" % len(visits_sorted))

    log.info("Setting up regions...")
    all_apertures = _get_apertures(grid)

    log.info("Found %i regions" % grid.n_regions)

    # Create dataframes which will contain results, having as indexes the same index as the grid
    # (i.e., the ids of the regions)

    df = pd.DataFrame(columns=obsid_sorted, index=grid.region_ids)
    df_errors = pd.DataFrame(columns=obsid_sorted, index=grid.region_ids)

    headers_prim = []

    n_failed = 0

    pool = multiprocessing.Pool(processes=4,
                                initializer=process_setup,
                                initargs=[all_apertures])

    # Loop through the visits

    try:

        for i, results in loop_with_progress(pool.imap(worker, feeder(visits_sorted, obsid_sorted)),
                                             len(visits_sorted), 1, log.info, with_enumerate=True):

            obsid, fluxes, fluxes_errors, header_prim = results

            log.info("Processed OBSID: %s" % (obsid))

            if fluxes is None:

                n_failed += 1

                fluxes = PHOTOMETRY_NULL_VALUE
                fluxes_errors = PHOTOMETRY_NULL_VALUE

            df[obsid] = fluxes
            df_errors[obsid] = fluxes_errors

            headers_prim.append(header_prim)

    except:

        raise

    finally:

        pool.close()

    return df, df_errors, headers_prim, obsid_sorted, n_failed


def feeder(visits_sorted, obsid_sorted):

    for (visit, obsid) in zip(visits_sorted, obsid_sorted):

        yield visit, obsid


all_apertures_global = None


def process_setup(all_apertures_local):

    global all_apertures_global

    all_apertures_global = all_apertures_local


def worker(args):

    visit_file, obsid = args

    with pyfits.open(visit_file) as f:
        headers_prim = f[0].header

        header_nobkgd = f[1].header

        data_nobkgd = np.array(f[1].data, dtype='float32')
        data_orig = np.array(f[3].data, dtype='float32')
        data_mask = np.array(f[2].data, dtype='float32')

    # Call helper methods to fill in the fluxes and conditions for these visits
    try:

        fluxes, fluxes_errors = _get_region_fluxes(header_nobkgd, data_nobkgd, data_orig, data_mask,
                                                   all_apertures_global)

    except ValueError:

        fluxes = None
        fluxes_errors = None

    return obsid, fluxes, fluxes_errors, headers_prim