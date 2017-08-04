import pandas as pd
import numpy as np
import os
import re
import functools
import multiprocessing
import ctypes

import astropy.units as u
import astropy.io.fits as pyfits

from circular_region import CircularRegion
from astropy import wcs
from oversample_image import oversample
from utils import database_io
from utils.chuncker import chunker
from utils.logging_system import get_logger
from utils.loop_with_progress import loop_with_progress

log = get_logger(os.path.basename(__file__))

# These two globals will be used to store a shared array between the processes when computing the fluxes
# Unfortunately having global variables is the only way to share memory between processes
data_t = [None, None, None]
orig_t = [None, None, None]


def shifter(ds9_string, w, mask_x, mask_y, mask_i_nonzero, mask_j_nonzero):

    this_region = CircularRegion.from_ds9_region(w, ds9_string)

    # Now shift the mask

    dx = float(this_region.x - mask_x)
    dy = float(this_region.y - mask_y)

    this_mask_i = mask_i_nonzero + int(dy)  # type: np.ndarray
    this_mask_j = mask_j_nonzero + int(dx)  # type: np.ndarray

    return this_mask_i, this_mask_j


def worker(ds9_string, w, mask_x, mask_y, mask_i_nonzero, mask_j_nonzero, bkgd_uncertainty, maxi, maxj):

    # Gather the data and the orig array from the shared multiprocessing array
    # and transform them back to numpy (note that this does not copy the data, it only re-interprets them)
    data = np.frombuffer(data_t[0].get_obj(), dtype=data_t[1]).reshape(data_t[2])
    orig = np.frombuffer(orig_t[0].get_obj(), dtype=orig_t[1]).reshape(orig_t[2])

    this_mask_i, this_mask_j = shifter(ds9_string, w, mask_x, mask_y, mask_i_nonzero, mask_j_nonzero)

    # Make sure that we do not "spill out"

    idx = (this_mask_i >= 0) & (this_mask_i <= maxi) & (this_mask_j >= 0) & (this_mask_j <= maxj)

    if np.sum(idx) == 0:

        # empty region
        flux = np.nan
        flux_error = np.nan

    else:

        this_mask_i = this_mask_i[idx]
        this_mask_j = this_mask_j[idx]

        # Add the fluxes of the pixels within the bounding box

        flux = np.sum(data[this_mask_i, this_mask_j])

        # Now compute the error
        counts_sum = np.sum(orig[this_mask_i, this_mask_j])

        flux_error = np.sqrt(counts_sum + bkgd_uncertainty ** 2)

    return flux, flux_error


class DataDatabase(object):
    '''
    The DataDatabase class is used to create and handle a _database and write to different tables: a region table,
    a flux table for each region, and a conditions table.
    '''

    def __init__(self, dbname, first=True):
        '''
        Initializes the _database and opens a connection to the _database.

        :param dbname: The name for the _database
        :param first: If the _database has not yet been created, True, False otherwise
        '''

        if first == False:
            assert os.path.exists(dbname), "Database %s does not exist" % dbname

        #self.db = database_io.SqliteDatabase(dbname)
        self.db = database_io.HDF5Database(dbname)

        self.db.connect()

    def fill_reg(self, reg_list_wcs, shape, angular_distance_arcsec, rotation_angle):
        '''
        Fills the _database with the string for each region as seen in DS9.

        :param reg_list_wcs: List of region centers (a list of tuples) in sky coordinates

        :return The number of regions
        '''

        if shape == "square":

            ds9_representation = map(lambda (ra, dec):'box(%f, %f, %f", %f", %f)' % (ra, dec,
                                                                                     angular_distance_arcsec,
                                                                                     angular_distance_arcsec,
                                                                                     360 - rotation_angle),
                                     reg_list_wcs)

        elif shape == "circle":

            ds9_representation = map(lambda (ra, dec): 'circle(%f, %f, %f")' % (ra, dec, angular_distance_arcsec),
                                     reg_list_wcs)

        else:

            raise ValueError("Shape %s is not valid. Only 'circle' and 'square' are supported. " % shape)

        # Fill the index array now that we know how many regions there were
        indices = range(1, len(ds9_representation)+1)

        series = pd.Series(ds9_representation, index=indices)

        reg_dataframe = pd.DataFrame.from_dict({'ds9_info': series})
        self.db.insert_dataframe(reg_dataframe, 'reg_dataframe')

        return len(indices)

    def _get_c_type(self, data):

        np_type = data.dtype.name

        if np_type == 'float64':

            c_type = ctypes.c_double

        elif np_type == 'float32':

            c_type = ctypes.c_float

        else:

            raise TypeError("Type %s not supported" % np_type)

        return c_type, np_type

    def _get_fluxes(self, reg, data, orig, header, bkgd_uncertainty, pool=None):
        '''
        Gets the fluxes for all of the regions in the image.

        :param reg: The region dataframe created with lsst_grid_generator_shapes.py
        :param data: The data array of the background-subtracted image (fluxes)
        :param orig: The data array of the image before background-subtraction
        :param header: The header of the visit file to be examined
        :param bkgd_uncertainty: The error value for the background

        :return fluxes: An array of the fluxes and flux errors for each region in the image
        '''

        global data_t
        global orig_t

        # Get the number of regions and use this number to initialize an array that will store the fluxes for each region
        num_regs = len(reg.index)
        fluxes = np.zeros(num_regs)
        fluxes_errors = np.zeros(num_regs)

        # Add the fluxes within each region by calling _sum_flux
        log.info("Measuring flux for each region\n")

        # We need the WCS

        w = wcs.WCS(header)

        # NOTE: regions are indexed starting from 1, not 0

        # First we consider the region in the center, and we compute its mask. Then we will just shift that mask
        # of the appropriate distance to get the mask for all other regions (much faster than recomputing it all the
        # time)

        # Get size of the first region
        region_diameter_pixels = CircularRegion.from_ds9_region(w, reg.get_value(1, "ds9_info")).d

        # Get the mask from the central region
        central_region = CircularRegion(np.ceil(data.shape[0] / 2.0),
                                        np.ceil(data.shape[1] / 2.0),
                                        region_diameter_pixels)

        mask = central_region.compute_mask(data)

        # Now get its center in pixel (so later we can compute by how much we need to move the mask)
        mask_x = central_region.x
        mask_y = central_region.y

        # Divide the mask into its i and j coordinates
        nz = mask.nonzero()
        mask_i_nonzero = nz[0]  # type: np.ndarray
        mask_j_nonzero = nz[1] # type: np.ndarray

        maxi = data.shape[0] - 1
        maxj = data.shape[1] - 1

        # Make a copy of the data in a form that multiprocessing can share
        data_c_type, data_np_type = self._get_c_type(data)
        scaled_data_nobkgd_shared = multiprocessing.Array(data_c_type, data.size)
        scaled_data_nobkgd_shared[:] = data.flatten()

        orig_c_type, orig_np_type = self._get_c_type(orig)
        orig_shared = multiprocessing.Array(orig_c_type, orig.size)
        orig_shared[:] = orig.flatten()

        data_t[:] = (scaled_data_nobkgd_shared, data_np_type, data.shape)
        orig_t[:] = (orig_shared, orig_np_type, orig.shape)

        partial_worker = functools.partial(worker, w=w, mask_x=mask_x, mask_y=mask_y,
                                           mask_i_nonzero=mask_i_nonzero, mask_j_nonzero=mask_j_nonzero,
                                           bkgd_uncertainty=bkgd_uncertainty,
                                           maxi=maxi, maxj=maxj)

        #print map(partial_worker, reg['ds9_info'].values[0:10])

        chunksize = 1000

        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

        try:

            for i, (fl, fl_error) in loop_with_progress(pool.imap(partial_worker, reg['ds9_info'], chunksize=chunksize),
                                                            reg['ds9_info'].shape[0], chunksize, log.info,
                                                            with_enumerate=True):

                fluxes[i] = fl
                fluxes_errors[i] = fl_error

        except:

            raise

        finally:

            pool.close()

        # for i in loop_with_progress(range(1, _n_regions+1), _n_regions, 10000, my_print):
        #
        #     # Get the ds9 definition
        #     ds9_string = reg.get_value(i, "ds9_info")
        #
        #     this_region = _get_region(w, ds9_string)
        #
        #     # Now shift the mask
        #     # NOTE that by definition these are integers because the regions are built by shifting them by integer
        #     # amount of pixels
        #     dx = int(this_region.x - mask_x)
        #     dy = int(this_region.y - mask_y)
        #
        #     this_mask_i = mask_x_nonzero + dx
        #     this_mask_j = mask_y_nonzero + dy
        #
        #     # Make sure that we do not "spill out"
        #
        #     idx = numexpr.evaluate("(this_mask_i >= 0) & (this_mask_i <= maxx) & "
        #                            "(this_mask_j >= 0) & (this_mask_j <= maxy)")
        #
        #     if np.sum(idx) == 0:
        #
        #         # empty region
        #         flux = np.nan
        #         flux_error = np.nan
        #
        #         n_zeros += 1
        #
        #     else:
        #
        #         this_mask_i = this_mask_i[idx]
        #         this_mask_j = this_mask_j[idx]
        #
        #         # Add the fluxes of the pixels within the bounding box
        #
        #         flux = np.sum(data[this_mask_i, this_mask_j])
        #
        #         # Now compute the error
        #         counts_sum = np.sum(orig[this_mask_i, this_mask_j])
        #
        #         flux_error = np.sqrt(counts_sum + bkgd_uncertainty ** 2)
        #
        #         n_non_zeros += 1
        #
        #     fluxes[i-1] = flux
        #     fluxes_errors[i-1] = flux_error

        log.info("Median value for the regions: %.3f" % np.nanmedian(fluxes))
        log.info("Median error for the regions: %.3f" % np.nanmedian(fluxes_errors))

        return fluxes, fluxes_errors


    def _get_background_uncertainty(self, oversampled_bkgsub_image, oversampled_counts_image, mask_data):
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


    def _get_data(self, dtype, *args, **kwargs):

        this_data = pyfits.getdata(*args, **kwargs)

        return np.array(this_data, dtype=dtype)


    def _fill_flux(self, headers_nobkgd, data_nobkgd, headers_orig, data_orig, headers_masks, data_masks):
        '''
        Fills the dataframe with the flux and the flux error for each region with the indices as the visit number.

        :param headers_nobkgd: An array of the headers from the background-subtracted images from each of the visits
        :param data_nobkgd: An array of the data from the background-subtracted images from each of the visits
        :param headers_orig: An array of the headers from the original images from each of the visits
        :param data_orig: An array of the data from the original images from each of the visits
        :param headers_masks: An array of the headers from the mask images from each of the visits
        :param data_masks: An array of the data from the mask images from each of the visits

        :return None
        '''

        # Arrays for all the visits that will store arrays of fluxes and flux errors for each region
        visit_fluxes = []
        visit_errs = []

        # Access the regions table in the _database in order to find the number of regions
        reg = self.db.get_table_as_dataframe('reg_dataframe')
        num_regs = len(reg.index)

        # Loop through the visit files
        n_visits_in_this_chunk = len(headers_nobkgd)

        for i in range(0, n_visits_in_this_chunk):

            log.info("Processing file %i of %i..." % ((i+1), n_visits_in_this_chunk))

            log.debug("oversampling")

            # Oversample the background-subtracted, the original images, and the mask images.
            scale_factor = 4

            scaled_data_nobkgd, scaled_wcs_nobkgd = oversample(data_nobkgd[i], headers_nobkgd[i], scale_factor)

            scaled_data_orig, scaled_wcs_orig = oversample(data_orig[i], headers_orig[i], scale_factor)

            scaled_data_mask, scaled_wcs_mask = oversample(data_masks[i], headers_masks[i], scale_factor,
                                                           method='nearest')

            log.debug("Getting background uncertainty")
            # Get the background error (which is the same for all regions)
            background_uncertainty = self._get_background_uncertainty(scaled_data_nobkgd, scaled_data_orig,
                                                                      scaled_data_mask)
            # Get the fluxes for each region for the scaled images
            log.debug("Getting fluxes for each region")

            fluxes_nobkgd, fluxes_errors = self._get_fluxes(reg, scaled_data_nobkgd, scaled_data_orig,
                                                            scaled_wcs_nobkgd, background_uncertainty)

            # Normalize the scaled images
            norm_factor = scale_factor **2

            fluxes_nobkgd /= norm_factor

            # Arrays store the fluxes and flux errors for each region for each visit.
            visit_fluxes.append(fluxes_nobkgd)
            visit_errs.append(fluxes_errors)

        # Switch the order of the arrays in order to store the data
        region_fluxes = np.swapaxes(visit_fluxes, 0, 1)
        region_errs = np.swapaxes(visit_errs, 0, 1)

        # Write the fluxes to the _database
        log.info("Writing to _database...\n")

        for r in loop_with_progress(range(num_regs), num_regs, 10000, log.info, False):

            flux_dataframe = pd.DataFrame(columns=['flux', 'err'])
            flux_dataframe['flux'] = region_fluxes[r]
            flux_dataframe['err'] = region_errs[r]

            self.db.insert_dataframe(flux_dataframe, 'flux_table_%i' % (r + 1))

        return None


    def _fill_cond(self, headers):
        '''
        Fills the dataframe with the conditions for each visit (seeing, duration, and date). Seeing at 5000 angstrom (sigma)

        :param headers: An array of the primary headers from each of the visits

        :return None
        '''

        seeings = []
        durations = []
        dates = []
        visit_index = range(1, len(headers) + 1)

        # Loop through the headers in order to read the seeing, duration, and date for each visit
        for header in headers:
            durations.append(float(header['EXPTIME']))
            seeings.append(1.00)
            dates.append(float(header['MJD-OBS']))

        series1 = pd.Series(durations, index=visit_index)
        series2 = pd.Series(seeings, index=visit_index)
        series3 = pd.Series(dates, index=visit_index)

        # Write the seeings and durations to the dataframe
        cond_dataframe = pd.DataFrame.from_dict(
            {'duration (s)': series1, 'seeing (")': series2, 'date (modified Julian)': series3})

        self.db.insert_dataframe(cond_dataframe, 'cond_table')

        return None

    def fill_visits(self, path, filter, flux, conditions, chunk_size=10):
        '''
        Fills the dataframes that are indexed by visits. It first fills the flux tables and then fills the conditions table.

        :param path: The path to folder with all visit files
        :param flux: True if the flux tables should be filled, False otherwise
        :param conditions: True if the conditions table should be filled, False otherwise
        :param filter: the filter to select
        :param chunk_size: The number of visit files to loop through at a time, default=10

        :return None
        '''

        log.info("Collecting headers and data from the visit files...\n")

        # Collect all the visit files from the directory
        files_set = []
        for root, dirs, files in os.walk(path):

            for name in files:

                if 'bkgd' not in name and '.fits' in name:

                    filename = os.path.join(root, name)

                    # Get filter name

                    try:

                        this_filter = pyfits.getval(filename, "FILTER", ext=0)

                    except:

                        # Could not read FILTER, probably not the image we are looking for

                        continue

                    else:

                        # Add if it matches the filter we are looking for

                        if this_filter.replace(" ", "") == filter:
                            files_set.append(filename)

        # Sort the visit files based on Modified Julian Time
        times_with_visits = []

        for f in files_set:
            times_with_visits.append((pyfits.getheader(f, 0)['MJD-OBS'], f))

        # Sort in tuples works by using the first element in each tuple
        times_with_visits.sort()

        visits_sorted = [files_set for times, files_set in times_with_visits]

        # Collect all the necessary headers from each visit file

        # Loop through in chunks of visits
        for i, chunk in enumerate(chunker(visits_sorted, chunk_size)):

            log.info("Processing chunk %i of %i..." % (i + 1, len(visits_sorted) // chunk_size + 1))

            # Arrays of headers and data. There will be one from each visit in the directory.
            headers_prim = []
            headers_nobkgd = []
            headers_orig = []
            headers_masks = []
            data_masks = []
            data_nobkgd = []
            data_orig = []

            for filename in chunk:
                with pyfits.open(filename) as f:
                    headers_prim.append(f[0].header)
                    headers_nobkgd.append(f[1].header)
                    data_nobkgd.append(f[1].data)
                    data_masks.append(f[2].data)
                    headers_masks.append(f[2].header)
                    headers_orig.append(f[3].header)
                    data_orig.append(f[3].data)

            # Call helper methods to fill in the fluxes and conditions for these visits

            if flux:
                self._fill_flux(headers_nobkgd, data_nobkgd, headers_orig, data_orig, headers_masks, data_masks)

            # Conditions are seeing, date, and so on, for each visit

            if conditions:
                self._fill_cond(headers_prim)

        return None

    def get_flux_list(self):
        '''

        :return:
        '''
        reg = self.db.get_table_as_dataframe('reg_dataframe')
        num_regs = len(reg.index)

        fluxes = {}

        for i in range(1,num_regs+1):

            df = self.db.get_table_as_dataframe('flux_table_%i' % i)

            fluxes_in_this_region = df['flux']

            fluxes["region_%i" % i] = fluxes_in_this_region

        return fluxes


    def close(self):
        '''
        Closes the connection to the _database.

        :return None
        '''

        # Disconnecting the _database writes any uncommitted information to the _database
        self.db.disconnect()

        log.info("Database closed")

        return None
