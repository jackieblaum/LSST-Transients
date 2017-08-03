import pandas as pd
import numpy as np
import os
import re
import photutils
import ctypes

import astropy.units as u
import astropy.io.fits as pyfits
from astropy.coordinates import SkyCoord
from astropy import wcs

from oversample_image import oversample
from circular_region import CircularRegion
from utils import database_io
from utils.chuncker import chunker
from utils.logging_system import get_logger
from utils.loop_with_progress import loop_with_progress

log = get_logger(os.path.basename(__file__))


class DataDatabase(object):
    '''
    The DataDatabase class is used to create and handle a database and write to different tables: a region table,
    a flux table for each region, and a conditions table.
    '''

    def __init__(self, dbname, first=True):
        '''
        Initializes the database and opens a connection to the database.

        :param dbname: The name for the database
        :param first: If the database has not yet been created, True, False otherwise
        '''

        if first == False:
            assert os.path.exists(dbname), "Database %s does not exist" % dbname

        #self.db = database_io.SqliteDatabase(dbname)
        self.db = database_io.HDF5Database(dbname)

        self.db.connect()

    def fill_reg(self, reg_list_wcs, shape, angular_distance_arcsec, rotation_angle):
        '''
        Fills the database with the string for each region as seen in DS9.

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

    def _get_fluxes(self, reg, data, orig, header, bkgd_uncertainty, apertures):
        '''
        Gets the fluxes for all of the regions in the image.

        :param reg: The region dataframe created with lsst_grid_generator_shapes.py
        :param data: The data array of the background-subtracted image (fluxes)
        :param orig: The data array of the image before background-subtraction
        :param header: The header of the visit file to be examined
        :param bkgd_uncertainty: The error value for the background

        :return fluxes: An array of the fluxes and flux errors for each region in the image
        '''

        # Add the fluxes within each region by calling _sum_flux
        log.info("Measuring flux for each region...")

        # We need the WCS

        w = wcs.WCS(header)

        fluxes_table = photutils.aperture_photometry(data, apertures, wcs=w)

        log.info("done")

        fluxes = fluxes_table['aperture_sum'].data

        log.info("Computing errors...")

        counts_sum = photutils.aperture_photometry(orig, apertures, wcs=w)

        fluxes_errors = np.sqrt(counts_sum['aperture_sum'].data + bkgd_uncertainty**2)

        log.info("done")

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


    def _get_region_fluxes(self, header_nobkgd, data_nobkgd, data_orig, data_masks, apertures):
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
        background_uncertainty = self._get_background_uncertainty(data_nobkgd, data_orig, data_masks)

        # Add the fluxes within each region by calling _sum_flux
        log.info("Measuring flux for each region...")

        # We need the WCS

        w = wcs.WCS(header_nobkgd)

        apertures_pix = apertures.to_pixel(w)

        fluxes_table = photutils.aperture_photometry(data_nobkgd, apertures_pix, error=np.sqrt(data_orig),
                                                     method='exact')

        log.info("done")

        fluxes = fluxes_table['aperture_sum'].data

        log.info("Computing errors...")

        # Adding the background uncertainty to the photometric uncertainty

        fluxes_errors = np.sqrt(fluxes_table['aperture_sum_err'].data ** 2 + background_uncertainty ** 2)

        log.info("done")

        log.info("Median value for the regions: %.3f" % np.nanmedian(fluxes))
        log.info("Median error for the regions: %.3f" % np.nanmedian(fluxes_errors))

        return fluxes, fluxes_errors

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

    def _find_visits_files(self, path, selected_filter):
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

        visits_sorted = map(lambda x:x[1], times_with_visits[:3])
        obsid_sorted = map(lambda x:x[2], times_with_visits[:3])

        return visits_sorted, obsid_sorted

    def _get_apertures(self):

        # Access the regions table in the database in order to find the number of regions
        reg = self.db.get_table_as_dataframe('reg_dataframe')
        num_regs = len(reg.index)

        # Make the aperture regions
        radiuses = np.zeros(num_regs)
        ras = np.zeros(num_regs)
        decs = np.zeros(num_regs)

        for i in range(num_regs):
            ds9_string = reg.get_value(i + 1, "ds9_info")

            split = re.split("[, (\")]+", ds9_string)
            shape = split[0]

            assert shape == "circle", "Only circles are supported at this point"

            ra = float(split[1])
            dec = float(split[2])
            radius = float(split[3])

            ras[i] = ra
            decs[i] = dec
            radiuses[i] = radius

        # NOTE: we assume all apertures have the same radius
        r = np.unique(radiuses)
        assert r.shape[0] == 1, "Regions have different radiuses!"

        apertures = photutils.SkyCircularAperture(SkyCoord(ra=ras * u.deg, dec=decs * u.deg, frame='icrs'),
                                                  r[0] * u.arcsec)

        return apertures, len(ras)

    def fill_visits(self, path, selected_filter, flux, conditions, chunk_size=10):
        '''
        Fills the dataframes that are indexed by visits. It first fills the flux tables and then fills the conditions table.

        :param path: The path to folder with all visit files
        :param flux: True if the flux tables should be filled, False otherwise
        :param conditions: True if the conditions table should be filled, False otherwise
        :param selected_filter: the filter to select
        :param chunk_size: The number of visit files to loop through at a time, default=10

        :return None
        '''

        log.info("Selecting visits...")

        visits_sorted, obsid_sorted = self._find_visits_files(path, selected_filter)

        log.info("Found %i visits" % len(visits_sorted))

        log.info("Setting up regions...")
        all_apertures, n_regions = self._get_apertures()

        log.info("Found %i regions" % n_regions)

        # Loop through the visits

        df = pd.DataFrame(columns=obsid_sorted)
        df_errors = pd.DataFrame(columns=obsid_sorted)

        headers_prim = []

        for i, (visit_file, obsid) in loop_with_progress(zip(visits_sorted, obsid_sorted),
                                                         len(visits_sorted), 1, log.info, with_enumerate=True):

            log.info("Processing visit %i of %i..." % (i + 1, len(visits_sorted)))

            with pyfits.open(visit_file) as f:

                headers_prim.append(f[0].header)

                header_nobkgd = f[1].header
                data_nobkgd = f[1].data

                data_orig = f[3].data
                data_mask = f[2].data

            # Call helper methods to fill in the fluxes and conditions for these visits

            fluxes, fluxes_errors = self._get_region_fluxes(header_nobkgd, data_nobkgd, data_orig, data_mask,
                                                            all_apertures)

            df[obsid] = fluxes
            df_errors[obsid] = fluxes_errors
            # Conditions are seeing, date, and so on, for each visit

        self._fill_cond(headers_prim)

        self.db.insert_dataframe(df, "fluxes")
        self.db.insert_dataframe(df_errors, "fluxes_errors")

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
        Closes the connection to the database.

        :return None
        '''

        # Disconnecting the database writes any uncommitted information to the database
        self.db.disconnect()

        log.info("Database closed")

        return None
