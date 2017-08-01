import pandas as pd
import numpy as np
import os
import re
import logging

import astropy.units as u
import astropy.io.fits as pyfits

from region import Region
from astropy import wcs
from astropy.wcs.utils import proj_plane_pixel_scales
from oversample_image import oversample
from utils.cartesian_product import cartesian_product
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

        self.db = database_io.SqliteDatabase(dbname)

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


    def init_flux_tables(self, num_regs):
        '''
        Inserts empty flux tables into the database.

        :param num_regs: Number of regions, and thus the number of flux tables

        :return None
        '''
        # Write the fluxes to the database


        # with database_io.bulk_operation(self.db):
        #
        #     # Create first table
        #
        #     flux_dataframe = pd.DataFrame(columns=['flux', 'err'], dtype=float)
        #     self.db.insert_dataframe(flux_dataframe, 'flux_table_1', commit=True)
        #
        #     # Now duplicate it (much faster than just creating new ones)
        #
        #     for r in loop_with_progress(range(1, num_regs), num_regs - 1, 1000, my_printer):
        #
        #         self.db.duplicate_table('flux_table_1', 'flux_table_%i' % (r + 1), commit=False)

        with database_io.bulk_operation(self.db):

            # Insert many tables
            class hack(object):

                def __init__(self, db):

                    self._giant_query = []
                    self._db = db

                def __call__(self, *args, **kwargs):

                    if self._giant_query:
                        log.info("Writing %s tables" % len(self._giant_query))
                        self._db.cursor.executescript(";".join(self._giant_query))
                        self._giant_query = []

                    log.info(*args, **kwargs)


            create_table_cmd = 'CREATE TABLE "flux_table_%i" ("flux" REAL, "err" REAL)'

            my_printer = hack(self.db)

            for r in loop_with_progress(range(num_regs), num_regs, 50000, my_printer):

                my_printer._giant_query.append(create_table_cmd % (r+1))

        self.db.connection.commit()

        return None


    def init_cond_table(self):
        '''
        Inserts an empty conditions table into the database.

        :return: None
        '''

        cond_dataframe = pd.DataFrame(columns=['date (modified Julian)', 'duration (s)', 'seeing (")'], dtype=float)
        self.db.insert_dataframe(cond_dataframe, 'cond_table')

        return None


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


    def _find_pixels_in_region(self, w, ds9_string, max_coords):
        """
        Returns a mask which selects all the pixels within the provided region

        :param w: The WCS from the header
        :param ds9_string: The information for a region in a string format
        :param max_coords: The maximum pixel coordinates of the image

        :return mask array, corners of bounding box

        """

        # Get info about the region from the DS9 string
        split = re.split("[, (\")]+", ds9_string)
        shape = split[0]
        ra = float(split[1])
        dec = float(split[2])
        diameter = float(split[3]) * 2

        # Get the diameter in pixels assuming that the diameter is the same in WCS for all regions
        pixel_scale = proj_plane_pixel_scales(w)
        assert np.isclose(pixel_scale[0], pixel_scale[1],
                          rtol=1e-2), "Pixel scale is different between X and Y direction"
        pixel_scale_with_units = pixel_scale[0] * u.Unit(w.wcs.cunit[0])
        pixel_scale_arcsec = pixel_scale_with_units.to("arcsec").value

        diameter_pix = np.ceil(diameter / pixel_scale_arcsec)

        if shape == 'box':
            rotation = split[5]

        # Find the smallest square around the region
        x_and_y = w.wcs_world2pix([[ra, dec]], 0)
        x_and_y = str(x_and_y).replace(']', '').replace('[', '').split()
        x = np.floor(float(x_and_y[0]))
        y = np.floor(float(x_and_y[1]))

        # Make a region object given the information found previously
        reg = Region(x, y, diameter_pix, shape)

        # Maximum pixel coordinates
        max_coords = str(max_coords[0]).replace('(', '').replace(')', '').split()
        max_x = max_coords[0].replace(',', '')
        max_y = max_coords[1]

        # Find the bounding box around the region
        corner1, corner2, corner3, corner4 = reg.get_boundingbox(max_x, max_y)

        # Get the pixel coordinates of the pixels within the bounding box
        in_bounding_box_x = np.arange(corner1[0], corner3[0])
        in_bounding_box_y = np.arange(corner1[1], corner2[1])
        cart_prod = cartesian_product([in_bounding_box_x, in_bounding_box_y]).T
        all_in_bdb_x = cart_prod[0]
        all_in_bdb_y = cart_prod[1]
        distances = self._get_distances(reg.x, reg.y, all_in_bdb_x, all_in_bdb_y)

        in_region = distances <= (diameter_pix / 2.0)

        return in_region, corner1, corner2, corner3, corner4


    def _get_data_in_region(self, data, in_region, corner1, corner2, corner3, corner4):
        '''

        :param data: Background-subtracted data from the image
        :param in_region: Mask to check whether the data is within the region
        :param corner1: Bottom left corner of the bounding box
        :param corner2: Top left corner of the bounding box
        :param corner3: Bottom right corner of the bounding box
        :param corner4: Top right corner of the bounding box

        :return: The masked data for the region
        '''

        # Select first elements within the bounding box
        reg_data = data[corner1[1]:corner2[1], corner1[0]:corner3[0]]

        # Now flatten the remaining oversampled_bkgsub_image and apply the in_region mask, selecting the pixels within the
        # region
        reg_data = reg_data.flatten()

        return reg_data[in_region]


    def _sum_flux(self, oversampled_bkgsub_image, in_region, corner1, corner2, corner3, corner4,
                  oversampled_counts_image, bkgd_uncertainty):
        '''
        Adds the flux from each pixel within the region in order to get the total flux for the region.

        :param oversampled_bkgsub_image: The data array of the background-subtracted image to be examined (the fluxes)
        :param in_region: Mask to check whether the data is within the region
        :param corner1: Bottom left corner of the bounding box
        :param corner2: Top left corner of the bounding box
        :param corner3: Bottom right corner of the bounding box
        :param corner4: Top right corner of the bounding box
        :param oversampled_counts_image: The data array of the image before background-subtraction
        :param bkgd_uncertainty: The error value for the background

        :return Total flux and flux error for this region
        '''

        # Add the fluxes of the pixels within the bounding box
        reg_data = self._get_data_in_region(oversampled_bkgsub_image,
                                            in_region,
                                            corner1, corner2, corner3, corner4)  # type: np.ndarray
        sum_flux = np.sum(reg_data)

        reg_data_counts = self._get_data_in_region(oversampled_counts_image,
                                                   in_region,
                                                   corner1, corner2, corner3, corner4)  # type: np.ndarray

        flux_error = np.sqrt(np.sum(reg_data_counts) + bkgd_uncertainty ** 2)

        return sum_flux, flux_error


    def _get_fluxes(self, reg, data, orig, header, bkgd_uncertainty):
        '''
        Gets the fluxes for all of the regions in the image.

        :param reg: The region dataframe created with lsst_grid_generator_shapes.py
        :param data: The data array of the background-subtracted image (fluxes)
        :param orig: The data array of the image before background-subtraction
        :param header: The header of the visit file to be examined
        :param bkgd_uncertainty: The error value for the background

        :return fluxes: An array of the fluxes and flux errors for each region in the image
        '''

        # Get the number of regions and use this number to initialize an array that will store the fluxes for each region
        num_regs = len(reg.index)
        fluxes = np.zeros(num_regs)
        fluxes_errors = np.zeros(num_regs)

        # Add the fluxes within each region by calling _sum_flux
        log.info("Measuring flux for each region\n")

        w = wcs.WCS(header)
        max_coords = [(header['NAXIS1'], header['NAXIS2'])]

        for i in range(num_regs):

            # Get the ds9 definition
            ds9_string = reg.get_value(i, "ds9_info")

            in_region, corner1, corner2, corner3, corner4 = self._find_pixels_in_region(w,
                                                                                        ds9_string,
                                                                                        max_coords)

            if (i + 1) % 1000 == 0:
                log.info("Processed region %i of %i" % (i + 1, num_regs))

            # Call the helper method to get the sum of the flux from all the pixels in the region
            flux, error = self._sum_flux(data, in_region, corner1, corner2, corner3, corner4, orig, bkgd_uncertainty)
            fluxes[i] = flux
            fluxes_errors[i] = error

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

        # Access the regions table in the database in order to find the number of regions
        reg = self.db.get_table_as_dataframe('reg_dataframe')
        num_regs = len(reg.index)

        # Loop through the visit files
        n_visits_in_this_chunk = len(headers_nobkgd)

        for i in range(0, n_visits_in_this_chunk):

            log.info("Processing file %i of %i..." % ((i+1), n_visits_in_this_chunk))

            log.debug("oversampling")

            # Oversample the background-subtracted, the original images, and the mask images.
            scale_factor = 2

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
            norm_factor = scale_factor * 2

            fluxes_nobkgd /= norm_factor

            # Arrays store the fluxes and flux errors for each region for each visit.
            visit_fluxes.append(fluxes_nobkgd)
            visit_errs.append(fluxes_errors)

        # Switch the order of the arrays in order to store the data
        region_fluxes = np.swapaxes(visit_fluxes, 0, 1)
        region_errs = np.swapaxes(visit_errs, 0, 1)

        # Write the fluxes to the database
        log.info("Writing to database...\n")
        dataframes = []

        for r in range(0, num_regs):

            if (r + 1) % 100 == 0:
                log.info("Processed region %i of %i" % (r + 1, num_regs))

            flux_dataframe = pd.DataFrame(columns=['flux', 'err'])
            flux_dataframe['flux'] = region_fluxes[r]
            flux_dataframe['err'] = region_errs[r]
            dataframes.append(flux_dataframe)

        # Insert many tables that will be committed when the database is disconnected.

        with database_io.bulk_operation(self.db):

            for r in range(0, num_regs):

                if r % 100 == 0:
                    log.info("Inserted table %i of %i" % (r + 1, num_regs))

                self.db.append_dataframe_to_table(dataframes[r], 'flux_table_%i' % (r + 1), commit=False)

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
        self.db.append_dataframe_to_table(cond_dataframe, 'cond_table')

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

        fluxes = []
        for i in range(1,num_regs+1):
            df = self.db.get_table_as_dataframe('flux_table_%i' % i)
            for j in range(len(df['flux'])):
                fluxes.append(df['flux'][j])
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
