import pandas as pd
import numpy as np
import os
import re
import logging

import astropy.units as u
import astropy.io.fits as pyfits

from lsst_transients.circular_region import CircularRegion
from astropy import wcs
from astropy.wcs.utils import proj_plane_pixel_scales
from lsst_transients.oversample_image import oversample
from lsst_transients.utils.cartesian_product import cartesian_product
from lsst_transients.utils import database_io
from lsst_transients.utils.chuncker import chunker

logging.basicConfig(format="%(asctime)s %(message)s")
log = logging.getLogger(os.path.basename(__file__))
log.setLevel(logging.DEBUG)


class Data_Database(object):
    '''
    The DataDatabase class is used to create a database and write to different tables: a region table, a flux table for
    each region, and a conditions table.
    '''
    
    
    def __init__(self, dbname, first=True):
        '''
        Initializes the database and opens a connection to the database.
        
        :param dbname: The name for the database
        :param first: If the database has not yet been created, True, False otherwise
        '''

        if first==False:

            assert os.path.exists(dbname), "Database %s does not exist" % dbname

        self.db = database_io.SqliteDatabase(dbname)

        self.db.connect()
        
        
    def fill_reg(self, regfile):
        '''
        Fills the database with the string for each region as seen in DS9.
        
        :param regfile: CircularRegion file created using lsst_grid_generator_shapes.py

        :return The number of regions
        '''
        
        with open(regfile) as f:
            
            # An array of strings where we will store the string of each region as it appears in the DS9 file.
            strs = []
            
            for i, line in enumerate(f):
                
                if i==0:

                    pass
                
                else:

                    #Store the string with the corresponding region in the array after trimming the "\n" off the end.
                    trimmed = line[0:len(line)-1]
                    strs.append(trimmed)
            
            # Fill the index array now that we know how many regions there were
            indices = range(1, i+1)
            
            series = pd.Series(strs, index=indices)
            
            reg_dataframe = pd.DataFrame.from_dict({'ds9_info': series})
            self.db.insert_dataframe(reg_dataframe, 'reg_dataframe')
            
            print(reg_dataframe)
            
            return len(indices)


    def init_flux_tables(self, num_regs):
        '''
        Inserts empty flux tables into the database.

        :param num_regs: Number of regions, and thus the number of flux tables

        :return None
        '''
        # Write the fluxes to the database

        # Insert many tables

        with database_io.bulk_operation(self.db):

            for r in range(0, num_regs):

                if r % 1000 == 0:
                    print("Initialized flux table %i of %i" % (r+1, num_regs))

                flux_dataframe = pd.DataFrame(columns=['flux', 'err'], dtype=float)
                self.db.insert_dataframe(flux_dataframe, 'flux_table_%i' % (r + 1), commit=False)

        return None


    def init_cond_table(self):
        '''
        Inserts an empty conditions table into the database.

        :return: None
        '''

        cond_dataframe = pd.DataFrame(columns=['date (modified Julian)','duration (s)','seeing (")'], dtype=float)
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
        
        distances = np.sqrt((x - x2)**2 + (y - y2)**2)
        
        return distances
        
        
    def _sum_flux(self, w, ds9_string, data, max_coords):
        '''
        Adds the flux from each pixel within the region in order to get the total flux for the region.
        
        :param w: The WCS from the header
        :param ds9_string: The information for a region in a string format 
        :param data: The data array of the visit file to be examined (the fluxes)
        :param max_coords: The maximum pixel coordinates of the image
        
        :return sum_flux: The total flux for this region
        '''

        # Get info about the region from the DS9 string
        split = re.split("[, (\")]+", ds9_string)
        shape = split[0]
        ra = float(split[1])
        dec = float(split[2])
        diameter = float(split[3]) * 2
        
        # Get the diameter in pixels assuming that the diameter is the same in WCS for all regions
        pixel_scale = proj_plane_pixel_scales(w)
        assert np.isclose(pixel_scale[0], pixel_scale[1], rtol=1e-2), "Pixel scale is different between X and Y direction"
        pixel_scale_with_units = pixel_scale[0] * u.Unit(w.wcs.cunit[0])
        pixel_scale_arcsec = pixel_scale_with_units.to("arcsec").value
        
        diameter_pix = np.ceil(diameter/pixel_scale_arcsec)
  
        if shape == 'box':
                        
                rotation = split[5]
     
        # Find the smallest square around the region
        x_and_y = w.wcs_world2pix([[ra, dec]], 0)
        x_and_y = str(x_and_y).replace(']','').replace('[','').split()
        x = np.floor(float(x_and_y[0]))
        y = np.floor(float(x_and_y[1]))

        # Make a region object given the information found previously
        reg = CircularRegion(x, y, diameter_pix, shape)
        
        # Maximum pixel coordinates
        max_coords = str(max_coords[0]).replace('(','').replace(')','').split()
        max_x = max_coords[0].replace(',','')
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
                                  
        in_region = distances <= (diameter_pix/2.0)
        reg_data = data.swapaxes(0,1)[corner1[0]:corner3[0], corner1[1]:corner2[1]]
        reg_data = reg_data.reshape(reg_data.size,)

        # Add the fluxes of the pixels within the bounding box
        sum_flux = float(np.sum(reg_data[in_region]))

        return sum_flux
    
        
    def _get_fluxes(self, reg, data, header):
        '''
        Gets the fluxes for all of the regions in the image.
        
        :param reg: The region dataframe created with lsst_grid_generator_shapes.py
        :param data: The data array of the visit file to be examined
        :param header: The header of the visit file to be examined
        
        :return fluxes: An array of the fluxes for each region in the image
        '''

        # Get the number of regions and use this number to initialize an array that will store the fluxes for each region
        num_regs = len(reg.index)
        fluxes = np.zeros(num_regs)

        # Add the fluxes within each region by calling _sum_flux
        log.info("Measuring flux for each region\n")
        
        w = wcs.WCS(header)
        max_coords = [(header['NAXIS1'], header['NAXIS2'])]       

        for i in range(num_regs):
                
            if (i+1) % 1000 == 0:
                    
                log.info("Processed region %i of %i" %(i+1, num_regs))
                 
            # Call the helper method to get the sum of the flux from all the pixels in the region  
            add = self._sum_flux(w, reg.get_value(i, "ds9_info"), data, max_coords)
            fluxes[i] = add

        return fluxes
    
    
    def _get_flux_errors(self, nobkgd_fluxes, orig_fluxes, mask_data):
        '''
        Computes the error of the flux for each region of the image.
        
        :param nobkgd_fluxes: A list of fluxes for each region for the background subtracted image
        :param orig_fluxes: A list of fluxes for each region for the original image (before subtraction)
        :param mask_data: A list of fluxes for each region for the mask image
        
        :return flux_errors: An array of flux errors for the background subtracted image, one error for each region
        '''
        
        # Find the error on the background
        mask = (np.array(mask_data) == 0)
        bkgd_level = np.median(orig_fluxes[mask])
        bkgd_level_check = np.median(orig_fluxes[mask]-nobkgd_fluxes[mask])
        assert np.isclose(bkgd_level, bkgd_level_check, rtol=0.05), "Background level is too different:" \
                                                                    "\nOriginal with mask: %f\n " \
                                                                    "Original with mask minus background-subtracted: %f" \
                                                                    % (bkgd_level, bkgd_level_check)

        bkgd_errors = np.std(orig_fluxes[mask])

        # Propagate the error to find the flux errors
        flux_errors = np.sqrt(orig_fluxes + bkgd_errors**2)

        return flux_errors

    
    def _get_data(self, dtype, *args, **kwargs):

        this_data = pyfits.getdata(*args, **kwargs)

        return np.array(this_data, dtype=dtype)


    def _fill_flux(self, headers_nobkgd, data_nobkgd, headers_orig, data_orig, headers_masks, data_masks):
        '''
        Fills the dataframe with the flux and the flux error for each region with the indices as the visit number.

        :param headers_nobkgd: An array of the headers from the background-subtracted images from each of the visits
        :param data_nobkgd: An array of the flux data from the background-subtracted images from each of the visits
        :param headers_orig: An array of the headers from the original images from each of the visits
        :param data_orig: An array of the flux data from the original images from each of the visits
        :param headers_masks: An array of the headers from the mask images from each of the visits
        :param data_masks: An array of the flux data from the mask images from each of the visits

        :return None
        '''
    
        # Arrays for all the visits that will store arrays of fluxes and flux errors for each region
        visit_fluxes = []
        visit_errs = []
        
        # Access the regions table in the database in order to find the number of regions
        reg = self.db.get_table_as_dataframe('reg_dataframe')
        num_regs = len(reg.index)

        # Loop through the visit files
        for i in range(0, len(headers_nobkgd)):
            
            print("File %i:\n\nOversampling...\n" % (i+1))
            
            # Oversample the background-subtracted, the original images, and the mask images.
            scale_factor = 2

            scaled_data_nobkgd, scaled_wcs_nobkgd = oversample(data_nobkgd[i], headers_nobkgd[i], scale_factor)

            scaled_data_orig, scaled_wcs_orig = oversample(data_orig[i], headers_orig[i], scale_factor)

            scaled_data_masks, scaled_wcs_masks = oversample(data_masks[i], headers_masks[i], scale_factor)
            
            # Get the fluxes for each region for the scaled images
            print("Scaled background-subtracted image\n")

            fluxes_nobkgd = self._get_fluxes(reg, scaled_data_nobkgd,
                                             scaled_wcs_nobkgd)
            
            print("\nScaled original image\n")
            fluxes_orig = self._get_fluxes(reg, scaled_data_orig,
                                           scaled_wcs_nobkgd)

            print("\nScaled mask image\n")
            fluxes_mask = self._get_fluxes(reg, scaled_data_masks,
                                           scaled_wcs_nobkgd)

            # Normalize the scaled images
            norm_factor = scale_factor * 2
            print("\nnormalization factor (background-subtracted): %f\n" % norm_factor)
            
            fluxes_nobkgd /= norm_factor
            
            fluxes_orig /= norm_factor

            fluxes_mask /= norm_factor
            
            # Use the helper method to find the errors for the flux in each region
            flux_errors_nobkgd = self._get_flux_errors(fluxes_nobkgd, fluxes_orig, fluxes_mask)
            
            # Arrays store the fluxes and flux errors for each region for each visit.
            visit_fluxes.append(fluxes_nobkgd)
            visit_errs.append(flux_errors_nobkgd)
        
        # Switch the order of the arrays in order to store the data
        region_fluxes = np.swapaxes(visit_fluxes, 0, 1)
        region_errs = np.swapaxes(visit_errs, 0, 1)
        
        # Write the fluxes to the database
        print("Writing to database...\n")
        dataframes = []

        for r in range(0,num_regs):
            
            if (r+1) % 100 == 0:
                    
                log.info("Processed region %i of %i" %(r+1, num_regs))

            flux_dataframe = pd.DataFrame(columns=['flux', 'err'])
            flux_dataframe['flux'] = region_fluxes[r]
            flux_dataframe['err'] = region_errs[r]
            dataframes.append(flux_dataframe)

        # Insert many tables that will be committed when the database is disconnected.

        with database_io.bulk_operation(self.db):

            for r in range(0,num_regs):

                if r % 100 == 0:
                    print("Inserted table %i of %i" % (r+1, num_regs))

                self.db.append_dataframe_to_table(dataframes[r], 'flux_table_%i' % (r+1), commit=False)
        
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
        visit_index = range(1, len(headers)+1)
        
        # Loop through the headers in order to read the seeing, duration, and date for each visit
        for header in headers:
            
            durations.append(float(header['EXPTIME']))
            seeings.append(1.00)
            dates.append(float(header['MJD-OBS']))
            
        series1 = pd.Series(durations, index=visit_index)
        series2 = pd.Series(seeings, index=visit_index)
        series3 = pd.Series(dates, index=visit_index)
        
        # Write the seeings and durations to the dataframe
        cond_dataframe = pd.DataFrame.from_dict({'duration (s)': series1, 'seeing (")': series2, 'date (modified Julian)': series3})
        self.db.append_dataframe_to_table(cond_dataframe, 'cond_table')
            
        return None

        
    def fill_visits(self, path, flux, conditions, chunk_size=10):
        '''
        Fills the dataframes that are indexed by visits. It first fills the flux tables and then fills the conditions table.

        :param path: The path to folder with all visit files
        :param flux: True if the flux tables should be filled, False otherwise
        :param conditions: True if the conditions table should be filled, False otherwise
        :param chunk_size: The number of visit files to loop through at a time, default=10

        :return None
        '''

        print("Collecting headers and data from the visit files...\n")

        # Collect all the visit files from the directory
        files_set = []
        for root, dirs, files in os.walk(path):

            for name in files:
                if 'bkgd' not in name and '.fits' in name:
                    files_set.append(os.path.join(root, name))

            for name in dirs:
                if 'bkgd' not in name and '.fits' in name:
                    files_set.append(os.path.join(root, name))

        # Sort the visit files based on Modified Julian Time
        times = []
        for f in files_set:
            times.append(pyfits.getheader(f,0)['MJD-OBS'])

        times_with_visits = zip(times, files_set)
        times_with_visits.sort()

        visits_sorted = [files_set for times, files_set in times_with_visits]

        # Collect all the necessary headers from each visit file

        # Loop through in chunks of visits
        for i, chunk in enumerate(chunker(visits_sorted, chunk_size)):

            print("Processing chunk %i of %i..." % (i+1, len(visits_sorted)//chunk_size+1))

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
            if flux == True:
                self._fill_flux(headers_nobkgd, data_nobkgd, headers_orig, data_orig, headers_masks, data_masks)
            if conditions == True:
                self._fill_cond(headers_prim)
        
        return None
        
        
    def close(self):
        '''
        Closes the connection to the database.

        :return None
        '''

        # Disconnecting the database writes any uncommitted information to the database
        self.db.disconnect()

        print("Closed successfully")
        
        return None
