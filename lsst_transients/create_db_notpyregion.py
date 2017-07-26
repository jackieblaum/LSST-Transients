import pandas as pd
import astropy.io.fits as pyfits
import numpy as np
import os
import re
import astropy.units as u

from region import Region
from astropy import wcs
from astropy.wcs.utils import proj_plane_pixel_scales
from oversample_image import oversample
from utils.cartesian_product import cartesian_product
from utils import database_io

import logging

logging.basicConfig(format="%(asctime)s %(message)s")
log = logging.getLogger(os.path.basename(__file__))
log.setLevel(logging.DEBUG)


class Data_Database(object):
    '''
    The DataDatabase classed is used to create a database and write to different tables: a region table, a flux table for each region, and a conditions table.
    '''
    
    
    def __init__(self, dbname):
        '''
        Initializes the database engine and opens a connection to the database.
        
        :param dbname: The name for the database (add .db at the end)
        '''

        assert os.path.exists(dbname), "Database %s does not exist" % dbname

        self.db = database_io.SqliteDatabase(dbname)

        self.db.connect()
        
        
    def fill_reg(self, regfile):
        '''
        #Fills the database with the string for each region as seen in DS9 and regID as the indices.
        
        :param regfile: Region file created using lsst_grid_generator_shapes.py
        '''
        
        with open(regfile) as f:
            
            # A counter for the number of lines in the file. It also allows us to skip the first line of the region file, which does not specify a region.
            i = 0
            
            # An array of strings where we will store the string of each region as it appears in the DS9 file.
            strs = []
            
            for line in f:
                
                if i==0:
                    pass
                
                else:
                    #Store the string with the corresponding region in the array after trimming the "\n" off the end.
                    trimmed = line[0:len(line)-1]
                    strs.append(trimmed)
                    
                # Move on to the next region
                i += 1
            
            # Fill the index array now that we know how many regions there were
            indices = range(1, i)
            
            series = pd.Series(strs, index=indices)
            
            reg_dataframe = pd.DataFrame.from_dict({'ds9_info': series})
            self.db.insert_dataframe(reg_dataframe, 'reg_dataframe', commit=False)
            
            print(reg_dataframe)
            
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
        
        reg = Region(x, y, diameter_pix, shape)
        
        # Maximum pixel coordinates
        max_coords = str(max_coords[0]).replace('(','').replace(')','').split()
        max_x = max_coords[0].replace(',','')
        max_y = max_coords[1]
        
        corner1, corner2, corner3, corner4 = reg.get_boundingbox(max_x, max_y)

        #in_bounding_box_x, in_bounding_box_y = np.where(data_wcs[corner1[0]:corner3[0], corner1[1]:corner2[1]])
        in_bounding_box_x = np.arange(corner1[0], corner3[0])
        in_bounding_box_y = np.arange(corner1[1], corner2[1])
        cart_prod = cartesian_product([in_bounding_box_x, in_bounding_box_y]).T
        all_in_bdb_x = cart_prod[0]
        all_in_bdb_y = cart_prod[1]
        distances = self._get_distances(reg.x, reg.y, all_in_bdb_x, all_in_bdb_y)
                                  
        in_region = distances <= (diameter_pix/2.0)
        reg_data = data.swapaxes(0,1)[corner1[0]:corner3[0], corner1[1]:corner2[1]]
        reg_data = reg_data.reshape(reg_data.size,)

        sum_flux = np.sum(reg_data[in_region])
        
        return sum_flux
    
        
    def _get_fluxes(self, reg, data, header):
        '''
        Gets the fluxes for all of the regions in the image.
        
        :param reg: The region dataframe created with lsst_grid_generator_shapes.py
        :param data: The data array of the visit file to be examined
        :param header: The header of the visit file to be examined
        
        :return fluxes: An array of the fluxes for each region in the image
        '''
        
        fluxes = []
        num_regs = len(reg.index)
        
        log.info("Measuring flux for each region\n")
        
        w = wcs.WCS(header)
        max_coords = [(header['NAXIS1'], header['NAXIS2'])]       

        for i in range(0, num_regs):
                
            if (i+1) % 1000 == 0:
                    
                log.info("Processed region %i of %i" %(i+1, num_regs))
                 
            # Call the helper method to get the sum of the flux from all the pixels in the region  
            add = self._sum_flux(w, reg.get_value(i, "ds9_info"), data, max_coords)
            fluxes.append(add)
        #sys.exit(-1)
        return fluxes
    
    
    def _get_flux_errors(self, nobkgd_fluxes, orig_fluxes, mask_data):
        '''
        Computes the error of the flux for each region of the image.
        
        :param nobkgd_fluxes: A list of fluxes for each region for the background subtracted image
        :param orig_fluxes: A list of fluxes for each region for the original image (before subtraction)
        
        :return flux_errors: An array of flux errors for the background subtracted image, one error for each region
        '''
        
        # The error for the flux of each region is given by the square root of the flux from the original image
        mask = (np.array(mask_data) == 0)
        bkgd_level = np.median(orig_fluxes[mask])
        bkgd_level_check = np.median(orig_fluxes[mask]-nobkgd_fluxes[mask])
        assert np.isclose(bkgd_level, bkgd_level_check, rtol=0.05), "Background level is too different:" \
                                                                    "\nOriginal with mask: %f\n " \
                                                                    "Original with mask minus background-subtracted: %f" \
                                                                    % (bkgd_level, bkgd_level_check)

        bkgd_errors = np.std(orig_fluxes[mask])

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
        '''
    
        # An array for all the visits that will store arrays of fluxes and flux errors for each region
        visit_fluxes = []
        visit_errs = []
        
        # Access the regions table in the database in order to find the number of regions
        reg = self.db.get_table_as_dataframe('reg_dataframe')

        num_regs = len(reg.index)
        
        # Loop through the visit files
        for i in range(0, len(headers_nobkgd)):
            
            print("File %i:\n\nOversampling...\n" % i)
            
            # Oversample both the background-subtracted and the original images.
            scale_factor = 2
            scaled_data_nobkgd, scaled_wcs_nobkgd = oversample(data_nobkgd[i], headers_nobkgd[i], scale_factor)
            
            scaled_data_orig, scaled_wcs_orig = oversample(data_orig[i], headers_orig[i], scale_factor)

            scaled_data_masks, scaled_wcs_masks = oversample(data_masks[i], headers_masks[i], scale_factor)
            
            # Write the oversampled data to a FITS file
            pyfits.writeto("nobkgd%i.fits" % i, scaled_data_nobkgd, header=scaled_wcs_nobkgd, clobber=True)
            pyfits.writeto("orig%i.fits" % i, scaled_data_orig, header=scaled_wcs_orig, clobber=True)
            pyfits.writeto("mask%i.fits" % i, scaled_data_masks, header=scaled_wcs_masks, clobber=True)
            
            # Get the fluxes for each region for the scaled images
            print("Scaled background-subtracted image\n")

            fluxes_nobkgd = self._get_fluxes(reg, self._get_data(float, "nobkgd%i.fits" % i),
                                             pyfits.getheader("nobkgd%i.fits" % i,0))
            
            print("\nScaled original image\n")
            fluxes_orig = self._get_fluxes(reg, self._get_data(float, "orig%i.fits" % i),
                                           pyfits.getheader("orig%i.fits" % i,0))

            print("\nScaled mask image\n")
            fluxes_mask = self._get_fluxes(reg, self._get_data(float, "mask%i.fits" % i),
                                           pyfits.getheader("mask%i.fits" % i,0))

            # Normalize the scaled images
            norm_factor = scale_factor * 2
            print("\nnormalization factor (background-subtracted): %f" % norm_factor)
            
            fluxes_nobkgd /= norm_factor
            
            fluxes_orig /= norm_factor


            fluxes_mask /= norm_factor
            
            # Use the helper method to find the errors for the flux in each region
            flux_errors_nobkgd = self._get_flux_errors(fluxes_nobkgd, fluxes_orig, fluxes_mask)
            
            # An array that stores the fluxes and flux errors for each region for each visit.
            visit_fluxes.append(fluxes_nobkgd)
            visit_errs.append(flux_errors_nobkgd)
        
        # Switch the order of the arrays in order to store the data
        region_fluxes = np.swapaxes(visit_fluxes, 0, 1)
        region_errs = np.swapaxes(visit_errs, 0, 1)
        
        # Write the fluxes to the database
        print("Writing to database...\n")
        dataframes = []
        
        #import pdb;pdb.set_trace()
        
        for r in range(0,num_regs):
            
            if (r+1) % 100 == 0:
                    
                log.info("Processed region %i of %i" %(r+1, num_regs))
                
            
            flux_dataframe = pd.DataFrame(index=range(1,len(headers_nobkgd)+1), columns=['flux', 'err'])
            flux_dataframe['flux'] = region_fluxes[r]
            flux_dataframe['err'] = region_errs[r]
            dataframes.append(flux_dataframe)

        # Insert many tables

        with database_io.bulk_operation(self.db):

            for r in range(0,num_regs):

                if r % 100 == 0:
                    print("Inserted table %i of %i" % (r+1, num_regs))

                self.db.insert_dataframe(dataframes[r], 'flux_table_%i' % (r+1), commit=False)
        
        return None
        
        
    def _fill_cond(self, headers):
        '''
        Fills the dataframe with the conditions for each visit (seeing, duration, etc.). Seeing at 5000 angstrom (sigma)

        :param headers: An array of the primary headers from each of the visits
        '''
        
        seeings = []
        durations = []
        dates = []
        visit_index = range(1, len(headers)+1)
        
        # Loop through the headers in order to read the seeing and duration for each visit
        for header in headers:
            
            durations.append(header['EXPTIME'])
            seeings.append(1)
            dates.append(header['MJD-OBS'])
            
        series1 = pd.Series(durations, index=visit_index)
        series2 = pd.Series(seeings, index=visit_index)
        series3 = pd.Series(dates, index=visit_index)
        
        # Write the seeings and durations to the dataframe
        cond_dataframe = pd.DataFrame.from_dict({'duration (s)': series1, 'seeing (")': series2,
                                                 'date (modified Julian)': series3})
        self.db.insert_dataframe(cond_dataframe, 'cond_table')
        
        print(cond_dataframe)
            
        return None

        
    def fill_visits(self, path, flux, conditions):
        '''
        Fills the two dataframes that are indexed by visits. It first fills the flux table and then fills the conditions table.

        :param path: The path to folder with all visit files
        '''

        print("Collecting headers and data from the visit files...\n")
        # Arrays of headers and data. There will be one from each visit in the directory.
        headers_prim = []
        headers_nobkgd = []
        headers_orig = []
        headers_masks = []
        data_masks = []
        data_nobkgd = []
        data_orig = []

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
        for f in visits_sorted:

            headers_prim.append(pyfits.getheader(f, 0))
            headers_nobkgd.append(pyfits.getheader(f, 1))
            data_nobkgd.append(pyfits.getdata(f, 1))
            data_masks.append(pyfits.getdata(f, 2))
            headers_masks.append(pyfits.getheader(f, 2))
            headers_orig.append(pyfits.getheader(f, 3))
            data_orig.append(pyfits.getdata(f, 3))
        
        # Call helper methods to fill in the fluxes and conditions for these visits
        if flux == "True":
            self._fill_flux(headers_nobkgd, data_nobkgd, headers_orig, data_orig, headers_masks, data_masks)
        if conditions == "True":
            self._fill_cond(headers_prim)
        
        return None
        
        
    def close(self):
        '''
        Closes the connection to the database.
        '''

        self.db.disconnect()

        print("Closed successfully")
        
        return None
