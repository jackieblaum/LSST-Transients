import sqlite3
import pandas as pd
import sqlalchemy.engine
import astropy.io.fits as pyfits
import numpy as np
import os, sys
import re
import pyregion
import pyregion._region_filter as filter

from lsst_transients.oversample_image import oversample
from lsst_transients.utils.cartesian_product import cartesian_product
from lsst_transients.utils.chuncker import chunker
from lsst_transients.utils.logging_system import get_logger
from lsst_transients.utils.loop_with_progress import loop_with_progress

import logging

logging.basicConfig(format="%(asctime)s %(message)s")
log = logging.getLogger(os.path.basename(__file__))
log.setLevel(logging.DEBUG)


class Database(object):
    '''
    The Database classed is used to create a database and write to three different tables: a region table, a flux table, and a conditions table.
    '''
    
    
    def __init__(self, dbname, first=True):
        '''
        Initializes the database engine and opens a connection to the database.
        
        :param dbname: The name for the database (no need to add .db at the end)
        '''

        if first == False:
            assert os.path.exists(dbname), "Database %s does not exist" % dbname

        self.engine = sqlalchemy.engine.create_engine('sqlite:///%s' % dbname)
        self.connection = self.engine.connect()
        
        
    def fill_reg(self, regfile):
        '''
        #Fills the database with the string for each region as seen in DS9 and regID as the indices.
        
        :param regfile: Region file created using grid_generator_shapes.py
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
            indices = [range(1, i)]
            
            series = pd.Series(strs, index=indices)
            
            reg_dataframe = pd.DataFrame.from_dict({'ds9_info': series})
            reg_dataframe.to_sql("reg_dataframe", self.connection, if_exists='replace', index_label='regID')
            
            print(reg_dataframe)
            
            return None
                
                
    def _sum_flux(self, ds9_string, data, header):
        '''
        Adds the flux from each pixel within the region in order to get the total flux for the region.
        
        :param ds9_string: The information for a region in a string format 
        :param data: The data array of the visit file to be examined
        :param header: The header of the visit file to be examined
        
        :return sum_flux: The total flux for this region
        '''
        
        reg_definition = 'icrs;%s' % ds9_string.replace(" ","")
        reg = pyregion.parse(reg_definition).as_imagecoord(header)
        filt = reg.get_filter()
        
        mask = filt.mask(data)
    
        sum_flux = np.sum(data[mask])
        
        return sum_flux
    
        
    def _get_fluxes(self, reg, data, header):
        '''
        Gets the fluxes for all of the regions in the image.
        
        :param reg: The region dataframe created with grid_generator_shapes.py
        :param data: The data array of the visit file to be examined
        :param header: The header of the visit file to be examined
        
        :return fluxes: An array of the fluxes for each region in the image
        '''

        num_regs = len(reg.index)
        fluxes = np.zeros(num_regs)

        
        log.info("Measuring flux for each region")
            
        for i in range(0, num_regs):
                
            if (i+1) % 10 == 0:
                    
                log.info("Processed region %i of %i" %(i+1, num_regs))
                 
                # Call the helper method to get the sum of the flux from all the pixels in the region    
            fluxes[i] = self._sum_flux(reg.get_value(i, "ds9_info"), data, header)
        
        #sys.exit(-1)
        
        return fluxes
    
        
    def _fill_flux(self, headers_nobkgd, data_nobkgd, headers_orig, data_orig):
        '''
        Fills the dataframe with the flux and the flux error for each region with the indices as the visit number.

        :param headers_nobkgd: An array of the headers from the background-subtracted images from each of the visits
        :param headers_orig: An array of the headers from the original images from each of the visits
        '''
    
        # An array for all the visits that will store arrays of fluxes and flux errors for each region
        visit_fluxes = []
        
        # Access the regions table in the database
        reg = pd.read_sql_table("reg_dataframe", self.connection)
        
        num_regs = len(reg.index)
        
        for i in range(0, len(headers_nobkgd)):

            scale_factor = 2
            
            # Oversample both the background-subtracted and the original images.            
            scaled_data_nobkgd, scaled_wcs_nobkgd = oversample(data_nobkgd[i], headers_nobkgd[i], scale_factor)
            
            scaled_data_orig, scaled_wcs_orig = oversample(data_orig[i], headers_orig[i], scale_factor)
            
            pyfits.writeto("nobkgd%i.fits" % i, scaled_data_nobkgd, header=scaled_wcs_nobkgd, overwrite=True)
            pyfits.writeto("orig%i.fits" % i, scaled_data_orig, header=scaled_wcs_orig, overwrite=True)
            
            
            fluxes_nobkgd = self._get_fluxes(reg, pyfits.getdata("nobkgd" + str(i) + ".fits"), pyfits.getheader("nobkgd%i.fits" % i,0))
            fluxes_orig = self._get_fluxes(reg, pyfits.getdata("orig" + str(i) + ".fits"), pyfits.getheader("orig%i.fits" % i,0))
            
            norm_factor = 2 * scale_factor
            
            fluxes_nobkgd /= norm_factor

            fluxes_orig /= norm_factor
            
            # An array that stores the fluxes and flux errors for each region for each visit.
            visit_fluxes.append(fluxes_nobkgd)

        # Switch the order of the arrays in order to store the data
        region_fluxes = np.swapaxes(visit_fluxes, 0, 1)
        
        cols = []
        
        # Add two columns per region, one for flux and the other for the flux error
        for r in range(0,num_regs):
            
            cols.append('flux%s' % str(r+1))
            
        flux_dataframe = pd.DataFrame(index=range(1,len(headers_nobkgd)+1), columns=cols)
            
        for r in range(0,num_regs):

            flux_dataframe['flux%s' % str(r+1)] = region_fluxes[r]
        print(flux_dataframe)
                

        flux_dataframe.to_sql("flux_table", self.connection, if_exists='replace', index_label='visitID')
        
        return None

        
        
    def _fill_cond(self, headers):
        '''
        Fills the dataframe with the conditions for each visit (seeing, duration, etc.). Seeing at 5000 angstrom (sigma)

        :param headers: An array of the primary headers from each of the visits
        '''
        
        seeings = []
        durations = []
        times = []
        visit_index = range(1, len(headers)+1)
        
        for header in headers:
            
            durations.append(header['EXPTIME'])
            seeings.append(1.0)
            times.append(header['MJD-OBS'])
            
        series1 = pd.Series(durations, index=visit_index)
        series2 = pd.Series(seeings, index=visit_index)
        series3 = pd.Series(times, index=visit_index)
    
        cond_dataframe = pd.DataFrame.from_dict({'duration (s)': series1, 'seeing (")': series2, 'date (modified Julian)': series3})
        cond_dataframe.to_sql("cond_table", self.connection, if_exists='replace', index_label='visitID')
        
        print(cond_dataframe)
            
        return None


        
        
    def fill_visits(self, path, filter, flux, conditions, chunk_size=10):
        '''
        Fills the two dataframes that are indexed by visits. It first fills the flux table and then fills the conditions table.

        :param path: The path to folder with all visit files
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
                self._fill_flux(headers_nobkgd, data_nobkgd, headers_orig, data_orig)

            # Conditions are seeing, date, and so on, for each visit

            if conditions:
                self._fill_cond(headers_prim)

        return None

    def get_flux_list(self):

        flux_table = pd.read_sql_table("flux_table", self.connection)

        fluxes = {}

        for i in range(1, len(flux_table.columns)):

            this_region_fluxes = flux_table['flux%i' % i]

            fluxes["region_%i" % i] = this_region_fluxes

        # for j in range(len(flux_table.index)):
        #
        #     this_region_fluxes = []
        #
        #     for i in range(1,len(flux_table.columns)):
        #
        #         this_region_fluxes.append(flux_table['flux%i'%i][j])
        #
        #     fluxes[] = this_region_fluxes)

        return fluxes
        
        
    def close(self):
        '''
        Closes the connection to the database.
        '''
        
        self.connection.close()
        print("Closed successfully")
        
        return None



