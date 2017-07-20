import sqlite3
import pandas as pd
import sqlalchemy.engine
import astropy.io.fits as pyfits
import numpy as np
import os
import re

from oversample_image import oversample

class Database(object):
    '''
    The DataDatabase classed is used to create a database and write to three different tables: a region table, a flux table, and a conditions table.
    '''
    
    
    def __init__(self, dbname):
        '''
        Initializes the database engine and opens a connection to the database.
        
        :param dbname: The name for the database (no need to add .db at the end)
        '''
        
        self.engine = sqlalchemy.engine.create_engine('sqlite:///%s.db' % dbname)
        self.connection = self.engine.connect()
        
        
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
            indices = [range(1, i)]
            
            series = pd.Series(strs, index=indices)
            
            reg_dataframe = pd.DataFrame.from_dict({'ds9_info': series})
            reg_dataframe.to_sql("reg_table", self.connection, if_exists='replace', index_label='regID')
            
            print(reg_dataframe)
            
            return None
                
                
    def _sum_flux(self, x, y, radius, shape):
        '''
        Adds the flux from each pixel within the region in order to get the total flux for the region.
        
        :param x: The x-coordinate for the center of this region
        :param y: The y-coordinate for the center of this region
        :param radius: The radius or side length of this region
        :param shape: The shape of the region, either box or circle
        
        :return sum_flux: The total flux for this region
        TODO: Get the flux from each pixel
        '''
        
    
        sum_flux = 0
        
        return sum_flux
    
        
    def _get_fluxes(self, regfile, header):
        '''
        Gets the fluxes for all of the regions in the image.
        
        :param regfile: The name of the region file created with lsst_grid_generator_shapes.py
        :param header: The header of the visit file to be examined
        
        :return fluxes: An array of the fluxes for each region in the image
        '''
        
        fluxes = []
        
        with open(regfile) as f:
            
            # A counter for the number of lines in the file. It also allows us to skip the first line of the region file, which does not specify a region.
            i = 0
            
            for line in f:
                
                if i==0:
                    
                    # Do nothing
                    pass
                
                else:
                    
                    # Get the position, radius/side length, and shape of the region by splitting each line of the file
                    
                    # Lines have one of the following formats:
                    # box(x, y, radius", radius", rotation)
                    # circle(x, y, radius")
                    
                    split = re.split('(|, |"|)', line)
                    shape = split[0]
                    x = split[1]
                    y = split[2]
                    radius = split[3]
                    
                    if shape == 'box':
                        
                        rotation = split[5]
                        
                    # Call the helper method to get the sum of the flux from all the pixels in the region    
                    fluxes.append(_sum_flux(x, y, radius, shape))
                    
                # Move on to the next region
                i += 1
                
        return fluxes
    
    
    def _get_flux_errors(self, regfile, nobkgd_fluxes, orig_fluxes):
        '''
        Computes the error of the flux for each region of the background-subtracted image.
        
        :param regfile: The name of the region file created with lsst_grid_generator_shapes.py
        :param nobkgd_fluxes: A list of fluxes for each region for the background subtracted image
        :param orig_fluxes: A list of fluxes for each region for the original image (before subtraction)
        
        :return flux_errors: An array of flux errors for the background subtracted image, one error for each region
        '''
        # This array will store the flux error for each region of the background-subtracted image
        flux_errors = []
        
        for i in len(nobkgd_fluxes):
            
            # The error for the flux of each region is given by the square root of the flux from the original image, minus the background, which is given by the difference of the orignal flux and the background-subtracted flux.
            flux_errors[i] = np.sqrt(orig_fluxes[i]) - (orig_fluxes[i] - nobkgd_fluxes[i])
        
        return flux_errors
    
        
    def _fill_flux(self, indices, regfile, headers_nobkgd, headers_orig):
        '''
        Fills the dataframe with the flux and the flux error for each region with the indices as the visit number.
        
        :param indices: An array of the indices for the dataframe, based on the number of visits
        :param headers_nobkgd: An array of the headers from the background-subtracted images from each of the visits
        :param headers_orig: An array of the headers from the original images from each of the visits
        '''
    
        # An array for all the visits that will store arrays of fluxes and flux errors for each region
        visit_info = []
        
        for i in len(headers_nobkgd):
            
            
            # Oversample both the background-subtracted and the original images.
            # TODO: Fill in inputs for oversample
            
            scaled_data_nobkgd, scaled_wcs_nobkgd = oversample(imageData, headers_nobkgd[i], scaleFactor, threshold)
            
            scaled_data_orig, scaled_wcs_orig = oversample(imageData, headers_orig[i], scaleFactor, threshold)

            
            fluxes_nobkgd = _get_fluxes(regfile, scaled_wcs_nobkgd)
            fluxes_orig = _get_fluxes(regfile, scaled_wcs_orig)
        
            flux_errors_nobkgd = _get_flux_errors(regfile, fluxes_nobkgd, fluxes_orig)
        
            # An array that stores the fluxes and flux errors for each region for each visit.
            visit_fluxes[i] = fluxes_nobkgd
            visit_errs[i] = flux_errors_nobkgd
        
        flux_dataframe = pd.Dataframe(index=[range(1,len(headers_nobkgd))])
        
        # Check to see how many regions we have (subtract 1 to get rid of the first line, which is not a region)
        with open(regfile) as f:
            regions = sum(1 for line in f) - 1
        
        # Switch the order of the arrays in order to store the data
        region_fluxes = []
        region_errs = []
        
        for r in regions:    
            for v in len(visit_fluxes):
                
                region_fluxes[r][v] = visit_fluxes[v][r]
                region_errs[r][v] = visit_errs[v][r]
        
        # Add two columns per region, one for flux and the other for the flux error
        for r in regions:
            
            flux_dataframe['flux' + str(r)] = pd.Series(region_fluxes[r], index=flux_dataframe.index)
            flux_dataframe['err' + str(r)] = pd.Series(region_errs[r], index=flux_dataframe.index)
                

        flux_dataframe.to_sql("flux_table", self.connection, if_exists='replace', index_label='visitID')
        
        print(flux_dataframe)
        
        return None

        
        
    def _fill_cond(self, indices, headers):
        '''
        Fills the dataframe with the conditions for each visit (seeing, duration, etc.). Seeing at 5000 angstrom (sigma)
        
        :param indices: An array of the indices for the dataframe, based on the number of visits
        :param headers: An array of the primary headers from each of the visits
        '''
        
        seeings = []
        durations = []
        visit_index = range(1, len(headers)+1)
        
        for header in headers:
            
            durations.append(header['EXPTIME'])
            seeings.append(1)
            
        series1 = pd.Series(durations, index=visit_index)
        series2 = pd.Series(seeings, index=visit_index)
    
        cond_dataframe = pd.DataFrame.from_dict({'duration (s)': series1, 'seeing (")': series2})
        cond_dataframe.to_sql("cond_table", self.connection, if_exists='replace', index_label='visitID')
        
        print(cond_dataframe)
            
        return None


        
        
    def fill_visits(self, regfile, path):
        '''
        Fills the two dataframes that are indexed by visits. It first fills the flux table and then fills the conditions table.
        
        :param regfile: The .reg file that stores the regions created by lsst_grid_generator_shapes.py
        :param path: The path to folder with all visit files
        '''
        
        # An array of headers. There will be one from each visit in the directory.
        headers_prim = []
        headers_nobkgd = []
        headers_orig = []
        
        for f in os.listdir(path):
            
            headers_prim.append(pyfits.getheader(path + f, 0))
            headers_nobkgd.append(pyfits.getheader(path + f, 1))
            headers_orig.append(pyfits.getheader(path + f, 3))
        
        # Fill the index array now that we know how many regions there were
        indices = [range(1, len(headers_prim))]
                         
        #self._fill_flux(indices, regfile, headers_nobkgd, headers_orig)
        self._fill_cond(indices, headers_prim)
        
        return None
        
        
    def close(self):
        '''
        Closes the connection to the database.
        '''
        
        self.connection.close()
        print("Closed successfully")
        
        return None



