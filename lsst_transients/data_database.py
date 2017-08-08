import pandas as pd
import os


from utils import database_io
from utils.logging_system import get_logger
from aperture_photometry import bulk_aperture_photometry

log = get_logger(os.path.basename(__file__))


class DataDatabase(object):
    '''
    The DataDatabase class is used to create and handle a _database and write to different tables: a region table,
    a flux table for each region, and a conditions table.
    '''

    def __init__(self, dbname):
        '''
        Initializes the _database and opens a connection to the _database.

        :param dbname: The name for the _database
        '''

        self.db = database_io.HDF5Database(dbname)

    def __enter__(self):

        self.open()

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):

        self.close()

    def open(self):

        self.db.connect()

    def _fill_cond(self, headers, obsids):
        '''
        Fills the dataframe with the conditions for each visit (seeing, duration, and date). Seeing at 5000 angstrom (sigma)

        :param headers: An array of the primary headers from each of the visits

        :return None
        '''

        seeings = []
        durations = []
        dates = []
        # visit_index = range(1, len(headers) + 1)

        # Loop through the headers in order to read the seeing, duration, and date for each visit
        for header in headers:
            durations.append(float(header['EXPTIME']))
            seeings.append(1.00)
            dates.append(float(header['MJD-OBS']))

        series1 = pd.Series(durations, index=obsids)
        series2 = pd.Series(seeings, index=obsids)
        series3 = pd.Series(dates, index=obsids)

        # Write the seeings and durations to the dataframe
        cond_dataframe = pd.DataFrame.from_dict(
            {'duration (s)': series1, 'seeing (")': series2, 'date (modified Julian)': series3})

        self.db.insert_dataframe(cond_dataframe, 'cond_table')

        return None

    def fill_visits(self, path, selected_filter, grid):
        '''
        Fills the dataframes that are indexed by visits. It first fills the flux tables and then fills the conditions table.

        :param path: The path to folder with all visit files
        :param selected_filter: The filter to be used
        :param grid: a Grid object

        :return None
        '''

        df, df_errors, headers_prim, obsid_sorted, n_failed = bulk_aperture_photometry(path, selected_filter, grid)

        # Conditions are seeing, date, and so on, for each visit

        self._fill_cond(headers_prim, obsid_sorted)

        self.db.insert_dataframe(df, "fluxes")
        self.db.insert_dataframe(df_errors, "fluxes_errors")

        log.info("Failed visits: %i" % n_failed)

    def close(self):
        '''
        Closes the connection to the _database.

        :return None
        '''

        # Disconnecting the _database writes any uncommitted information to the _database
        self.db.disconnect()

        log.info("Database closed")
