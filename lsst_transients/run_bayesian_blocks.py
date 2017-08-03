import matplotlib.pyplot as plt
import numpy as np
import yaml

from astropy.stats import bayesian_blocks
from lsst_transients.data_database import DataDatabase
from lsst_transients.utils.loop_with_progress import loop_with_progress
from lsst_transients.utils.logging_system import get_logger

log = get_logger(__name__)


class BayesianBlocks(object):
    '''
    Runs the Bayesian Block algorithm on the data.
    '''

    def __init__(self, dbname, filename, min_blocks):

        self.database = DataDatabase("%s.db" % dbname)

        # Get the number of regions
        reg = self.database.db.get_table_as_dataframe('reg_dataframe')
        self.num_regs = len(reg.index)
        self.file = '%s%s' % (filename, '.yml')
        self.min_blocks = min_blocks


    def _sort_visits(self):
        '''
        Sort the visits by modified Julian time.

        :return: None
        '''

        # Sort the visits in the conditions table based on time
        df = self.database.db.get_table_as_dataframe('cond_table')
        indices = df['date (modified Julian)'].argsort()
        df = df.iloc[indices]

        # Replace the data frame in the database with the sorted version
        self.database.db.remove_table('cond_table')
        self.database.db.insert_dataframe(df, 'cond_table')
        print(df)

        # Use the same ordering for all of the flux tables

        # Loop through the flux tables and sort the visits
        for i in range(1, self.num_regs+1):

            table_name = 'flux_table_%i' % i
            df = self.database.db.get_table_as_dataframe(table_name)
            df = df.iloc[indices]

            self.database.db.remove_table(table_name)
            self.database.db.insert_dataframe(df, table_name)
            print(df)

        return None


    def _append_edges(self, d):
        '''

        :param d:
        :return:
        '''

        with open(self.file, "a+") as f:
            yaml.dump(d, f)


    def run_algorithm(self, sort=True):
        '''
        Sort the visits if needed, then run the Bayesian Block algorithm on the data

        :param sort: True if the visits still need to be sorted by time, false otherwise

        :return: None
        '''

        if sort:

            print("Sorting the visits by time...")
            self._sort_visits()

        # Get the array of times for the visits from the database
        conditions = self.database.db.get_table_as_dataframe('cond_table')

        conditions['date (modified Julian)'] = np.array(conditions['date (modified Julian)'], float)
        times = conditions['date (modified Julian)'].values

        # Loop over all the regions and run Bayesian_Blocks on each of them
        print("Running the Bayesian Block algorithm on each region...")

        # Erase any existing file with the same name, create new file
        with open(self.file, "w") as f:
            pass

        df = self.database.db.get_table_as_dataframe('fluxes')
        df_errors = self.database.db.get_table_as_dataframe('fluxes_errors')

        for i, regid in loop_with_progress(df.index, len(df.index), len(df.index) // 50,
                                           log.info, with_enumerate=True):

            # flux_tables.append(df)

            edges = bayesian_blocks(t=times, x=df.loc[regid].values, sigma=df_errors.loc[regid].values,
                                    fitness='measures', p0=1e-3)
            d = {}

            if len(edges) >= self.min_blocks:
                d['reg%i' % i] = edges.tolist()
                self._append_edges(d)

        self.database.close()

        return None