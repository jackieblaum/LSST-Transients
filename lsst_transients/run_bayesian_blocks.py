from astropy.stats import bayesian_blocks
from lsst_transients.create_db_notpyregion import Data_Database

class BayesianBlocks(object):

    def __init__(self, dbname):
        self.db = Data_Database("%s.db" % dbname)

        # Get the number of regions
        reg = self.db.get_table_as_dataframe('reg_dataframe')
        self.num_regs = len(reg.index)

    def _sort_visits(self):

        # Sort the visits in the conditions table based on time
        df = self.db.get_table_as_dataframe('cond_table')
        indices = df['date (modified Julian)'].argsort()
        df = df.iloc[indices]

        # Replace the data frame in the database with the sorted version
        self.db.remove_table('cond_table')
        self.db.insert_dataframe(df, 'cond_table')

        # Use the same ordering for all of the flux tables

        # Loop through the flux tables and sort the visits
        for i in range(1, self.num_regs+1):

            table_name = 'flux_table_%i' % i
            df = self.db.get_table_as_dataframe(table_name)
            df = df.iloc[indices]

            self.db.remove_table(table_name)
            self.db.insert_dataframe(df, table_name)

    def run_algorithm(self, sort):

        if sort == "True":
            self.sort_visits()

        # Loop over all the regions and run Bayesian_Blocks on each of them
        for i in range(1, self.num_regs+1):

            df = self.db.get_table_as_dataframe('flux_table_%i' % i)
            edges = bayesian_blocks(x=df['flux'].values, sigma=df['err'].values, fitness='measures')
            print edges

        self.db.close()