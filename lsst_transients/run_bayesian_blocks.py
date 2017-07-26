import matplotlib.pyplot as plt

from astropy.stats import bayesian_blocks
from lsst_transients.create_db_notpyregion import Data_Database

class BayesianBlocks(object):

    def __init__(self, dbname):
        self.database = Data_Database("%s.db" % dbname)

        # Get the number of regions
        reg = self.database.db.get_table_as_dataframe('reg_dataframe')
        self.num_regs = len(reg.index)

    def _sort_visits(self):

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

    def run_algorithm(self, sort):

        if sort == "True":
            print("Sorting the visits by time...")
            self._sort_visits()

        # Get the array of times for the visits from the database
        conditions = self.database.db.get_table_as_dataframe('cond_table')
        times = conditions['date (modified Julian)'].values

        # Get the array of flux dataframes from the database
        print("Getting the fluxes from each data frame...")
        flux_tables = []
        for i in range(1, self.num_regs+1):
            if i%100 == 0:
                print("Processed table %i of %i" % (i, self.num_regs))
            flux_tables.append(self.database.db.get_table_as_dataframe('flux_table_%i' % i))

        # We will plot later
        fig, subs = plt.subplots(1,self.num_regs,figsize=(75,10))
        # Loop over all the regions and run Bayesian_Blocks on each of them
        print("Running the Bayesian Block algorithm on each region...")
        for i in range(1, self.num_regs+1):
            edges = bayesian_blocks(t=times, x=flux_tables[i-1]['flux'].values, sigma=flux_tables[i-1]['err'].values, fitness='measures')
            print("Completed region %i of %i" % (i, self.num_regs))
            print edges

            # Plot and save
            subs[i-1].errorbar(times, flux_tables[i-1]['flux'].values, yerr=flux_tables[i-1]['err'].values, fmt='.')
            subs[i-1].set_title('reg%i' % i)
            for j in range(0,len(edges)):
                subs[i-1].axvline(edges.item(j), linestyle='--')

        fig.savefig("plot.png")

        self.database.close()