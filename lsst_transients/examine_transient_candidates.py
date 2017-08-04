import yaml
import re
import matplotlib.pyplot as plt


from data_database import DataDatabase

def get_regs_and_blk_edges(yaml_file):

    reg_ids = []
    block_edges = []

    with open(yaml_file, "r+") as f:
        file_contents = yaml.load(f)

        # Iterate over the dictionary with region IDs as keys and block edge arrays as values
        for key, value in file_contents.items():

            reg_ids.append(key)
            block_edges.append(value)

    return reg_ids, block_edges


def write_ds9_region_files(reg_ids, db):
    '''
    Writes all transient candidate regions to their own DS9 region files.

    :param reg_ids: The IDs of the transient candidate regions as found in the yaml file (ie. reg1)
    :param db: The name of the _database
    :return: None
    '''

    df = db.db.get_table_as_dataframe('reg_dataframe')

    # Loop over the regions that contain transient candidates
    for i, region in enumerate(reg_ids):

        # Region indices start at 0, so subtract one from the region ID
        reg_index = int(re.compile(r'(\d+)$').search(region).group(1)) - 1

        # Get the DS9 region string from the _database and write it to a new region file
        with open('%s.reg' % region, "w+") as f:

            ds9_string = df['ds9_info'][reg_index]

            f.write('icrs\n%s' % ds9_string)

def make_lightcurves(reg_ids, block_edges, db):

    # Loop over the transient candidate regions
    for i in range(len(reg_ids)):

        reg_database_index = re.compile(r'(\d+)$').search(reg_ids[i]).group(1)

        # Get the dataframe that corresponds with the region
        df = db.db.get_table_as_dataframe('flux_table_%s' % reg_database_index)
        cond = db.db.get_table_as_dataframe('cond_table')
        times = cond['date (modified Julian)']

        # Plot and save

        plt.xlabel('Time (MJD)')
        plt.ylabel('Flux')

        plt.errorbar(times, df['flux'].values, yerr=list(df['err'].values), fmt='.')
        plt.title('Region %i Lightcurve' % i)

        plt.yscale("symlog")

        for j in range(len(block_edges[i])):
            edge_lines = plt.axvline(block_edges[i].item(j), linestyle='--')
            plt.setp(edge_lines, color='r', linewidth=2.0)

        plt.savefig("%s.png" % reg_ids[i])
        plt.show()



def examine_transient_candidates(database, block_edges_file):

    db = DataDatabase("%s.db" % database, first=False)

    # Get the region IDs and corresponding block edges of transient candidates from the block edges file
    reg_ids, block_edges = get_regs_and_blk_edges(block_edges_file)

    # Make a new DS9 region file for each individual transient candidate region
    write_ds9_region_files(reg_ids, db)

    make_lightcurves(reg_ids, block_edges, db)

    db.close()
