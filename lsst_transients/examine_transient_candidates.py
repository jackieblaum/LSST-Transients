import yaml
import re
import os
import matplotlib.image as mpimg
import numpy as np
import matplotlib.animation as animation

from matplotlib import colors, cm, pyplot as plt
from astropy import wcs

from astropy.io import fits as pyfits
from data_database import DataDatabase
from circular_region import CircularRegion

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


def write_ds9_region_file(region, df, directory):
    '''
    Writes all transient candidate regions to their own DS9 region files.

    :param reg_ids: The IDs of the transient candidate regions as found in the yaml file (ie. reg1)
    :param db: The name of the _database
    :return: None
    '''

    reg_index = int(re.compile(r'(\d+)$').search(region).group(1))+1

    # Get the DS9 region string from the database and write it to a new region file
    with open('%s/%s.reg' % (directory, region), "w+") as f:

        ds9_string = df['ds9_info'][reg_index]

        f.write('icrs\n%s' % ds9_string)

    return ds9_string


def make_lightcurve(region, block_edges, df_fluxes, df_errors, cond, times, directory):


    reg_database_index = re.compile(r'(\d+)$').search(region).group(1)

    # Plot and save

    plt.xlabel('Time (MJD)')
    plt.ylabel('Flux')

    plt.errorbar(times, df_fluxes.iloc[int(reg_database_index)].values,
                 yerr=list(df_errors.iloc[int(reg_database_index)].values), fmt='.')
    plt.title('Region %s Lightcurve' % reg_database_index)

    for j in range(len(block_edges)):
        edge_lines = plt.axvline(block_edges[j], linestyle='--')
        plt.setp(edge_lines, color='r', linewidth=2.0)

    plt.savefig("%s/%s.png" % (directory, region))


def make_image(region_str, multiply, visit):

    # Get needed info from FITS file
    with pyfits.open(visit) as f:
        all_data = f[1].data
        header = f[1].header

    w = wcs.WCS(header)

    reg = CircularRegion.from_ds9_region(w, region_str)
    xcenter = reg.x
    ycenter = reg.y

    # Convert the diameter of the region to pixel coordinates
    diameter_pix = reg.d

    # Calculate the starting and ending coordinates that we want from the FITS file
    xstart = xcenter - (diameter_pix * multiply)
    ystart = ycenter - (diameter_pix * multiply)

    istart = int(ystart)
    jstart = int(xstart)

    xend = xcenter + (diameter_pix * multiply)
    yend = ycenter + (diameter_pix * multiply)

    iend = int(yend)
    jend = int(xend)

    selected_data = all_data[istart:iend, jstart:jend]
    normalized_data = (selected_data - selected_data.min() + 0.01) / (selected_data.max() - selected_data.min() + 0.01)
    norm = colors.LogNorm(0.01, normalized_data.max())
    image = plt.imshow(normalized_data, cmap=cm.gray, norm=norm, origin="lower",
                       animated=True)

    return image


def make_movie(region_str, region, directory, multiply, visits):

    fig = plt.figure()

    images = []

    for visit in visits:
        images.append(make_image(region_str, multiply, visit))

    ani = animation.FuncAnimation(fig, images, interval=50, blit=True, repeat_delay=1000)
    ani.save('%s/%s.mp4' % (directory, region))
    plt.show()


def examine_transient_candidates(database, block_edges_file, multiply, visits, selected_filter):

    db = DataDatabase("%s.db" % database, first=False)

    # Get the region IDs and corresponding block edges of transient candidates from the block edges file
    reg_ids, block_edges = get_regs_and_blk_edges(block_edges_file)

    # Make a directory to store the region files and lightcurves
    directory = 'transient_candidates'
    if not os.path.exists(directory):
        os.makedirs(directory)

    df_regions = db.db.get_table_as_dataframe('reg_dataframe')

    df_fluxes = db.db.get_table_as_dataframe('fluxes')
    df_errors = db.db.get_table_as_dataframe('fluxes_errors')

    # Get the dataframe that corresponds with the region
    df_cond = db.db.get_table_as_dataframe('cond_table')
    times = df_cond['date (modified Julian)']

    visits, obsids = db.find_visits_files(visits, selected_filter)

    for i, region in enumerate(reg_ids):
        # Make a new DS9 region file for each individual transient candidate region
        region_str = write_ds9_region_file(region, df_regions, directory)

        make_lightcurve(region, block_edges[i], df_fluxes, df_errors, df_cond, times, directory)

        make_movie(region_str, region, directory, multiply, visits)

    db.close()
