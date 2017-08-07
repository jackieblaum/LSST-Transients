import json
import yaml
import re
import os
import pyregion

import matplotlib
matplotlib.use("Agg")

import matplotlib.image as mpimg
import numpy as np
import matplotlib.animation as animation

from matplotlib import colors, cm, pyplot as plt
from astropy import wcs
from astropy.nddata.utils import Cutout2D
from astropy import units as u

from astropy.io import fits as pyfits
from data_database import DataDatabase
from circular_region import CircularRegion

def get_blk_edges(regid, json_file):

    with open(json_file, "r+") as f:
        file_contents = yaml.load(f)

        # Iterate over the dictionary with region IDs as keys and block edge arrays as values
        for key, value in file_contents.items():
            if key == regid:
                block_edges = value
                return block_edges

    return None


def write_ds9_region_file(region, df, directory):
    '''
    Writes all transient candidate regions to their own DS9 region files.

    :param reg_ids: The IDs of the transient candidate regions as found in the json file (ie. reg1)
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


def make_image(center, size, visit):

    # Get needed info from FITS file
    with pyfits.open(visit) as f:
        all_data = f[1].data
        header = f[1].header

    w = wcs.WCS(header)

    #filt = region.as_imagecoord(header).get_filter()
    #mask = filt.mask(all_data)
    selected_data = Cutout2D(all_data, center, size, wcs=w)

    normalized_data = (selected_data.data - selected_data.data.min() + 0.1) / (selected_data.data.max() - selected_data.data.min() + 0.1)
    norm = colors.LogNorm(normalized_data.min(), normalized_data.max())
    fig = plt.figure()
    sub = fig.add_subplot(111, projection=w)

    image = sub.imshow(normalized_data, cmap=cm.gray, norm=norm, origin="lower",
                       animated=True)

    return fig, sub, image


def make_movie(region_str, region_name, directory, multiply, visits):

    fig = plt.figure()

    print("Making animation...")

    images = []

    # Change the region to a square
    split = re.split("[, (\")]+", region_str)
    ra = float(split[1])
    dec = float(split[2])
    side = float(split[3])*2*multiply
    size = u.Quantity((side, side), u.arcsec)
    center = (ra,dec)

    for visit in visits:
        this_figure, _, image = make_image(center, size, visit)
        if image:
            images.append([image])
            #plt.close(this_figure)

    ani = animation.ArtistAnimation(fig, images, interval=500, blit=True, repeat_delay=1000)

    ani.save('%s/%s.mp4' % (directory, region_name))


def examine_transient_candidates(database, regid, block_edges_file, multiply, visits, selected_filter):

    db = DataDatabase("%s.db" % database, first=False)

    # Get the corresponding block edges of the specified transient candidate
    block_edges = get_blk_edges(regid, block_edges_file)

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

    # Make a new DS9 region file
    region_str = write_ds9_region_file(str(regid), df_regions, directory)

    make_lightcurve(str(regid), block_edges, df_fluxes, df_errors, df_cond, times, directory)

    make_movie(region_str, str(regid), directory, multiply, visits)

    db.close()
