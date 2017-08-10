import json
import re
import os
import numpy as np
import shutil
import multiprocessing
import functools

import matplotlib
matplotlib.use("Agg")

import reproject

import imageio

from matplotlib import colors, cm, pyplot as plt
from matplotlib.patches import Circle
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.nddata.utils import Cutout2D
from astropy import units as u

from astropy.io import fits as pyfits
from astropy.wcs.utils import proj_plane_pixel_scales
from data_database import DataDatabase
from aperture_photometry import find_visits_files
from utils.logging_system import get_logger
from utils.loop_with_progress import loop_with_progress

log = get_logger(__name__)


def make_lightcurve(region_id, block_dict, df_fluxes, df_errors, times, directory):


    #reg_database_index = re.compile(r'(\d+)$').search(region_id).group(1)

    # Plot and save

    plt.xlabel('Time (MJD)')
    plt.ylabel('Flux')

    plt.errorbar(times, df_fluxes.loc[region_id].values,
                 yerr=list(df_errors.loc[region_id].values), fmt='.')
    plt.title('Region %s Lightcurve' % region_id)

    edges = block_dict['tstart']
    edges.append(block_dict['tstop'][-1])

    for j in range(len(edges)):

        plt.axvline(edges[j], linestyle='--', color='r', lw=2.0)

    plt.ylim([0, np.nanmax(df_fluxes.loc[region_id].values)])

    plt.savefig("%s/%s.png" % (directory, region_id))

    return np.nanmax(df_fluxes.loc[region_id].values)


def reproject_onto_original_wcs(center, size, visit, original_wcs):

    # Reproject new image onto WCS from original image

    with pyfits.open(visit) as f:

        # Compute background level
        original_image = f[3].data
        mask_data = f[2].data
        mask = (np.array(mask_data) == 0)

        bkgd_level = np.median(original_image[mask])

        this_hdu = f[3]

        this_w = wcs.WCS(this_hdu.header)

        selected_data = Cutout2D(this_hdu.data, center, size, wcs=this_w)

        # Get the MJD
        mjd = f[0].header['MJD-OBS']

    # Reproject into the original WCS
    reprojected_data, _ = reproject.reproject_interp((selected_data.data, selected_data.wcs), original_wcs,
                                                     shape_out=selected_data.shape)

    return reprojected_data, bkgd_level, mjd


def make_plot(reprojected_data, bkg_level, filename, ra, dec, radius, orig_wcs, max_flux, autonorm=True):

    if autonorm:

        # Circle area
        pixel_scale = proj_plane_pixel_scales(orig_wcs)

        # We take the geometric average of the pixel scales in the X and Y direction, as done
        # in pyregion

        pixel_scale_with_units = (pixel_scale[0] * pixel_scale[1]) ** 0.5 * u.Unit(orig_wcs.wcs.cunit[0])

        pixel_scale_deg = pixel_scale_with_units.to("deg").value
        radius_pixel = radius / pixel_scale_deg
        area = np.pi * radius_pixel**2

        norm = colors.LogNorm(vmin=bkg_level, vmax=(max_flux/area + bkg_level)*1.3)

    else:

        norm = colors.LogNorm(vmin=np.nanmin(reprojected_data), vmax=np.nanmax(reprojected_data))

    idx = np.isnan(reprojected_data)
    reprojected_data[idx] = bkg_level

    fig = plt.figure()
    sub = fig.add_subplot(111, projection=orig_wcs)

    sub.imshow(reprojected_data, norm=norm, origin="lower", cmap='jet', interpolation='none')

    c = Circle((ra, dec), radius, edgecolor='white', facecolor='none',
               transform=sub.get_transform('fk5'), linewidth=2.0)

    sub.add_patch(c)

    fig.savefig(filename)

    plt.close(fig)


def worker(args, center, size, wcs, temp_dir, ra, dec, radius, max_flux):

    i, visit = args

    this_data, bkg_level, this_mjd = reproject_onto_original_wcs(center, size, visit, wcs)

    this_filename = os.path.join(temp_dir, "frame%010i.png" % i)

    make_plot(this_data, bkg_level, this_filename, ra, dec, radius, wcs, max_flux)

    return this_filename, this_data, bkg_level, this_mjd


def make_movie(region_str, region_name, directory, multiply, visits, max_flux, block_edges):

    # Change the region to a square
    split = re.split("[, (\")]+", region_str)
    ra = float(split[1])
    dec = float(split[2])
    radius = float(split[3])
    side = radius*2*multiply
    size = u.Quantity((side, side), u.arcsec)
    center = SkyCoord(ra=ra, dec=dec, frame='fk5', unit='deg')

    log.info("Generating animation...")

    # Get WCS from first visit
    with pyfits.open(visits[0]) as f:

        original_header = f[1].header
        original_data = f[1].data
        original_wcs = wcs.WCS(original_header)

    selected_data = Cutout2D(original_data, center, size, wcs=original_wcs)

    temp_dir = "__frame__"

    try:

        shutil.rmtree(temp_dir)

    except:

        pass

    os.makedirs(temp_dir)

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

    worker_partial = functools.partial(worker,
                                       center=center, size=size, wcs=selected_data.wcs, temp_dir=temp_dir,
                                       ra=ra, dec=dec, radius=(radius * u.arcsec).to(u.deg).value,
                                       max_flux=max_flux)

    frames = []
    bkg_levels = np.zeros(len(visits))
    mjds = np.zeros(len(visits))

    with imageio.get_writer(os.path.join(directory, "%s.gif" % region_name), mode='I', duration=0.3) as writer:

        try:

            for i, (this_file, this_data, this_bkg, this_mjd) in loop_with_progress(pool.imap(worker_partial,
                                                                                    zip(range(len(visits)), visits)),
                                                                                    len(visits), 10, log.info,
                                                                                    with_enumerate=True):

                image = imageio.imread(this_file)
                writer.append_data(image)

                frames.append(this_data)
                mjds[i] = this_mjd
                bkg_levels[i] = this_bkg

        except:

            raise

        finally:

            pool.close()

    # Now coadd the frames to generate the animation of the blocks

    with imageio.get_writer(os.path.join(directory, "%s_coadd.gif" % region_name), mode='I', duration=1.0) as writer:

        for t1, t2 in zip(block_edges['tstart'], block_edges['tstop']):

            # Find frames in this block
            # NOTE that the bayesian blocks are half-points in-between bins except the first and the last, that's
            # why we need the >= and the <= on both sides
            idx = (mjds >= t1) & (mjds <= t2)

            assert np.sum(idx) > 0

            # Co-add them
            elements = idx.nonzero()[0]

            print elements

            this_coadd = frames[elements[0]]
            this_bkg_level = bkg_levels[elements[0]]

            for el in elements[1:]:

                this_coadd += frames[el]
                this_bkg_level += bkg_levels[el]

            this_file = os.path.join(temp_dir, "%s_%s.png" % (t1, t2))

            make_plot(this_coadd, this_bkg_level, this_file, ra, dec, radius, selected_data.wcs, max_flux, False)

            image = imageio.imread(this_file)
            writer.append_data(image)

    shutil.rmtree(temp_dir)

    log.info("done")


def examine_transient_candidates(database, grid, regid, block_edges_file, multiply, visits, selected_filter):

    # Get the corresponding block edges of the specified transient candidate

    with open(block_edges_file, "r") as f:

        file_contents = json.load(f)

    block_edges = file_contents[regid]

    # Make a directory to store the region files and lightcurves
    directory = 'transient_candidates'

    if not os.path.exists(directory):

        os.makedirs(directory)

    with DataDatabase("%s" % database) as db:

        df_fluxes = db.db.get_table_as_dataframe('fluxes')
        df_errors = db.db.get_table_as_dataframe('fluxes_errors')

        # Get the dataframe that corresponds with the region
        df_cond = db.db.get_table_as_dataframe('cond_table')

    times = df_cond['date (modified Julian)']

    # Note: this list is time ordered (find_visits_files orders by time)

    visits, obsids = find_visits_files(visits, selected_filter)

    # Make a new DS9 region file
    #region_str = write_ds9_region_file(str(regid), grid, directory)

    region_str = grid.get_ds9_region(regid)

    max_flux = make_lightcurve(str(regid), block_edges, df_fluxes, df_errors, times, directory)

    log.info("Maximum net flux: %.3f" % max_flux)

    make_movie(region_str, str(regid), directory, multiply, visits, max_flux, block_edges)
