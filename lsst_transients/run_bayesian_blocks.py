import numpy as np
import pandas as pd
import numexpr
import json
import multiprocessing

#from astropy.stats import bayesian_blocks
from lsst_transients.data_database import DataDatabase
from lsst_transients.utils.loop_with_progress import loop_with_progress
from lsst_transients.utils.logging_system import get_logger
from lsst_transients.aperture_photometry import PHOTOMETRY_NULL_VALUE
from lsst_transients.results import Results

log = get_logger(__name__)

###############################
# numexpr optimization
##############################

# Speed tricks

# HACK: susbstitute the getArguments of numexpr to be faster in our particular case,
# where we always supply the local dictionary (see the bayesian block loop)
def getArguments(names, local_dict, *args, **kwargs):
    """Get the arguments based on the names."""

    return map(local_dict.__getitem__, names)
# Monkey patch numexpr
numexpr.necompiler.getArguments = getArguments

# Set numexpr precision to low (more than enough for us), which is
# faster than high
oldaccuracy = numexpr.set_vml_accuracy_mode('low')
numexpr.set_num_threads(1)
numexpr.set_vml_num_threads(1)

###############################
# end of numexpr optimization
##############################


def simulate(n_points, n_sims, p0, min_blocks):

    n_false_positives = 0

    for i in loop_with_progress(range(n_sims), n_sims, int(np.ceil(n_sims/100)), log.info):

        x = np.random.poisson(500, size=n_points) - 500

        t = np.arange(100, 10000, n_points)

        sigma = np.sqrt(x+500)

        edges = bayesian_blocks(t, x, sigma, p0)

        if len(edges) >= min_blocks:

            n_false_positives += len(edges) - min_blocks + 1

    return n_false_positives / float(n_sims)


def bayesian_blocks(t, x, sigma, p0):

    """Fit the Bayesian Blocks model given the specified fitness function.

            Parameters
            ----------
            t : array_like
                data times (one dimensional, length N)
            x : array_like (optional)
                data values
            sigma : array_like or float (optional)
                data errors

            Returns
            -------
            edges : ndarray
                array containing the (M+1) edges defining the M optimal bins
            """

    # compute values needed for computation, below
    s2 = sigma**2
    ak_raw = np.ones_like(x) / s2
    bk_raw = x / s2

    # create length-(N + 1) array of cell edges
    edges = np.concatenate([t[:1],
                            0.5 * (t[1:] + t[:-1]),
                            t[-1:]])

    # arrays to store the best configuration
    N = len(t)
    best = np.zeros(N, dtype=float)
    last = np.zeros(N, dtype=int)

    ncp_prior = 4 - np.log(73.53 * p0 * (N ** -0.478))

    # Speed tricks: resolve once for all the functions which will be used
    # in the loop
    numexpr_evaluate = numexpr.evaluate
    numexpr_re_evaluate = numexpr.re_evaluate

    cumsum = np.cumsum

    # ----------------------------------------------------------------
    # Start with first data cell; add one cell at each iteration
    # ----------------------------------------------------------------
    for R in range(N):

        # a_k: eq. 31
        a_k = 0.5 * cumsum(ak_raw[:R + 1][::-1])[::-1]

        # b_k: eq. 32
        b_k = - cumsum(bk_raw[:R + 1][::-1])[::-1]

        # evaluate fitness function
        if R == 0:

            fit_vec = numexpr_evaluate("(b_k * b_k) / (4 * a_k)",
                                       optimization='aggressive', local_dict={'b_k': b_k, 'a_k': a_k})

        else:

            fit_vec = numexpr_re_evaluate(local_dict={'b_k': b_k, 'a_k': a_k})

        A_R = fit_vec - ncp_prior
        A_R[1:] += best[:R]

        i_max = np.argmax(A_R)
        last[R] = i_max
        best[R] = A_R[i_max]

    # ----------------------------------------------------------------
    # Now find changepoints by iteratively peeling off the last block
    # ----------------------------------------------------------------
    change_points = np.zeros(N, dtype=int)
    i_cp = N
    ind = N
    while True:
        i_cp -= 1
        change_points[i_cp] = ind
        if ind == 0:
            break
        ind = last[ind - 1]
    change_points = change_points[i_cp:]

    final_blocks = edges[change_points]

    # Return also flux and errors
    n_intervals = len(final_blocks) - 1
    fluxes = np.zeros(n_intervals)
    errors = np.zeros_like(fluxes)

    for i in range(n_intervals):

        tstart = final_blocks[i]
        tstop = final_blocks[i+1]

        idx = (t >= tstart) & (t < tstop)

        n = float(np.sum(idx))

        fluxes[i] = np.sum(x[idx]) / n
        errors[i] = np.sqrt(np.sum(sigma[idx]**2)) / n

    return final_blocks, fluxes, errors


class BayesianBlocks(object):
    '''
    Runs the Bayesian Block algorithm on the data.
    '''

    def __init__(self, dbname, filename, min_blocks):

        self._database = dbname

        self._out_filename = filename

        self._minimum_number_of_blocks = min_blocks

    def run_algorithm(self, p0):
        '''
        Sort the visits if needed, then run the Bayesian Block algorithm on the data

        :param p0: Null hypothesis for bayesian blocks

        :return: None
        '''

        with DataDatabase(self._database) as data:

            # Get the array of times for the visits from the _database
            conditions = data.db.get_table_as_dataframe('cond_table')

            conditions['date (modified Julian)'] = np.array(conditions['date (modified Julian)'], float)
            times = conditions['date (modified Julian)'].values

            df = data.db.get_table_as_dataframe('fluxes')
            df_errors = data.db.get_table_as_dataframe('fluxes_errors')

        # Analyze in parallel
        def feeder():

            for regid in df.index:

                yield times, df.loc[regid].values, df_errors.loc[regid].values, p0, regid

        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

        results = Results()

        # Loop over all the regions and run Bayesian_Blocks on each of them
        print("Running the Bayesian Block algorithm on each region...")

        try:

            for i, (regid, edges, fluxes, errors) in loop_with_progress(pool.imap_unordered(worker, feeder(),
                                                                                            chunksize=1000),
                                                                        len(df.index),
                                                                        1000 * multiprocessing.cpu_count(),
                                                                        log.info, with_enumerate=True):

                if len(edges) >= self._minimum_number_of_blocks:

                    results.add_candidate(regid, edges, fluxes, errors)

        except:

            raise

        finally:

            pool.close()

        log.info("Writing results...")

        # Write file
        results.write_to(self._out_filename)


def worker(args):

    t, x, sigma, p0, regid = args

    x = np.array(x, dtype=float)
    sigma = np.array(sigma, dtype=float)
    t = np.array(t, dtype=float)

    idx = (x != PHOTOMETRY_NULL_VALUE)

    edges, fluxes, errors = bayesian_blocks(t=t[idx], x=x[idx], sigma=sigma[idx], p0=p0)

    return regid, edges, fluxes, errors