import pandas as pd
import json
import os
import matplotlib.pyplot as plt


class Results(object):

    def __init__(self):

        self._candidates = {}

    def add_candidate(self, region_id, edges, fluxes, errors):

        self._candidates[region_id] = {'tstart': list(edges[:-1]),
                                       'tstop': list(edges[1:]),
                                       'flux': list(fluxes),
                                       'error': list(errors)}
    @property
    def region_ids(self):

        return self._candidates.keys()

    def has(self, region_id):

        return region_id in self._candidates

    def get_candidate(self, region_id):

        assert self.has(region_id), "Region %s not known" % region_id

        return pd.DataFrame.from_dict(self._candidates[region_id])

    def write_to(self, filename):

        with open(filename, "w+") as f:

            json.dump(self._candidates, f)

    @classmethod
    def read_from(cls, filename):

        instance = cls()

        assert os.path.exists(filename), "Filename %s does not exist" % filename

        with open(filename, "r") as f:

            instance._candidates = json.load(f)

        return instance

    def plot_lightcurve(self, region_id):

        assert self.has(region_id), "Region %s not known" % region_id

        fig, sub = plt.subplots(1,1)

        candidate = self.get_candidate(region_id)

        tc = (candidate['tstart'] + candidate['tstop']) / 2.0
        dt = (candidate['tstop'] - candidate['tstart']) / 2.0

        sub.errorbar(tc, candidate['flux'], xerr=dt, yerr=candidate['error'], fmt='.')

        return fig

