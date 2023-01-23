import numpy as np
import pandas as pd
from astropy.table import Table
from sn_tools.sn_calcFast import CalcSN
from sn_fit.sn_utils import Selection


class Fit_LC(Selection):
    """
    class to fit simulated light curves

    Parameters
    ---------------
    model: str, opt
      model to use for the fit (default: None)
    version: float, opt
       model version (default: -1.0)
    display: bool, opt
      to display the fit online (default: False)
    bands: str, opt
      bands to consider (default: ugrizy)
    snrmin: float, opt
      min SNR for considered LC points for the fit (default: 5)

    """

    def __init__(self, model=None, version=-1.0,
                 display=False, bands='ugrizy', snrmin=1., **kwargs):
        super().__init__(snrmin)

        self.display = display
        self.bands = bands

        # name parameters
        self.parNames = dict(zip(['z', 't0', 'x0', 'x1', 'c'], [
                             'z', 't0', 'x0', 'x1', 'color']))

    def __call__(self, lc):
        """
        call method: this is where the fit is effectively performed
        as we are here with input from a fast simulator
        there are no real fit but rather a transformation of
        Fisher matrix elements to sigmas'

        Parameters
        ---------------
        lc: astropy table
           lc points

        Returns
        -----------

        """

        select = self.select(lc)

        if select is None:
            return Table()

        sn = CalcSN(select, nBef=0, nAft=0,
                    nPhamin=0, nPhamax=0,
                    params=['x0', 'x1', 'daymax', 'color'])

        # Make a dict of the fitted result (plus metadata)
        meta = lc.meta
        resa = self._transform(meta, sn.sn, sn.fitstatus)

        """
        resb = self._get_infos(
            select.meta['z'], select.meta['daymax'], select)

        resa.update(resb)
        """
        output = Table(rows=[list(resa.values())], names=list(resa.keys()))

        return output

    def _transform(self, meta, sn, fitstatus):
        """
        Method to transform input data to a coherent dictionary

        Parameters
        ---------------
        meta: dict
          LC metadata
        par_names: list(str)
          parameter list
        params: list(float)
          parameter values
        vpar_names: list(str)
          covariance parameter list
        covmat: matrix
           covariance matrix
        mbfit: float
          mbfit value
        fitstatus: str
           status of the fit

        Returns
        ----------
        dict with the following keys
        - SN parameters(from metadata):
           Dec, Ra, SNID, color, dL, daymax, epsilon_color, epsilon_daymax,
           epsilon_x0, epsilon_x1,
           index_hdf5, pixDec, pixID, pixRa, season, survey_area, x0, x1, z
        - result of the fit:
          z_fit, t0_fit, x0_fit, x1_fit, color_fit, Cov_t0t0, Cov_t0x0,
          Cov_t0x1, Cov_t0color, Cov_x0x0,
         Cov_x0x1, Cov_x0color, Cov_x1x1, Cov_x1color, Cov_colorcolor,
         mbfit, fitstatus
        """
        snres = Table(sn)

        snres['Cov_daymaxdaymax'].name = 'Cov_t0t0'
        res = {}
        for vv in ['z_fit', 't0_fit', 'x0_fit', 'x1_fit',
                   'color_fit', 'mbfit', 'Cov_t0x0', 'Cov_t0x1',
                   'Cov_t0color', 'Cov_x0x1', 'Cov_x0color',  'Cov_x1color']:
            res[vv] = -99.
        for vv in ['Cov_t0t0', 'Cov_x0x0', 'Cov_x1x1',  'Cov_colorcolor']:
            res[vv] = snres[vv].data.item()

        for key, value in meta.items():
            res[key] = value

        res['fitstatus'] = fitstatus
        return res

    def _get_infos(self, z, T0, lc):
        """
        Method to estimate LC features

        Parameters
        ---------------
        z: float
          redshift value
        T0:  float
          daymax value
        lc: astropy table
            LC points

        Returns
        -----------
        dict with the following keys:
        phase_min, phase_max, N_bef, N_aft, N_bef_u, N_aft_u,
        SNR_u, N_bef_g, N_aft_g, SNR_g, N_bef_r, N_aft_r, SNR_r,
        N_bef_i, N_aft_i, SNR_i, N_bef_z, N_aft_z, SNR_z, N_bef_y,
        N_aft_y, SNR_y

        """

        lcdf = lc.to_pandas()

        lcdf['phase'] = (lcdf['time']-T0)/(1.+z)

        lcdf['N_bef'] = lcdf['phase'] <= 0
        lcdf['N_aft'] = lcdf['phase'] > 0

        resdf = lcdf.groupby(['band']).apply(
            lambda x: self.calc(x)).reset_index()

        resfi = pd.DataFrame([lcdf['phase'].min()], columns=['phase_min'])
        resfi['phase_max'] = lcdf['phase'].max()
        resfi['N_bef'] = lcdf['N_bef'].sum()
        resfi['N_aft'] = lcdf['N_aft'].sum()

        for b in self.bands:
            idx = resdf['band'] == 'LSST::'+b
            sel = resdf[idx]
            if len(sel) > 0:
                for var in ['N_bef', 'N_aft', 'SNR']:
                    resfi['{}_{}'.format(var, b)] = sel[var].values
            else:
                for var in ['N_bef', 'N_aft', 'SNR']:
                    resfi['{}_{}'.format(var, b)] = 0

        return resfi.to_dict(orient='records')[0]

    def SNR(self, lc):
        """
        Method to estimate SNR of lc

        Parameters
        ---------------
        lc: pandas df
           data to process

        Returns
        -----------
        Signal-to-Noise Ratio (SNR) defined by
        SNR = sqrt(sum(flux/fluxerr)**2)

        """
        idx = lc['flux']/lc['fluxerr_photo'] > 0.
        sel = lc[idx]

        return np.sqrt(np.sum((sel['flux']/sel['fluxerr_photo'])**2))

    def calc(self, grp):
        """
        Method to estimate a set of infos on groups

        Parameters
        --------------
        grp: pandas df group
          data to process

        Returns
        ----------

        """
        return pd.DataFrame({'N_bef': [grp['N_bef'].sum()],
                             'N_aft': [grp['N_aft'].sum()],
                             'SNR': [self.SNR(grp)]})
