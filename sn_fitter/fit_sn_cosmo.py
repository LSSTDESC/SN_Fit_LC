import sncosmo
import numpy as np
from astropy import (cosmology, units as u, constants as const)
import pandas as pd
from astropy.table import Table


class Fit_LC:
    """
    class to fit simulated light curves

    Parameters
    ---------------
    model: str, opt
      model to use for the fit (default: salt2-extended)
    version: float, opt
       model version (default: 1.0)
    telescope: sn_tools.telescope,opt
      telescope model (default: None)
    display: bool, opt
      to display the fit online (default: False)
    bands: str, opt
      bands to consider (default: ugrizy)
    snrmin: float, opt
      min SNR for considered LC points for the fit (default: 5)
    nbef: int, opt
      minimal number of LC points before max (default: 4)
    naft: int, opt
      minimal number of LC points after max (default: 5) 
    
    """

    def __init__(self, model='salt2-extended', version=1.0, telescope=None, display=False, bands='ugrizy',snrmin=5.,nbef=4,naft=5):

        self.display = display
        self.bands = bands
        self.snrmin = snrmin
        self.nbef = nbef
        self.naft = naft
        
        # get the bands for sncosmo registration - with and without airmass

        for band in bands:
            if telescope.airmass > 0:
                band = sncosmo.Bandpass(
                    telescope.atmosphere[band].wavelen, telescope.atmosphere[band].sb, name='LSST::'+band, wave_unit=u.nm)
            else:
                band = sncosmo.Bandpass(
                    telescope.system[band].wavelen, telescope.system[band].sb, name='LSST::'+band, wave_unit=u.nm)
            sncosmo.registry.register(band, force=True)

        # get the source
        source = sncosmo.get_source(model, version=str(version))
        # get the dust
        dust = sncosmo.OD94Dust()

        # sn_fit_model instance
        self.SN_fit_model = sncosmo.Model(source=source)

        # name parameters
        self.parNames = dict(zip(['z', 't0', 'x0', 'x1', 'c'], [
                             'z', 't0', 'x0', 'x1', 'color']))

    def __call__(self, lc):
        """
        call method: this is where the fit is effectively performed

        Parameters
        ---------------
        lc: astropy table
           lc points

        Returns
        -----------

        """

        if lc is None:
            return None, None

        res_param_names = ['z', 't0', 'x0', 'x1', 'c']
        res_params_values = np.zeros((5,), dtype=float)
        vparam_names = ['t0', 'x0', 'x1', 'c']
        covariance = np.zeros((4, 4,), dtype=float)
        mbfit = -1
        z = -1
        fitstatus = 'nofit'
        meta = None
        if lc is not None:
            # get metadata
            meta = lc.meta

            # set redshift for the fit
            z = meta['z']
            self.SN_fit_model.set(z=z)

            # LC points selection: make sure enough points with minimal snr

            idx = lc['flux'] > 0.
            idx &= lc['fluxerr'] > 0.

            """
            selecta = lc[idx]
            idx = selecta['flux']/selecta['fluxerr'] >= snr_min
            selecta = selecta[idx][['flux', 'fluxerr',
                                    'band', 'zp', 'zpsys', 'time']]
            """
            
            
            select = lc[idx]
            select['snr'] = select['flux']/select['fluxerr_photo']
            idx = select['snr'] >=self.snrmin
            select = select[idx]
            select['diff_time'] = select.meta['daymax']-select['time']
            nlc_bef = len(select[select['diff_time']>=0])
            nlc_aft = len(select[select['diff_time']<0])
            # check the total number of LC points here
            assert((nlc_bef+nlc_aft)==len(select))
            #select = select[['flux', 'fluxerr', 'band', 'zp', 'zpsys', 'time']]
            #print('Selection', len(selecta))
            #if len(select) >= 5:
            
            if nlc_bef >= self.nbef and nlc_aft >= self.naft:
                try:
                    # fit here
                    res, fitted_model = sncosmo.fit_lc(select, self.SN_fit_model, [
                        't0', 'x0', 'x1', 'c'], bounds={'z': (z-0.01, z+0.01)}, minsnr=0.0)
                    # get parameters
                    if res['success']:
                        mbfit = fitted_model._source.peakmag(
                            'bessellb', 'vega')
                        res_param_names = res['param_names']
                        res_params_values = res['parameters']
                        vparam_names = res['vparam_names']
                        covariance = res['covariance']
                        fitstatus = 'fitok'
                    else:
                        fitstatus = 'badfit'
                except (RuntimeError, TypeError, NameError):
                    fitstatus = 'crash'
            else:
                fitstatus = 'nodat'

       # Make a dict of the fitted result (plus metadata)
        resa = self._transform(meta, res_param_names, list(
            res_params_values), vparam_names, covariance, mbfit, fitstatus)

        resb = self._get_infos(
            z, res_params_values[res_param_names.index('t0')], select)

        if self.display and len(select) >= 5:
            print(select['band','time','flux','fluxerr'])
            import pylab as plt
            sncosmo.plot_lc(select, model=fitted_model,
                            color='r', pulls=False, errors=res.errors)
            plt.show()

        resa.update(resb)

        output =  Table(rows=[list(resa.values())], names=list(resa.keys()))

        return output

    def _transform(self, meta, par_names, params, vpar_names, covmat, mbfit, fitstatus):
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
        - SN parameters (from metadata):
           Dec, Ra, SNID, color, dL, daymax, epsilon_color, epsilon_daymax, epsilon_x0, epsilon_x1,
           index_hdf5, pixDec, pixID, pixRa, season, survey_area, x0, x1, z
        - result of the fit:
          z_fit, t0_fit, x0_fit, x1_fit, color_fit, Cov_t0t0, Cov_t0x0, Cov_t0x1, Cov_t0color, Cov_x0x0,
         Cov_x0x1, Cov_x0color, Cov_x1x1, Cov_x1color, Cov_colorcolor, mbfit, fitstatus
        """
        res = {}
        for key, value in meta.items():
            res[key] = value

        for i in range(len(par_names)):
            res[self.parNames[par_names[i]]+'_fit'] = params[i]

            for i, name in enumerate(vpar_names):
                for j, nameb in enumerate(vpar_names):
                    if j >= i:
                        if covmat is not None:
                            res['Cov_'+self.parNames[name] +
                                self.parNames[nameb]] = covmat[i, j]
                        else:
                            res['Cov_'+self.parNames[name] +
                                self.parNames[nameb]] = 0.

        res['mbfit'] = mbfit
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
        idx = lc['flux']/lc['fluxerr'] > 0.
        sel = lc[idx]

        return np.sqrt(np.sum((sel['flux']/sel['fluxerr'])**2))

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
    
    def plotLC(self, table, time_display=10):
        """ Light curve plot using sncosmo methods
        
        Parameters
        ---------------
        table: astropy table
         table with LS informations (flux, ...)
       time_display: float
         duration of the window display
        """

        import matplotlib.pyplot as plt

        if 'x1' in table.meta.keys():
            x1 = table.meta['x1']
            color = table.meta['color']
            x0 = table.meta['x0']
        else:
            x1 = 0.
            color = 0.
            x0 = 0.
        daymax = table.meta['daymax']

        print('plotting')
        model = sncosmo.Model('salt2')
        model.set(z=z,
                  c=color,
                  t0=daymax,
                  #x0=x0,
                  x1=x1)

        sncosmo.plot_lc(data=table)

        plt.draw()
        plt.pause(time_display)
        plt.close()
