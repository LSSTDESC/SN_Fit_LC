# import sncosmo
import numpy as np
import pandas as pd
from astropy.table import Table
from sn_fit.sn_utils import Selection
import warnings
warnings.filterwarnings("ignore", category=UserWarning)


class Fit_LC(Selection):
    """
    class to fit simulated light curves

    Parameters
    ---------------
    model: str, opt
      model to use for the fit (default: salt2-extended)
    version: float, opt
       model version (default: 1.0)
    display: bool, opt
      to display the fit online (default: False)
    bands: str, opt
      bands to consider (default: ugrizy)
    snrmin: float, opt
      min SNR for considered LC points for the fit (default: 1)

    """

    def __init__(self, sncosmo_emul, model='salt2-extended', version=1.0,
                 bands='ugrizy',
                 snrmin=1, fit_selected=0,
                 vparam_names=['t0', 'x0', 'x1', 'c'],
                 outType='astropyTable', telescope=None,
                 sigmaz=1.e-5):
        # airmassType='const', airmass=1.2,
        # pwv=4.0, ozone=300., aerosol=0.0):
        super().__init__(snrmin)

        self.sncosmo = sncosmo_emul
        self.bands = bands
        self.fit_selected = fit_selected
        self.outType = outType

        # get the source
        source = self.sncosmo.get_source(model, version=str(version))
        # get the dust
        dustmap = self.sncosmo.OD94Dust()

        if model == 'salt3':
            source._wave[0] = 1500.  # used to be 1700
            source._wave[-1] = 24990.

        # sn_fit_model instance
        # self.SN_fit_model = sncosmo.Model(source=source)
        self.SN_fit_model = self.sncosmo.Model(source=source,
                                               effects=[dustmap, dustmap],
                                               effect_names=['host', 'mw'],
                                               effect_frames=['rest', 'obs'])

        # parameters to fit
        self.vparam_names = vparam_names
        # name parameters
        self.parNames = dict(zip(['z', 't0', 'x0', 'x1', 'c',
                                  'hostebv', 'hostr_v', 'mwebv', 'mwr_v'],
                                 ['z', 't0', 'x0', 'x1', 'color', 'hostebv',
                                  'hostr_v', 'mwebv', 'mwr_v']))

        self.telescope = telescope

        self.sigmaz = sigmaz
        """
        self.airmassType = airmassType
        self.airmass = airmass
        self.pwv = pwv
        self.ozone = ozone
        self.aerosol = aerosol
        """
        # self.register_bands()

    def register_bands_deprecated(self):
        """
        Method to register bands in sncosmo

        Returns
        -------
        None.

        """
        from sn_tools.sn_utils import register_bands_sncosmo_old
        airmass = [self.airmass]
        pwvs = [self.pwv]
        ozs = [self.ozone]
        aerosols = [self.aerosol]

        # values in pandas df
        cols = ['airmass', 'pwv', 'ozone', 'aerosol']
        vals = [airmass, pwvs, ozs, aerosols]
        df = pd.DataFrame.from_dict(dict(zip(cols, vals)))

        """
        print(df_dict)

        df = df_dict['airmass']
        for col in cols[1:]:
            df = df.merge(df_dict[col], how='cross')
        """
        for i, row in df.iterrows():
            airmass = row['airmass']
            aerosol = row['aerosol']
            pwv = row['pwv']
            ozone = row['ozone']
            register_bands_sncosmo_old(sncosmo,
                                       self.telescope,
                                       airmass, aerosol, pwv, ozone)

    def register_bands_params_deprecated(self, airmass, aerosol, pwv, ozone):
        """
        Method to register bands in sncosmo

        Parameters
        ----------
        airmass : float
            airmass value.
        aerosol : float
            aerosol value.
        pwv : float
            precipitable water vapor value.
        ozone : float
            ozone value.

        Returns
        -------
        None.

        """

        # band registery in sncosmo
        from astropy import units as u

        self.telescope.new_atmosphere(site_name=self.telescope.site_name,
                                      airmass=airmass,
                                      aerosol=aerosol,
                                      pwv=pwv, ozone=ozone)
        for band in 'grizy':
            name = '{}::{}_{}'.format(
                self.telescope.site_name, band, int(10*airmass))
            throughput = self.telescope.throughputs[band]
            bandcosmo = sncosmo.Bandpass(throughput.wavelen,
                                         throughput.sb,
                                         name=name,
                                         wave_unit=u.nm)
            sncosmo.registry.register(bandcosmo, force=True)

    def register_bands_on_the_fly_deprecated(self, telescope, data):
        """
        Method to register bands on sncosmo

        Parameters
        ----------
        telescope: Telescope class
            telescope to use
        data: pandas df
            data to register

        Returns
        -------
        None.

        """

        for i, row in data.iterrows():
            bandname = row['band_cosmo']
            band = row['filter']
            airmass = row['airmass']
            pwv = row['pwv']
            ozone = row['ozone']
            aerosol = row['aerosol']
            register_bands_sncosmo(sncosmo, telescope,
                                   bandname, band,
                                   airmass, pwv, ozone, aerosol)

    def register_bands_deprecated(self):
        """
        Method to register bands in sncosmo

        Returns
        -------
        None.

        """
        # band registery in sncosmo
        from astropy import units as u
        print('boooo', self.telescope.atmosDir)

        if self.airmassType != 'const':
            # for various airmass
            for airmass in range(10, 31, 1):
                for band in 'grizy':
                    name = '{}::{}_{}'.format(
                        self.telescope.name, band, airmass)
                    self.telescope.load_atmosphere(airmass/10, 'aerosol')
                    throughput = self.telescope.lsst_atmos_aerosol[band]
                    bandcosmo = sncosmo.Bandpass(
                        throughput.wavelen,
                        throughput.sb, name=name,
                        wave_unit=u.nm)
                    sncosmo.registry.register(bandcosmo, force=True)

        else:
            for band in 'grizy':
                name = '{}::{}'.format(self.telescope.name, band)
                self.telescope.load_atmosphere(1.2, 'aerosol')
                throughput = self.telescope.lsst_atmos_aerosol[band]
                bandcosmo = sncosmo.Bandpass(
                    throughput.wavelen,
                    throughput.sb, name=name,
                    wave_unit=u.nm)
                sncosmo.registry.register(bandcosmo, force=True)

    def __call__(self, lc, remove_sat=False):
        """
        call method: this is where the fit is effectively performed

        Parameters
        ---------------
        lc: astropy table
           lc points
       remove_sat : bool, optional
            To remove saturated fluxes. The default is False.

        Returns
        -----------

        """
        dict_res = {}

        dict_res['res_param_names'] = ['z', 't0', 'x0', 'x1', 'c']
        dict_res['res_params_values'] = np.zeros((5, 1), dtype=float)
        # vparam_names = ['t0', 'x0', 'x1', 'c']
        dict_res['vparam_names'] = self.vparam_names
        nc = len(dict_res['vparam_names'])
        dict_res['covariance'] = np.zeros((nc, nc), dtype=float)
        dict_res['mbfit'] = -1.
        dict_res['z'] = -1.
        dict_res['fitstatus'] = 'nofit'
        dict_res['chisq'] = 99999.
        dict_res['ndof'] = -1
        errors = dict(zip(self.vparam_names, [0.]*len(self.vparam_names)))
        dict_res['errors'] = errors
        meta = None

        if lc is None:
            return None

        # get metadata
        meta = lc.meta

        fit_please = True

        if len(lc) == 0:
            fit_please = False

        # remove saturated flux here
        if len(lc) > 0:
            if remove_sat:
                idx = lc['sat'] == 0
                lc = lc[idx]

            if len(lc) == 0:
                fit_please = False

        if 'selected' in lc.meta.keys():
            if self.fit_selected and not lc.meta['selected']:
                fit_please = False

        if fit_please:
            dict_res = self.fitIt(meta,
                                  # dict_res['vparam_names'],
                                  lc)
        else:
            fitstatus = 'nodat'
            dict_res['fitstatus'] = fitstatus

        # Make a dict of the fitted result (plus metadata)
        res = dict_res
        if self.outType == 'astropyTable':
            resa = self._transform(meta,
                                   dict_res['res_param_names'],
                                   list(dict_res['res_params_values']),
                                   dict_res['vparam_names'],
                                   dict_res['covariance'],
                                   dict_res['mbfit'],
                                   dict_res['fitstatus'],
                                   dict_res['chisq'],
                                   dict_res['ndof'],
                                   dict_res['errors'])

            output = Table(rows=[list(resa.values())], names=list(resa.keys()))
            res = output

        return res

    def fitIt(self, meta, lc):
        """
        Method to (try) to perform LC fit

        Parameters
        ----------
        meta : dict
            LC metedata.
        lc : astropy table
            Light curve.
        Returns
        -------
        res_param_names : TYPE
            DESCRIPTION.
        res_params_values : TYPE
            DESCRIPTION.
        vparam_names : TYPE
            DESCRIPTION.
        covariance : TYPE
            DESCRIPTION.
        mbfit : TYPE
            DESCRIPTION.
        fitstatus : TYPE
            DESCRIPTION.
        chisq : TYPE
            DESCRIPTION.
        ndof : TYPE
            DESCRIPTION.

        """

        res_param_names = ['z', 't0', 'x0', 'x1', 'c']
        res_params_values = np.zeros((5, 1), dtype=float)
        # vparam_names = ['t0', 'x0', 'x1', 'c']
        vparam_names = self.vparam_names
        nc = len(self.vparam_names)
        errors = dict(zip(vparam_names, [0.]*len(vparam_names)))
        covariance = np.zeros((nc, nc), dtype=float)
        mbfit = -1.
        z = -1.
        fitstatus = 'nofit'
        chisq = 99999.
        ndof = -1
        fitted_model = None
        res = None

        """
        if 'filter' in lc.columns and 'band' in lc.columns:
            del lc['filter']
        """

        # set redshift for the fit
        z = meta['z']
        # daymax = meta['daymax']

        zmin = z-self.sigmaz*(1+z)
        zmax = z+self.sigmaz*(1+z)
        bounds = {'z': (zmin, zmax),
                  'x1': (-5.0, 5.0), 'c': (-0.5, 0.5)}

        if 'z' not in self.vparam_names:
            self.SN_fit_model.set(z=meta['zmeas'])
            bounds = {'x1': (-5.0, 5.0), 'c': (-0.5, 0.5)}

        # apply extinction here
        self.SN_fit_model.set(mwebv=meta['ebvofMW'])

        select_lc = self.select(lc)

        if select_lc is not None:
            try:
                """
                # register bands in sn_cosmo here
                the_data = select_lc[['band_cosmo',
                                      'band',
                                      'airmass',
                                      'pwv', 'ozone', 'aerosol', 'filter']].to_pandas()
                import time
                time_ref = time.time()
                self.register_bands_on_the_fly(self.telescope, the_data)
                print('after registry', time.time()-time_ref)
                """
                # fit here
                selfit = select_lc.copy()

                if 'filter' in selfit.columns and 'band' in selfit.columns:
                    del selfit['filter']

                idx = selfit['fluxerr'] > 0.
                idx &= selfit['flux'] >= 0.
                selfit = selfit[idx]

                res, fitted_model = self.sncosmo.fit_lc(
                    selfit, model=self.SN_fit_model,
                    vparam_names=self.vparam_names,
                    bounds=bounds, minsnr=self.snrmin)

                # get parameters
                if res['success']:
                    mbfit = fitted_model._source.peakmag(
                        'bessellb', 'vega')
                    res_param_names = res['param_names']
                    res_params_values = res['parameters']
                    vparam_names = res['vparam_names']
                    covariance = res['covariance']
                    errors = res['errors']
                    fitstatus = 'fitok'
                    chisq = res.chisq
                    ndof = res.ndof

                else:
                    # print('badfit',res)
                    fitstatus = 'bafit'
            except (RuntimeError, TypeError, NameError) as err:
                print('crash fit', err)
                fitstatus = 'crash'
                # set the simulation values here
                if meta['sn_type'] == 'SN_Ia':
                    res_params_values = np.array(
                        [meta['z'], meta['daymax'], meta['x0'],
                         meta['x1'], meta['color']])
                else:
                    res_params_values = np.array(
                        [meta['z'], meta['daymax'], -1.0, -1.0, -1.0])
        else:
            fitstatus = 'nodat'

        dict_res = {}
        dict_res['res_param_names'] = res_param_names
        dict_res['res_params_values'] = res_params_values
        dict_res['vparam_names'] = vparam_names
        dict_res['covariance'] = covariance
        dict_res['mbfit'] = mbfit
        dict_res['fitstatus'] = fitstatus
        dict_res['chisq'] = chisq
        dict_res['ndof'] = ndof
        dict_res['fitted_model'] = fitted_model
        dict_res['errors'] = errors
        if res is not None:
            dict_res['res_errors'] = res.errors
        else:
            dict_res['res_errors'] = None
        dict_res['fitstatus'] = fitstatus
        dict_res['lc'] = lc

        return dict_res

    def plotIt(self, select, fitted_model, errors,
               fitstatus,
               figtext=''):
        """
        Method to plot LC

        Parameters
        ----------
        select : astropy table
            lc to plot.
        fitted_model : sn_cosmo model
            fitted model.
        errors : array
            fit params error.
        fitstatus : str
            fit status.
        block: bool, opt
         for remnent display. Default: False

        Returns
        -------
        None.

        """

        if len(select) >= 1:
            if fitstatus == 'fitok':
                fig = self.sncosmo.plot_lc(select, model=fitted_model,
                                           errors=errors, xfigsize=10, pulls=False,
                                           figtext=figtext)
            else:
                fig = self.sncosmo.plot_lc(
                    select, xfigsize=10, figtext=figtext)
            return fig
        return None

    def _transform(self, meta, par_names, params, vpar_names, covmat,
                   mbfit, fitstatus, chisq, ndof, errors):
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
        chisq: float
          chi2 value of the fit
        ndof: int
          number of degree of freedom of the fit
        errors: dict
         dict of parameter errors

        Returns
        ----------
        dict with the following keys
        - SN parameters (from metadata):
           Dec, Ra, SNID, color, dL, daymax, epsilon_color, epsilon_daymax,
           epsilon_x0, epsilon_x1,
           index_hdf5, pixDec, pixID, pixRa, season, survey_area, x0, x1, z
        - result of the fit:
          z_fit, t0_fit, x0_fit, x1_fit, color_fit, Cov_t0t0, Cov_t0x0,
          Cov_t0x1, Cov_t0color, Cov_x0x0,
         Cov_x0x1, Cov_x0color, Cov_x1x1, Cov_x1color, Cov_colorcolor,
         mbfit, fitstatus,chisq,ndof
        """
        res = {}
        for key, value in meta.items():
            res[key] = value

        for i in range(len(par_names)):
            res[self.parNames[par_names[i]]+'_fit'] = params[i].item()

        for i, name in enumerate(vpar_names):
            for j, nameb in enumerate(vpar_names):
                if j >= i:
                    if covmat is not None:
                        res['Cov_'+self.parNames[name] +
                            self.parNames[nameb]] = covmat[i, j]
                    else:
                        res['Cov_'+self.parNames[name] +
                            self.parNames[nameb]] = 0.

        for name in vpar_names:
            res['sigma_{}'.format(name)] = errors[name]

        res['mbfit'] = mbfit
        res['fitstatus'] = fitstatus
        res['chisq'] = chisq
        res['ndof'] = ndof

        vvs = ['hostebv_fit', 'hostr_v_fit', 'mwebv_fit', 'mwr_v_fit']

        for vv in vvs:
            if vv not in res.keys():
                res[vv] = -1.0

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
        z = table.meta['z']

        print('plotting')
        model = self.sncosmo.Model('salt2')
        model.set(z=z,
                  c=color,
                  t0=daymax,
                  # x0=x0,
                  x1=x1)

        self.sncosmo.plot_lc(data=table)

        plt.draw()
        plt.pause(time_display)
        plt.close()
