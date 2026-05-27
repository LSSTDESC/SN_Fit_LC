# from sn_tools.sn_utils import multiproc
from importlib import import_module
import os
import h5py
from astropy.table import Table, vstack, unique
import sncosmo as sncosmo_emul
from sn_tools.sn_utils import register_bands_sncosmo
import numpy as np
import pandas as pd


class Fitting:
    """
    class to perform fits

    Parameters
    ----------------
    fitter_config: dict
      dict of parameters

    """

    def __init__(self, fitter_config, covmb=None):

        # load instrument
        from sn_telmodel.sn_throughputs import load_throughputs_from_config
        self.telescope = load_throughputs_from_config(
            fitter_config['InstrumentFit'])
        """
        tel_par = fitter_config['InstrumentFit']
        telescope = Telescope(name=tel_par['name'],
                              throughput_dir=tel_par['throughputDir'],
                              atmos_dir=tel_par['atmosDir'],
                              atmos=tel_par['atmos'],
                              aerosol=tel_par['aerosol'],
                              airmass=tel_par['airmass'])
        """
        self.mbcalc = fitter_config['mbcov']['estimate']
        self.covmb = covmb
        display_lc = fitter_config['Display']
        LC_sel = fitter_config['LCSelection']

        module = import_module(fitter_config['Fitter']['name'])

        # fit instance
        par_names = fitter_config['Fitter']['parnames'].split(',')
        sigmaz = fitter_config['Fitter']['sigmaz']
        self.snrmin = LC_sel['snrmin']
        fit_selected = fitter_config['fit']['selected']
        # config_inst = fitter_config['InstrumentFit']
        self.fit_coadded = fitter_config['fit']['coadded']

        self.fitter = module.Fit_LC(sncosmo_emul,
                                    model=fitter_config['Fitter']['model'],
                                    version=fitter_config['Fitter']['version'],
                                    snrmin=self.snrmin,
                                    fit_selected=fit_selected,
                                    vparam_names=par_names,
                                    telescope=self.telescope,
                                    sigmaz=sigmaz)

        if fitter_config['OutputFit']['save']:
            self.prepareSave(
                fitter_config['OutputFit']['directory'],
                fitter_config['ProductionIDFit'])

        self.alpha = 0.14
        self.beta = 3.1

        self.val = []

    def fit_lc(self, lc, params=None, j=-1, output_q=None):
        """
        fit_lc method: this is where the fit LC is performed

        Parameters
        ---------------
        lc: astropy table
          data to fit (LC points)
        params: dict
          dict of parameters
        j: int, opt
           internal parameter for multi processing (default: -1)
        output_q: multiprocessing.queue, opt
           queue for multiprocessing (default: None)


        Returns
        -----------
        astropytable with fitted values

        """
        remove_sat = False
        if params is not None:
            remove_sat = params['remove_sat']
        """
        if 'filter' in lc.columns:
            lc.remove_columns(['filter'])
        """
        # LC fit here
        resfit = self.fitter(lc, remove_sat=remove_sat)

        # estimate mbcov if requested
        if self.mbcalc and resfit:
            idx = resfit['fitstatus'] == 'fitok'
            if len(resfit[idx]) > 0:
                covDict = self.mbcovCalc(resfit)
                for key in covDict.keys():
                    resfit[key] = [covDict[key]]
            else:
                pp = ['Cov_x0mb', 'Cov_x1mb', 'Cov_colormb',
                      'Cov_mbmb', 'mb_recalc', 'sigma_mu']
                for key in pp:
                    resfit[key] = [-1.]

        if output_q is not None:
            output_q.put({j: resfit})
        else:
            return resfit

    def prepareSave(self, outdir, prodid):
        """
        Method to prepare to save output results on disk

        Parameters
        ---------------
        outdir: str
          output directory
        prodid: str
          outpu directory name (Fit_prodid.hdf5)

        """
        if not os.path.exists(outdir):
            print('Creating output directory', outdir)
            os.makedirs(outdir)

        self.fit_out = outdir+'/Fit_'+prodid+'.hdf5'
        if os.path.exists(self.fit_out):
            os.remove(self.fit_out)

    def dump(self, tab, inum):
        """
        Method to dump the results in a hdf5 file

        Parameters
        ---------------
        tab: astropy table
           data to dump
        inum: int
          internal parameter - key in the hdf5 file

        """

        # res = Table(np.rec.fromrecords(val,names = names))
        # tab = Table(rows=val, names=names)
        tab['fitstatus'] = tab['fitstatus'].astype(
            h5py.special_dtype(vlen=str))
        tab.write(self.fit_out, 'fit_lc_{}'.format(
            inum), append=True, compression=True)

    def mbcovCalc(self, vals):
        """
        Method to estimate mb covariance data

        Parameters
        ---------------
        covmb: Mbcov class
           class to estimate mb covariance data
        vals: astropy table
            fitted parameters

        Returns
        ----------
        dict with the following keys:
        Cov_x0mb,Cov_x1mb,Cov_colormb,Cov_mbmb,mb_recalc,sigma_mu

        """
        import numpy as np

        cov = np.ndarray(shape=(3, 3), dtype=float, order='F')
        cov[0, 0] = vals['Cov_x0x0'].data
        cov[1, 1] = vals['Cov_x1x1'].data
        cov[2, 2] = vals['Cov_colorcolor'].data
        cov[0, 1] = vals['Cov_x0x1'].data
        cov[0, 2] = vals['Cov_x0color'].data
        cov[1, 2] = vals['Cov_x1color'].data
        cov[2, 1] = cov[1, 2]
        cov[1, 0] = cov[0, 1]
        cov[2, 0] = cov[0, 2]

        x0f = 'x0_fit'
        x1f = 'x1_fit'
        colorf = 'color_fit'

        """
        if vals['x0_fit'] <= -90.:
            # this is probably fast fitter -> take simulated values as input
            x0f = 'x0'
            x1f = 'x1'
            colorf = 'color'
        """
        params = dict(zip(['x0', 'x1', 'c'], [vals[x0f].data,
                                              vals[x1f].data,
                                              vals[colorf].data]))

        resu = self.covmb.mbCovar(params, cov, ['x0', 'x1', 'c'])
        sigmu_sq = resu['Cov_mbmb']
        sigmu_sq += self.alpha**2 * vals['Cov_x1x1'].data + \
            self.beta**2 * vals['Cov_colorcolor'].data
        sigmu_sq += 2.*self.alpha*resu['Cov_x1mb']
        sigmu_sq += -2.*self.alpha*self.beta*vals['Cov_x1color'].data
        sigmu_sq += -2.*self.beta*resu['Cov_colormb']
        sigmu = np.array([0.])
        if sigmu_sq >= 0.:
            sigmu = np.sqrt(sigmu_sq)

        resu['sigma_mu'] = sigmu.item()
        resu['alpha'] = self.alpha
        resu['beta'] = self.beta
        return resu

    def fit_multiproc(self, lc_list, remove_sat=False, nproc=8):
        """
        Method to fit light curves

        Parameters
        ----------
        lc_list : list(astropy table)
            LC to fit
        remove_sat : bool, optional
            To remove saturated fluxes. The default is False.

        Returns
        -------
        None.

        """
        """
     from astropy.table import Table, vstack
     res = Table()
     for lc in lc_list:
         lc.convert_bytestring_to_unicode()
         resfit = self.fit(lc)
         if resfit is not None:
             res = vstack([res, resfit])

     return res
     """

        from sn_tools.sn_utils import multiproc
        params = {}
        params['remove_sat'] = remove_sat

        res = multiproc(lc_list, params, self.fit_lcs, nproc)

        return res

    def fit_lcs(self, lc_list, params, j=0, output_q=None):
        """
        Method to fit LCs

        Parameters
        ----------
        lc_list : list(astropy table)
            light-curves to fit.
        params : dict
            parameters.
        j : int, optional
            Tag for multiprocessing. The default is 0.
        output_q : multiprocessing queue, optional
            queue managing multiprocessing run. The default is None.

        Returns
        -------
        astropytable
            Result of the fit.

        """

        from astropy.table import Table, vstack

        # coadd if requested
        """
        import time
        time_ref = time.time()
        """

        lc_list = self.prepare_for_fit(lc_list)

        """
        if self.fit_coadded:
            lc_list = self.coadd_lcs(lc_list)

        print('coadd', time.time()-time_ref)
        # register bands in sn_cosmo here (gain time)

        time_ref = time.time()
        self.register_bands(lc_list)
        print('registry', time.time()-time_ref)
        # loop on lc_list and fit
        """
        res = Table()

        import time
        for io, lc in enumerate(lc_list):
            # print('fitting', j, io, len(lc))
            # register
            self.register_band(lc)
            #time_ref = time.time()
            # self.register_band(lc)
            # 'fitting', lc[['band_cosmo', 'airmass', 'pwv', 'ozone', 'aerosol']])
            lc.convert_bytestring_to_unicode()
            resfit = self.fit_lc(lc, params)
            # print('after fit', j, len(lc), time.time()-time_ref)
            if resfit is not None:
                resfit = self.check_correct(resfit)
                res = vstack([res, resfit])

        if output_q is not None:
            return output_q.put({j: res})
        else:
            return res

    def prepare_for_fit(self, lc_list):
        """
        Method to be ready for the fit:
            - LC points with SNR >= snrmin
            - possibility to coadd LC points
            - register bands in sn_cosmo

        Parameters
        ----------
        lc_list : list(light curves)
            LC list to process.

        Returns
        -------
        lc_res : list(LC)
            output LC list.

        """

        import time

        ccols = ['night', 'airmass', 'ozone', 'aerosol', 'mean_wave', 'band',
                 'pwv', 'zp', 'time', 'band_cosmo', 'zpsys', 'flux', 'fluxerr',
                 'snr_m5', 'snr', 'filter', 'sat', 'phase']

        mycols = ['flux', 'fluxerr', 'airmass',
                  'pwv', 'ozone', 'aerosol', 'filter', 'zp']

        lc_res = []
        for lc in lc_list:
            if len(lc) == 0:
                continue
            # SNR selection
            idx = lc['snr'] >= self.snrmin
            sel = lc[idx]

            """
            print('before coadd', len(sel))
            print(sel[mycols])
            rra = sel[mycols].to_pandas()
            """

            if self.fit_coadded:
                from sn_tools.sn_lcana import coadd_lc
                sel = coadd_lc(sel)

                """
                print('after coadd', len(sel))
                print(sel[mycols])
                rrb = sel[mycols].to_pandas()
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots()
                ax.hist(rra['flux'], histtype='step', bins=20)
                ax.hist(rrb['flux'], histtype='step', bins=20)
                figb, axb = plt.subplots()
                axb.hist(rra['zp'], histtype='step', bins=40)
                axb.hist(rrb['zp'], histtype='step', bins=40)
                plt.show()
                """
            lc_res.append(sel)

        # time_ref = time.time()

        """
        if len(lc_res) > 0:
            self.register_bands(lc_res)
        """
        # print('registry', time.time()-time_ref)

        return lc_res

    def coadd_lcs_deprecated(self, lc_list):
        """
        Method to coadd list of lcs

        Parameters
        ----------
        lc_list : list(lc)
            LC list.

        Returns
        -------
        res: list(LC)
          list of coadded lcs

        """

        ccols = ['night', 'airmass', 'ozone', 'aerosol', 'mean_wave', 'band',
                 'pwv', 'zp', 'time', 'band_cosmo', 'zpsys', 'flux', 'fluxerr',
                 'snr_m5', 'snr', 'filter', 'sat', 'phase']

        lc_res = []
        for lc in lc_list:
            # SNR selection
            idx = lc['snr'] >= self.snrmin
            sel = lc[idx]
            df = sel[ccols].to_pandas()
            dfb = df.groupby(['filter', 'night']).apply(
                lambda x: self.coadd_lc(x)).reset_index()
            dfb['band_cosmo'] = self.telescope.site_name+'::' + \
                dfb['filter']+'_' + \
                dfb['airmass'].astype(str)+'_' + \
                dfb['pwv'].astype(str)+'_' + \
                dfb['ozone'].astype(str)+'_' +\
                dfb['aerosol'].astype(str)
            dfb['band'] = dfb['band_cosmo']
            rr = Table.from_pandas(dfb)
            rr.meta = lc.meta
            lc_res.append(rr)

        return lc_res

    def coadd_lc_deprecated(self, grp,
                            col_means_weighted=[('flux', 'fluxerr')],
                            col_means=['airmass', 'pwv', 'ozone',
                                       'aerosol', 'mean_wave', 'zp', 'time'],
                            col_round=['airmass', 'pwv', 'ozone',
                                       'aerosol'],
                            round_vals=[1, 1, 1, 1],
                            col_unique=['zpsys']):
        """
        Method to coadd light-curve points per night/filter

        Parameters
        ----------
        grp : pandas df
            Data to process.
        col_means_weighted : list(str), optional
            list of cols for weighted mean estimation. 
            The default is [('flux','fluxerr')].
        col_means : list(str), optional
            list of cols for mean estimation. 
            The default is ['airmass','pwv','ozone','aerosol',
                            'mean_wave','zp','time'].
        col_round : list(str), optional
            list of cols to round. 
            The default is ['airmass','pwv','ozone',
                            'aerosol','zp','mean_wave'].
        round_vals : list(int), optional
            list of rounding values corresponding to col_round. 
            The default is [2,1,1,1,2,2].
        col_unique : list(str), optional
            list of cols with unique value. The default is ['zpsys'].

        Returns
        -------
        astropy table
        output value

        """

        """
        print('in coadd', len(grp))
        print(grp[['flux', 'fluxerr']])
        """

        grp['weight_flux'] = 1./grp['fluxerr']**2

        dictout = {}
        for vv in col_means_weighted:
            pp = vv[0]
            pp_weight = 'weight_{}'.format(pp)
            weight_sum = np.sum(grp[pp_weight])
            mean_weighted = np.sum(grp[pp]*grp[pp_weight])/weight_sum
            dictout[pp] = [mean_weighted]
            dictout[vv[1]] = [1./np.sqrt(weight_sum)]

        for vv in col_means:
            val = grp[vv].mean()
            if vv in col_round:
                idx = col_round.index(vv)
                val = np.round(val, round_vals[idx])

            dictout[vv] = [val]

        for vv in col_unique:
            dictout[vv] = grp[vv].unique().tolist()

        res_df = pd.DataFrame.from_dict(dictout)
        res_df['snr'] = res_df['flux']/res_df['fluxerr']

        """
        print('finally')
        print(res_df[['flux', 'fluxerr']])
        """
        return res_df

    def register_bands(self, lc_list):
        """
        Method to register bands in sncosmo

        Parameters
        ----------
        lc_list : list(LC)
            List of light curves to fit.

        Returns
        -------
        None.

        """

        tt = Table()
        ccols = ['band_cosmo', 'band', 'airmass',
                 'pwv', 'ozone', 'aerosol', 'filter']
        for lc in lc_list:
            tt = vstack([tt, lc[ccols]], metadata_conflicts='silent')

        if len(tt) > 1:
            tt = unique(tt)

        self.register_bands_on_the_fly(tt.to_pandas())

    def register_band(self, lc):
        """
        Method to register a lc in sncosmo

        Parameters
        ----------
        lc : astropy table
            LC to register.

        Returns
        -------
        None.

        """

        ccols = ['band_cosmo', 'band', 'airmass',
                 'pwv', 'ozone', 'aerosol', 'filter']

        tt = lc[ccols]

        if len(tt) > 1:
            tt = unique(tt)

        self.register_bands_on_the_fly(tt.to_pandas())

    def check_correct(self, sn):
        """
        Method to correct for Cov_xy col names

        Parameters
        ----------
        sn : astropy Table
            Data to process.

        Returns
        -------
        sn : astropy Table
            Processed data

        """

        varlist = ['z', 't0', 'x0', 'x1', 'color']

        if 'Cov_zz' not in sn.columns:
            varlist = ['t0', 'x0', 'x1', 'color']

        for i, namea in enumerate(varlist):
            for j, nameb in enumerate(varlist):
                if j >= i:
                    vva = 'Cov_{}{}'.format(namea, nameb)
                    vvb = 'Cov_{}{}'.format(nameb, namea)
                    if vva not in sn.columns:
                        sn.rename_column(vvb, vva)

        return sn

    def register_bands_on_the_fly(self, data):
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

            register_bands_sncosmo(sncosmo_emul, self.telescope,
                                   bandname, band,
                                   airmass, pwv, ozone, aerosol)
