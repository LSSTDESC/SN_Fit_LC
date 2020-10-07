from sn_tools.sn_telescope import Telescope
from astropy.table import Table, vstack
from importlib import import_module
import os
import numpy as np
import h5py


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
        tel_par = fitter_config['Instrument']
        telescope = Telescope(name=tel_par['name'],
                              throughput_dir=tel_par['throughputDir'],
                              atmos_dir=tel_par['atmosDir'],
                              atmos=tel_par['atmos'],
                              aerosol=tel_par['aerosol'],
                              airmass=tel_par['airmass'])

        self.mbcalc = fitter_config['Fitter']['covmb']
        self.covmb = covmb
        display_lc = fitter_config['Display']
        LC_sel = fitter_config['LCSelection']

        module = import_module(fitter_config['Fitter']['name'])
        # fit instance
        self.fitter = module.Fit_LC(
            model=fitter_config['Fitter']['model'],
            version=fitter_config['Fitter']['version'],
            telescope=telescope, display=display_lc,
            snrmin=LC_sel['snrmin'],
            nbef=LC_sel['nbef'],
            naft=LC_sel['naft'])

        if fitter_config['Output']['save']:
            self.prepareSave(
                fitter_config['Output']['directory'], fitter_config['ProductionID'])

        self.alpha = 0.14
        self.beta = 3.1

        self.val = []

    def __call__(self, lc, j=-1, output_q=None):
        """
        call method: this is where the fit LC is performed

        Parameters
        ---------------
        lc: astropy table
          data to fit (LC points)
        j: int, opt
           internal parameter for multi processing (default: -1)
        output_q: multiprocessing.queue, opt
           queue for multiprocessing (default: None)


        Returns
        -----------
        astropytable with fitted values

        """

        # LC fit here
        resfit = self.fitter(lc)

        # estimate mbcov if requested
        if self.mbcalc:
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

        params = dict(zip(['x0', 'x1', 'c'], [vals['x0_fit'].data,
                                              vals['x1_fit'].data, vals['color_fit'].data]))

        resu = self.covmb.mbCovar(params, cov, ['x0', 'x1', 'c'])
        sigmu_sq = resu['Cov_mbmb']
        sigmu_sq += self.alpha**2 * vals['Cov_x1x1'].data + \
            self.beta**2 * vals['Cov_colorcolor'].data
        sigmu_sq += 2.*self.alpha*resu['Cov_x1mb']
        sigmu_sq += -2.*self.alpha*self.beta*vals['Cov_x1color'].data
        sigmu_sq += -2.*self.beta*resu['Cov_colormb']
        sigmu = 0.
        if sigmu_sq >= 0.:
            sigmu = np.sqrt(sigmu_sq)

        resu['sigma_mu'] = sigmu.item()

        return resu
