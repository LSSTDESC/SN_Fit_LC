from astropy.table import Table, vstack
import numpy as np


class Selection:
    """
    class to select a ligh curve according to input criterias

    Parameters
    --------------
    snrmin: float
       min SNR of LC points to be considered
    nbef: float
       number of LC points (SNR>=snrmin) before max
    naft: float
       number of LC points (SNR>=snrmin) after max
    nbands: int
      number of bands with at least two points with SNR>=5.
    phase_min: float
      min phase for nphase_min estimator
    phase_max: float
      max phase for nphase_max estimator
    nphase_min: int
      min number of points with phase<=phase_min
    nphase_max: int
      min number of points with phase>=phase_max
    errmodrel: float,opt
      max errormodel relative error (default: -1.)
    include_errmodel_in_lcerror: bool, opt
      to include the error model in lc point errors
    """

    def __init__(self, snrmin, nbef, naft, nbands, phase_min, phase_max, nphase_min, nphase_max, errmodrel=-1., include_errmodel_in_lcerror=False):

        self.snrmin = snrmin
        self.nbef = nbef
        self.naft = naft
        self.nbands = nbands
        self.phase_min = phase_min
        self.phase_max = phase_max
        self.nphase_min = nphase_min
        self.nphase_max = nphase_max
        self.errmodrel = errmodrel
        self.include_errmodel_in_lcerror = include_errmodel_in_lcerror

    def select(self, lc):
        """
        Method to select LC according to criteria

        Parameters
        --------------
        lc : astropy table
          light curve to process

        Returns
        ----------
        selected astropy table

        """

        if lc is None or lc.meta['status'] != 1:
            return None

        if not self.include_errmodel_in_lcerror:
            lc['fluxerr'] = lc['fluxerr_photo']

        # some cleaning here
        idx = lc['flux'] > 0.
        idx &= lc['fluxerr'] > 0.

        selecta = lc[idx]
        # add snr
        selecta['snr'] = selecta['flux']/selecta['fluxerr']

        # remove points with too high errormodel
        if self.errmodrel > 0.:
            selecta = self.select_error_model(selecta)

        # select LC points according to SNRmin
        idx = selecta['snr'] >= self.snrmin
        selecta = selecta[idx]

        if len(selecta) == 0:
            return None

        # number of points before/after
        if 'phase' not in selecta.columns:
            selecta['phase'] = (
                selecta['time']-selecta.meta['daymax'])/(1.+selecta.meta['z'])
        nlc_bef = len(selecta[selecta['phase'] <= 0])
        nlc_aft = len(selecta[selecta['phase'] > 0])

        # check the total number of LC points here
        assert((nlc_bef+nlc_aft) == len(selecta))

        if nlc_bef < self.nbef or nlc_aft < self.naft:
            return None

        # phase min and phase max sel
        idd = selecta['phase'] <= self.phase_min
        nph_min = len(selecta[idd])
        idd = selecta['phase'] >= self.phase_max
        nph_max = len(selecta[idd])

        if nph_min < self.nphase_min or nph_max < self.nphase_max:
            return None

        if self.nbands > 0:
            # select the number of bands to be used
            selb = Table()
            for b in np.unique(selecta['band']):
                io = selecta['band'] == b
                selo = selecta[io]
                ibo = selo['snr'] >= 5.
                # print(b,len(selo[ibo]))
                if len(selo[ibo]) >= 2.:
                    selb = vstack([selb, selo])

            selecta = Table(selb)

            nbands = 0
            if len(selecta) > 0.:
                nbands = len(np.unique(selecta['band']))

            if nbands < self.nbands:
                return None

        return selecta

    def select_error_model(self, lc):
        """
        function to select LCs

        Parameters
        ---------------
        lc: astropy table
          lc to consider

        Returns
        ----------
        lc with filtered values
       """

        if self.errmodrel < 0.:
            return lc

        # first: select iyz bands

        bands_to_keep = []

        lc_sel = Table()
        for b in 'izy':
            bands_to_keep.append('LSST::{}'.format(b))
            idx = lc['band'] == 'LSST::{}'.format(b)
            lc_sel = vstack([lc_sel, lc[idx]])

        # now apply selection on g band for z>=0.25
        sel_g = self.sel_band(lc, 'g', 0.25)

        # now apply selection on r band for z>=0.6
        sel_r = self.sel_band(lc, 'r', 0.6)

        lc_sel = vstack([lc_sel, sel_g])
        lc_sel = vstack([lc_sel, sel_r])

        return lc_sel

    def sel_band(self, tab, b, zref):
        """
        Method to performe selections depending on the band and z

        Parameters
        ---------------
        tab: astropy table
          lc to process
        b: str
          band to consider
        zref: float
           redshift below wiwh the cut wwill be applied

        Returns
        ----------
        selected lc
        """

        idx = tab['band'] == 'LSST::{}'.format(b)
        sel = tab[idx]
        if len(sel) == 0:
            return Table()

        if sel.meta['z'] >= zref:
            idb = sel['fluxerr_model']/sel['flux'] <= self.errmodrel
            selb = sel[idb]
            return selb

        return sel
