from astropy.table import Table, vstack
import numpy as np


class Selection:
    """
    class to select a ligh curve according to input criterias

    Parameters
    --------------
    snrmin: float
       min SNR of LC points to be considered
   
    """

    def __init__(self, snrmin):

        self.snrmin = snrmin

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

        """
        if not self.include_errmodel_in_lcerror:
            lc['fluxerr'] = lc['fluxerr_photo']
        """
        # some cleaning here
        idx = lc['flux'] > 0.
        idx &= lc['fluxerr'] > 0.

        selecta = lc[idx]
        # add snr
        selecta['snr'] = selecta['flux']/selecta['fluxerr']


        # select LC points according to SNRmin
        idx = selecta['snr'] >= self.snrmin
        selecta = selecta[idx]

        if len(selecta) == 0:
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
           redshift below which the cut will be applied

        Returns
        ----------
        selected lc
        """

        idx = tab['band'] == 'LSST::{}'.format(b)
        sel = tab[idx]
        if len(sel) == 0:
            return Table()

        if sel.meta['z'] >= zref:
            idb = sel['fluxerr_model'] >= 0.
            idb &= sel['fluxerr_model']/sel['flux'] <= self.errmodrel
            selb = sel[idb]
            return selb

        return sel
