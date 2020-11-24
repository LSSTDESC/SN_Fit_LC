from astropy.table import Table,vstack
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

    """
    def __init__(self,snrmin, nbef, naft,nbands):
        self.snrmin = snrmin
        self.nbef = nbef
        self.aft = naft
        self.nbands = nbands
        

    def select(self,lc):
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

        # some cleaning here
        idx = lc['flux'] > 0.
        idx &= lc['fluxerr'] > 0.
      
        selecta = lc[idx]
        # add snr
        selecta['snr'] = selecta['flux']/selecta['fluxerr']

        # select LC points according to SNRmin
        idx = selecta['snr']>= self.snrmin
        selecta = selecta[idx]

        if len(selecta) == 0:
            return None
        
        # number of points before/after
        selecta['diff_time'] = selecta.meta['daymax']-selecta['time']
        nlc_bef = len(selecta[selecta['diff_time']>=0])
        nlc_aft = len(selecta[selecta['diff_time']<0])

        # check the total number of LC points here
        assert((nlc_bef+nlc_aft)==len(selecta))

        if nlc_bef < self.nbef or nlc_aft < self.naft:
            return None
        # select the number of bands to be used
        selb = Table()
        for b in np.unique(selecta['band']):
            io = selecta['band']==b
            selo = selecta[io]
            ibo = selo['snr']>=5.
            #print(b,len(selo[ibo]))
            if len(selo[ibo]) >= 2.:
                selb = vstack([selb,selo])

        select = Table(selb)

        nbands = 0
        if len(select) > 0.:
            nbands = len(np.unique(select['band']))
        
        if nbands < self.nbands:
            return None

        return select
