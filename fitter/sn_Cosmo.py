import sncosmo
import numpy as np
from astropy import (cosmology, units as u, constants as const)

class Fit_LC:
    def __init__(self,model='salt2-extended',version=1.0,telescope=None,display=False,bands='ugrizy'):
      
        self.display = display
        self.bands = bands

        for band in bands:
            if telescope.airmass > 0:
                band=sncosmo.Bandpass(telescope.atmosphere[band].wavelen,telescope.atmosphere[band].sb, name='LSST::'+band,wave_unit=u.nm)
            else:
                band=sncosmo.Bandpass(telescope.system[band].wavelen,telescope.system[band].sb, name='LSST::'+band,wave_unit=u.nm) 
            sncosmo.registry.register(band, force=True)

        source=sncosmo.get_source(model,version=str(version))
        dust = sncosmo.OD94Dust()

        self.SN_fit_model=sncosmo.Model(source=source)
        self.parNames = dict(zip(['z', 't0', 'x0', 'x1', 'c'],['z', 't0', 'x0', 'x1', 'color']))


    def __call__(self,lc):

        if lc is None:
            return None, None
        
        res_param_names = ['z', 't0', 'x0', 'x1', 'c']
        res_params_values = np.zeros((5,), dtype=float)
        vparam_names = ['t0', 'x0', 'x1', 'c']
        covariance = np.zeros((4,4,), dtype=float)
        mbfit = -1
        z = -1
        fitstatus = 'nofit'
        meta = None
        if lc is not None:
            meta = lc.meta
            z = meta['z']
            self.SN_fit_model.set(z=z)
        
            
            select=lc[np.where(np.logical_and(lc['flux']/lc['fluxerr']>5.,lc['flux']>0.))]
            select = select[['flux','fluxerr','band','zp','zpsys','time']]
            select = lc[['flux','fluxerr','band','zp','zpsys','time']]
            if len(select) >= 5:
                try:
                    
                    res, fitted_model = sncosmo.fit_lc(select, self.SN_fit_model,['t0', 'x0', 'x1', 'c'],bounds={'z':(z-0.001, z+0.001)},min_snr=0.)
                    mbfit=fitted_model._source.peakmag('bessellb','vega')
                    res_param_names = res['param_names']
                    res_params_values = res['parameters']
                    vparam_names = res['vparam_names']
                    covariance = res['covariance']
                    fitstatus = 'fitok'
                except (RuntimeError, TypeError, NameError):
                    fitstatus = 'crash'
            else:
                fitstatus = 'nodat'
                
        resa = self._transform(meta,res_param_names, list(res_params_values), vparam_names,covariance,mbfit,fitstatus)
        resb = self._get_infos(z,res_params_values[res_param_names.index('t0')],select)

        if self.display and len(select) >= 5:
            import pylab as plt
            sncosmo.plot_lc(select, model=fitted_model,color='r',pulls=False,errors=res.errors) 
            plt.show()
            
        resa.update(resb)
        
        return resa.keys(), resa.values()
    
    def _transform(self,meta,par_names,params,vpar_names,covmat,mbfit,fitstatus):
        
        res = {}
        for key, value in meta.items(): 
            res[key]=value

        for i in range(len(par_names)):
            res[self.parNames[par_names[i]]+'_fit']=params[i]

        for i, name in enumerate(vpar_names):
            for j, nameb in enumerate(vpar_names):
                if j >= i:
                    res['Cov_'+self.parNames[name]+self.parNames[nameb]] = covmat[i,j]
                    
        res['mbfit'] = mbfit
        res['fitstatus'] = fitstatus
        return res

    def _get_infos(self, z, T0, lc):

        res = {}
    
        res['phase_min'] = 0.
        res['phase_max'] = 0.
        for band in self.bands:
            res['N_bef_'+band] = 0
            res['N_aft_'+band] = 0
            res['SNR_'+band] = 0.
            
        res['N_bef_all'] = 0
        res['N_aft_all'] = 0
        
        if len(lc) == 0:
            return res
        
        phases = (lc['time']-T0)/(1.+z)
        res['phase_min'] = np.min(phases)
        res['phase_max'] = np.max(phases)

        nbef = 0
        naft = 0
        for band in self.bands:
            idx = lc['band'] == 'LSST::'+band
            sel = lc[idx]
            if len(sel) > 0:
                diff = sel['time']-T0
                idb = diff >=0
                selb = sel[idb]
                res['N_bef_'+band] = len(sel)-len(selb)
                res['N_aft_'+band] = len(selb)
                nbef +=  len(sel)-len(selb)
                naft += len(selb)
                res['SNR_'+band] = self._calcSNR(sel)
            else:
                res['N_bef_'+band] = 0
                res['N_aft_'+band] = 0
                res['SNR_'+band] = 0.
        res['N_bef_all'] = nbef
        res['N_aft_all'] = naft
        
        return res

    def _calcSNR(self,lc):
        
        idx = lc['flux']/lc['fluxerr'] > 0.
        sel = lc[idx]
     
        sum_flux = np.sum(sel['flux'])
        rms_flux = np.sum(sel['fluxerr']**2)

        return sum_flux/np.sqrt(rms_flux)
