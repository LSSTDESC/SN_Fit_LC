import sncosmo
import numpy as np
from astropy import (cosmology, units as u, constants as const)

class Fit_LC:
    def __init__(self,model='salt2-extended',version=1.0,telescope=None,display=False,bands='ugrizy'):
      
        self.display = display

        for band in bands:
            if telescope.airmass > 0:
                band=sncosmo.Bandpass(telescope.atmosphere[band].wavelen,telescope.atmosphere[band].sb, name='LSST::'+band,wave_unit=u.nm)
            else:
                band=sncosmo.Bandpass(telescope.system[band].wavelen,telescope.system[band].sb, name='LSST::'+band,wave_unit=u.nm) 
            sncosmo.registry.register(band, force=True)

        source=sncosmo.get_source(model,version=str(version))
        dust = sncosmo.OD94Dust()

        self.SN_fit_model=sncosmo.Model(source=source)
        
    def __call__(self,lc):
        
        z = lc.meta['z']
        self.SN_fit_model.set(z=z)

        select=lc[np.where(np.logical_and(lc['flux']/lc['fluxerr']>5.,lc['flux']>0.))]
        select = select[['flux','fluxerr','band','zp','zpsys','time']]
        print(select)
        res, fitted_model = sncosmo.fit_lc(select, self.SN_fit_model,['t0', 'x0', 'x1', 'c'],bounds={'z':(z-0.001, z+0.001)})
        mbfit=fitted_model._source.peakmag('bessellb','vega')

        #print('hello',self._transform(lc.meta,res['param_names'],list(res['parameters']),res['vparam_names'],res['covariance'],mbfit))
        """
        print('res',res)
        print('fitted',fitted_model)
        print('mbfit',mbfit)
        """
        
        if self.display:
            import pylab as plt
            sncosmo.plot_lc(select, model=fitted_model,color='r',pulls=False,errors=res.errors) 
            plt.show()

        return self._transform(lc.meta,res['param_names'],list(res['parameters']),res['vparam_names'],res['covariance'],mbfit)
    
    def _transform(self,meta,par_names,params,vpar_names,covmat,mbfit):

        rv = []
        names = []
        for key, value in meta.items(): 
            names.append(key)
            rv.append(value)

        
        names += [par_names[i]+'_fit' for i in range(len(par_names))]
        rv += params

        for i, name in enumerate(vpar_names):
            for j, nameb in enumerate(vpar_names):
                if j >= i:
                    names.append('Cov_'+name+nameb)
                    rv.append(covmat[i,j])
        names.append('mbfit')
        rv.append(mbfit)
        #print(names,rv)
        return names, rv
