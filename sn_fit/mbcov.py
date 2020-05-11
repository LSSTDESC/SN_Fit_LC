import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
import os


class MbCov:
    def __init__(self, salt2Dir, paramNames=dict(zip(['x0', 'x1', 'color'], ['x0', 'x1', 'color'])), interp=True):
        """
        Class to estimate covariance matrix with mb
        Parameters
        ----------
        salt2Dir : str
         director where SALT2 reference files are to be found
        """

        self.load(salt2Dir)
        self.paramNames = paramNames
        # self.transNames = dict(
        #    zip(['t0', 'x0', 'x1', 'c'], ['t0', 'x0', 'x1', 'color']))
        self.transNames = dict(map(reversed, paramNames.items()))
        self.interp = interp
        if self.interp:
            # check whether outputdir is ready
            self.ratName = '{}/RatInt_for_mb.npy'.format(salt2Dir)
            if not os.path.exists(self.ratName):
                print('You would like to use the fast calculation for mB')
                print('But you need a file that does not exist')
                print('Will create it for you (it will take ~25 min)')
                self.genRat_int()
                print('The file has been generated.')
            self.ratio_Int = np.load(self.ratName)

    def load(self, salt2Dir):
        """
        Load a set of SALT2 files requested for mb cov estimations
        """
        # from F. Mondon 2017/10/20
        # wavelength limits for salt2 model
        wl_min_sal = 3000
        wl_max_sal = 7000

        # interpolation of TB and Trest
        # filt2 = np.genfromtxt('{}/snfit_data/Instruments/SNLS3-Landolt-model/sb-shifted.dat'.format(salt2Dir))
        filt2 = np.genfromtxt(
            '{}/Instruments/SNLS3-Landolt-model/sb-shifted.dat'.format(salt2Dir))
        filt2 = np.genfromtxt(
            '{}/Instruments/Landolt/sb_-41A.dat'.format(salt2Dir))
        wlen = filt2[:, 0]
        tran = filt2[:, 1]
        self.splB = Spline1d(wlen, tran, k=1, ext=1)

        # interpolation of ref spectrum
        # data = np.genfromtxt(thedir+'/snfit_data/MagSys/bd_17d4708_stisnic_002.ascii')
        data = np.genfromtxt(
            '{}/MagSys/bd_17d4708_stisnic_002.ascii'.format(salt2Dir))
        dispersion = data[:, 0]
        flux_density = data[:, 1]
        self.splref = Spline1d(dispersion, flux_density, k=1, ext=1)

        # interpolation of the spectrum model
        template_0 = np.genfromtxt(
            '{}/snfit_data/salt2-4/salt2_template_0.dat'.format(salt2Dir))
        template_1 = np.genfromtxt(
            '{}/snfit_data/salt2-4/salt2_template_1.dat'.format(salt2Dir))

        wlM0 = []
        M0 = []
        for i in range(len(template_0[:, 0])):
            if template_0[:, 0][i] == 0.0:
                wlM0.append(template_0[:, 1][i])
                M0.append(template_0[:, 2][i])
        self.splM0 = Spline1d(wlM0, M0, k=1, ext=1)

        wlM1 = []
        M1 = []
        for i in range(len(template_1[:, 0])):
            if template_1[:, 0][i] == 0.0:
                wlM1.append(template_1[:, 1][i])
                M1.append(template_1[:, 2][i])
        self.splM1 = Spline1d(wlM1, M1, k=1, ext=1)

        # computation of the integral
        dt = 100000
        self.xs = np.linspace(float(wl_min_sal), float(wl_max_sal), dt)
        self.dxs = (float(wl_max_sal-wl_min_sal)/(dt-1))

        self.I2 = np.sum(self.splref(self.xs)*self.xs *
                         self.splB(self.xs)*self.dxs)

    def genRat_int(self):
        """
        Estimate set of ratios
        """
        x1 = np.arange(-3.0, 3.0, 0.1)
        color = np.arange(-0.3, 0.3, 0.01)
        x1_all = np.repeat(x1, len(color))
        color_all = np.tile(color, len(x1))

        r = []
        mref = 9.907
        for (x1v, colorv) in zip(x1_all, color_all):
            r.append((x1v, colorv, self.ratInt(x1v, colorv), mref))

        tab = np.rec.fromrecords(r, names=['x1', 'color', 'ratioInt', 'mref'])

        np.save(self.ratName, tab)

    def ratInt(self, x1, color):
        """
        Estimate a ratio of two sums requested to estimated mb
        Parameters
        ---------------
        x1: float
         x1 of the supernova
        color: float
         color of the supernova
        Returns
        -----------
        float
         ratio value
        """
        I1 = np.sum((self.splM0(self.xs)*10**-12+x1*self.splM1(self.xs)*10**-12)*(
            10**(-0.4*self.color_law_salt2(self.xs)*color))*self.xs*self.splB(self.xs)*self.dxs)

        return I1/self.I2

    def mB_interp(self, x0, x1, color):
        """
        Estimate mB interpolation for supernovae
        Parameters
        ----------------
        params: dict
         dict of parameters: x0, x1, color
        Returns
        -----------
        mb : float
         mb value
        """

        rat = griddata((self.ratio_Int['x1'], self.ratio_Int['color']),
                       self.ratio_Int['ratioInt'], (x1, color), method='cubic')
        mb = -2.5*np.log10(x0*rat)+np.mean(self.ratio_Int['mref'])

        return mb

    def mB(self, params):
        """
        Estimate mB for supernovae
        Parameters
        ---------------
        params: dict
         dict of parameters: x0, x1, color
        Returns
        -----------
        mb : float
         mb value
        """

        rat = self.ratInt(params[self.paramNames['x1']],
                          params[self.paramNames['color']])
        # computation of mb
        mref = 9.907
        mb = -2.5*np.log10(params[self.paramNames['x0']]*rat)+mref

        return mb

    def mB_old(self, params):
        """ Estimate mB for supernovae
        Parameters
        ----------
        params: dict
         dict of parameters: x0, x1, color
        Returns
        -------
        mb : float
         mb value
        """

        #    I1=np.sum((splM0(xs)*10**-12+res.parameters[3]*splM1(xs)*10**-12)*(10**(-0.4*salt2source.colorlaw(xs)*res.parameters[4]))*xs*splB(xs)*dxs)
        I1 = np.sum((self.splM0(self.xs)*10**-12+params[self.paramNames['x1']]*self.splM1(self.xs)*10**-12)*(
            10**(-0.4*self.color_law_salt2(self.xs)*params[self.paramNames['color']]))*self.xs*self.splB(self.xs)*self.dxs)
        I2 = np.sum(self.splref(self.xs)*self.xs*self.splB(self.xs)*self.dxs)
        # print(I1, I2,params['x1'],params['c'])

        # computation of mb
        mref = 9.907
        mb = -2.5*np.log10(params[self.paramNames['x0']]*(I1/I2))+mref

        return mb
    """
    def calcInteg(self, bandpass, signal,wavelen):
        fa = interpolate.interp1d(bandpass.wavelen,bandpass.sb)
        fb = interpolate.interp1d(wavelen,signal)
        min_wave=np.max([np.min(bandpass.wavelen),np.min(wavelen)])
        max_wave=np.min([np.max(bandpass.wavelen),np.max(wavelen)])
        # print 'before integrand',min_wave,max_wave
        wavelength_integration_step=5
        waves=np.arange(min_wave,max_wave,wavelength_integration_step)
        integrand=fa(waves) *fb(waves)
        # print 'rr',len(f2(wavelen)),len(wavelen),len(integrand)
        range_inf=min_wave
        range_sup=max_wave
        n_steps = int((range_sup-range_inf) / wavelength_integration_step)
        x = np.core.function_base.linspace(range_inf, range_sup, n_steps)
        # print len(waves),len(x)
        return integrate.simps(integrand,x=waves)
    def Get_Mag(self,filename,name,band):
        sfile=open(filename,'rb')
        spectrum_file='unknown'
        for line in sfile.readlines():
            if 'SPECTRUM' in line:
                spectrum_file=line.split(' ')[1].strip()
            if name in line and band in line:
                return float(line.split(' ')[2]),spectrum_file
        sfile.close()
    """

    def color_law_salt2(self, wl):
        """ Color law for SALT2
        """
        B_wl = 4302.57
        V_wl = 5428.55
        l = (wl-B_wl)/(V_wl-B_wl)
        l_lo = (2800.-B_wl)/(V_wl-B_wl)
        l_hi = (7000.-B_wl)/(V_wl-B_wl)
        a = -0.504294
        b = 0.787691
        c = -0.461715
        d = 0.0815619
        cst = 1-(a+b+c+d)
        cl = []
        for i in range(len(l)):
            if l[i] > l_hi:
                cl.append(-(cst*l_hi+l_hi**2*a+l_hi**3*b+l_hi**4*c+l_hi**5*d +
                            (cst+2*l_hi*a+3*l_hi**2*b+4*l_hi**3*c+5*l_hi**4*d)*(l[i]-l_hi)))
            if l[i] < l_lo:
                cl.append(-(cst*l_lo+l_lo**2*a+l_lo**3*b+l_lo**4*c+l_lo**5*d +
                            (cst+2*l_lo*a+3*l_lo**2*b+4*l_lo**3*c+5*l_lo**4*d)*(l[i]-l_lo)))
            if l[i] >= l_lo and l[i] <= l_hi:
                cl.append(-(cst*l[i]+l[i]**2*a+l[i]**3*b+l[i]**4*c+l[i]**5*d))
        return np.array(cl)

    def mbCovar(self, params, covar, vparam_names):
        """ mb covariance matrix wrt fit parameters
        Parameters
        ----------
        params: dict
         parameter values
        covar: matrix
         covariance matrix of the parameters
        vparam_names: list
         names of the parameters
        Returns
        -------
        res: dict
         final covariance dict
        """

        res = {}
        h_ref = 1.e-8
        Der = np.zeros(shape=(len(vparam_names), 1))

        par_var = params.copy()
        ider = -1
        for i, key in enumerate(vparam_names):
            h = h_ref
            if np.abs(par_var[key]) < 1.e-5:
                h = 1.e-10

            par_var[key] += h
            ider += 1
            Der[ider] = (self.mB(par_var)-self.mB(params))/h

            par_var[key] -= h

        Prod = np.dot(covar, Der)

        for i, key in enumerate(vparam_names):
            res['Cov_{}mb'.format(self.transNames[key])] = Prod[i, 0]
            """
            if key != 'c':
                res['Cov_'+key.upper()+'mb']=Prod[i,0]
            else:
               res['Cov_Colormb']=Prod[i,0]
            """
        res['Cov_mbmb'] = np.dot(Der.T, Prod).item()
        res['mb_recalc'] = self.mB(par_var).item()

        return res

    def mbCovar_int(self, params, covar, vparam_names):
        """ mb covariance matrix wrt fit parameters
            uses mb_int (griddata from template of mb)
        Parameters
        ----------
        params: dict
         parameter values
        covar: matrix
         covariance matrix of the parameters
        vparam_names: list
         names of the parameters
        Returns
        -------
        res: dict
         final covariance dict
        """

        res = {}
        h_ref = 1.e-8
        Der = np.zeros(shape=(len(vparam_names), 1))

        rt = []
        r = []
        for i, key in enumerate(vparam_names):
            r.append(params[key])
        r.append(0.0)
        rt.append(tuple(r))

        for i, key in enumerate(vparam_names):
            rot = list(rt[0])
            h = h_ref
            if np.abs(rot[i]) < 1.e-5:
                h = 1.e-10
            rot[i] += h
            rot[-1] = h
            rt.append(tuple(rot))

        tabDiff = np.rec.fromrecords(rt, names=vparam_names+['h'])
        mbVals = self.mB_interp(
            tabDiff[self.paramNames['x0']], tabDiff[self.paramNames['x1']], tabDiff[self.paramNames['color']])
        tabDiff = rf.append_fields(tabDiff, 'mB', mbVals)

        ider = -1
        for i, key in enumerate(vparam_names):
            ider += 1
            Der[ider] = (tabDiff['mB'][i+1]-tabDiff['mB'][0])/tabDiff['h'][i+1]

        Prod = np.dot(covar, Der)

        for i, key in enumerate(vparam_names):
            res['Cov_{}mb'.format(self.transNames[key])] = Prod[i, 0]

        res['Cov_mbmb'] = np.asscalar(np.dot(Der.T, Prod))
        res['mb_recalc'] = self.mB_interp(
            params[self.paramNames['x0']], params[self.paramNames['x1']], params[self.paramNames['color']]).item()

        return res
    """
    def mbDeriv(self,params,vparam_names):
        res={}
        h=1.e-6
        # Der=np.zeros(shape=(len(vparam_names),1))
        Der={}
        # print params
        par_var=params.copy()
        ider=-1
        for i,key in enumerate(vparam_names):
            par_var[key]+=h
            ider+=1
            Der[key]=(self.mB(par_var)-self.mB(params))/h
            par_var[key]-=h
        return Der
    """

    def test(self):
        """ Test function
        To test whether this class is usable or not
        """

        """
        Salt2Model
        BEGIN_OF_FITPARAMS Salt2Model
        DayMax 53690.0336018 0.105513809169
        Redshift 0.1178 0 F
        Color -0.0664131339433 0.0234330339301
        X0 0.00030732251016 8.89813428854e-06
        X1 -0.0208012409076 0.160846457522
        CovColorColor 0.00054910707917 -1
        CovColorDayMax 0.00040528682468 -1
        CovColorX0 -1.68238293879e-07 -1
        CovColorX1 0.00114702847231 -1
        CovDayMaxDayMax 0.0111331639253 -1
        CovDayMaxX0 -2.94345317778e-07 -1
        CovDayMaxX1 0.0131008809199 -1
        CovX0X0 7.91767938168e-11 -1
        CovX0X1 -7.23852420336e-07 -1
        CovX1X1 0.0258715828973
        """

        salt2_res = {}
        salt2_res['DayMax'] = 53690.0336018
        salt2_res['Color'] = -0.0664131339433
        salt2_res['X0'] = 0.00030732251016
        salt2_res['X1'] = -0.0208012409076
        salt2_res['Color'] = 0.0
        salt2_res['X1'] = 0.0
        salt2_res['CovColorColor'] = 0.00054910707917
        salt2_res['CovColorDayMax'] = 0.00040528682468
        salt2_res['CovColorX0'] = -1.68238293879e-07
        salt2_res['CovColorX1'] = 0.00114702847231
        salt2_res['CovDayMaxDayMax'] = 0.0111331639253
        salt2_res['CovDayMaxX0'] = -2.94345317778e-07
        salt2_res['CovDayMaxX1'] = 0.0131008809199
        salt2_res['CovX0X0'] = 7.91767938168e-11
        salt2_res['CovX0X1'] = -7.23852420336e-07
        salt2_res['CovX1X1'] = 0.0258715828973
        # salt2_res['']=
        vparam_names = [self.paramNames['t0'], self.paramNames['color'],
                        self.paramNames['x0'], self.paramNames['x1']]
        covar = np.zeros(shape=(len(vparam_names), len(vparam_names)))

        covar[0, 1] = salt2_res['CovColorDayMax']
        covar[0, 2] = salt2_res['CovDayMaxX0']
        covar[0, 3] = salt2_res['CovDayMaxX1']

        covar[1, 2] = salt2_res['CovColorX0']
        covar[1, 3] = salt2_res['CovColorX1']

        covar[2, 3] = salt2_res['CovX0X1']

        covar = covar+covar.T

        covar[0, 0] = salt2_res['CovDayMaxDayMax']
        covar[1, 1] = salt2_res['CovColorColor']
        covar[2, 2] = salt2_res['CovX0X0']
        covar[3, 3] = salt2_res['CovX1X1']

        # print covar

        params = {}
        params[self.paramNames['t0']] = salt2_res['DayMax']
        params[self.paramNames['color']] = salt2_res['Color']
        params[self.paramNames['x0']] = salt2_res['X0']
        params[self.paramNames['x1']] = salt2_res['X1']

        cov = self.mbCovar(params, covar, vparam_names)
        print(cov)
        if self.interp:
            cov_int = self.mbCovar_int(params, covar, vparam_names)
            print(cov_int)
