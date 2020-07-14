from builtins import zip
import numpy as np
import unittest
import lsst.utils.tests
from sn_fitter.fit_sn_cosmo import Fit_LC
from sn_tools.sn_telescope import Telescope
from sn_fit.process_fit import Fitting
from sn_fit.mbcov import MbCov
import os
import h5py
from astropy.table import Table, vstack

main_repo = 'https://me.lsst.eu/gris/DESC_SN_pipeline'


def getFile(refdir, fname):
    fullname = '{}/{}/{}'.format(main_repo, refdir, fname)

    # check whether the file is available; if not-> get it!
    if not os.path.isfile(fname):
        print('wget path:', fullname)
        cmd = 'wget --no-clobber --no-verbose {}'.format(fullname)
        os.system(cmd)


def getRefDir(dirname):
    fullname = '{}/{}'.format(main_repo, dirname)

    if not os.path.exists(dirname):
        print('wget path:', fullname)
        cmd = 'wget --no-verbose --recursive {} --directory-prefix={} --no-clobber --no-parent -nH --cut-dirs=3 -R \'index.html*\''.format(
            fullname+'/', dirname)
        os.system(cmd)


def getconfig(prodid, simu_name, lc_name, saveData=False):

    conf = {}
    conf['ProductionID'] = prodid

    conf['Instrument'] = {}
    conf['Instrument']['name'] = 'LSST'
    conf['Instrument']['throughput_dir'] = 'LSST_THROUGHPUTS_BASELINE'
    conf['Instrument']['atmos_dir'] = 'THROUGHPUTS_DIR'
    conf['Instrument']['airmass'] = 1.1
    conf['Instrument']['atmos'] = True
    conf['Instrument']['aerosol'] = False

    conf['Simulations'] = {}
    conf['Simulations']['filename'] = simu_name
    conf['Simulations']['dirname'] = '.'

    conf['Fitter'] = {}
    conf['Fitter']['name'] = 'sn_fitter.fit_sn_cosmo'
    conf['Fitter']['model'] = 'salt2-extended'
    conf['Fitter']['version'] = 1.0
    conf['Fitter']['covmb'] = 1

    conf['Display'] = 0

    conf['Output'] = {}
    conf['Output']['directory'] = '.'
    conf['Output']['save'] = saveData

    conf['Multiprocessing'] = {}
    conf['Multiprocessing']['nproc'] = 1

    return conf


class TestSNFitcosmo(unittest.TestCase):

    def testFit(self):

        # fit instance
        telescope = Telescope(airmass=1.2)
        _fit = Fit_LC(telescope=telescope)
        print('test')

        # grab LCs
        simuLC = 'LC_sncosmo_Fake_Fake_DESC_seas_-1_-2.0_0.2.hdf5'

        # get this file from server if necessary
        getFile('unittests', simuLC)

        # open LC file from simulation
        fFile = h5py.File(simuLC, 'r')
        keys = list(fFile.keys())

        # fit LC and store results in data - an astropy Table
        data = Table()
        for key in keys:
            # get LC
            lcsimu = Table.read(fFile, path=key)
            # print(lcsimu.meta)
            # fit LC
            resfit = _fit(lcsimu)
            data = vstack([data, resfit])

        # select a set of results to make the test
        idx = data['z'] < 0.9
        sel = data[idx]
        sel['z'] = np.round(sel['z'], 2)
        ll = np.round(list(np.arange(0.1, 0.9, 0.02)), 2)
        selb = sel[np.in1d(sel['z'], ll)]

        """
        for col in ['z', 'Cov_colorcolor', 'mbfit']:
            print('dictRef[\'{}\']='.format(col), selb[col].tolist())
        """
        dictRef = {}
        dictRef['z'] = [0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46,
                        0.48, 0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88]

        dictRef['Cov_colorcolor'] = [7.847411468307449e-07, 1.2926429633578325e-06, 2.182689313770883e-06, 3.4064682540206656e-06, 5.094755132641675e-06, 7.418516965058399e-06, 1.0449342672175223e-05, 1.4530838906902605e-05, 1.9838442223282924e-05, 5.488615265833384e-05, 6.79687434393572e-05, 8.154211056708287e-05, 9.585755607011092e-05, 0.00011249971780490519, 0.0001332468082708582, 0.0001583244598719656, 0.00018974040864400726, 0.00023307069810571821, 0.0002929759524889111,
                                     0.00036750366827063667, 0.0004505526493535665, 0.0005418476331991809, 0.0006556090670946579, 0.0007993634180891086, 0.0009706867250971499, 0.0011684978183072229, 0.0013961226214366633, 0.0054156739068237205, 0.006306687006515156, 0.007213958623759652, 0.00819136400938683, 0.00939627974924506, 0.010844715095289384, 0.01246028538030565, 0.01424202236460317, 0.016338703851064878, 0.018937841331963138, 0.022182703639333848, 0.025822864238779813, 0.0]
        dictRef['mbfit'] = [20.006435534930837, 20.430453337465664, 20.79137803592443, 21.106957852340685, 21.388860253425335, 21.643240843203596, 21.875275946771513, 22.088544601198045, 22.28621053426997, 22.470159114763483, 22.642592920968404, 22.804801863419183, 22.958205505175076, 23.103367585657537, 23.241190478068184, 23.372313980135566, 23.497345309887262, 23.616965168330655, 23.73171763384528,
                            23.842071730391293, 23.948421905233754, 24.051095990149463, 24.150012532855033, 24.245404432237017, 24.337458733240545, 24.426662484018227, 24.513165417734825, 24.597371243206602, 24.679404190398955, 24.758758483003128, 24.83545473841958, 24.909948797457137, 24.982528069838313, 25.053561604365488, 25.12335255171979, 25.191183726268534, 25.257131587610104, 25.32151390441974, 25.384753330779187, -1.0]

        # check results here
        for key, vals in dictRef.items():
            assert(np.isclose(vals, selb[key].tolist()).all())

        # cleaning directory
        if os.path.isfile(simuLC):
            os.system('rm {}'.format(simuLC))


class TestSNFit(unittest.TestCase):
    def testSNFit(self):

        # get data
        # two files requested:
        # Simu_*: astropy table with a list of SN simulation parameters
        # LC_*: astropy tables with LC

        # grab LCs
        prodid = 'sncosmo_Fake_Fake_DESC_seas_-1_-2.0_0.2.hdf5'
        simu_name = 'Simu_{}'.format(prodid)
        lc_name = 'LC_{}'.format(prodid)

        # get these files from server if necessary
        getFile('unittests', simu_name)
        getFile('unittests', lc_name)
        getRefDir('SALT2_Files')

        # build configuration file
        conf = getconfig(prodid, simu_name, lc_name)

        # covmb
        covmb = MbCov('SALT2_Files', paramNames=dict(
            zip(['x0', 'x1', 'color'], ['x0', 'x1', 'c'])))

        # fit instance
        fit = Fitting(conf, covmb=covmb)

        # getting the simu file
        f = h5py.File(simu_name, 'r')
        print(f.keys())
        # reading the simu file
        for i, key in enumerate(f.keys()):
            simul = Table.read(simu_name, path=key)

        simul['z'] = np.round(simul['z'], 2)
        ll = np.round(list(np.arange(0.1, 0.9, 0.1)), 2)
        simul = simul[np.in1d(simul['z'], ll)]

        res = Table()
        for simu in simul:
            lc = None
            lc = Table.read(lc_name, path='lc_{}'.format(simu['id_hdf5']))
            resfit = fit(lc)
            res = vstack([res, resfit])

        idx = res['z'] < 0.8
        sel = res[idx]

        keychecks = ['z', 'Cov_colorcolor', 'mbfit', 'mb_recalc', 'sigma_mu']

        """
        for key in keychecks:
            print('dictRef[\'{}\']='.format(key), res[key].tolist())
        """
        dictRef = {}

        dictRef['z'] = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
        dictRef['Cov_colorcolor'] = [7.84149720548657e-07, 7.411367293053516e-06, 6.793756102580619e-05,
                                     0.00015816525343875836, 0.0004502363398936092, 0.001167662290074013, 0.008186441880889383, 0.01633751529759251]
        dictRef['mbfit'] = [20.005795111491146, 21.64206956980012, 22.64256913193122, 23.372165710405174,
                            23.948145769722828, 24.42636698333764, 24.835436823080904, 25.19122934728612]
        dictRef['mb_recalc'] = [20.02731344585068, 21.663516325567763, 22.66404477842765,
                                23.39363924396809, 23.969613792664823, 24.4478277942898, 24.85691319019894, 25.212677718902857]
        dictRef['sigma_mu'] = [0.0025744440149275205, 0.008148971048222674, 0.024267233165841652,
                               0.03934245428396615, 0.06908538342541443, 0.10804914411477154, 0.288358020033676, 0.4019877529973541]

        for key in keychecks:
            assert(np.isclose(dictRef[key], res[key].tolist()).all())

        # cleaning directory
        for fi in [simu_name, lc_name]:
            if os.path.isfile(fi):
                os.system('rm {}'.format(fi))
        for ddir in ['SALT2_Files']:
            if os.path.isdir(ddir):
                os.system('rm -rf {}'.format(ddir))


fit_snfit = TestSNFit
# fit_sncosmo = TestSNFitcosmo

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main(verbosity=5)


# unittest.main(verbosity=5)
