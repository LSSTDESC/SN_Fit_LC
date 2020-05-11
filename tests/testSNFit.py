from builtins import zip
import numpy as np
import unittest
import lsst.utils.tests
from sn_fitter.fit_sncosmo import Fit_LC
from sn_tools.sn_telescope import Telescope
from sn_fit.sn_fit import Fitting
import os
import h5py
from astropy.table import Table, vstack

main_repo = 'https://me.lsst.eu/gris/Reference_Files'


def getFile(refdir, fname):
    fullname = '{}/{}/{}'.format(main_repo, refdir, fname)

    # check whether the file is available; if not-> get it!
    if not os.path.isfile(fname):
        print('wget path:', fullname)
        cmd = 'wget --no-clobber --no-verbose {}'.format(fullname)
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
    conf['Fitter']['name'] = 'sn_fitter.fit_sncosmo'
    conf['Fitter']['model'] = 'salt2-extended'
    conf['Fitter']['version'] = 1.0
    conf['Fitter']['covmb'] = 0

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

        # get this file from server if necessary
        getFile('unittests', simu_name)
        getFile('unittests', lc_name)

        # build configuration file
        conf = getconfig(prodid, simu_name, lc_name)

        # fit instance
        fit = Fitting(conf)

        # getting the simu file
        f = h5py.File(simu_name, 'r')
        print(f.keys())
        # reading the simu file
        for i, key in enumerate(f.keys()):
            simul = Table.read(simu_name, path=key)

        res = Table()
        for simu in simul:
            lc = None
            lc = Table.read(lc_name, path='lc_{}'.format(simu['id_hdf5']))
            resfit = fit(lc)
            res = vstack([res, resfit])

        idx = res['z'] < 0.8
        sel = res[idx]

        keychecks = ['z', 'Cov_colorcolor']
        """
        for key in keychecks:
            print('dictRef[\'{}\']='.format(key), res[key].tolist())
        """

        dictRef = {}
        dictRef['z'] = [0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54,
                        0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7000000000000001, 0.72, 0.74, 0.76, 0.78, 0.8, 0.8200000000000001, 0.84, 0.86, 0.88, 0.9, 0.92, 0.9400000000000001, 0.96, 0.98, 1.0]

        dictRef['Cov_colorcolor'] = [5.0188024838012865e-09, 2.079459024964159e-08, 9.572515257546943e-08, 2.5579974791635596e-07, 5.381968114366316e-07, 7.84149720548657e-07, 1.2904117589079427e-06, 2.1835747872003875e-06, 3.4031247974842205e-06, 5.089387982701995e-06, 7.411367293053516e-06, 1.0440058224771452e-05, 1.4516057013468932e-05, 1.9820846580639528e-05, 5.4848824107982635e-05, 6.793756102580619e-05, 8.14878997971749e-05, 9.578321673721968e-05, 0.00011241858051157339, 0.00013313850190210713, 0.00015816525343875836, 0.00018953338140741282,
                                     0.00023279260332887108, 0.0002926926540771169, 0.00036724081727145983, 0.0004502363398936092, 0.0005414746852036682, 0.0006552500455930502, 0.0007989447815061136, 0.0009701086108708995, 0.001167662290074013, 0.0013948261230415425, 0.005411630796879028, 0.006302585916700493, 0.007209322714217017, 0.008186441880889383, 0.009392125822218945, 0.010841303440276704, 0.012457675161195515, 0.014240225551167742, 0.01633751529759251, 0.018936075819973055, 0.02217601051431927, 0.02581103860697213, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        for key in keychecks:
            assert(np.isclose(dictRef[key], res[key].tolist()).all())

        # cleaning directory
        for fi in [simu_name, lc_name]:
            if os.path.isfile(fi):
                os.system('rm {}'.format(fi))


fit_snfit = TestSNFit
# fit_sncosmo = TestSNFitcosmo

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main(verbosity=5)


# unittest.main(verbosity=5)
