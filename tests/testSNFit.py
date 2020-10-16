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
import yaml

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


class TestSNFit(unittest.TestCase):
    def testSNFit(self):

        # get data
        # two files requested:
        # Simu_*: astropy table with a list of SN simulation parameters
        # LC_*: astropy tables with LC

        # grab LCs
        #prodid = 'sncosmo_Fake_Fake_DESC_seas_-1_-2.0_0.2.hdf5'
        prodid = 'DD_descddf_v1.4_10yrs_COSMOS_0'
        simu_name = 'Simu_{}.hdf5'.format(prodid)
        lc_name = 'LC_{}.hdf5'.format(prodid)

        # get these files from server if necessary
        getFile('unittests', simu_name)
        getFile('unittests', lc_name)
        getRefDir('SALT2_Files')

        # get configuration file
        with open('param_fit_test.yaml') as file:
            # The FullLoader parameter handles the conversion from YAML
            # scalar values to Python the dictionary format
            conf= yaml.load(file, Loader=yaml.FullLoader)

        conf['Simulations']['dirname']='.'
        conf['Simulations']['prodid']= prodid
        conf['mbcov']['estimate'] = 1
        print(conf)
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

        print('number of LC to fit',len(simul))
        simul['z'] = np.round(simul['z'], 2)
        ll = np.round(list(np.arange(0.1, 0.9, 0.1)), 2)
        simul = simul[np.in1d(simul['z'], ll)]

        res = Table()
        for simu in simul:
            lc = None
            lc = Table.read(lc_name, path='lc_{}'.format(simu['index_hdf5']))
            lc.convert_bytestring_to_unicode()
            #print('here lc',lc)
            resfit = fit(lc)
            res = vstack([res, resfit])

        idx = res['z'] < 0.8
        sel = res[idx]

        
        ref_fit =  'data_tests/ref_fit_{}.hdf5'.format(prodid)

        # save the reference file if it does not exist
        if not os.path.exists(ref_fit):
            sel.write(ref_fit,'lc_fits', append=True, compression=True)

        # load fit reference values
        fit_lc_ref = Table.read(ref_fit,path='lc_fits')
        
        keychecks = ['z', 'Cov_colorcolor', 'mbfit', 'mb_recalc', 'sigma_mu']

        for key in keychecks:
            assert(np.isclose(fit_lc_ref[key].tolist(), sel[key].tolist()).all())

        """
        print(res.columns)
        for key in keychecks:
            print('dictRef[\'{}\']='.format(key), res[key].tolist())
        print(testoo)
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
        """
        # cleaning directory
        for fi in [simu_name, lc_name]:
            if os.path.isfile(fi):
                os.system('rm {}'.format(fi))
        for ddir in ['SALT2_Files','Output_Fit']:
            if os.path.isdir(ddir):
                os.system('rm -rf {}'.format(ddir))

#fit_snfit = TestSNFit
# fit_sncosmo = TestSNFitcosmo

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main(verbosity=5)


# unittest.main(verbosity=5)
