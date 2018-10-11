import yaml
import argparse
import time
import h5py
from astropy.table import Table
from importlib import import_module
from SN_Telescope import Telescope
import os
import numpy as np

class Fit_All:
    def __init__(self,telescope,output_config,display_lc,fitter_config,prodid):

        module = import_module(fitter_config['name'])
        self.fitter = module.Fit_LC(model = fitter_config['model'],version=fitter_config['version'],telescope=telescope,display=display_lc)

        if output_config['save']:
            self._prepare_Save(output_config['directory'], prodid)

        self.val =[]
    
    def __call__(self,lc):
        
        names,val = self.fitter(lc)

        self.names = names
        self.val.append(tuple(val))

    def _prepare_Save(self, outdir, prodid):

        if not os.path.exists(outdir):
            print('Creating output directory', outdir)
            os.makedirs(outdir)
        
        self.fit_out = outdir+'/Fit_'+prodid+'.hdf5'
        if os.path.exists(self.fit_out):
            os.remove(self.fit_out)
            
    def Finish(self):
        
        # dump the results in a hdf5 file
        Table(rows = self.val, names=self.names).write(self.fit_out, 'fit_lc', compression=True)

        
parser = argparse.ArgumentParser(
    description='Run a LC fitter from a configuration file')
parser.add_argument('config_filename',
                    help='Configuration file in YAML format.')


def run(config_filename):
    # YAML input file.
    config = yaml.load(open(config_filename))
    print(config)
    
    # load telescope
    tel_par = config['Instrument']

    # this is for output
    save_status = config['Output']['save']
    outdir = config['Output']['directory']
    prodid = config['ProductionID']

    simu_name = config['Simulations']['dirname']+'/Simu_'+prodid+'.hdf5'
    lc_name = config['Simulations']['dirname']+'/LC_'+prodid+'.hdf5'

    telescope=Telescope(atmos=True,airmass=1.2)
    
    Fit = Fit_All(telescope,config['Output'],config['Display'],config['Fitter'],prodid)

    # Loop on the simu_file to grab simulated LCs
    
    f = h5py.File(simu_name, 'r')
    print(f.keys())
    for i, key in enumerate(f.keys()):
        simu = Table.read(simu_name, path=key)
        print(len(simu), simu.dtype)
        for ko,info in enumerate(simu):
            print(info)
            if info['n_lc_points'] > 0:
                print(info['n_lc_points'],info['id_hdf5'])
                lc =  Table.read(lc_name, path='lc_'+str(info['id_hdf5']))
                Fit(lc)

            if ko > 4:
                break

    Fit.Finish()
    
def main(args):
    print('running')
    time_ref = time.time()
    run(args.config_filename)
    print('Time', time.time()-time_ref)


if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
