import yaml
import argparse
import time
import h5py
from astropy.table import Table
from importlib import import_module
from SN_Telescope import Telescope
import os
import numpy as np
import multiprocessing

def Multiproc(simu_name='', lc_name='',Fit=None,nproc=1):

    f = h5py.File(simu_name, 'r')
    print(f.keys())
    for i, key in enumerate(f.keys()):
        simu = Table.read(simu_name, path=key)

    nlc = len(simu)
    print('total number of LC',nlc)
    nlc = 32
    batch = range(0,nlc,nproc)
    batch = np.append(batch,nlc)
    
    names = None
    val =[]
    for i in range(len(batch)-1):
        result_queue = multiprocessing.Queue()
        
        ida = batch[i]
        idb = batch[i+1]
    
        for j in range(ida,idb):
            lc =  Table.read(lc_name, path='lc_'+str(simu[ida]['id_hdf5']))
            p=multiprocessing.Process(name='Subprocess-'+str(j),target=Fit,args=(lc,j,result_queue))
            p.start()

        
        resultdict = {}
        for j in range(ida,idb):
            resultdict.update(result_queue.get())
        
        for p in multiprocessing.active_children():
            p.join()

        for j in range(ida,idb):
            names = resultdict[j][0]
            val.append(tuple(resultdict[j][1]))
            
    Fit.Finish(names,val)

class Fit_All:
    def __init__(self,telescope,output_config,display_lc,fitter_config,prodid):

        module = import_module(fitter_config['name'])
        self.fitter = module.Fit_LC(model = fitter_config['model'],version=fitter_config['version'],telescope=telescope,display=display_lc)

        if output_config['save']:
            self._prepare_Save(output_config['directory'], prodid)

        self.val =[]
    
    def __call__(self,lc,j=-1,output_q=None):
        
        names,val = self.fitter(lc)

        if output_q is not None:
            output_q.put({j : (list(names),list(val))})
            
    def _prepare_Save(self, outdir, prodid):

        if not os.path.exists(outdir):
            print('Creating output directory', outdir)
            os.makedirs(outdir)
        
        self.fit_out = outdir+'/Fit_'+prodid+'.hdf5'
        if os.path.exists(self.fit_out):
            os.remove(self.fit_out)
            
    def Finish(self,names,val):
        
        # dump the results in a hdf5 file
        Table(rows = val, names=names).write(self.fit_out, 'fit_lc', compression=True)
        
        
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
    Multiproc(simu_name=simu_name,lc_name=lc_name,Fit=Fit,nproc=config['Multiprocessing']['nproc'])
  
    
def main(args):
    print('running')
    time_ref = time.time()
    run(args.config_filename)
    print('Time', time.time()-time_ref)


if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
