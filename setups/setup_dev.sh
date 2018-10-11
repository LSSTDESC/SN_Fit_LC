#!/bin/bash

#setup lsst_sims

#source /global/common/software/lsst/cori-haswell-gcc/stack/setup_w_2018_13-sims_2_7_0.sh
source /cvmfs/sw.lsst.eu/linux-x86_64/lsst_sims/sims_2_8_0//loadLSST.bash
setup lsst_sims

export PYTHONPATH=${PWD}/SN_Fit_LC/Fit_SNCosmo:$PYTHONPATH
export PYTHONPATH=${PWD}/SN_Utils/Utils:$PYTHONPATH
export SN_UTILS_DIR=${PWD}/SN_Utils
export SALT2_DIR=${PWD}/SN_Utils/SALT2_Files

#checking whether hdf5 is accessible localy or not
lib='h5py'
thedir=${PWD}/lib/python3.6/site-packages/
echo $thedir
if [ -d ${thedir}$lib ]
then
    echo $lib 'already installed -> updating PYTHONPATH'
else
    echo $lib 'not installed -> installing with pip'
    pip install --prefix=${PWD} ${lib}==2.7.1
fi

export PYTHONPATH=${thedir}:$PYTHONPATH
