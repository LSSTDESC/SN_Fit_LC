#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 14:21:15 2025

@author: philippe.gris@clermont.in2p3.fr
"""
from sn_tools.sn_utils import load_config


class FitWrapper:
    def __init__(self, config):
        """
        Class to fit a set of light curves

        Parameters
        ----------
        config_fit : dict
            parameters fot

        Returns
        -------
        None.

        """
        from sn_fit.process_fit import Fitting

        # Fit instance
        #config = load_config(yaml_config_fit)

        self.fit = Fitting(config)
        self.nproc = config['MultiprocessingFit']['nproc']

        self.saveData = config['OutputFit']['save']

        self.outDir = config['OutputFit']['directory']

        self.prodid = config['Simulations']['prodid']

        self.ccolref = []

        if self.saveData:
            from sn_tools.sn_io import checkDir
            checkDir(self.outDir)
            outFile = 'SN_{}.hdf5'.format(self.prodid)
            self.outName = '{}/{}'.format(self.outDir, outFile)
            # check wether this file already exist and remove it
            import os
            if os.path.isfile(self.outName):
                os.system('rm {}'.format(self.outName))

    def __call__(self, lc_list, remove_sat=False):
        """
        Main fit method using multiprocessing

        Parameters
        ----------
        lc_list : list(lc)
            List of light curves to fit.
        remove_sat : bool, optional
            To remove saturated points. The default is False.

        Returns
        -------
        res : pandas df
            output results.

        """

        if self.nproc > 1:
            res = self.fit.fit_multiproc(lc_list, remove_sat, self.nproc)
        else:
            params = {}
            params['remove_sat'] = remove_sat
            res = self.fit.fit_lcs(lc_list,params)

        return res

    def __call__deprecated(self, lc_list, remove_sat=False):
        """
        Method to fit light curves

        Parameters
        ----------
        lc_list : list(astropy table)
            LC to fit
        remove_sat : bool, optional
            To remove saturated fluxes. The default is False.

        Returns
        -------
        None.

        """
        """
        from astropy.table import Table, vstack
        res = Table()
        for lc in lc_list:
            lc.convert_bytestring_to_unicode()
            resfit = self.fit(lc)
            if resfit is not None:
                res = vstack([res, resfit])

        return res
        """
        # from sn_tools.sn_utils import multiproc
        params = {}
        params['remove_sat'] = remove_sat

        res = self.multiproc(lc_list, params, self.fit_lcs, self.nproc)

        return res

    def fit_lcs_deprecated(self, lc_list, params, j=0, output_q=None):
        """
        Method to fit LCs

        Parameters
        ----------
        lc_list : list(astropy table)
            light-curves to fit.
        params : dict
            parameters.
        j : int, optional
            Tag for multiprocessing. The default is 0.
        output_q : multiprocessing queue, optional
            queue managing multiprocessing run. The default is None.

        Returns
        -------
        astropytable
            Result of the fit.

        """

        from astropy.table import Table, vstack
        res = Table()
        # print('processing fit', j)

        for lc in lc_list:
            lc.convert_bytestring_to_unicode()
            resfit = self.fit(lc, params)
            if resfit is not None:
                resfit = self.check_correct(resfit)
                res = vstack([res, resfit])

        if output_q is not None:
            return output_q.put({j: res})
        else:
            return res

    def check_correct(self, sn):
        """
        Method to correct for Cov_xy col names

        Parameters
        ----------
        sn : astropy Table
            Data to process.

        Returns
        -------
        sn : astropy Table
            Processed data

        """

        varlist = ['z', 't0', 'x0', 'x1', 'color']

        if 'Cov_zz' not in sn.columns:
            varlist = ['t0', 'x0', 'x1', 'color']

        for i, namea in enumerate(varlist):
            for j, nameb in enumerate(varlist):
                if j >= i:
                    vva = 'Cov_{}{}'.format(namea, nameb)
                    vvb = 'Cov_{}{}'.format(nameb, namea)
                    if vva not in sn.columns:
                        sn.rename_column(vvb, vva)

        return sn

    def dump(self, fitlc):
        """


        Parameters
        ----------
        fitlc : pandas df
            data to dump

        Returns
        -------
        None.

        """
        """
        if self.outName != '':
            keyhdf = '{}'.format(int(sn['healpixID'].mean()))
            sn.write(self.outName, keyhdf, append=True, compression=True)
        """
        import pandas as pd
        if self.saveData:
            fitlc.convert_bytestring_to_unicode()
            df = pd.DataFrame(fitlc.to_pandas())

            if 'selected' in df.columns:
                df = df.drop(columns=['selected'])

            if not self.ccolref:
                self.ccolref = df.columns.to_list()
            else:
                df = df.reindex(columns=self.ccolref)

            """
            print('chisq', df['chisq'])
            for vv in df.columns:
                print(vv, df[vv].dtype)
            """
            """
            for vv in self.ccolref:
                print(vv, df[vv].dtype)
            """
            """
            cols = ['sn_type', 'sn_model', 'sn_version', 'fitstatus']

            print(df['fitstatus'].unique(), df['SNID'].unique())
            """

            df.to_hdf(self.outName, key='SN', append=True)
