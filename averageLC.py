#!/usr/bin/env python
'''
@author: S. Rest
'''

from SNloop import SNloopclass
from statistics import median
from copy import deepcopy
import sigmacut
import sys
import numpy as np
import pandas as pd

class averagelcclass(SNloopclass):
    def __init__(self):
        SNloopclass.__init__(self)

    def calcaveragelc(self,SNindex,controlindex=0,MJDbinsize=1):
        Nclip_max = self.cfg.params['averageLC']['Nclip_max']
        Ngood_min = self.cfg.params['averageLC']['Ngood_min']
        X2norm_max = self.cfg.params['averageLC']['X2norm_max']
        if self.verbose>1:print('Max Nskipped: %d, min Nused: %d, max X2norm: %f' % (Nclip_max,Ngood_min,X2norm_max))

        MJD = int(np.amin(self.lc.t['MJD']))
        MJDmax = int(np.amax(self.lc.t['MJD']))+1

        # clear averagelc table and set dtypes
        self.averagelc.t = self.averagelc.t.iloc[0:0]
        self.averagelc.t['ControlID'] = pd.Series([], dtype=np.int32)
        self.averagelc.t['MJD'] = pd.Series([], dtype=np.float64)
        self.averagelc.t['MJDbin'] = pd.Series([], dtype=np.float64)
        self.averagelc.t['uJy'] = pd.Series([], dtype=np.float64)
        self.averagelc.t['duJy'] = pd.Series([], dtype=np.float64)
        self.averagelc.t['stdev'] = pd.Series([], dtype=np.float64)
        self.averagelc.t['X2norm'] = pd.Series([], dtype=np.float64)
        self.averagelc.t['Nclip'] = pd.Series([], dtype=np.int64)
        self.averagelc.t['Ngood'] = pd.Series([], dtype=np.int64)
        self.averagelc.t['Nexcluded'] = pd.Series([], dtype=np.int64)
        self.averagelc.t['Mask'] = pd.Series([], dtype=np.int32)

        while MJD <= MJDmax:
            # get 4 measurements
            if self.verbose>1: 
                print('MJD range: ',MJD,' to ',MJD+MJDbinsize)
            ix1 = self.lc.ix_inrange(colnames=['MJD'],lowlim=MJD,uplim=MJD+MJDbinsize,exclude_uplim=True)
            if len(ix1)==0:
                if self.verbose>1: 
                    print('Length of MJD range = 0, skipping MJD range...')
                MJD += MJDbinsize
                continue
            
            # get good and ok measurements out of 4
            ix2 = self.lc.ix_unmasked('Mask',maskval=self.flag_c2_bad|self.flag_c0_X2norm|self.flag_c0_uncertainty,indices=ix1)

            # add row to averagelc table
            df = {'ControlID':controlindex,'MJDbin':MJD+0.5*MJDbinsize,'Nclip':0,'Ngood':0,'Nexcluded':len(ix1)-len(ix2),'Mask':0}
            lcaverageindex = self.averagelc.newrow(df)

            # No good measurements?
            if len(ix2)==0:
                flag_array = np.full(len(ix1),self.flag_day_bad)
                self.lc.t.loc[ix1,'Mask'] = np.bitwise_or(self.lc.t.loc[ix1,'Mask'],flag_array)
                self.averagelc.t.loc[lcaverageindex,'Mask'] = int(self.averagelc.t.loc[lcaverageindex,'Mask']) |  self.flag_day_bad
                if self.verbose>1: 
                    print('Length of good and ok measurements = 0, no average flux values')

                MJD += MJDbinsize
                continue
            else:
                if self.verbose>2: 
                    print('Good and ok indices (ix2): ',ix2)

            # sigmacut good and ok measurements
            self.lc.calcaverage_sigmacutloop('uJy',noisecol='duJy',indices=ix2,verbose=1,Nsigma=3.0,median_firstiteration=True)
            fluxstatparams = deepcopy(self.lc.statparams)
            if self.verbose>1: 
                print('Nclip: {}, Ngood: {}, X2norm: {}'.format(fluxstatparams['Nclip'],fluxstatparams['Ngood'],fluxstatparams['X2norm']))
            if fluxstatparams['mean'] is None or len(fluxstatparams['ix_good'])<1:
                if self.verbose>1: 
                    print('Mean uJy is None OR length index good < 1, flagging bad day and skipping MJD range...')
                flag_array = np.full(len(ix1),self.flag_day_bad)
                self.lc.t.loc[ix1,'Mask'] = np.bitwise_or(self.lc.t.loc[ix1,'Mask'],flag_array)
                self.averagelc.t.loc[lcaverageindex,'Mask'] = int(self.averagelc.t.loc[lcaverageindex,'Mask']) | self.flag_day_bad

                MJD += MJDbinsize
                continue

            # get average mjd
            self.lc.calcaverage_sigmacutloop('MJD',noisecol='duJy',indices=fluxstatparams['ix_good'],
                                            verbose=1,Nsigma=0,median_firstiteration=False)
            averagemjd = self.lc.statparams['mean']
            
            # add row to averagelc table
            df = {'MJD':averagemjd,
             self.flux_colname:fluxstatparams['mean'],self.dflux_colname:fluxstatparams['mean_err'],
             'stdev':fluxstatparams['stdev'],'X2norm':fluxstatparams['X2norm'],
             'Nclip':fluxstatparams['Nclip'],'Ngood':fluxstatparams['Ngood'],
             'Mask':0}
            self.averagelc.add2row(lcaverageindex,df)

            # flag clipped sigmacut measurements in original lc
            ix_bad = fluxstatparams['ix_clip']
            if self.verbose>2:
                print('Sigmacut clipped indices (ix_bad): ',ix_bad)
            if len(ix_bad) > 0:
                flag_array = np.full(len(ix_bad),self.flag_daysigma)
                self.lc.t.loc[ix_bad,'Mask'] = np.bitwise_or(self.lc.t.loc[ix_bad,'Mask'],flag_array)

            badflag = 0
            if len(ix2)<=2:
                # flag in original and averaged lc
                flag_array = np.full(len(ix2),self.flag_daysmallnumber)
                self.lc.t.loc[ix2,'Mask'] = np.bitwise_or(self.lc.t.loc[ix2,'Mask'],flag_array)
                self.averagelc.t.loc[lcaverageindex,'Mask'] = int(self.averagelc.t.loc[lcaverageindex,'Mask']) | self.flag_daysmallnumber
            else:
                # check sigmacut stats and, if badflag, flag as bad day in original and averaged lc
                if fluxstatparams['Ngood'] < Ngood_min: 
                    badflag = 1
                if fluxstatparams['Nclip'] > Nclip_max: 
                    badflag = 1
                if not(fluxstatparams['X2norm'] is None) and (fluxstatparams['X2norm'] > X2norm_max): 
                    badflag = 1
            if badflag == 1:
                if self.verbose>1:
                    print('# Flagged as bad day!')
                flag_array = np.full(len(ix1),self.flag_day_bad)
                self.lc.t.loc[ix1,'Mask'] = np.bitwise_or(self.lc.t.loc[ix1,'Mask'],flag_array)
                self.averagelc.t.loc[lcaverageindex,'Mask'] = int(self.averagelc.t.loc[lcaverageindex,'Mask']) | self.flag_day_bad
            else:
                if self.verbose>1:
                    print('# Flagged as good day!')

            MJD += MJDbinsize

        for col in ['ControlID','Nclip','Ngood','Nexcluded','Mask']: 
            self.averagelc.t[col] = self.averagelc.t[col].astype(np.int32)

        avglcfilename = self.lcbasename(SNindex=SNindex,filt=self.filt,controlindex=controlindex,MJDbinsize=MJDbinsize)+'.txt'
        if self.verbose>1: self.averagelc.write()
        self.averagelc.write(avglcfilename,overwrite=True,verbose=True)

    def averagelcloop(self,SNindex):
        print('###################################\nAveraging LCs...\n###################################')
        
        # loop through SN and control lcs
        for controlindex in range(0,len(self.RADECtable.t)):    
            # stop loop if only SN should be done
            if (not self.cfg.params['averageLC']['apply2offsets']) and (controlindex>0):
                break
            # load control lc
            if self.verbose: print('Averaging lc for controlID %d' % self.RADECtable.t['ControlID'][controlindex])
            self.load_lc(SNindex,controlindex=controlindex,filt=self.filt)
            if len(self.lc.t)==0: 
                print('WARNING!!! no data in lc')
                return(1)
            # get sigmacut info and flag for 4-day bins
            if not(args.MJDbinsize is None):
                MJDbinsize = args.MJDbinsize
                self.calcaveragelc(SNindex,controlindex=controlindex,MJDbinsize=MJDbinsize)
            else:
                self.calcaveragelc(SNindex,controlindex=controlindex)
            # save lc
            self.save_lc(SNindex=SNindex,controlindex=controlindex,filt=self.filt,overwrite=True)
        print('### Averaging LCs done')

if __name__ == '__main__':

    averagelc = averagelcclass()
    parser = averagelc.define_options()
    args = parser.parse_args()

    SNindexlist = averagelc.initialize(args)

    for SNindex in SNindexlist:
        averagelc.loadRADEClist(SNindex,filt=averagelc.filt)
        averagelc.averagelcloop(SNindex)
