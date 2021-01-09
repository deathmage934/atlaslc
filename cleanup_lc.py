#!/usr/bin/env python
'''
@author: S. Rest
'''

from SNloop import SNloopclass
from statistics import median
from astropy.table import Table, Column, MaskedColumn
from copy import deepcopy
import sigmacut
import sys,re
import numpy as np
import pandas as pd
from pdastro import pdastrostatsclass

class cleanuplcclass(SNloopclass):
    def __init__(self):
        SNloopclass.__init__(self)

    def c0_PSF_uncertainty_cut(self,N_dflux_max):
        if N_dflux_max is None:
            if self.verbose>1: print('No cut on the detection uncertainties...')
            return(0)
        
        dflux_max = N_dflux_max * median(self.lc.t[self.dflux_colname])
        if self.verbose: print('Flagging all measurements with %s bigger than %f...' % (self.dflux_colname, dflux_max))

        # get indices
        a_indices = self.lc.ix_inrange(self.dflux_colname,dflux_max,None)
        #a_indices = np.where(self.lc.t[self.dflux_colname]>a)
        #a_indices = list(a_indices[0])
        #if self.verbose>1: print('Indices: ',a_indices)
        if len(a_indices)>0:
            #if self.verbose: print('# %s above %i: %i/%i' % (self.dflux_colname, a, len(self.lc.t.loc[a_indices,self.dflux_colname]),len(self.lc.t[self.dflux_colname])))
            if self.verbose: print('# %s above %f: %i/%i' % (self.dflux_colname,dflux_max,len(a_indices),len(self.lc.getindices())))
        else:
            if self.verbose: print('# No measurements flagged!')
            return(0)

        # update 'Mask' column
        flag_c0_uncertainty = np.full(self.lc.t.loc[a_indices,'Mask'].shape, self.flag_c0_uncertainty|self.flag_c0_bad)
        self.lc.t.loc[a_indices,'Mask'] = np.bitwise_or(self.lc.t.loc[a_indices,'Mask'],flag_c0_uncertainty) 
        return(0)
        
    def c0_PSF_X2norm_cut(self,X2norm_max):
        if X2norm_max is None:
            if self.verbose>1: print('No cut on the X2norm...')
            return(0)
        
        if self.verbose: print('Flagging all measurements with chi/N bigger than %f...' % X2norm_max)
        
        # get indices
        a_indices = self.lc.ix_inrange('chi/N',X2norm_max,None)
        #a_indices = np.where(self.lc.t['chi/N']>c0max_X2norm)
        #a_indices = list(a_indices[0])
        #if self.verbose>1: print('Indices: ',a_indices) 
        
        if len(a_indices)>0:
            if self.verbose: print('# chi/N above %f: %i/%i' % (X2norm_max,len(a_indices),len(self.lc.getindices())))
        else:
            if self.verbose: print('# No measurements flagged!')
            return(0)

        # update 'Mask' column
        flag_c0_X2norm = np.full(self.lc.t.loc[a_indices,'Mask'].shape, self.flag_c0_X2norm|self.flag_c0_bad)
        
        ### CHANGE ME BACK!!!!!
        #flag_c0_X2norm = np.full(self.lc.t.loc[a_indices,'Mask'].shape, self.flag_c0_X2norm|self.flag_daysigma)
        
        self.lc.t.loc[a_indices,'Mask'] = np.bitwise_or(self.lc.t.loc[a_indices,'Mask'], flag_c0_X2norm)
        return(0)
    
    def make_c0_cuts(self, SNindex, prepare_c1c2_cuts=False):
        if prepare_c1c2_cuts:
            # load main SN LC
            self.load_lc(SNindex, filt=self.filt, controlindex=0, MJDbinsize=None)
            
            # construct arrays for control lc data
            MJD_SN = self.lc.t['MJD']

            N_lc = len(self.RADECtable.t)
            N_MJD = len(self.lc.t['MJD'])

            uJy = np.full((N_lc,N_MJD),np.nan)
            duJy = np.full((N_lc,N_MJD),np.nan)
            Mask = np.full((N_lc,N_MJD),0,dtype=np.int32)
        else:
            MJD_SN = uJy = duJy = Mask = None
            
        for controlindex in self.RADECtable.getindices():
            # load lc
            self.load_lc(SNindex, filt=self.filt, controlindex=controlindex, MJDbinsize=None)
            
            if len(self.lc.t) == 0:
                raise RuntimeError('No data in control lc with index {}'.format(controlindex))
                
            # set Mask to 0
            self.lc.t['Mask'] = 0
            
            # get rid of all old statistic columns
            dropcols=[]
            if 'Noffsetlc' in self.lc.t.columns: dropcols.append('Noffsetlc')
            for col in self.lc.t.columns:
                if re.search('^c\d_',col): dropcols.append(col)
                if re.search('^o\d_',col): dropcols.append(col)
            if len(dropcols)>0: self.lc.t.drop(columns=dropcols,inplace=True)
            
            # make c0 cuts
            self.c0_PSF_uncertainty_cut(self.cfg.params['cleanlc']['cut0']['N_dflux_max'])
            self.c0_PSF_X2norm_cut(self.cfg.params['cleanlc']['cut0']['PSF_X2norm_max'])
            
            #if self.cfg.params['cleanlc']['cut0']['PSF_uncertainty']['apply'] is True:
            #    self.c0_PSF_uncertainty_cut()
            #if self.cfg.params['cleanlc']['cut0']['PSF_X2norm']['apply'] is True:
            #    self.c0_PSF_X2norm_cut()

            # for later c1/c2 cuts: prepare the flux,dflux,mask arrays of all light curves
            if prepare_c1c2_cuts and controlindex!=0:
                # make sure MJD_SN is the same as self.lc.t['MJD'], then fill array of control lc uJy, duJy, Mask
                if (len(self.lc.t) != N_MJD) or (np.array_equal(MJD_SN, self.lc.t['MJD']) is False):
                    print('ERROR!!!')
                    print('SN MJD:',MJD_SN)
                    print('Control LC MJD:',self.lc.t['MJD'])
                    raise RuntimeError('SN lc not equal to control lc for controlindex {}! Please run verifyMJD.py'.format(controlindex))
                else:
                    uJy[controlindex,:] = self.lc.t[self.flux_colname]
                    duJy[controlindex,:] = self.lc.t[self.dflux_colname]
                    Mask[controlindex,:] = self.lc.t['Mask']

            # Save the light curve with cuts in Mask
            self.save_lc(SNindex=SNindex,filt=self.filt,overwrite=True,controlindex=controlindex)

        return(MJD_SN,uJy,duJy,Mask)
        

    def calc_c1c2_stats(self,uJy,duJy,mask):
        # load main lc
        self.load_lc(SNindex,controlindex=0)

        N_MJD = uJy.shape[-1]

        c1_param2columnmapping = self.lc.intializecols4statparams(prefix='c1_',format4outvals='{:.2f}',parammapping={'Ngood':'Nvalid'},skipparams=['converged','i','Nclip','Nmask'])
        c2_param2columnmapping = self.lc.intializecols4statparams(prefix='c2_',format4outvals='{:.2f}',skipparams=['converged','i'])
 
        for index in range(N_MJD):
            pda4MJD = pdastrostatsclass()
            pda4MJD.t['uJy']=uJy[1:,index]
            pda4MJD.t['duJy']=duJy[1:,index]
            pda4MJD.t['Mask']=np.bitwise_and(mask[1:,index],self.flag_c0_uncertainty|self.flag_c0_X2norm)
            
            # c1 stats ...
            pda4MJD.calcaverage_sigmacutloop('uJy',noisecol='duJy',verbose=1,Nsigma=0.0,median_firstiteration=False)
            # ... and save them into the table
            self.lc.statresults2table(pda4MJD,c1_param2columnmapping,destindex=index)

            # c2 stats ...
            pda4MJD.calcaverage_sigmacutloop('uJy',noisecol='duJy',maskcol='Mask',maskval=(self.flag_c0_uncertainty|self.flag_c0_X2norm),verbose=1,Nsigma=3.0,median_firstiteration=True)
            # ... and save them into the table

            self.lc.statresults2table(pda4MJD,c2_param2columnmapping,destindex=index)            

        # Save the light curve with cuts in Mask
        self.save_lc(SNindex=SNindex,overwrite=True,controlindex=0)
           
        return(0)

    def make_c1c2_cuts(self):
        #c1_X2norm_max = self.cfg.params['cleanlc']['cut1']['c1_X2norm_max']
        #c1_absmeannorm_max = self.cfg.params['cleanlc']['cut1']['c1_absmeannorm_max']
        #c2_X2norm_max = self.cfg.params['cleanlc']['cut2']['c2_X2norm_max']
        #c2_Nclipped_max = self.cfg.params['cleanlc']['cut2']['c2_Nclipped_max']
        #c2_Nused_max = self.cfg.params['cleanlc']['cut2']['c2_Nused_max']

        # load main lc
        self.load_lc(SNindex,controlindex=0)

        for index in self.lc.getindices():
            # c1 cut!
            mask = 0
            if (self.lc.t.loc[index,'c1_X2norm']>=self.cfg.params['cleanlc']['cut1']['c1_X2norm_max']): 
                mask |= self.flag_c1_X2norm
            if (np.fabs(self.lc.t.loc[index,'c1_mean']/self.lc.t.loc[index,'c1_mean_err'])>=self.cfg.params['cleanlc']['cut1']['c1_absmeannorm_max']):  
                mask |= self.flag_c1_absnormmean
                
            if mask == 0:
                self.lc.t.loc[index,'Mask'] |= self.flag_c1_good
                # if the c1 cut indicates a good measurement, skip c2 cut
                continue                
            else:
                self.lc.t.loc[index,'Mask'] = np.bitwise_or(self.lc.t.loc[index,'Mask'],mask)
                
                
            # c2 cut!
            mask = 0
            if (self.lc.t.loc[index,'c2_X2norm']>self.cfg.params['cleanlc']['cut2']['c2_X2norm_max']):
                mask |= self.flag_c2_X2norm
            if (np.fabs(self.lc.t.loc[index,'c2_mean']/self.lc.t.loc[index,'c2_mean_err'])>self.cfg.params['cleanlc']['cut2']['c2_absmeannorm_max']):
                mask |= self.flag_c2_absmeanerr
            if (self.lc.t.loc[index,'c2_Nclip']>self.cfg.params['cleanlc']['cut2']['c2_Nclipped_max']):
                mask |= self.flag_c2_Nclip
            if (self.lc.t.loc[index,'c2_Ngood']<self.cfg.params['cleanlc']['cut2']['c2_Ngood_min']):
                mask |= self.flag_c2_Nused

            if mask == 0 and self.lc.t.loc[index,'c2_Nclip']==0:
                self.lc.t.loc[index,'Mask'] = np.bitwise_or(self.lc.t.loc[index,'Mask'],self.flag_c2_good)
                # all good!!
                continue         
            elif mask == 0:
                self.lc.t.loc[index,'Mask'] = np.bitwise_or(self.lc.t.loc[index,'Mask'],self.flag_c2_ok)
                # all good!!
                continue         
            else:
                self.lc.t.loc[index,'Mask'] = np.bitwise_or(self.lc.t.loc[index,'Mask'],mask|self.flag_c2_bad)
                
        if self.verbose>2:
            self.lc.write()
        # save main lc
        self.save_lc(SNindex=SNindex,controlindex=0,filt=self.filt,overwrite=True)
        
        print('Copying over SN c2 mask column to control lc mask column...')
        # get c1 and c2 masks from SN lc and copy to control lc mask column
        flags_array = np.full(self.lc.t['Mask'].shape,self.flags_c1c2)
        
        
        omask = np.bitwise_and(self.lc.t['Mask'],flags_array)
        # loop through control lcs
        for controlindex in range(1,len(self.RADECtable.t)):
            # load control lc
            self.load_lc(SNindex,controlindex=controlindex,filt=self.filt)
            if self.verbose>1: print('Length of self.lc.t: ',len(self.lc.t))
            if len(self.lc.t)==0: return(1)
            # copy over SN c1 and c2 masks to control lc mask column
            self.lc.t['Mask'] = np.bitwise_or(self.lc.t['Mask'],omask)
            # save lc
            print('... done, saving lc')
            self.save_lc(SNindex=SNindex,controlindex=controlindex,filt=self.filt,overwrite=True)

    def cleanuplcloop(self,args,SNindex):

        # cut0 - mask lcs based on PSF X2norm and uncertainty
        (MJD_SN,uJy,duJy,Mask) = self.make_c0_cuts(SNindex,prepare_c1c2_cuts = self.cfg.params['cleanlc']['apply_c1c2'])
        
        if not self.cfg.params['cleanlc']['apply_c1c2']:
            print('Skipping c1c2 cuts since apply_c1c2 is ',self.cfg.params['cleanlc']['apply_c1c2'])
            return(0)
        
        if self.verbose>2:
            print(self.flux_colname,': ',uJy)
            print(self.dflux_colname,': ',duJy)
            print('Mask: ',Mask)
            
        # calculate control lc stats
        print('Calculating control LCs statistics...')
        self.calc_c1c2_stats(uJy,duJy,Mask)    

        
        # make cuts on good, bad, and ok measurements
        print('Making cuts based on control LCs statistics...')
        self.make_c1c2_cuts()

if __name__ == '__main__':

    cleanuplc = cleanuplcclass()
    parser = cleanuplc.define_options()
    args = parser.parse_args()

    SNindexlist = cleanuplc.initialize(args)

    for SNindex in SNindexlist:
        cleanuplc.loadRADEClist(SNindex=SNindex,filt=cleanuplc.filt)
        cleanuplc.cleanuplcloop(args,SNindex)

    print('\n')
