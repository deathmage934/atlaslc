#!/usr/bin/env python
"""
Created on Wed Dec 30 12:32:14 2020

@author: arest
"""

from SNloop import SNloopclass
import sys
import numpy as np
from pdastro import AnotB

class verifyMJDclass(SNloopclass):
    def __init__(self):
        SNloopclass.__init__(self)
    
    def verifyMJD(self,SNindex,skip_savelc=False):
        """
        Parameters
        ----------
        SNindex : interger
            index of SN in SNlist.


        Returns
        -------
        None.

        """
  
        # load main SN lc
        self.load_lc(SNindex,controlindex=0,filt=self.filt)
        if self.verbose:
            print('Length of self.lc.t: ',len(self.lc.t))
        if len(self.lc.t) == 0:
            return(1)
        
        # sort the lc by MJD ...
        ix_sorted = self.lc.ix_sort_by_cols('MJD')
        # ... and save it sorted
        if not(skip_savelc):
            if self.verbose: print('Saving SN LC sorted by MJD')
            self.save_lc(SNindex,controlindex=0,filt=self.filt,indices=ix_sorted)
        # ... and get the sorted MJD array
        MJD_SN = self.lc.t.loc[ix_sorted,'MJD'].to_numpy()
        if self.verbose>2:
            print('SN LC sorted by MJD')
            self.lc.write(indices=ix_sorted)
         
        # loop through the control LCs
        for controlindex in range(1,len(self.RADECtable.t)):
            if self.verbose: print('Control LC %d' % controlindex)

            # load offset lc
            self.load_lc(SNindex,controlindex=controlindex,filt=self.filt)
            
            # sort it by MJD...
            ix_sorted = self.lc.ix_sort_by_cols('MJD')
            # ... and get the sorted MJD array
            MJD_controllc = self.lc.t.loc[ix_sorted,'MJD'].to_numpy()
            
            # Compare the sorted MJD arrays of the SN and the control LC. 
            # If not in agreement, fix it!
            if (len(MJD_SN) != len(MJD_controllc)) or (np.array_equal(MJD_SN, MJD_controllc) is False):
                if self.verbose: print('Inconsistent MJD array, fixing...')

                # get the MJDs that are only in the SN LC
                MJDs_onlySN = AnotB(MJD_SN,MJD_controllc)
                # get the MJDs that are only in the control LC
                MJDs_onlycontrol = AnotB(MJD_controllc,MJD_SN)

                # For the MJDs only in the SN LC: add a line with that MJD 
                # to control LC, with all values of other columns NaN
                if len(MJDs_onlySN)>0:
                    if self.verbose: print('MJDs only in SN LC:',MJDs_onlySN)
                    for MJD in MJDs_onlySN:
                        self.lc.newrow({'MJD':MJD,'Mask':0})

                # remove the indices of rows in control LC for which there is 
                # no MJD in the SN LC: not needed!
                if len(MJDs_onlycontrol)>0:
                    if self.verbose: print('MJDs only in control LC:',MJDs_onlycontrol)
                    indices2skip=[]
                    for MJD in MJDs_onlycontrol:
                        ix = self.lc.ix_equal('MJD',MJD)
                        if len(ix)!=1:
                            raise RuntimeError('Couldn\'t find MJD={} in column MJD, but should be there!!!'.format(MJD))
                        indices2skip.extend(ix)
                    #self.lc.t.drop(index=indices2skip)
                    indices = AnotB(self.lc.getindices(),indices2skip)
                else:
                    indices = self.lc.getindices()
                
                ix_sorted = self.lc.ix_sort_by_cols('MJD',indices=indices)
                if not(skip_savelc):
                    if self.verbose: print('Saving control LC sorted by MJD')
                    self.lc.t['Mask'] = self.lc.t['Mask'].astype(np.int32)
                    self.save_lc(SNindex,controlindex=controlindex,filt=self.filt,indices=ix_sorted)
                if self.verbose>2:
                    print('control LC sorted by MJD')
                    self.lc.write(indices=ix_sorted)

            else:
                if self.verbose: print('MJDs fit!!')

if __name__ == '__main__':

    verifyMJD = verifyMJDclass()
    parser = verifyMJD.define_options()
    args = parser.parse_args()

    SNindexlist = verifyMJD.initialize(args)

    for SNindex in SNindexlist:
        verifyMJD.loadRADEClist(SNindex,filt=verifyMJD.filt)
        verifyMJD.verifyMJD(SNindex)
 