#!/usr/bin/env python

import astropy.table as at
import sys,socket,os,re
from astropy.io import ascii
from datetime import datetime as dt
import argparse
from astrotable import astrotableclass
from tools import yamlcfgclass
from tools import rmfile
from download_atlas_lc import download_atlas_lc_class
from SNloop import SNloopclass
import numpy as np

class downloadloopclass(SNloopclass):
    def __init__(self):
        SNloopclass.__init__(self)
        self.download_atlas_lc = download_atlas_lc_class()        

    def getlc4SN(self, SNindex,lookbacktime_days=None, savelc=False, overwrite=False, fileformat=None, offsetindex=None):
        if offsetindex is None:
           RA = self.t[SNindex]['ra']
           Dec = self.t[SNindex]['dec']
        else:
            RA = self.RADECtable.t['RaNew'][offsetindex]
            Dec = self.RADECtable.t['DecNew'][offsetindex]

        self.download_atlas_lc.verbose = self.verbose
        self.download_atlas_lc.debug = self.debug
        self.download_atlas_lc.get_lc(RA,Dec,lookbacktime_days=lookbacktime_days)

        # read the lc into an ascii table
        self.lc.load_spacesep(self.download_atlas_lc.lcraw, formatMapping={'MJD':'%.6f','chi/N':'%.2f','x':'%.2f','y':'%.2f','m':'%.3f','dm':'%.3f'})

        # save the lc file with the automatic output filename
        if savelc:
            self.save_lc(SNindex,overwrite=overwrite,fileformat=fileformat, offsetindex=offsetindex)
            for filt in ['c','o']:
                filename = self.lcbasename(SNindex,filt=filt, offsetindex=offsetindex)+'.txt'
                if fileformat is None: fileformat = self.cfg.params['output']['fileformat']
                detections4filt=np.where(self.lc.t['F']==filt)
                if len(detections4filt[0]) is 0:
                    print('Removing %s because nothing to save...' % filename)
                    rmfile(filename)
                else: 
                    self.lc.write(filename, indeces=detections4filt, format=fileformat, overwrite=overwrite, verbose=(self.verbose>0))
 
if __name__ == '__main__':

    downloadATLASlc = downloadloopclass()
    parser = downloadATLASlc.download_atlas_lc.define_optional_args()
    parser = downloadATLASlc.define_options(parser=parser)
    args = parser.parse_args()

    # load config files
    downloadATLASlc.cfg.loadcfgfiles(args.cfgfile,
                                    extracfgfiles=args.extracfgfile,
                                    params=args.params,
                                    params4all=args.pall,
                                    params4sections=args.pp,
                                    verbose=args.verbose)

    downloadATLASlc.setoutdir(outrootdir=args.outrootdir, outsubdir=args.outsubdir)
    downloadATLASlc.verbose = args.verbose
    downloadATLASlc.debug = args.debug

    downloadATLASlc.download_atlas_lc.connect(args.atlasmachine,'sofia','Starswirl1410!@')

    downloadATLASlc.load_spacesep(args.snlistfilename)
    print(downloadATLASlc.t)
    print(args.SNlist)
    indexlist = downloadATLASlc.getSNlist(args.SNlist)

    for i in indexlist:
        downloadATLASlc.getlc4SN(i,
                                lookbacktime_days=args.lookbacktime_days,
                                savelc=args.savelc,
                                overwrite=args.overwrite,
                                fileformat=args.fileformat)
