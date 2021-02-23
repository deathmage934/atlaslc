#!/usr/bin/env python
'''
@author: S. Rest
'''

import numpy as np
import math
import sys,socket,os,re
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import Angle
import pandas as pd
import argparse
from tools import yamlcfgclass
from tools import rmfile
from tools import RaInDeg
from tools import DecInDeg
#from astrotable import astrotableclass
from download_atlas_lc import download_atlas_lc_class
import sigmacut
from pdastro import pdastroclass, pdastrostatsclass, AnotB
import mastcasjobs
import pylab
import json
import getpass

class SNloopclass(pdastroclass):
    def __init__(self):
        pdastroclass.__init__(self)

        self.lctype = None # can be 'avg' or 'single'

        # config file
        self.cfg = yamlcfgclass()

        # miscellaneous vars
        self.verbose = 0
        self.debug = False
        self.outrootdir = None
        self.filt = None
        
        # tables
        self.lc = pdastrostatsclass(hexcols=['Mask'])
        self.RADECtable = pdastroclass()
        self.averagelc = pdastrostatsclass(hexcols=['Mask'])

        # FLAGS
        self.flag_c0_X2norm      = 0x1 
        self.flag_c0_uncertainty = 0x2
        
        #self.flag_c1_good = 0x20 
        self.flag_c1_X2norm      = 0x10
        self.flag_c1_absnormmean = 0x20
        
        self.flag_c2_X2norm      = 0x100
        self.flag_c2_absnormmean = 0x200
        self.flag_c2_Nclip       = 0x400
        self.flag_c2_Nused       = 0x800

        self.flag_daysigma       = 0x1000
        self.flag_daysmallnumber = 0x2000

        #self.flag_c0_good = 0x10000
        #self.flag_c1_good        = 0x20000
        #self.flag_c2_good        = 0x40000
        self.flag_c2_ok          = 0x80000

        self.flag_c0_bad         = 0x100000
        #self.flag_c1_bad    = 0x200000
        self.flag_c2_bad         = 0x400000
        self.flag_day_bad        = 0x800000

        self.flags={'flag_c0_uncertainty':self.flag_c0_uncertainty,
                    'flag_c0_X2norm':self.flag_c0_X2norm,
                    'flag_c1_X2norm':self.flag_c1_X2norm,
                    'flag_c1_absnormmean':self.flag_c1_absnormmean,
                    'flag_c2_X2norm':self.flag_c2_X2norm,
                    'flag_c2_absnormmean':self.flag_c2_absnormmean,
                    'flag_c2_Nclip':self.flag_c2_Nclip,
                    'flag_c2_Nused':self.flag_c2_Nused,
                    'flag_daysigma':self.flag_daysigma,
                    'flag_daysmallnumber':self.flag_daysmallnumber,
                    'flag_c0_bad':self.flag_c0_bad,
                    'flag_c2_bad':self.flag_c2_bad,
                    'flag_day_bad':self.flag_day_bad}
        
        self.flags_c1c2 = self.flag_c1_X2norm|self.flag_c1_absnormmean|self.flag_c2_X2norm|self.flag_c2_absnormmean|self.flag_c2_Nclip|self.flag_c2_Nused|self.flag_c2_bad|self.flag_c2_ok#|self.flag_c1_good|self.flag_c2_good

    def define_options(self, parser=None, usage=None, conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)

        # options for config file
        if 'ATLASLC_SOURCEDIR' in os.environ and os.environ['ATLASLC_SOURCEDIR'] != '':
            cfgfile = '%s/precursor.cfg' % os.environ['ATLASLC_SOURCEDIR']
            outrootdir = os.environ['ATLASLC_DATA']
        else:
            cfgfile = None
            outrootdir = None

        # can be a sn name from snlist.txt or 'all' (loops through all sn names in snlist.txt)
        parser.add_argument('SNlist', nargs='+')
        parser.add_argument('-f','--filt', default=None, choices=['c','o'], help=('specify default filter'))
        parser.add_argument('-m','--MJDbinsize', default=1.0, help=('specify MJD bin size for averaging lcs'),type=float)
        parser.add_argument('--forcedphot_offset', default=False, help=("download offsets (settings in config file)"))
        parser.add_argument('--api', default=False, help=('use API instead of SSH to get light curves from ATLAS'))
        parser.add_argument('--plot', default=False, help=('plot lcs'))
        parser.add_argument('--plot_avg', default=False, help=('plot average lcs'))
        parser.add_argument('--xlim_lower', default=None, type=float, help=('set lower x limit when plotting'))
        parser.add_argument('--xlim_upper', default=None, type=float, help=('set upper x limit when plotting'))
        parser.add_argument('--ylim_lower', default=None, type=float, help=('set lower y limit when plotting'))
        parser.add_argument('--ylim_upper', default=None, type=float, help=('set upper y limit when plotting'))
        parser.add_argument('--averagelc', default=False, help=('average lcs'))
        parser.add_argument('--detectbumps', default=False, help=('detect bumps in lcs'))
        parser.add_argument('-v','--verbose', default=0, action='count')
        parser.add_argument('-d', '--debug', action='count', help="debug")
        parser.add_argument('--snlistfilename', default=None, help=('filename of SN list (default=%(default)s)'))
        parser.add_argument('-s','--savelc', default=False, action="store_true", help=("save lc"))
        parser.add_argument('--outrootdir', default=outrootdir, help=('output root directory.''(default=%(default)s)'))
        parser.add_argument('--outsubdir', default=None, help=('subdir added to the output root directory (and filename) ''(default=%(default)s)'))
        parser.add_argument('-c', '--cfgfile', default=cfgfile, help='main config file. (default=%(default)s)')
        parser.add_argument('-e', '--extracfgfile', default=None, action='append', help=('additional config file. These cfg files do not need to have all ''parameters. They overwrite the parameters in the main cfg file.'))
        parser.add_argument('-p', '--params', default=None, action='append', nargs=2, help=('"param val": change parameter in config file (not in section, only ''main part) (default=%(default)s)'))
        parser.add_argument('--pall', action='append', default=None, nargs=2, help=('"param val". change parameter in all sections of config file ''(section independent) (default=%(default)s)'))
        parser.add_argument('--pp', action='append', default=None, nargs=3, help=('"section param val". change parameters in given section of ''config file (default=%(default)s)'))
        return parser

    def setoutdir(self,args,outrootdir=None,outsubdir=None):
        if outrootdir is None:
            basedir = self.cfg.params['output']['outrootdir']
        else:
            basedir = outrootdir

        # add outsubdir
        outsubdir = self.cfg.params['output']['outsubdir']
        if not(outsubdir is None):
            outsubdir = outsubdir
        if not(outsubdir is None):
            basedir += '/'+outsubdir
        self.outrootdir = basedir

    def loadcfgfiles(self, *pargs, filt=None, **kwargs):
        cfgfiles=self.cfg.loadcfgfiles(*pargs, **kwargs)
        # set filter to filt in precursor.cfg, then check if args.filt set
        self.filt = self.cfg.params['filter']
        if not(filt is None):
            self.filt=filt
        return(cfgfiles)

    def lcbasename(self, SNindex=None, yse=False, TNSname=None, controlindex=None, filt=None, MJDbinsize=None,addsuffix=None):
        # define address and file name of the data table
        if yse is True:
            SNID = TNSname
        else:
            SNID = self.t['tnsname'][SNindex]
        
        if filt is None:
            filt=self.filt
        
        basename = '%s/%s/%s' % (self.outrootdir,SNID,SNID)
        
        if not(controlindex is None):
            basename += '_i%03d' % self.RADECtable.t['ControlID'][controlindex]
        
        self.filt = self.cfg.params['filter']
        
        if not(filt is None):
            self.filt=filt
            basename += '.%s' % filt
        
        if not(MJDbinsize is None):
            if int(MJDbinsize)==MJDbinsize:
                basename += '.%ddays' % int(MJDbinsize)
            else:
                basename += '.%.2fdays' % int(MJDbinsize)
        
        basename += '.lc'
        if not(addsuffix is None):
            basename += addsuffix
        return(basename)

    def getSNlist(self, SNlist):
        # if --snlist all, get data for all SN in snlist.txt; otherwise, only for listed SN in cmd
        if len(SNlist)==1 and SNlist[0]=='all':
            SNindexlist = range(len(self.t))
        else:
            SNindexlist = []
            for index in range(0,len(self.t)):
                if self.t.at[index,'tnsname'] in SNlist:
                    SNindexlist.append(index)
        return(SNindexlist)

    def load_lc(self, SNindex, filt=None, controlindex=None, MJDbinsize=None,addsuffix=None,hexcols=None):
        # get lc from already existing file
        if MJDbinsize is None:
            self.lctype = 'og'
        else:
            self.lctype = 'avg'
        filename = self.lcbasename(SNindex=SNindex,filt=filt,controlindex=controlindex,MJDbinsize=MJDbinsize,addsuffix=addsuffix)+'.txt' 
        self.lc.load_spacesep(filename, delim_whitespace=True, hexcols=hexcols,verbose=(self.verbose>1))
        if self.lc.default_formatters is None: self.lc.default_formatters = {}
        self.lc.default_formatters['Mask']='0x{:06x}'.format
            
        return(0)

    def save_lc(self, SNindex=None, yse=False, TNSname=None, indices=None, filt=None, overwrite=False, controlindex=None, MJDbinsize=None,addsuffix=None):
        # write table and save lc as file
        filename = self.lcbasename(SNindex=SNindex, yse=False, filt=filt, controlindex=controlindex, MJDbinsize=MJDbinsize,addsuffix=addsuffix)+'.txt'
        self.lc.write(filename, indices=indices,overwrite=True,verbose=self.verbose)
        #self.lc.write(filename,format=fileformat, overwrite=overwrite,verbose=(self.verbose>0))
        return(0)

    def saveRADEClist(self, SNindex, filt=None):
        RADEClistfilename = self.lcbasename(SNindex=SNindex,filt=filt)+'.RADEClist.txt'
        self.RADECtable.write(RADEClistfilename,overwrite=True,verbose=True)
        return(0)

    def loadRADEClist(self, SNindex, filt=None):
        # get RADEClist from alreadt existing file
        RADEClistfilename = self.lcbasename(SNindex=SNindex,filt=filt)+'.RADEClist.txt'
        if self.verbose>1: print('Loading RADEClist %s' % RADEClistfilename)
        self.RADECtable.load_spacesep(RADEClistfilename, delim_whitespace=True)
        
        # temporary fix to move OffsetID to ControlID
        if 'OffsetID' in self.RADECtable.t.columns:
            if 'ControlID' in self.RADECtable.t.columns: self.RADECtable.t.drop(columns=['ControlID'])
            self.RADECtable.t = self.RADECtable.t.rename(columns={'OffsetID':'ControlID'})
        
        if self.verbose>2:
            self.RADECtable.write()
        return(0)

    def makecuts_indices(self,SNindex,controlindex,procedure1):
        # use when cleaning up data in plot_lc.py or average_lc.py; makes cuts based on mask column created in cleanup_lc.py
        # set flags in precursor.cfg to control what data to cut based on uncertainties and/or chi/N
        if procedure1 == 'plotlc': 
            print('Using cfg flags: ',self.cfg.params['plotlc']['flags2apply'])
            maskval = 0
            for key in self.cfg.params['plotlc']['flags2apply']:
                if not (key in self.flags):
                    raise RuntimeError('Bad flag name: %s' % key)
                maskval |= self.flags[key]
            flags = maskval
        elif procedure1 =='upltoyse':
            flags = self.cfg.params['upltoyse']['flags']
        else:
            raise RuntimeError('Procedure %s must be plotlc or upltoyse!' % procedure1)
        print('Setting indices using flag values: 0x%x' % flags)
        
        mask=np.bitwise_and(self.lc.t['Mask'], flags)
        cuts_indices = np.where(mask==0)
        cuts_indices = cuts_indices[0]
        #cuts_indices = self.lc.ix_unmasked('Mask',maskval=flags)
        #bad_data = self.lc.AnotB(self.lc.getindices(),cuts_indices)
        bad_data = np.where(mask!=0)
        bad_data = bad_data[0]
        datacut = len(self.lc.t)-len(self.lc.t.loc[cuts_indices])
        print('Length original lc: ',len(self.lc.t),', length cleaned lc: ',len(self.lc.t.loc[cuts_indices]),', data points cut: ',datacut)

        lc_uJy = self.lc.t.loc[cuts_indices, self.flux_colname]
        lc_duJy = self.lc.t.loc[cuts_indices, self.dflux_colname]
        lc_MJD = self.lc.t.loc[cuts_indices, 'MJD']
        return(lc_uJy, lc_duJy, lc_MJD, cuts_indices, bad_data)

    # get all measurements that have not been flagged as bad or questionable
    def getgoodindices(self,indices=None):
        if self.lctype is None:
            if 'ControlID' in self.lc.columns:
                self.lctype = 'avg'
            else:
                self.lctype = 'og'
        masks = 0 #self.flag_c2_bad|self.flag_c0_X2norm|self.flag_c0_uncertainty
        if self.lctype == 'og':
            flagslist = self.cfg.params['cleanlc']['questionable_flags_og']
            flagslist = flagslist.append(self.cfg.params['cleanlc']['exclude_flags_og'])
            print('Excluding flags: ',flagslist)
            for key in flagslist:
                if not(key in self.flags):
                    raise RuntimeError('Bad flag name: %s' % key)
                masks |= self.flags[key]
        elif self.lctype == 'avg':
            flagslist = self.cfg.params['cleanlc']['questionable_flags_avg']
            flagslist = flagslist.append(self.cfg.params['cleanlc']['exclude_flags_avg'])
            print('Excluding flags: ',flagslist)
            for key in flagslist:
                if not(key in self.flags):
                    raise RuntimeError('Bad flag name: %s' % key)
                masks |= self.flags[key]
        else:
            raise RuntimeError('lc type must be og or avg!')
        indices = self.lc.ix_unmasked('Mask',maskval=masks)
        return(indices)

    # get all measurements that have been flagged as questionable
    def getquestionableindices(self):
        if self.lctype is None:
            if 'ControlID' in self.lc.columns:
                self.lctype = 'avg'
            else:
                self.lctype = 'og'
        masks = 0 #self.flag_c2_bad|self.flag_c0_X2norm|self.flag_c0_uncertainty
        if self.lctype == 'og':
            for key in self.cfg.params['cleanlc']['questionable_flags_og']:
                if not(key in self.flags):
                    raise RuntimeError('Bad flag name: %s' % key)
                masks |= self.flags[key]
        elif self.lctype == 'avg':
            for key in self.cfg.params['cleanlc']['questionable_flags_avg']:
                if not(key in self.flags):
                    raise RuntimeError('Bad flag name: %s' % key)
                masks |= self.flags[key]
        else:
            raise RuntimeError('lc type must be og or avg!')
        indices = self.lc.ix_masked('Mask',maskval=masks)
        return(indices)

    # get all measurements that have not been flagged as bad
    def getusableindices(self):
        if self.lctype is None:
            if 'ControlID' in self.lc.columns:
                self.lctype = 'avg'
            else:
                self.lctype = 'og'
        masks = 0 #self.flag_c2_bad|self.flag_c0_X2norm|self.flag_c0_uncertainty
        if self.lctype == 'og':
            for key in self.cfg.params['cleanlc']['exclude_flags_og']:
                if not(key in self.flags):
                    raise RuntimeError('Bad flag name: %s' % key)
                masks |= self.flags[key]
        elif self.lctype == 'avg':
            for key in self.cfg.params['cleanlc']['exclude_flags_avg']:
                if not(key in self.flags):
                    raise RuntimeError('Bad flag name: %s' % key)
                masks |= self.flags[key]
        else:
            raise RuntimeError('lc type must be og or avg!')
        indices = self.lc.ix_unmasked('Mask',maskval=masks)
        return(indices)

    # get all measurements that have been flagged as bad
    def getbadindices(self):
        if self.lctype is None:
            if 'ControlID' in self.lc.columns:
                self.lctype = 'avg'
            else:
                self.lctype = 'og'
        masks = 0 #self.flag_c2_bad|self.flag_c0_X2norm|self.flag_c0_uncertainty
        if self.lctype == 'og':
            for key in self.cfg.params['cleanlc']['exclude_flags_og']:
                if not(key in self.flags):
                    raise RuntimeError('Bad flag name: %s' % key)
                masks |= self.flags[key]
        elif self.lctype == 'avg':
            for key in self.cfg.params['cleanlc']['exclude_flags_avg']:
                if not(key in self.flags):
                    raise RuntimeError('Bad flag name: %s' % key)
                masks |= self.flags[key]
        else:
            raise RuntimeError('lc type must be og or avg!')
        indices = self.lc.ix_masked('Mask',maskval=masks)
        return(indices)

    def addnanrows(self,indices=None):
        # make new lc
        lc2 = pdastroclass(hexcols='Mask')
        lc2.t = pd.DataFrame(columns=self.lc.t.columns)
        print(lc2.t)
        # get indices
        if indices is None:
            indices = self.lc.getindices()
        print('Indices: ',indices)
        # get MJD and mjd range
        MJD = int(np.amin(self.lc.t['MJD']))
        print('MJD range: ',MJD,' to ',MJD+1)
        
        MJDbin = self.lc.t.loc[0,'MJDbin'] # get first MJD bin
        index = min(indices)
        while index < max(indices):
            if (self.lc.t.loc[index,'MJDbin']==MJDbin): # check to see if MJD bin matches
                # add row with MJDbin and other columns
                lc2.t = lc2.t.append(self.lc.t.loc[index],ignore_index=True)
                index += 1
            else:
                # add row with MJDbin and nans
                row = pd.DataFrame([[self.lc.t.loc[index,'ControlID'],np.nan,MJDbin,np.nan,np.nan,np.nan,np.nan,0,0,0,0x800000]],columns=self.lc.t.columns)
                lc2.t = lc2.t.append(row,ignore_index=True)
            MJDbin += 1 
        lc2.t['Mask'] = lc2.t['Mask'].astype(np.int32)
        return(lc2.t)

    def autosearch(self, ra, dec, search_size):
        os.environ['CASJOBS_WSID'] = str(self.cfg.params['casjobs_wsid'])
        print('Casjobs WSID set to %s in precursor.cfg...' % self.cfg.params['casjobs_wsid'])
        os.environ['CASJOBS_PW'] = getpass.getpass('Enter Casjobs password:')

        query = """select o.ObjID, o.raMean, o.decMean, o.nDetections
        from fGetNearbyObjEq("""+str(ra)+','+str(dec)+","+str(search_size/60)+""") nb
        join MeanObjectView o on o.ObjID=nb.ObjID
        where o.nDetections > 5
        and o.rmeankronmag < 18
        """
        jobs = mastcasjobs.MastCasJobs(context="PanSTARRS_DR2")
        results = jobs.quick(query, task_name="python cone search")
        return(results)

    def initialize(self,args):
        # load config files
        self.loadcfgfiles(args.cfgfile,
                            filt=args.filt,
                            extracfgfiles=args.extracfgfile,
                            params=args.params,
                            params4all=args.pall,
                            params4sections=args.pp,
                            verbose=args.verbose)

        snlistfilename = self.cfg.params['snlistfilename']
        if not(args.snlistfilename is None):
            snlistfilename = self.outrootdir+args.snlistfilename

        if not(os.path.isfile(snlistfilename)):
            raise RuntimeError("SN list file %s does not exist, exiting!!" % snlistfilename)
        
        self.setoutdir(args,outrootdir=args.outrootdir,outsubdir=args.outsubdir)
        self.verbose = args.verbose
        self.debug = args.debug
        self.flux_colname = self.cfg.params['flux_colname']
        self.dflux_colname = self.cfg.params['dflux_colname']

        self.RADECtable = pdastroclass(columns=['ControlID','PatternID','Ra','Dec','RaOffset','DecOffset','Radius','Ndet','Ndet_c','Ndet_o'])
        self.RADECtable.default_formatters = {'ControlID':'{:3d}'.format,
                                              'PatternID':'{:2d}'.format,
                                              'Ra':'{:.8f}'.format,
                                              'Dec':'{:.8f}'.format,
                                              'RaOffset':'{:.2f}'.format,
                                              'DecOffset':'{:.2f}'.format,
                                              'Radius':'{:.2f}'.format,
                                              'Ndet':'{:4d}'.format,
                                              'Ndet_c':'{:4d}'.format,
                                              'Ndet_o':'{:4d}'.format}


        self.averagelc = pdastroclass(columns=['ControlID','MJD','MJDbin',self.flux_colname,self.dflux_colname,'stdev','X2norm','Nclip','Ngood','Nexcluded','Mask'],hexcols='Mask')
        self.averagelc.default_formatters = {'ControlID':'{:3d}'.format,
                                            'MJD':'{:.6f}'.format,
                                            'MJDbin':'{:.2f}'.format,
                                            self.flux_colname:'{:.3f}'.format,
                                            self.dflux_colname:'{:.3f}'.format,
                                            'stdev':'{:.2f}'.format,
                                            'X2norm':'{:.3f}'.format,
                                            'Nclip':'{:4d}'.format,
                                            'Ngood':'{:4d}'.format,
                                            'Nexcluded':'{:4d}'.format,
                                            'Mask':'0x{:06x}'.format}

        self.load_spacesep(snlistfilename)
        if self.verbose>1:
            print(self.t)
        #print(args.SNlist)

        SNindexlist = self.getSNlist(args.SNlist)
        return(SNindexlist)
