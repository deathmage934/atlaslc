#!/usr/bin/env python

import numpy as np
import math
import sys,socket,os,re
#import astropy.table as at
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import Angle
#import matplotlib.pyplot as plt
#import pylab as matlib
#from datetime import datetime as dt
import argparse
from tools import yamlcfgclass
from tools import rmfile
from tools import RaInDeg
from tools import DecInDeg
from astrotable import astrotableclass
from download_atlas_lc import download_atlas_lc_class
import sigmacut

class SNloopclass(astrotableclass):
	def __init__(self):
		astrotableclass.__init__(self)

		# config file
		self.cfg = yamlcfgclass()

		self.verbose = 0
		self.debug = False
		self.outrootdir = None
		self.filt=None
		self.lc = astrotableclass()
		self.RADECtable = astrotableclass(names=('OffsetID','Ra','Dec','RaNew','DecNew','RaDistance','RaOffset','DecOffset','Radius','Ndet','Ndet_c','Ndet_o'), dtype=(np.int32,np.float64,np.float64,np.float64,np.float64,np.float32,np.float32,np.float32,np.int32,np.int32,np.int32,np.int32))
		self.averagelctable = astrotableclass(names=('OffsetID','MJD','uJy','duJy','stdev','X2norm','Nused','Nclipped','MJDNused','MJDNskipped'), dtype=(np.int32,np.float64,np.float64,np.float64,np.float64,np.float64,np.int64,np.int64,np.int64,np.int64))

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

		parser.add_argument('SNlist', nargs='+')
		parser.add_argument('--verbose', '-v', default=0, action='count')
		parser.add_argument('-d', '--debug', help="debug", action='count')
		parser.add_argument('--snlistfilename', default=None, help=('filename of SN list (default=%(default)s)'))
		parser.add_argument('-s','--savelc', help=("save lc"), action="store_true", default=False)
		parser.add_argument('--outrootdir', default=outrootdir, help=('output root directory.''(default=%(default)s)'))
		parser.add_argument('--outsubdir', default=None, help=('subdir added to the output root directory (and filename) ''(default=%(default)s)'))
		parser.add_argument('-c', '--cfgfile', default=cfgfile, help='main config file. (default=%(default)s)')
		parser.add_argument('-e', '--extracfgfile', action='append', default=None, help=('additional config file. These cfg files do not need to have all ''parameters. They overwrite the parameters in the main cfg file.'))
		parser.add_argument('-p', '--params', action='append', default=None, nargs=2, help=('"param val": change parameter in config file (not in section, only ''main part) (default=%(default)s)'))
		parser.add_argument('--pall', action='append', default=None, nargs=2, help=('"param val". change parameter in all sections of config file ''(section independent) (default=%(default)s)'))
		parser.add_argument('--pp', action='append', default=None, nargs=3, help=('"section param val". change parameters in given section of ''config file (default=%(default)s)'))
		parser.add_argument('-f','--filt', default=None, help=('specify filter'), choices=['c','o'])
		parser.add_argument('--MJDbinsize', default=None, help=('specify MJD bin size'),type=float)
		parser.add_argument('--forcedphot_offset', default=False)
		parser.add_argument('--pattern', default='circular',help=('offset pattern, defined in the config file, (default=%(default)s)'))
		parser.add_argument('--plot', default=False, help=('plot lcs'))
		parser.add_argument('--averagelc',default=False,help=('average lcs'))
		parser.add_argument('--skip_uncert',default=False,help=('skip cleanup lcs using uncertainties'))
		parser.add_argument('--skip_chi',default=False,help=('skip cleanup lcs using chi/N'))

		return parser

	def setoutdir(self,outrootdir=None, outsubdir=None):
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
		result=self.cfg.loadcfgfiles(*pargs, **kwargs)
		# set filter to filt in precursor.cfg, then check if args.filt set
		self.filt = self.cfg.params['filter']
		if not(filt is None):
			self.filt=filt
		return(result)

	def lcbasenameID(self, SNID, offsetindex=None, filt=None, MJDbinsize=None):
		basename = '%s/%s/%s' % (self.outrootdir,SNID,SNID)
		if not(offsetindex is None):
			basename += '_i%03d' % self.RADECtable.t['OffsetID'][offsetindex]

		self.filt = self.cfg.params['filter']
		if not(filt is None):
			self.filt=filt
			basename += '.%s' % filt

		if not(MJDbinsize is None):
			basename += '.%ddays' % int(MJDbinsize)
		basename += '.lc'

		print(basename)
		return(basename)

	def lcbasename(self, SNindex, offsetindex=None, filt=None, MJDbinsize=None):	
		print(self.t)
		SNID = self.t[SNindex]['tnsname']
		if filt==None:
			filt=self.filt
		return(self.lcbasenameID(SNID,filt=filt,offsetindex=offsetindex,MJDbinsize=MJDbinsize))

	def getSNlist(self, SNlist):
		if len(SNlist)==1 and SNlist[0]=='all':
			SNindexlist = range(len(self.t))
		else:
			SNindexlist = []
			for i in range(len(self.t)):
				if self.t[i]['tnsname'] in SNlist:
					SNindexlist.append(i) 
		return(SNindexlist)

	def load_lc(self, SNindex, filt=None, fileformat=None, offsetindex=None, MJDbinsize=None):
		filename = self.lcbasename(SNindex,filt=filt,offsetindex=offsetindex,MJDbinsize=MJDbinsize)+'.txt'
		if self.verbose: print('Loading lc %s' % filename)
		self.lc.load_generic(filename)
		if fileformat is None: fileformat = self.cfg.params['output']['fileformat']
		return(0)

	def save_lc(self, SNindex, filt=None, overwrite=False, fileformat=None, offsetindex=None, MJDbinsize=None):
		filename = self.lcbasename(SNindex, filt=filt, offsetindex=offsetindex, MJDbinsize=MJDbinsize)+'.txt'
		if fileformat is None: fileformat = self.cfg.params['output']['fileformat']
		self.lc.write(filename,format=fileformat, overwrite=overwrite,verbose=(self.verbose>0))
		return(0)

	def saveRADEClist(self, SNindex, filt=None):
		RADEClistfilename = self.lcbasename(SNindex,filt=filt)+'.RADEClist.txt'
		if self.verbose: print('Saving RADEClist %s' % RADEClistfilename)
		self.RADECtable.write(RADEClistfilename,format='fixed_width_two_line')
		return(0)

	def loadRADEClist(self, SNindex, filt=None):
		RADEClistfilename = self.lcbasename(SNindex,filt=filt)+'.RADEClist.txt'
		if self.verbose: print('Loading RADEClist %s' % RADEClistfilename)
		self.RADECtable.load_generic(RADEClistfilename)
		print(self.RADECtable.t)
		return(0)

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
			snlistfilename=args.snlistfilename

		if not(os.path.isfile(snlistfilename)):
			raise RuntimeError("SN list file %s does not exist, exiting!!" % snlistfilename)
			
		self.setoutdir(outrootdir=args.outrootdir, outsubdir=args.outsubdir)
		self.verbose = args.verbose
		self.debug = args.debug

		self.load_spacesep(snlistfilename)
		print(self.t)
		print(args.SNlist)

		SNindexlist = self.getSNlist(args.SNlist)
		return(SNindexlist)
