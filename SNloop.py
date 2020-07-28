#!/usr/bin/env python

import astropy.table as at
import sys,socket,os,re
from astropy.io import ascii
from datetime import datetime as dt
import argparse
from astrotable import astrotableclass
from tools import yamlcfgclass
from download_atlas_lc import download_atlas_lc_class
import numpy as np

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
			snlistfilename = '%s/atlastest.txt' % os.environ['ATLASLC_DATA']
			outrootdir = os.environ['ATLASLC_DATA']
		else:
			cfgfile = None
			snlistfilename = None
			outrootdir = None

		parser.add_argument('SNlist', nargs='+')
		parser.add_argument('--verbose', '-v', default=0, action='count')
		parser.add_argument('-d', '--debug', help="debug", action='count')
		parser.add_argument('--snlistfilename', default=snlistfilename, help=('filename of SN list''(default=%(default)s)'))
		parser.add_argument('-s','--savelc', help=("save lc"), action="store_true", default=False)
		parser.add_argument('--outrootdir', default=outrootdir, help=('output root directory.''(default=%(default)s)'))
		parser.add_argument('--outsubdir', default=None,
							help=('subdir added to the output root directory (and filename) ''(default=%(default)s)'))
		parser.add_argument('-c', '--cfgfile', default=cfgfile, help='main config file. (default=%(default)s)')
		parser.add_argument('-e', '--extracfgfile', action='append', default=None,
							help=('additional config file. These cfg files do not need to have all ''parameters. They overwrite the parameters in the main cfg file.'))
		parser.add_argument('-p', '--params', action='append', default=None, nargs=2,
							help=('"param val": change parameter in config file (not in section, only ''main part) (default=%(default)s)'))
		parser.add_argument('--pall', action='append', default=None, nargs=2,
							help=('"param val". change parameter in all sections of config file ''(section independent) (default=%(default)s)'))
		parser.add_argument('--pp', action='append', default=None, nargs=3,
							help=('"section param val". change parameters in given section of ''config file (default=%(default)s)'))
		parser.add_argument('-f','--filt', default=None, help=('specify filter'), choices=['c','o'])
		parser.add_argument('--MJDbinsize', default=None, help=('specify MJD bin size'),type=float)
		
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

	def loadcfgfiles(self, *pargs, filt=None,**kwargs):
		result=self.cfg.loadcfgfiles(*pargs, **kwargs)
		self.filt = self.cfg.params['filter']
		if not(filt is None):
			self.filt=filt
		return(result)

	def lcbasenameID(self, SNID, offsetindex=None, filt=None, MJDbinsize=None):
		basename = '%s/%s/%s' % (self.outrootdir,SNID,SNID)
		if not(offsetindex is None):
			basename += '_i%03d' % self.RADECtable.t['OffsetID'][offsetindex]
		if not(filt is None):
			basename += '.%s' % filt
		if not(MJDbinsize is None):
			basename += '.%ddays' % int(MJDbinsize)
		basename += '.lc'


		print(basename)

		return(basename)

	def lcbasename(self, SNindex, offsetindex=None, filt=None, MJDbinsize=None):
		SNID = self.t[SNindex]['tnsname']
		if filt==None:
			filt=self.filt
		return(self.lcbasenameID(SNID,filt=filt,offsetindex=offsetindex,MJDbinsize=MJDbinsize))

	def getSNlist(self, SNlist):
		if len(SNlist)==1 and SNlist[0]=='all':
			indexlist = range(len(self.t))
		else:
			indexlist = []
			for i in range(len(self.t)):
				if self.t[i]['tnsname'] in SNlist:
					indexlist.append(i) 
		return(indexlist)

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
			
if __name__ == '__main__':

	SNloop = SNloopclass()
	parser = SNloop.download_atlas_lc.define_optional_args()
	parser = SNloop.define_options(parser=parser)
	args = parser.parse_args()

	# load config files
	SNloop.loadcfgfiles(args.cfgfile,
							extracfgfiles=args.extracfgfile,
							params=args.params,
							params4all=args.pall,
							params4sections=args.pp,
							verbose=args.verbose)

	SNloop.setoutdir(outrootdir=args.outrootdir, outsubdir=args.outsubdir)
	SNloop.verbose = args.verbose
	SNloop.debug = args.debug
	
	SNloop.download_atlas_lc.connect(args.atlasmachine,'arest','Sofilena50%')

	SNloop.load(args.snlistfilename)
	print(SNloop.t)
	print(args.SNlist)
	indexlist = SNloop.getSNlist(args.SNlist)

	for i in indexlist:
		SNloop.downloadlc4SN(i,
						lookbacktime_days=args.lookbacktime_days,
						savelc=args.savelc,
						overwrite=args.overwrite,
						fileformat=args.fileformat)
