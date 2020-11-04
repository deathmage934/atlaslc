#!/usr/bin/env python

# S. Rest

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
from astrotable import astrotableclass
from download_atlas_lc import download_atlas_lc_class
import sigmacut
from pdastro import pdastroclass
import mastcasjobs
import pylab
import json
import getpass

class SNloopclass(pdastroclass):
	def __init__(self):
		pdastroclass.__init__(self)

		# config file
		self.cfg = yamlcfgclass()
		self.verbose = 0
		self.debug = False
		self.outrootdir = None
		self.filt = None
		self.lc = pdastroclass()
		self.RADECtable = pdastroclass()
		self.averagelctable = pdastroclass()

		self.flag_cut0_uncertainty = 0x1
		self.flag_cut0_X2norm_dynamic = 0x2
		self.flag_cut0_X2norm_static = 0x4
		self.flag_o1_X2norm = 0x20
		self.flag_o2_X2norm = 0x200
		self.flag_o1_meannorm = 0x40
		self.flag_o2_meannorm = 0x400

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

		# set default filter
		parser.add_argument('-f','--filt', default=None, choices=['c','o'], help=('specify default filter'))
		
		# set MJD bin size for averagelc.py
		parser.add_argument('-m','--MJDbinsize', default=None, help=('specify MJD bin size'),type=float)
		
		# intialize forced photometry offsets
		parser.add_argument('--forcedphot_offset', default=False, help=("download offsets (settings in config file)"))
		
		# set forced photometry offset pattern
		parser.add_argument('--pattern', choices=['circle','box','closebright'], help=('offset pattern, defined in the config file; options are circle, box, or closebright'))
		
		# initialize plotlc.py
		parser.add_argument('--plot', default=False, help=('plot lcs'))
		
		# initialize averagelc.py
		parser.add_argument('--averagelc', default=False, help=('average lcs'))
		
		# skip uncertainty cleanup when cleaning lcs
		parser.add_argument('--skip_uncert', default=False, help=('skip cleanup lcs using uncertainties'))
		
		# skip chi square cleanup when cleaning lcs
		parser.add_argument('--skip_chi', default=False, help=('skip cleanup lcs using chi/N'))
		
		# specify whether or not to make cuts using cleaned data when averaging data
		parser.add_argument('--avg_makecuts', default=None, choices=['True','False'], help=('skip cutting measurements using mask column when averaging'))
		
		# specify procedure used in offsetstats.py
		parser.add_argument('--procedure', default='mask_nan', choices=['mask_nan','mask4mjd'],help=('define offsetstats.py procedure type. can be mask1 or mask2'))
		
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
		#parser.add_argument('--skip_makecuts_offsetstats',default=False,help=('skip cutting measurements using mask colum when '))
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
		cfgfiles=self.cfg.loadcfgfiles(*pargs, **kwargs)
		# set filter to filt in precursor.cfg, then check if args.filt set
		self.filt = self.cfg.params['filter']
		if not(filt is None):
			self.filt=filt
		return(cfgfiles)

	def lcbasename(self, SNindex=None, yse=False, TNSname=None, offsetindex=None, filt=None, MJDbinsize=None):
		# define address and file name of the data table
		if yse is True:
			SNID = TNSname
		else:
			SNID = self.t['tnsname'][SNindex]
		
		if filt is None:
			filt=self.filt
		
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

	def load_lc(self, SNindex, filt=None, fileformat=None, offsetindex=None, MJDbinsize=None):
		# get lc from already existing file
		filename = self.lcbasename(SNindex=SNindex,filt=filt,offsetindex=offsetindex,MJDbinsize=MJDbinsize)+'.txt' 
		print('Loading lc: %s' % filename)
		self.lc.load_spacesep(filename, delim_whitespace=True)
		if fileformat is None: 
			fileformat = self.cfg.params['output']['fileformat']
		return(0)

	def save_lc(self, SNindex=None, yse=False, TNSname=None, indices=None, filt=None, overwrite=False, fileformat=None, offsetindex=None, MJDbinsize=None):
		# write table and save lc as file
		filename = self.lcbasename(SNindex=SNindex, yse=False, filt=filt, offsetindex=offsetindex, MJDbinsize=MJDbinsize)+'.txt'
		if fileformat is None: 
			fileformat = self.cfg.params['output']['fileformat']
		self.lc.write(filename,indices=indices,overwrite=True,verbose=True)
		#self.lc.write(filename,format=fileformat, overwrite=overwrite,verbose=(self.verbose>0))
		return(0)

	def saveRADEClist(self, SNindex, filt=None):
		RADEClistfilename = self.lcbasename(SNindex=SNindex,filt=filt)+'.RADEClist.txt'
		if self.verbose: 
			print('Saving RADEClist %s' % RADEClistfilename)
		self.RADECtable.write(RADEClistfilename,overwrite=True,verbose=True)
		return(0)

	def loadRADEClist(self, SNindex, filt=None):
		# get RADEClist from alreadt existing file
		RADEClistfilename = self.lcbasename(SNindex=SNindex,filt=filt)+'.RADEClist.txt'
		if self.verbose: 
			print('Loading RADEClist %s' % RADEClistfilename)
		self.RADECtable.load_spacesep(RADEClistfilename, delim_whitespace=True)
		if self.verbose>1:
			print(self.RADECtable.write())
		return(0)

	def makecuts_indices(self,SNindex,offsetindex,procedure1):

		# use when cleaning up data in plot_lc.py or average_lc.py; makes cuts based on mask column created in cleanup_lc.py
		# set flags in precursor.cfg to control what data to cut based on uncertainties and/or chi/N
		if procedure1 is 'averagelc': 
			flags = self.cfg.params['averagelc']['flags']
		elif procedure1 is 'plotlc': 
			flags = self.cfg.params['plotlc']['flags']
		elif procedure1 is 'offsetstats':
			flags = self.cfg.params['offsetstats']['flags']
		else:
			raise RuntimeError('Procedure %s must be averagelc or offsetstats!' % procedure1)
		print('Setting indices using flags: %x' % flags)
		
		mask=np.bitwise_and(self.lc.t['Mask'], flags)
		cuts_indices = np.where(mask==0)
		cuts_indices = cuts_indices[0]
		bad_data = np.where(mask!=0)
		bad_data = bad_data[0]
		datacut = len(self.lc.t)-len(self.lc.t.loc[cuts_indices])
		print('Length original lc: ',len(self.lc.t),', length cleaned lc: ',len(self.lc.t.loc[cuts_indices]),', data points cut: ',datacut)

		lc_uJy = self.lc.t.loc[cuts_indices, self.flux_colname]
		lc_duJy = self.lc.t.loc[cuts_indices, self.dflux_colname]
		lc_MJD = self.lc.t.loc[cuts_indices, 'MJD']
		return(lc_uJy, lc_duJy, lc_MJD, cuts_indices, bad_data)

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
			snlistfilename=args.snlistfilename

		if not(os.path.isfile(snlistfilename)):
			raise RuntimeError("SN list file %s does not exist, exiting!!" % snlistfilename)
			
		self.setoutdir(outrootdir=args.outrootdir, outsubdir=args.outsubdir)
		self.verbose = args.verbose
		self.debug = args.debug
		self.flux_colname = self.cfg.params['flux_colname']
		self.dflux_colname = self.cfg.params['dflux_colname']

		self.RADECtable = pdastroclass(columns=['OffsetID','PatternID','Ra','Dec','RaOffset','DecOffset','Radius','Ndet','Ndet_c','Ndet_o'])
		self.averagelctable = pdastroclass(columns=['OffsetID','MJD',self.flux_colname,self.dflux_colname,'stdev','X2norm','Nused','Nclipped','MJDNused','MJDNskipped'])
		self.RADECtable.default_formatters = {'OffsetID':'{:3d}'.format,
											  'PatternID':'{:2d}'.format,
											  'Ra':'{:.8f}'.format,
											  'Dec':'{:.8f}'.format,
											  'RaOffset':'{:.2f}'.format,
											  'DecOffset':'{:.2f}'.format,
											  'Radius':'{:.2f}'.format,
											  'Ndet':'{:4d}'.format,
											  'Ndet_c':'{:4d}'.format,
											  'Ndet_o':'{:4d}'.format}
		self.averagelctable.default_formatters = {'OffsetID':'{:3d}'.format,
												  'MJD':'{:.5f}'.format,
												  self.flux_colname:'{:.2f}'.format,
												  self.dflux_colname:'{:.2f}'.format,
												  'stdev':'{:.2f}'.format,
												  'X2norm':'{:.3f}'.format,
												  'Nused':'{:4d}'.format,
												  'Nclipped':'{:4d}'.format,
												  'MJDNused':'{:4d}'.format,
												  'MJDNskipped':'{:4d}'.format}

		self.load_spacesep(snlistfilename)
		print(self.t)
		print(args.SNlist)

		SNindexlist = self.getSNlist(args.SNlist)
		return(SNindexlist)
