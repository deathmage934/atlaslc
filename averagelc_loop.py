#!/usr/bin/env python

import numpy as np
import math
import sys,socket,os,re
import astropy.table as at
from astropy.io import ascii
from astropy.io.votable.tree import Values as v
from astropy import units as u
from astropy.coordinates import Angle
import matplotlib.pyplot as plt
import pylab as matlib
from datetime import datetime as dt
import argparse
from tools import yamlcfgclass,rmfile,RaInDeg,DecInDeg
import astrotable
from astrotable import astrotableclass
from download_atlas_lc import download_atlas_lc_class
from SNloop import SNloopclass
from download_atlas_lc_loop import downloadloopclass
import sigmacut

class averagelcloopclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)

	def calcaveragelc(self, SNindex, MJDbinsize=None, offsetindex=None):
		self.load_lc(SNindex, filt=self.filt, offsetindex=offsetindex, MJDbinsize=None)
		self.averagelctable.formattable(formatMapping={'OffsetID':'3d','MJD':'.5f','uJy':'.2f','duJy':'.2f','stdev':'.2f','X2norm':'.3f'})
		print('Length of self.lc.t: ',len(self.lc.t))
		if len(self.lc.t) == 0:
			return(1)

		MJD = int(np.amin(self.lc.t['MJD']))
		MJDmax = np.amax(self.lc.t['MJD'])
		print('MJD: ',MJD)

		while MJD <= MJDmax:
			windowindeces = np.logical_and(self.lc.t['MJD']>=MJD, self.lc.t['MJD']<(MJD+MJDbinsize))
			uJy_windowdata = self.lc.t['uJy'][windowindeces]
			print(uJy_windowdata)
			duJy_windowdata = self.lc.t['duJy'][windowindeces]
			print(duJy_windowdata)
			MJD_windowdata = self.lc.t['MJD'][windowindeces]
			print(MJD_windowdata)

			calcaverage=sigmacut.calcaverageclass()
			calcaverage.calcaverage_sigmacutloop(uJy_windowdata,noise=duJy_windowdata,verbose=2,Nsigma=3.0,median_firstiteration=True,saveused=True)

			if calcaverage.Nused<=0:
				print('No data values. Skipping...')
				MJD += MJDbinsize
				continue

			print("mean:%f (uncertainty:%f) " % (calcaverage.mean,calcaverage.mean_err))
			lcaverageindex = len(self.averagelctable.t)
			self.averagelctable.t.add_row({'OffsetID':self.RADECtable.t['OffsetID'][offsetindex],'uJy':calcaverage.mean,'duJy':calcaverage.mean_err, 'stdev':calcaverage.stdev, 'X2norm':calcaverage.X2norm, 'Nused':calcaverage.Nused, 'Nclipped':calcaverage.Nskipped})
			
			print(';bbbbbb',lcaverageindex,calcaverage.use)
			mask=np.logical_not(calcaverage.use)
			calcaverage.calcaverage_sigmacutloop(MJD_windowdata,noise=duJy_windowdata,mask=mask,verbose=2,Nsigma=0.0,median_firstiteration=True,saveused=True)
			self.averagelctable.t['MJD'][lcaverageindex] = calcaverage.mean
			self.averagelctable.t['MJDNused'][lcaverageindex] = calcaverage.Nused
			self.averagelctable.t['MJDNskipped'][lcaverageindex] = calcaverage.Nskipped

			MJD += MJDbinsize

		print(self.averagelctable.t)
		return(0)
	
	def averagelcs(self, SNindex, MJDbinsize=None,fileformat=None,overwrite=True):
		for offsetindex in range(len(self.RADECtable.t)):
			if offsetindex>0:
				self.averagelctable.t.remove_rows(slice(0,len(self.averagelctable.t)))
			result = self.calcaveragelc(SNindex,offsetindex=offsetindex,MJDbinsize=MJDbinsize)

			filename = self.lcbasename(SNindex,filt=self.filt,offsetindex=offsetindex,MJDbinsize=MJDbinsize)+'.txt'
			if fileformat is None: 
				fileformat = self.cfg.params['output']['fileformat']

			if result == 0: 
				self.averagelctable.write(filename,format=fileformat,overwrite=True,verbose=(self.verbose>0))
			else:
				print('Removing file because length of self.lc.t is 0...')
				rmfile(filename)

if __name__ == '__main__':

	averagelc = averagelcloopclass()
	parser = averagelc.define_options()
	args = parser.parse_args()

	averagelc.loadcfgfiles(args.cfgfile,
						filt=args.filt,
						extracfgfiles=args.extracfgfile,
						params=args.params,
						params4all=args.pall,
						params4sections=args.pp,
						verbose=args.verbose)

	averagelc.setoutdir(outrootdir=args.outrootdir, outsubdir=args.outsubdir)
	averagelc.verbose = args.verbose
	averagelc.debug = args.debug

	averagelc.load_spacesep(args.snlistfilename)
	print(args.SNlist)
	indexlist = averagelc.getSNlist(args.SNlist)

	MJDbinsize = averagelc.cfg.params['output']['MJDbinsize']
	if not(args.MJDbinsize is None):
		MJDbinsize = args.MJDbinsize

	for i in indexlist:
		averagelc.loadRADEClist(i)
		averagelc.averagelcs(i, MJDbinsize=MJDbinsize)
