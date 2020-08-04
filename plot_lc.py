#!/usr/bin/env python

import numpy as np
import math
import sys,socket,os,re
import astropy.table as at
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import Angle
import matplotlib.pyplot as plt
import pylab as matlib
from datetime import datetime as dt
import argparse
from tools import yamlcfgclass
from tools import rmfile
from tools import RaInDeg
from tools import DecInDeg
import astrotable
from astrotable import astrotableclass
from download_atlas_lc import download_atlas_lc_class
from SNloop import SNloopclass
from download_atlas_lc_loop import downloadloopclass
import sigmacut

def dataPlot(x, y, dx=None, dy=None, sp=None, fmt='bo', ecolor='k', elinewidth=None, barsabove = False, capsize=1, logx=False, logy=False):
	if sp == None:
		sp = matlib.subplot(111)
	if dx is None and dy is None:
		if logy:
			if logx:
				plot, = sp.loglog(x, y, fmt)
			else:
				 plot, = sp.semilogy(x, y, fmt)
		elif logx:
			plot, = sp.semilogx(x, y, fmt)
		else:
			if barsabove:
				plot, dplot,dummy = sp.errorbar(x, y, fmt=fmt, capsize=capsize, barsabove=barsabove)
			else:
				plot, = sp.plot(x, y, fmt)
		return sp, plot, None
	else:
		if logy:
			sp.set_yscale("log", nonposx='clip')
		if logx:
			sp.set_xscale("log", nonposx='clip')
		plot, dplot,dummy = sp.errorbar(x, y, xerr=dx, yerr=dy, fmt=fmt, ecolor=ecolor, elinewidth=elinewidth, capsize=capsize, barsabove=barsabove)
		return sp, plot, dplot

class plotlcclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)

	def plot_lc(self,SNindex,filt=None,sp=None,plotoffsetlcflag=False):
		print('Plotting SN lc and offsets...')

		MJD_offsetlc = []
		uJy_offsetlc = []
		duJy_offsetlc = []

		if sp is None:
			sp = matlib.subplot(111)

		self.loadRADEClist(SNindex)

		for i in range(len(self.RADECtable.t)-1,-1,-1):
			self.load_lc(SNindex, filt=self.filt, offsetindex=self.RADECtable.t['OffsetID'][i])
			print('Offset index: ',i)

			if i==0:
				(sp, plotSN, dplotSN)=dataPlot(self.lc.t['MJD'], self.lc.t['uJy'], dy=self.lc.t['duJy'],sp=sp)
				matlib.setp(plotSN,ms=5,color='r')
			else:
				MJD_offsetlc.extend(self.lc.t['MJD'])
				uJy_offsetlc.extend(self.lc.t['uJy'])
				duJy_offsetlc.extend(self.lc.t['duJy'])
				(sp, plotOffset, dplotOffset)=dataPlot(self.lc.t['MJD'],self.lc.t['uJy'],dy=self.lc.t['duJy'],sp=sp)
				matlib.setp(plotOffset,ms=5,color='b')

		if len(self.cfg.params['forcedphotpatterns']['circular']['radii'])==1:
			plt.legend(('SN', '%s %s" Offset' % (self.cfg.params['forcedphotpatterns']['circular']['n'], self.cfg.params['forcedphotpatterns']['circular']['radii'][0])))
		else:
			plt.legend(('SN', '%s %s" Offset and %s %s" Offset' % (self.cfg.params['forcedphotpatterns']['circular']['n'], self.cfg.params['forcedphotpatterns']['circular']['radii'][0],self.cfg.params['forcedphotpatterns']['circular']['n'], self.cfg.params['forcedphotpatterns']['circular']['radii'][1])))
		
		plt.axhline(linewidth=1,color='k')
		#plt.xlim(59000,59050)
		#plt.ylim(-5000,5000)
		plt.xlabel('MJD')
		plt.ylabel('uJy')

if __name__ == '__main__':

	plot_lc = plotlcclass()
	parser = plot_lc.define_options()
	args = parser.parse_args()

	indexlist = plot_lc.initialize(args)

	for i in indexlist:
		plot_lc.plot_lc(i)
		plotfilename = plot_lc.lcbasename(i)+'.png'
		print('Plot file name: ',plotfilename)
		plt.savefig(plotfilename)
