#!/usr/bin/env python

from SNloop import SNloopclass
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pylab as matlib
from pdastro import AnotB
from matplotlib.backends.backend_pdf import PdfPages

def dataPlot(x, y, dx=None, dy=None, sp=None, label=None, fmt='bo', ecolor='k', elinewidth=None, barsabove = False, capsize=1, logx=False, logy=False, zorder=None):
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
				plot, dplot,dummy = sp.errorbar(x, y, label=label, fmt=fmt, capsize=capsize, barsabove=barsabove)
			else:
				plot, = sp.plot(x, y, fmt)
		return sp, plot, None
	else:
		if logy:
			sp.set_yscale("log", nonposx='clip')
		if logx:
			sp.set_xscale("log", nonposx='clip')
		plot, dplot,dummy = sp.errorbar(x, y, xerr=dx, yerr=dy, label=label, fmt=fmt, ecolor=ecolor, elinewidth=elinewidth, capsize=capsize, barsabove=barsabove, zorder=zorder)
		return sp, plot, dplot

class lightlcclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)

	def savelightweight(self):
		"""
		load original lc
		new lc: only MJD, uJy, duJy, m, dm, and Mask cols
		new lc: only good/usable measurements with only MJD, uJy, duJy, m, dm, and Mask cols

		load avg lc
		new avg lc: only MJD, uJy, duJy, m, dm, Mask, Nused, and Nclip cols
		new avg lc: only good/usable measurements with only MJD, uJy, duJy, m, dm, Mask, Nused, and Nclip cols
		"""

		# original
		self.load_lc(SNindex,filt=self.filt,controlindex=0)
		self.lc.t = self.lc.t.drop(columns=['F','err','chi/N','RA','Dec','x','y','maj','min','phi','apfit','mag5sig','Sky','Obs','c1_mean','c1_mean_err','c1_stdev','c1_stdev_err','c1_X2norm','c1_Nvalid','c1_Nnan','c2_mean','c2_mean_err','c2_stdev','c2_stdev_err','c2_X2norm','c2_Ngood','c2_Nclip','c2_Nmask','c2_Nnan'])
		basename = '%s/%s/lightweight/%s.%s' % (self.outrootdir,self.t['tnsname'][SNindex],self.t['tnsname'][SNindex],self.filt)
		filename = basename + '.light.txt'
		print('Saving original light curve: ',filename)
		self.lc.write(filename,overwrite=True)

		# original clean
		indices = self.getusableindices()
		filename = basename + '.light.clean.txt'
		print('Saving cleaned original light curve: ',filename)
		self.lc.write(filename,indices=indices,overwrite=True)

		# average
		self.load_lc(SNindex,filt=self.filt,controlindex=0,MJDbinsize=args.MJDbinsize)
		self.lc.t = self.lc.t.drop(columns=['MJDbin','stdev','X2norm','Nexcluded'])
		if int(args.MJDbinsize) == args.MJDbinsize:
			basename += '.%ddays' % int(args.MJDbinsize)
		else:
			basename += '.%.2fdays' % args.MJDbinsize
		filename = basename + '.light.txt'
		print('Saving averaged light curve: ',filename)
		self.lc.write(filename,overwrite=True)

		# average clean
		indices = self.getusableindices()
		filename = basename + '.light.clean.txt'
		print('Saving cleaned averaged light curve: ',filename)
		self.lc.write(filename,indices=indices,overwrite=True)

	def addflagdescriptions(self):
		basename = '%s/%s/lightweight' % (self.outrootdir,self.t['tnsname'][SNindex])
		f = open(basename+"/flagdescriptions.txt","w+")
		f.write("We've put flags in the mask column that will tell you which measurements to use. Measurements were classified as bad, questionable, or good based on a variety of factors, including PSF statistics and statistics garnered from control light curve data.")
		f.write("\n\n- The measurements to exclude in the single-measurement original light curve will have one or more of the following flags: 0x100000 (based on PSF statistics) 0x400000 (based on control light curve statistics), and 0x800000 (based on statistics of a measurement's 1-day epoch). If the flag 0x80000 is set, measurements in its 1-day epoch were excluded in the daily averaging, so this measurement is classified as questionable but usable.")
		f.write("\n\n- The measurements to exclude in the averaged light curve will have the following flag: 0x800000 (based on statistics of a measurement's 1-day epoch). Flags for questionable measurements include 0x1000 and 0x2000, both of which also take into account statistics of a measurement's 1-day epoch.")
		f.write("\n\n- In sum, 0x100000, 0x400000, and 0x800000 flag bad measurements; 0x80000, 0x1000, and 0x2000 flag questionable measurements that can be used at the user's discretion.")
		f.close()

	def plotlightweight(self,SNindex,pdf,MJDbinsize=None):
		fig = plt.figure()
		sp = matlib.subplot(111)
		#plt.rcParams['font.family'] = 'serif'
		#plt.rcParams["font.serif"] = 'times'
		for controlindex in range(len(self.RADECtable.t)-1,-1,-1):
			self.load_lc(SNindex,filt=self.filt,controlindex=self.RADECtable.t.at[controlindex,'ControlID'],MJDbinsize=MJDbinsize)
			goodix = self.getgoodindices()
			#badix = self.getbadindices()
			allix = self.lc.getindices()
			badix = AnotB(allix,goodix)
			if controlindex == 0:
				# plot bad data with open red circles
				sp, plotbad, dplotbad = dataPlot(self.lc.t.loc[badix,'MJD'],self.lc.t.loc[badix,self.flux_colname],dy=self.lc.t.loc[badix,self.dflux_colname],sp=sp)
				matlib.setp(plotbad,mfc='white',ms=4,color='r')
				# plot good and usable data with closed red circless
				sp, plotSN, dplotSN = dataPlot(self.lc.t.loc[goodix,'MJD'],self.lc.t.loc[goodix,self.flux_colname],dy=self.lc.t.loc[goodix,self.dflux_colname],sp=sp)
				matlib.setp(plotSN,ms=4,color='r')
				maxlc = max(self.lc.t.loc[goodix,self.flux_colname])
				minlc = min(self.lc.t.loc[goodix,self.flux_colname])
			else:
				# plot bad data with open red circles
				sp, plotControlLCBad, dplotControlLCBad = dataPlot(self.lc.t.loc[badix,'MJD'],self.lc.t.loc[badix,self.flux_colname],dy=self.lc.t.loc[badix,self.dflux_colname],sp=sp)
				matlib.setp(plotControlLCBad,mfc='white',ms=4,color='b')
				# plot good data in closed blue circles
				sp, plotControlLC, dplotControlLC = dataPlot(self.lc.t.loc[goodix,'MJD'],self.lc.t.loc[goodix,self.flux_colname],dy=self.lc.t.loc[goodix,self.dflux_colname],sp=sp)
				matlib.setp(plotControlLC,ms=4,color='b')
		# get control lc legend label
		if len(self.RADECtable.t)>1:
			# control lc PatternID circle gets specific legend
			if max(self.RADECtable.t['PatternID']) == 1:
				if len(self.cfg.params['forcedphotpatterns']['circle']['radii'])==1:
					controlLClabel = '%s %s" Control LCs' % (self.cfg.params['forcedphotpatterns']['circle']['n'], self.cfg.params['forcedphotpatterns']['circle']['radii'][0])
				else:
					controlLClabel = '%s %s" Control LCs and %s %s" Control LCs' % (self.cfg.params['forcedphotpatterns']['circle']['n'], self.cfg.params['forcedphotpatterns']['circle']['radii'][0],self.cfg.params['forcedphotpatterns']['circle']['n'], self.cfg.params['forcedphotpatterns']['circle']['radii'][1])
			# greater/multiple different controlLC PatternIDs get simpler legend
			else:
				controlLCtotal = len(self.RADECtable.t)-1
				controlLClabel = '%d Total Control LCs' % controlLCtotal
		plt.legend((plotSN,plotbad,plotControlLC,plotControlLCBad),('SN %s Accurate & Usable Data' % self.t.at[SNindex,'tnsname'],'SN %s Inaccurate Data' % self.t.at[SNindex,'tnsname'],controlLClabel+' Accurate & Usable Data',controlLClabel+' Inaccurate Data'))
		title = 'SN %s ' % self.t.at[SNindex,'tnsname']
		if self.lctype == 'avg': 
			title += 'Averaged '
		title += 'All Detections'
		plt.title(title)
		plt.axhline(linewidth=1,color='k')
		plt.xlabel('MJD')
		plt.ylabel(self.flux_colname)
		plt.ylim(1.1*minlc,1.1*maxlc)
		basename = '%s/%s/lightweight/%s.%s.alldetections' % (self.outrootdir,self.t['tnsname'][SNindex],self.t.at[SNindex,'tnsname'],self.filt)
		if self.lctype == 'avg':
			if int(args.MJDbinsize) == args.MJDbinsize:
				basename += '.%ddays' % int(args.MJDbinsize)
			else:
				basename += '.%.2fdays' % args.MJDbinsize
		#print('Saving fig "SN %s All Detections" at %s.png' % (self.t.at[SNindex,'tnsname'],basename))
		pdf.savefig(fig)
		#plt.savefig(basename+".png",dpi=200)
		plt.clf()

		fig = plt.figure(1)
		sp = matlib.subplot(111)
		for controlindex in range(len(self.RADECtable.t)-1,-1,-1):
			self.load_lc(SNindex,filt=self.filt,controlindex=self.RADECtable.t.at[controlindex,'ControlID'],MJDbinsize=MJDbinsize)
			goodix = self.getgoodindices()
			if controlindex == 0:
				# plot good and usable data with closed red circless
				sp, plotSN, dplotSN = dataPlot(self.lc.t.loc[goodix,'MJD'],self.lc.t.loc[goodix,self.flux_colname],dy=self.lc.t.loc[goodix,self.dflux_colname],sp=sp)
				matlib.setp(plotSN,ms=4,color='r')
				maxlc = max(self.lc.t.loc[goodix,self.flux_colname])
				minlc = min(self.lc.t.loc[goodix,self.flux_colname])
			else:
				# plot good data in closed blue circles
				sp, plotControlLC, dplotControlLC = dataPlot(self.lc.t.loc[goodix,'MJD'],self.lc.t.loc[goodix,self.flux_colname],dy=self.lc.t.loc[goodix,self.dflux_colname],sp=sp)
				matlib.setp(plotControlLC,ms=4,color='b')
		# get control lc legend label
		if len(self.RADECtable.t)>1:
			# control lc PatternID circle gets specific legend
			if max(self.RADECtable.t['PatternID']) == 1:
				if len(self.cfg.params['forcedphotpatterns']['circle']['radii'])==1:
					controlLClabel = '%s %s" Control LCs' % (self.cfg.params['forcedphotpatterns']['circle']['n'], self.cfg.params['forcedphotpatterns']['circle']['radii'][0])
				else:
					controlLClabel = '%s %s" Control LCs and %s %s" Control LCs' % (self.cfg.params['forcedphotpatterns']['circle']['n'], self.cfg.params['forcedphotpatterns']['circle']['radii'][0],self.cfg.params['forcedphotpatterns']['circle']['n'], self.cfg.params['forcedphotpatterns']['circle']['radii'][1])
			# greater/multiple different controlLC PatternIDs get simpler legend
			else:
				controlLCtotal = len(self.RADECtable.t)-1
				controlLClabel = '%d Total Control LCs' % controlLCtotal
		plt.legend((plotSN,plotControlLC),('SN %s' % self.t.at[SNindex,'tnsname'],controlLClabel))
		title = 'SN %s ' % self.t.at[SNindex,'tnsname']
		if self.lctype == 'avg': 
			title += 'Averaged '
		title += 'Good Detections Only'
		plt.title(title)
		plt.axhline(linewidth=1,color='k')
		plt.xlabel('MJD')
		plt.ylabel(self.flux_colname)
		plt.ylim(1.1*minlc,1.1*maxlc)
		basename = '%s/%s/lightweight/%s.%s.usabledetections' % (self.outrootdir,self.t['tnsname'][SNindex],self.t.at[SNindex,'tnsname'],self.filt)
		if self.lctype == 'avg':
			if int(args.MJDbinsize) == args.MJDbinsize:
				basename += '.%ddays' % int(args.MJDbinsize)
			else:
				basename += '.%.2fdays' % args.MJDbinsize
		#print('Saving fig "SN %s Good Detections Only" at %s.png' % (self.t.at[SNindex,'tnsname'],basename))
		pdf.savefig(fig)
		#plt.savefig(basename+".png",dpi=200)
		plt.clf()

	def plotloop(self,SNindex):
		self.loadRADEClist(SNindex)
		basename = '%s/%s/lightweight/%s.%s.plots' % (self.outrootdir,self.t['tnsname'][SNindex],self.t.at[SNindex,'tnsname'],self.filt)
		pdf = PdfPages(basename+'.pdf')

		# first 2 plots: individal detections: good only, then all detections with bad in open circles
		self.plotlightweight(SNindex,pdf)
		# second 2 plots: averaged detections: good only, then all detections with bad in open circles
		self.plotlightweight(SNindex,pdf,MJDbinsize=args.MJDbinsize)
		
		pdf.close()

if __name__ == '__main__':

	lightlc = lightlcclass()
	parser = lightlc.define_options()
	args = parser.parse_args()

	SNindexlist = lightlc.initialize(args)

	for SNindex in SNindexlist:
		if args.filt is None:
			filtlist = ['o','c']
			print('Looping through c and o filters...')
		else:
			filtlist = [args.filt]
		for filt in filtlist:
			lightlc.filt = filt
			print('### FILT SET: ',lightlc.filt)
			lightlc.loadRADEClist(SNindex,filt=lightlc.filt)
			lightlc.savelightweight()
			lightlc.addflagdescriptions()
			lightlc.plotloop(SNindex)

		# pdf file
		# single and avg lc plots: 1. only good, 2. all detections with bad in open circles