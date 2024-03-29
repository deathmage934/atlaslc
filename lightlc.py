#!/usr/bin/env python

from SNloop import SNloopclass
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pylab as matlib
from pdastro import AnotB
from matplotlib.backends.backend_pdf import PdfPages
import argparse

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['legend.fontsize'] = 'small'

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
		# original
		self.load_lc(SNindex,filt=self.filt,controlindex=0)
		self.lc.t = self.lc.t.drop(columns=['F','err','chi/N','RA','Dec','x','y','maj','min','phi','apfit','mag5sig','Sky','Obs'])
		if '__tmp_SN' in self.lc.t.columns:
			self.lc.t = self.lc.t.drop(columns=['__tmp_SN'])
		if 'c1_mean' in self.lc.t.columns:
			self.lc.t = self.lc.t.drop(columns=['c1_mean','c1_mean_err','c1_stdev','c1_stdev_err','c1_X2norm','c1_Nvalid','c1_Nnan','c2_mean','c2_mean_err','c2_stdev','c2_stdev_err','c2_X2norm','c2_Ngood','c2_Nclip','c2_Nmask','c2_Nnan'])
		basename = '%s/%s/%s/%s.%s' % (self.outrootdir,self.t['tnsname'][SNindex],self.t['tnsname'][SNindex],self.t['tnsname'][SNindex],self.filt)
		filename = basename + '.txt'
		print('Saving original light curve: ',filename)
		self.lc.write(filename,overwrite=True)

		# original clean
		indices = self.getusableindices()
		filename = basename + '.clean.txt'
		print('Saving cleaned original light curve: ',filename)
		self.lc.write(filename,indices=indices,overwrite=True)

		# average
		self.load_lc(SNindex,filt=self.filt,controlindex=0,MJDbinsize=args.MJDbinsize)
		self.lc.t = self.lc.t.drop(columns=['MJDbin','stdev','X2norm','Nexcluded'])
		if int(args.MJDbinsize) == args.MJDbinsize:
			basename += '.%ddays' % int(args.MJDbinsize)
		else:
			basename += '.%.2fdays' % args.MJDbinsize
		filename = basename + '.txt'
		print('Saving averaged light curve: ',filename)
		self.lc.write(filename,overwrite=True)

		# average clean
		indices = self.getusableindices()
		filename = basename + '.clean.txt'
		print('Saving cleaned averaged light curve: ',filename)
		self.lc.write(filename,indices=indices,overwrite=True)

	def addflagdescriptions(self):
		basename = '%s/%s/%s' % (self.outrootdir,self.t['tnsname'][SNindex],self.t['tnsname'][SNindex])
		f = open(basename+"/flagdescriptions.txt","w+")
		f.write("We've put flags in the mask column that will tell you which measurements to use. Measurements were classified as bad, questionable, or good based on a variety of factors, including PSF statistics and statistics garnered from control light curve data.")
		f.write("\n\n- The measurements to exclude in the single-measurement original light curve will have one or more of the following flags: 0x100000 (based on PSF statistics) 0x400000 (based on control light curve statistics), and 0x800000 (based on statistics of a measurement's 1-day epoch). If the flag 0x80000 is set, measurements in its 1-day epoch were excluded in the daily averaging, so this measurement is classified as questionable but usable.")
		f.write("\n\n- The measurements to exclude in the averaged light curve will have the following flag: 0x800000 (based on statistics of a measurement's 1-day epoch). Flags for questionable measurements include 0x1000 and 0x2000, both of which also take into account statistics of a measurement's 1-day epoch.")
		f.write("\n\n- In sum, 0x100000, 0x400000, and 0x800000 flag bad measurements; 0x80000, 0x1000, and 0x2000 flag questionable measurements that can be used at the user's discretion.")
		f.close()

	def setplotlims(self,args,xlims=True,ylims=True,maxlc=None,minlc=None):
		# get limits from args; else, leave as is
		if xlims is True:
			xlim_lower, xlim_upper = plt.xlim()
			if not(args.xlim_lower is None): 
				xlim_lower = args.xlim_lower
			if not(args.xlim_upper is None):
				xlim_upper = args.xlim_upper
			plt.xlim(xlim_lower,xlim_upper)
			print('xlim lower: ',xlim_lower,'. xlim upper: ',xlim_upper,'. ')
		if ylims is True:
			ylim_lower, ylim_upper = plt.ylim()
			if not(maxlc is None):
				ylim_upper = 1.1*maxlc
			if not(minlc is None):
				ylim_lower = 1.1*minlc
			if not(args.ylim_lower is None): 
				ylim_lower = args.ylim_lower
			if not(args.xlim_upper is None):
				ylim_upper = args.ylim_upper
			plt.ylim(ylim_lower,ylim_upper)
			print('ylim lower: ',ylim_lower,'. ylim upper: ',ylim_upper,'. ')

	def plotindepth(self,args,SNindex,pdf,MJDbinsize=None):
		if self.filt == 'c':
			color = 'cyan'
		else:
			color = 'orange'

		fig = plt.figure()
		sp = matlib.subplot(111)
		self.loadRADEClist(args,SNindex)
		for controlindex in range(len(self.RADECtable.t)-1,-1,-1):
			self.load_lc(SNindex,filt=self.filt,controlindex=self.RADECtable.t.at[controlindex,'ControlID'],MJDbinsize=MJDbinsize)
			goodix = self.getgoodindices()
			allix = self.lc.getindices()
			badix = AnotB(allix,goodix)
			if controlindex == 0:
				# plot bad data with open red circles
				sp, plotbad, dplotbad = dataPlot(self.lc.t.loc[badix,'MJD'],self.lc.t.loc[badix,self.flux_colname],dy=self.lc.t.loc[badix,self.dflux_colname],sp=sp)
				matlib.setp(plotbad,mfc='white',ms=4,color=color)
				# plot good and usable data with closed red circless
				sp, plotSN, dplotSN = dataPlot(self.lc.t.loc[goodix,'MJD'],self.lc.t.loc[goodix,self.flux_colname],dy=self.lc.t.loc[goodix,self.dflux_colname],sp=sp)
				matlib.setp(plotSN,ms=4,color=color)
				maxlc = max(self.lc.t.loc[goodix,self.flux_colname])
				minlc = min(self.lc.t.loc[goodix,self.flux_colname])
				maxmjd = max(self.lc.t.loc[goodix,'MJD'])
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
				controlLClabel = '%s %s" Control LCs' % (self.cfg.params['forcedphotpatterns']['circle']['n'],self.cfg.params['forcedphotpatterns']['circle']['radii'][0])
				if not(len(self.cfg.params['forcedphotpatterns']['circle']['radii'])==1):
					controlLClabel += ' and %s %s" Control LCs' % (self.cfg.params['forcedphotpatterns']['circle']['n'],self.cfg.params['forcedphotpatterns']['circle']['radii'][1])
			# greater/multiple different controlLC PatternIDs get simpler legend
			else:
				controlLClabel = '%d Total Control LCs' % (len(self.RADECtable.t)-1)
		plt.legend((plotSN,plotbad,plotControlLC,plotControlLCBad),('SN %s %s-band Good Data' % (self.t.at[SNindex,'tnsname'],self.filt),'SN %s %s-band Bad Data' % (self.t.at[SNindex,'tnsname'],self.filt),controlLClabel+' %s-band Good Data' % self.filt,controlLClabel+' %s-band Bad Data' % self.filt))
		filename = '%s.%s' % (self.t['tnsname'][SNindex],self.filt)
		if not(MJDbinsize is None):
			if int(args.MJDbinsize) == args.MJDbinsize:
				filename += '.%ddays' % int(args.MJDbinsize)
			else:
				filename += '.%.2fdays' % args.MJDbinsize
		filename += '.txt'
		title = 'SN %s, %s-band' % (self.t.at[SNindex,'tnsname'], self.filt)
		if self.lctype == 'avg': 
			title += ', Averaged'
		title += '\nFilename: %s' % (filename)
		plt.title(title)
		print('Adding plot to PDF: "%s"' % title)
		plt.axhline(linewidth=1,color='k')
		plt.xlabel('MJD')
		plt.ylabel(self.flux_colname)
		if not(args.deltat is None):
			xlim_lower, xlim_upper = plt.xlim()
			print('xlim_lower set: ',maxmjd-args.deltat)
			plt.xlim(maxmjd-args.deltat,xlim_upper)
		else:
			self.setplotlims(args,ylims=False)
		plt.ylim(1.1*minlc,1.1*maxlc)
		pdf.savefig(fig)
		plt.clf()

		fig = plt.figure()
		sp = matlib.subplot(111)
		for controlindex in range(len(self.RADECtable.t)-1,-1,-1):
			self.load_lc(SNindex,filt=self.filt,controlindex=self.RADECtable.t.at[controlindex,'ControlID'],MJDbinsize=MJDbinsize)
			goodix = self.getgoodindices()
			if controlindex == 0:
				# plot good and usable data with closed red circless
				sp, plotSN, dplotSN = dataPlot(self.lc.t.loc[goodix,'MJD'],self.lc.t.loc[goodix,self.flux_colname],dy=self.lc.t.loc[goodix,self.dflux_colname],sp=sp)
				matlib.setp(plotSN,ms=4,color=color)
				maxlc = max(self.lc.t.loc[goodix,self.flux_colname])
				minlc = min(self.lc.t.loc[goodix,self.flux_colname])
				maxmjd = max(self.lc.t.loc[goodix,'MJD'])
			else:
				# plot good data in closed blue circles
				sp, plotControlLC, dplotControlLC = dataPlot(self.lc.t.loc[goodix,'MJD'],self.lc.t.loc[goodix,self.flux_colname],dy=self.lc.t.loc[goodix,self.dflux_colname],sp=sp)
				matlib.setp(plotControlLC,ms=4,color='b')
		# get control lc legend label
		if len(self.RADECtable.t)>1:
			# control lc PatternID circle gets specific legend
			if max(self.RADECtable.t['PatternID']) == 1:
				controlLClabel = '%s %s" Control LCs' % (self.cfg.params['forcedphotpatterns']['circle']['n'],self.cfg.params['forcedphotpatterns']['circle']['radii'][0])
				if not(len(self.cfg.params['forcedphotpatterns']['circle']['radii'])==1):
					controlLClabel += ' and %s %s" Control LCs' % (self.cfg.params['forcedphotpatterns']['circle']['n'],self.cfg.params['forcedphotpatterns']['circle']['radii'][1])
			# greater/multiple different controlLC PatternIDs get simpler legend
			else:
				controlLClabel = '%d Total Control LCs' % (len(self.RADECtable.t)-1)
		plt.legend((plotSN,plotControlLC),('SN %s %s-band'%(self.t.at[SNindex,'tnsname'],self.filt),controlLClabel+' %s-band'%self.filt))
		title = 'SN %s, %s-band' % (self.t.at[SNindex,'tnsname'], self.filt)
		if self.lctype == 'avg': 
			title += ', Averaged'
		filename = '%s.%s' % (self.t['tnsname'][SNindex],self.filt)
		if not(MJDbinsize is None):
			if int(args.MJDbinsize) == args.MJDbinsize:
				filename += '.%ddays' % int(args.MJDbinsize)
			else:
				filename += '.%.2fdays' % args.MJDbinsize
		filename += '.clean.txt'
		title += ', Good Detections Only\nFilename: %s' % filename
		plt.title(title)
		print('Adding plot to PDF: "%s"' % title)
		plt.axhline(linewidth=1,color='k')
		plt.xlabel('MJD')
		plt.ylabel(self.flux_colname)
		if not(args.deltat is None):
			xlim_lower, xlim_upper = plt.xlim()
			print('xlim_lower set: ',maxmjd-args.deltat)
			plt.xlim(maxmjd-args.deltat,xlim_upper)
		else:
			self.setplotlims(args,ylims=False)
		plt.ylim(1.1*minlc,1.1*maxlc)
		pdf.savefig(fig)
		plt.clf()

		return pdf

	def plotloop(self,args,SNindex):
		pdf = PdfPages('%s/%s/%s/%s_plots.pdf' % (self.outrootdir,self.t['tnsname'][SNindex],self.t['tnsname'][SNindex],self.t.at[SNindex,'tnsname']))

		# summary plot: not clean, both c and o
		fig = plt.figure()
		plt.clf()
		sp = matlib.subplot(111)
		self.loadRADEClist(args,SNindex)
		# load c band
		self.load_lc(SNindex,filt='c',controlindex=0,MJDbinsize=None)
		sp,plotc,dplot = dataPlot(self.lc.t['MJD'],self.lc.t[self.flux_colname],dy=self.lc.t[self.dflux_colname],sp=sp)
		matlib.setp(plotc,ms=4,color='cyan')
		# load o band
		self.load_lc(SNindex,filt='o',controlindex=0,MJDbinsize=None)
		sp,ploto,dplot = dataPlot(self.lc.t['MJD'],self.lc.t[self.flux_colname],dy=self.lc.t[self.dflux_colname],sp=sp)
		matlib.setp(ploto,ms=4,color='orange')
		# set x and y lims
		goodix = self.getgoodindices()
		maxlc = max(self.lc.t.loc[goodix,self.flux_colname])
		minlc = min(self.lc.t.loc[goodix,self.flux_colname])
		if max(self.lc.t.loc[goodix,self.flux_colname]) > maxlc:
			maxlc = max(self.lc.t.loc[goodix,self.flux_colname])
		if min(self.lc.t.loc[goodix,self.flux_colname]) < minlc:
			minlc = min(self.lc.t.loc[goodix,self.flux_colname])
		self.setplotlims(args,maxlc=maxlc,minlc=minlc)
		plt.ylim(1.1*minlc,1.1*maxlc)
		# rest of plot
		plt.title('SN %s, c- and o-band' % self.t.at[SNindex,'tnsname'])
		plt.axhline(linewidth=1,color='k')
		plt.xlabel('MJD')
		plt.ylabel(self.flux_colname)
		plt.legend((plotc,ploto),('c-band','o-band'))
		print('Adding summary plot to PDF: "SN %s, c- and o-band"' % self.t.at[SNindex,'tnsname'])
		pdf.savefig(fig)
		plt.clf()

		# summary plot 2: avg sn lc: o good, o bad, c good, c bad
		fig = plt.figure()
		sp = matlib.subplot(111)
		legend_vars = ()
		for filt in ['c','o']:
			self.load_lc(SNindex,filt=filt,controlindex=0,MJDbinsize=args.MJDbinsize)
			goodix = self.getgoodindices()
			badix = AnotB(self.lc.getindices(),goodix)
			if filt == 'c':
				color = 'cyan'
			else:
				color = 'orange'
			sp, plotbad, dplotbad = dataPlot(self.lc.t.loc[badix,'MJD'],self.lc.t.loc[badix,self.flux_colname],dy=self.lc.t.loc[badix,self.dflux_colname],sp=sp)
			matlib.setp(plotbad,mfc='white',ms=4,color=color)
			sp, plotgood, dplotgood = dataPlot(self.lc.t.loc[goodix,'MJD'],self.lc.t.loc[goodix,self.flux_colname],dy=self.lc.t.loc[goodix,self.dflux_colname],sp=sp)
			matlib.setp(plotgood,ms=4,color=color)
			sp,plot,dplot = dataPlot(self.lc.t['MJD'],self.lc.t[self.flux_colname],dy=self.lc.t[self.dflux_colname],sp=sp)
			matlib.setp(plot,ms=4,color=color)
			legend_vars += (plotbad, plotgood)
		# set x and y lims
		goodix = self.getgoodindices()
		maxlc = max(self.lc.t.loc[goodix,self.flux_colname])
		minlc = min(self.lc.t.loc[goodix,self.flux_colname])
		if max(self.lc.t.loc[goodix,self.flux_colname]) > maxlc:
			maxlc = max(self.lc.t.loc[goodix,self.flux_colname])
		if min(self.lc.t.loc[goodix,self.flux_colname]) < minlc:
			minlc = min(self.lc.t.loc[goodix,self.flux_colname])
		self.setplotlims(args,maxlc=maxlc,minlc=minlc)
		plt.ylim(1.1*minlc,1.1*maxlc)
		# rest of plot
		plt.title('SN %s, c- and o-band, Averaged' % self.t.at[SNindex,'tnsname'])
		plt.axhline(linewidth=1,color='k')
		plt.xlabel('MJD')
		plt.ylabel(self.flux_colname)
		plt.legend(legend_vars,('c-band Bad Data','c-band Good Data','o-band Bad Data','o-band Good Data'))
		print('Adding summary plot to PDF: "SN %s, c- and o-band"' % self.t.at[SNindex,'tnsname'])
		pdf.savefig(fig)
		plt.clf()

		for filt in ['c','o']:
			self.filt = filt
			print('### FILT SET: ',self.filt)
			if filt == 'c':
				color = 'cyan'
			else:
				color = 'orange'

			# first 2 plots: individal detections: good only, then all detections with bad in open circles
			pdf = self.plotindepth(args,SNindex,pdf)

			# second 2 plots: averaged detections: good only, then all detections with bad in open circles
			pdf = self.plotindepth(args,SNindex,pdf,MJDbinsize=args.MJDbinsize)

			# baseline plot 
			fig = plt.figure()
			sp = matlib.subplot(111)
			self.loadRADEClist(args,SNindex)
			for controlindex in range(len(self.RADECtable.t)-1,-1,-1):
				self.load_lc(SNindex,filt=self.filt,controlindex=self.RADECtable.t.at[controlindex,'ControlID'])
				goodix = self.getgoodindices()
				allix = self.lc.getindices()
				badix = AnotB(allix,goodix)
				if controlindex == 0:
					# plot bad data with open red circles
					sp, plotbad, dplotbad = dataPlot(self.lc.t.loc[badix,'MJD'],self.lc.t.loc[badix,self.flux_colname],dy=self.lc.t.loc[badix,self.dflux_colname],sp=sp)
					matlib.setp(plotbad,mfc='white',ms=4,color=color)
					# plot good and usable data with closed red circless
					sp, plotSN, dplotSN = dataPlot(self.lc.t.loc[goodix,'MJD'],self.lc.t.loc[goodix,self.flux_colname],dy=self.lc.t.loc[goodix,self.dflux_colname],sp=sp)
					matlib.setp(plotSN,ms=4,color=color)
					maxlc = max(self.lc.t.loc[goodix,self.flux_colname])
					minlc = min(self.lc.t.loc[goodix,self.flux_colname])
					maxmjd = max(self.lc.t.loc[goodix,'MJD'])
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
					controlLClabel = '%s %s" Control LCs' % (self.cfg.params['forcedphotpatterns']['circle']['n'],self.cfg.params['forcedphotpatterns']['circle']['radii'][0])
					if not(len(self.cfg.params['forcedphotpatterns']['circle']['radii'])==1):
						controlLClabel += ' and %s %s" Control LCs' % (self.cfg.params['forcedphotpatterns']['circle']['n'],self.cfg.params['forcedphotpatterns']['circle']['radii'][1])
				# greater/multiple different controlLC PatternIDs get simpler legend
				else:
					controlLClabel = '%d Total Control LCs' % (len(self.RADECtable.t)-1)
			plt.legend((plotSN,plotbad,plotControlLC,plotControlLCBad),('SN %s %s-band Good Data' % (self.t.at[SNindex,'tnsname'],self.filt),'SN %s %s-band Bad Data' % (self.t.at[SNindex,'tnsname'],self.filt),controlLClabel+' %s-band Good Data' % self.filt,controlLClabel+' %s-band Bad Data' % self.filt))
			plt.title('SN %s, %s-band, Zoomed to Baseline' % (self.t.at[SNindex,'tnsname'], self.filt))
			print('Adding plot to PDF: "SN %s, %s-band, Zoomed to Baseline"' % (self.t.at[SNindex,'tnsname'], self.filt))
			plt.axhline(linewidth=1,color='k')
			plt.xlabel('MJD')
			plt.ylabel(self.flux_colname)
			self.setplotlims(args,ylims=False)
			plt.ylim(-200,200)
			pdf.savefig(fig)
			plt.clf()


		pdf.close()

if __name__ == '__main__':

	lightlc = lightlcclass()

	parser = argparse.ArgumentParser()
	parser.add_argument('-t','--deltat', help="lookback time in days for filter-separated plots", default=None, type=float)
	#args = parser.parse_args()
	parser = lightlc.define_options(parser=parser)
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
			lightlc.loadRADEClist(args,SNindex,filt=lightlc.filt)
			lightlc.savelightweight()
			lightlc.addflagdescriptions()
		lightlc.plotloop(args,SNindex)

		# pdf file
		# single and avg lc plots: 1. only good, 2. all detections with bad in open circles