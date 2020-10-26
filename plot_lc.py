#!/usr/bin/env python

from SNloop import SNloopclass
import pylab as matlib
import matplotlib.pyplot as plt

def dataPlot(x, y, dx=None, dy=None, sp=None, label=None, fmt='bo', ecolor='k', elinewidth=None, barsabove = False, capsize=1, logx=False, logy=False):
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
		plot, dplot,dummy = sp.errorbar(x, y, xerr=dx, yerr=dy, label=label, fmt=fmt, ecolor=ecolor, elinewidth=elinewidth, capsize=capsize, barsabove=barsabove)
		return sp, plot, dplot

class plotlcclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)

		self.flag_uncertainty = 0x1
		self.flag_dynamic = 0x2
		self.flag_static = 0x4

	def plot_lc(self,args,SNindex,sp=None):
		print('Plotting SN lc and offsets...')
		plt.figure()
		if sp is None:
			sp = matlib.subplot(111)
		self.loadRADEClist(SNindex)
		totalcuts = 0

		# plot SN in red first, then loop through offsets to plot in blue
		for offsetindex in range(len(self.RADECtable.t)-1,-1,-1):
			self.load_lc(SNindex, filt=self.filt, offsetindex=self.RADECtable.t.at[offsetindex,'OffsetID'])
			if self.verbose>2: # FIX
				print('Offset index: ',offsetindex)

			# check if make cuts on data to plot
			makecuts_apply = self.cfg.params['plotlc']['makecuts']
			if not(args.avg_makecuts) is None:
				if args.avg_makecuts is True:
					makecuts_apply = True
				else:
					makecuts_apply = False
			if makecuts_apply == True:
				if not('Mask' in self.lc.t.columns):
					raise RuntimeError('No "Mask" column exists! Please run "cleanup_lc.py %s" beforehand.' % self.t.at[SNindex,'tnsname'])
				# make cuts on lc, save bad data to possibly plot
				lc_uJy, lc_duJy, lc_MJD, cuts_indices, bad_data = self.makecuts_indices(SNindex, offsetindex=offsetindex, procedure1='plotlc')
				lc_uJy_bad = self.lc.t.loc[bad_data, self.flux_colname]
				lc_duJy_bad = self.lc.t.loc[bad_data, self.dflux_colname]
				lc_MJD_bad = self.lc.t.loc[bad_data, 'MJD']
			else:
				print('Skipping makecuts using mask column...')
				lc_uJy = self.lc.t[self.flux_colname]
				lc_duJy = self.lc.t[self.dflux_colname]
				lc_MJD = self.lc.t['MJD']

			if offsetindex==0:
				if self.cfg.params['plotlc']['plot_bad_data'] is True:
					# plot bad data with open red circles
					sp, plotbad, dplotbad = dataPlot(lc_MJD_bad, lc_uJy_bad, dy=lc_duJy_bad)
					matlib.setp(plotbad,mfc='white',ms=4,color='r')
					# plot good data in closed red circles
					sp, plotSN, dplotSN = dataPlot(lc_MJD, lc_uJy, dy=lc_duJy)
					matlib.setp(plotSN,ms=4,color='r')

					# fix?
					plt.ylim(min(lc_uJy)*1.1,max(lc_uJy)*1.1)
				else:
					sp, plotSN, dplotSN = dataPlot(lc_MJD, lc_uJy, dy=lc_duJy)
					matlib.setp(plotSN,ms=4,color='r')
			elif len(self.RADECtable.t)==1:
				print('No offsets, skipping plotting...')
			else: 
				if self.cfg.params['plotlc']['plot_offset_data'] is True:
					sp, plotOffset, dplotOffset = dataPlot(lc_MJD,lc_uJy,dy=lc_duJy,sp=sp)
					matlib.setp(plotOffset,ms=4,color='b')

		# determine legend
		# check if offset data plotted
		if (len(self.RADECtable.t)>1) and (self.cfg.params['plotlc']['plot_offset_data'] is True):
			# offset PatternID circle gets specific legend
			if max(self.RADECtable.t['PatternID']) == 1:
				if len(self.cfg.params['forcedphotpatterns']['circle']['radii'])==1:
					offsetlabel = '%s %s" Offset' % (self.cfg.params['forcedphotpatterns']['circle']['n'], self.cfg.params['forcedphotpatterns']['circle']['radii'][0])
				else:
					offsetlabel = '%s %s" Offset and %s %s" Offset' % (self.cfg.params['forcedphotpatterns']['circle']['n'], self.cfg.params['forcedphotpatterns']['circle']['radii'][0],self.cfg.params['forcedphotpatterns']['circle']['n'], self.cfg.params['forcedphotpatterns']['circle']['radii'][1])
			# greater/multiple different offset PatternIDs get simpler legend
			else:
				offsettotal = len(self.RADECtable.t)-1
				offsetlabel = '%d Total Offsets' % offsettotal

			# check if bad data plotted
			if self.cfg.params['plotlc']['plot_bad_data'] is True:
				plt.legend((plotSN,plotOffset,plotbad),('SN %s cleaned' % self.t.at[SNindex,'tnsname'],offsetlabel,'SN %s cut data' % self.t.at[SNindex,'tnsname']))
			else:
				plt.legend((plotSN,plotOffset),('SN %s' % self.t.at[SNindex,'tnsname'],offsetlabel))
		else:
			# check if bad data plotted
			if self.cfg.params['plotlc']['plot_bad_data'] is True:
				plt.legend((plotSN,plotbad),('SN %s cleaned' % self.t.at[SNindex,'tnsname'],'SN %s cut data' % self.t.at[SNindex,'tnsname']))

		plt.title('SN %s' % self.t.at[SNindex,'tnsname'])
		plt.axhline(linewidth=1,color='k')
		plt.xlabel('MJD')
		plt.ylabel(self.flux_colname)
		#if not(len(lc_MJD)==0):
			#plt.ylim(min(lc_uJy)*1.1,max(lc_uJy)*1.1)

		plotfilename = self.lcbasename(SNindex=SNindex)+'.png'
		print('Plot file name: ',plotfilename)
		plt.savefig(plotfilename)

	def plot_lc_offsetstats(self,args,SNindex,sp=None,o1_flag=False,o2_flag=False):
		if sp is None:
			sp = matlib.subplot(111)
		self.load_lc(SNindex, filt=self.filt, offsetindex=0)

		makecuts_apply = self.cfg.params['plotlc']['makecuts']
		if not(args.avg_makecuts) is None:
			if args.avg_makecuts is True:
				makecuts_apply = True
			else:
				makecuts_apply = False
		if makecuts_apply is True:
			if not('Mask' in self.lc.t.columns):
				raise RuntimeError('No "Mask" column exists! Please run "cleanup_lc.py %s" beforehand.' % self.t.at[SNindex,'tnsname'])
			lc_uJy, lc_duJy, lc_MJD, cuts_indices, bad_data = self.makecuts_indices(SNindex, offsetindex=0, procedure1='plotlc')
			lc_uJy_bad = self.lc.t.loc[bad_data, self.flux_colname]
			lc_duJy_bad = self.lc.t.loc[bad_data, self.dflux_colname]
			lc_MJD_bad = self.lc.t.loc[bad_data, 'MJD']
		else:
			print('Skipping makecuts using mask column...')
			lc_uJy = self.lc.t[self.flux_colname]
			lc_duJy = self.lc.t[self.dflux_colname]
			lc_MJD = self.lc.t['MJD']

		# plot SN lc and o1 offsetstats
		if o1_flag is True:
			plt.figure(3)
			o1_stddev1 = self.lc.t['o1_mean']+self.lc.t['o1_stddev']
			o1_stddev2 = self.lc.t['o1_mean']-self.lc.t['o1_stddev']

			if makecuts_apply is True:
				sp, plot2, dplot2 = dataPlot(lc_MJD_bad, lc_uJy_bad, dy=lc_duJy_bad)
				matlib.setp(plot2,mfc='white',ms=4,color='r')
				sp, plot, dplot = dataPlot(lc_MJD, lc_uJy, dy=lc_duJy)
				matlib.setp(plot,ms=4,color='r')
			else:
				sp, plot, dplot = dataPlot(lc_MJD, lc_uJy, dy=lc_duJy)
				matlib.setp(plot,ms=4,color='r')
			plt.fill_between(self.lc.t['MJD'],o1_stddev1,o1_stddev2)
			# plot data points in addition to fill_between
			'''
			sp, o1_stddev1, dplot_o1_stddev1 = dataPlot(lc_MJD,o1_stddev1)
			matlib.setp(o1_stddev1,ms=4,color='k')
			sp, o1_stddev2, dplot_o1_stddev2 = dataPlot(lc_MJD,o1_stddev2)
			matlib.setp(o1_stddev2,ms=4,color='k')
			'''

			plt.title('%s with cleaned offset data' % self.t.at[SNindex,'tnsname'])
			plt.xlabel('MJD')
			plt.ylabel('uJy')
			if makecuts_apply is True:
				plot_legend = self.t.at[SNindex,'tnsname']+' cleaned'
				plot2_legend = self.t.at[SNindex,'tnsname']+' cut data'
				plt.legend((plot,plot2,o1_stddev1,o1_stddev2),(plot_legend,plot2_legend,'o1_stddev1','o1_stddev2'))
			else:
				plt.legend((plot,o1_stddev1,o1_stddev2),(self.t.at[SNindex,'tnsname'],'o1_stddev1','o1_stddev2'))
			plt.axhline(linewidth=1,color='k')

			plotfilename = self.lcbasename(SNindex=SNindex)+'.mask4mjd.png'
			print('Plot file name: ',plotfilename)
			plt.savefig(plotfilename)

		# plot SN lc and o2 offsetstats
		if o2_flag is True:
			plt.figure(4)
			o2_stddev1 = self.lc.t['o2_mean']+self.lc.t['o2_stddev']
			o2_stddev2 = self.lc.t['o2_mean']-self.lc.t['o2_stddev']

			if makecuts_apply is True:
				sp, plot2, dplot2 = dataPlot(lc_MJD_bad, lc_uJy_bad, dy=lc_duJy_bad)
				matlib.setp(plot2,mfc='white',ms=4,color='r')
				sp, plot, dplot = dataPlot(lc_MJD, lc_uJy, dy=lc_duJy)
				matlib.setp(plot,ms=4,color='r')
			else:
				sp, plot, dplot = dataPlot(lc_MJD, lc_uJy, dy=lc_duJy)
				matlib.setp(plot,ms=4,color='r')
			plt.fill_between(self.lc.t['MJD'],o2_stddev1,o2_stddev2)
			# plot data points in addition to fill_between
			'''
			sp, o2_stddev1, dplot_o2_stddev1 = dataPlot(lc_MJD,o2_stddev1)
			matlib.setp(o2_stddev1,ms=4,color='k')
			sp, o2_stddev2, dplot_o2_stddev2 = dataPlot(lc_MJD,o2_stddev2)
			matlib.setp(o2_stddev2,ms=4,color='k')
			'''

			plt.title('%s with offset nans cut' % self.t.at[SNindex,'tnsname'])
			plt.xlabel('MJD')
			plt.ylabel('uJy')
			if makecuts_apply is True:
				plot_legend = self.t.at[SNindex,'tnsname']+' cleaned'
				plot2_legend = self.t.at[SNindex,'tnsname']+' cut data'
				plt.legend((plot,plot2,o2_stddev1,o2_stddev2),(plot_legend,plot2_legend,'o2_stddev1','o2_stddev2'))
			else:
				plt.legend((plot,o2_stddev1,o2_stddev2),(self.t.at[SNindex,'tnsname'],'o2_stddev1','o2_stddev2'))
			plt.axhline(linewidth=1,color='k')

			plotfilename = self.lcbasename(SNindex=SNindex)+'.mask_nan.png'
			print('Plot file name: ',plotfilename)
			plt.savefig(plotfilename)

	def plotlcloop(self,args,SNindex):
		self.plot_lc(args,SNindex)
		
		o1_flag = False
		o2_flag = False
		if self.cfg.params['plotlc']['plot_mask4mjd'] is True:
			o1_flag = True
		if self.cfg.params['plotlc']['plot_mask_nan'] is True:
			o2_flag = True

		if (o1_flag is True) and (o2_flag is True):
			print('mask4mjd and mask_nan data detected, plotting mask4mjd and mask_nan offsetstats...')
			self.plot_lc_offsetstats(args,SNindex,o1_flag=True,o2_flag=True)
		elif o1_flag is True:
			print('mask4mjd data detected, plotting mask4mjd offsetstats...')
			self.plot_lc_offsetstats(args,SNindex,o1_flag=True)
		elif o2_flag is True:
			print('mask_nan data detected, plotting mask_nan offsetstats...')
			self.plot_lc_offsetstats(args,SNindex,o2_flag=True)
		else:
			print('No offsetstats data detected!! Please run offsetstats.py first.')

if __name__ == '__main__':

	plotlc = plotlcclass()
	parser = plotlc.define_options()
	args = parser.parse_args()

	SNindexlist = plotlc.initialize(args)

	for SNindex in SNindexlist:
		SNindex_counter = 0
		print(SNindex,plotlc.t.at[SNindex,'tnsname']) # delete me
		plotlc.plotlcloop(args,SNindex)

	print('\n')
