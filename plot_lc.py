#!/usr/bin/env python
'''
@author: S. Rest
'''

from SNloop import SNloopclass
import pylab as matlib
import matplotlib.pyplot as plt
import sys,os

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

class plotlcclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)

	def plot_lc(self,args,SNindex,sp=None):
		print('Plotting SN lc and control lc...')
		plt.figure()
		if sp is None: sp = matlib.subplot(111)
		self.loadRADEClist(SNindex)

		# plot SN in red first, then loop through control lcs to plot in blue
		for controlindex in range(len(self.RADECtable.t)-1,-1,-1):
			plt.rcParams['font.family'] = 'serif'
			plt.rcParams["font.serif"] = 'times'

			self.load_lc(SNindex, filt=self.filt, controlindex=self.RADECtable.t.at[controlindex,'ControlID'])
			if self.verbose>2: print('control index: ',controlindex)

			# check if plotting good data
			makecuts_apply = self.cfg.params['plotlc']['makecuts']
			# check if plotting bad SN lc data ONLY IF making cuts
			plot_bad_data = self.cfg.params['plotlc']['plot_bad_data']
			if makecuts_apply is False:
				print("WARNING: Cannot plot bad data without making cuts! Setting plot_bad_data to False.")
				plot_bad_data = False
			# check if plotting cotrol lc data (makecuts also applies)
			plot_controllc_data = self.cfg.params['plotlc']['plot_controllc_data']
			
			if makecuts_apply is True:
				if not('Mask' in self.lc.t.columns): raise RuntimeError('No "Mask" column exists! Please run "cleanup_lc.py %s" beforehand.' % self.t.at[SNindex,'tnsname'])
				# make cuts on lc, save bad data to possibly plot
				lc_uJy, lc_duJy, lc_MJD, cuts_indices, bad_data = self.makecuts_indices(SNindex, controlindex=controlindex, procedure1='plotlc')
				lc_uJy_bad = self.lc.t.loc[bad_data, self.flux_colname]
				lc_duJy_bad = self.lc.t.loc[bad_data, self.dflux_colname]
				lc_MJD_bad = self.lc.t.loc[bad_data, 'MJD']
			else:
				print('Skipping makecuts using mask column...')
				lc_uJy = self.lc.t[self.flux_colname]
				lc_duJy = self.lc.t[self.dflux_colname]
				lc_MJD = self.lc.t['MJD']

			if controlindex==0:
				
				if plot_bad_data is True:
					# plot bad data with open red circles
					sp, plotbad, dplotbad = dataPlot(lc_MJD_bad,lc_uJy_bad,dy=lc_duJy_bad,sp=sp)
					matlib.setp(plotbad,mfc='white',ms=4,color='r')
				# plot good data in closed red circles
				sp, plotSN, dplotSN = dataPlot(lc_MJD,lc_uJy,dy=lc_duJy,sp=sp)
				matlib.setp(plotSN,ms=4,color='r')
				maxlc = max(lc_uJy)
				minlc = min(lc_uJy)
				
			else: 
				if plot_controllc_data is True:
					if plot_bad_data is True:
						# plot bad data with open blue circles
						sp, plotControlLCBad, dplotControlLCBad = dataPlot(lc_MJD_bad,lc_uJy_bad,dy=lc_duJy_bad,sp=sp)
						matlib.setp(plotControlLCBad,mfc='white',ms=4,color='b')
					# plot good data in closed blue circles
					sp, plotControlLC, dplotControlLC = dataPlot(lc_MJD,lc_uJy,dy=lc_duJy,sp=sp)
					matlib.setp(plotControlLC,ms=4,color='b')
		
		# determine legend and check if control lc data plotted
		# if control lcs
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
			# check if bad data plotted
			if plot_bad_data is True:
				plt.legend((plotSN,plotbad,plotControlLC,plotControlLCBad),('SN %s Accurate & Usable Data' % self.t.at[SNindex,'tnsname'],'SN %s Inaccurate Data' % self.t.at[SNindex,'tnsname'],controlLClabel+' Accurate & Usable Data',controlLClabel+' Inaccurate Data'))
			else:
				plt.legend((plotSN,plotControlLC),('SN %s' % self.t.at[SNindex,'tnsname'],controlLClabel))
		# if only sn
		else:
			plt.legend(('SN %s cleaned' % self.t.at[SNindex,'tnsname']))
		
		plt.title('SN %s' % self.t.at[SNindex,'tnsname'])
		#plt.title('8 Control Light Curves around SN 2021pb, 17" Radius')
		plt.axhline(linewidth=1,color='k')
		plt.xlabel('MJD')
		plt.ylabel(self.flux_colname)

		# get x limits from args; else, leave as is
		xlim_lower, xlim_upper = plt.xlim()
		if not(args.xlim_lower is None): 
			xlim_lower = args.xlim_lower
		if not(args.xlim_upper is None):
			xlim_upper = args.xlim_upper
		
		# get y limits from args; else, set to min and max lc
		if not(args.ylim_lower is None): 
			ylim_lower = args.ylim_lower
		else: 
			ylim_lower = minlc*1.1
		if not(args.ylim_upper is None): 
			ylim_upper = args.ylim_upper
		else:
			ylim_upper = maxlc*1.1
		
		# set x and y limits
		plt.xlim(xlim_lower,xlim_upper)
		plt.ylim(ylim_lower,ylim_upper)
		print('xlim lower: ',xlim_lower,'. xlim upper: ',xlim_upper,'. ')
		print('ylim lower: ',ylim_lower,'. ylim upper: ',ylim_upper,'. ')

		# save plot
		plotfilename = self.lcbasename(SNindex=SNindex)+'.png'
		print('Plot file name: ',plotfilename)
		plt.savefig(plotfilename,dpi=200)

	def plot_lc_controlLC(self,args,SNindex,sp=None,c1_flag=False,c2_flag=False):
		plt.rcParams['font.family'] = 'serif'
		plt.rcParams["font.serif"] = 'times'

		if sp is None:
			sp = matlib.subplot(111)
		self.load_lc(SNindex, filt=self.filt, controlindex=0)

		makecuts_apply = self.cfg.params['plotlc']['makecuts']
		if not(args.avg_makecuts) is None:
			if args.avg_makecuts is True:
				makecuts_apply = True
			else:
				makecuts_apply = False
		if makecuts_apply is True:
			if not('Mask' in self.lc.t.columns):
				raise RuntimeError('No "Mask" column exists! Please run "cleanup_lc.py %s" beforehand.' % self.t.at[SNindex,'tnsname'])
			lc_uJy, lc_duJy, lc_MJD, cuts_indices, bad_data = self.makecuts_indices(SNindex, controlindex=0, procedure1='plotlc')
			lc_uJy_bad = self.lc.t.loc[bad_data, self.flux_colname]
			lc_duJy_bad = self.lc.t.loc[bad_data, self.dflux_colname]
			lc_MJD_bad = self.lc.t.loc[bad_data, 'MJD']
		else:
			print('Skipping makecuts using mask column...')
			lc_uJy = self.lc.t[self.flux_colname]
			lc_duJy = self.lc.t[self.dflux_colname]
			lc_MJD = self.lc.t['MJD']

		# plot SN lc and o1 controlLC stats
		if c1_flag is True:
			plt.figure(3)
			c1_stddev1 = self.lc.t['c1_mean']+self.lc.t['c1_stddev']
			c1_stddev2 = self.lc.t['c1_mean']-self.lc.t['c1_stddev']

			if makecuts_apply is True:
				sp, plot2, dplot2 = dataPlot(lc_MJD_bad, lc_uJy_bad, dy=lc_duJy_bad)
				matlib.setp(plot2,mfc='white',ms=4,color='r')
				sp, plot, dplot = dataPlot(lc_MJD, lc_uJy, dy=lc_duJy)
				matlib.setp(plot,ms=4,color='r')
			else:
				sp, plot, dplot = dataPlot(lc_MJD, lc_uJy, dy=lc_duJy)
				matlib.setp(plot,ms=4,color='r')
			plt.fill_between(self.lc.t['MJD'],c1_stddev1,c1_stddev2)
			# plot data points in addition to fill_between
			'''
			sp, c1_stddev1, dplot_c1_stddev1 = dataPlot(lc_MJD,c1_stddev1)
			matlib.setp(c1_stddev1,ms=4,color='k')
			sp, c1_stddev2, dplot_c1_stddev2 = dataPlot(lc_MJD,c1_stddev2)
			matlib.setp(c1_stddev2,ms=4,color='k')
			'''

			plt.title('%s with cleaned control lc data' % self.t.at[SNindex,'tnsname'])
			plt.xlabel('MJD')
			plt.ylabel('uJy')
			if makecuts_apply is True:
				plot_legend = self.t.at[SNindex,'tnsname']+' cleaned'
				plot2_legend = self.t.at[SNindex,'tnsname']+' cut data'
				plt.legend((plot,plot2),(plot_legend,plot2_legend))
			else:
				plt.legend((plot),(self.t.at[SNindex,'tnsname']))
			plt.axhline(linewidth=1,color='k')
			plt.xlim(58900,59250) # delete me
			plt.ylim(min(lc_uJy)*1.1,max(lc_uJy)*1.1) # delete me

			plotfilename = self.lcbasename(SNindex=SNindex)+'.mask4mjd.png'
			print('Plot file name: ',plotfilename)
			plt.savefig(plotfilename)

		# plot SN lc and o2 controlLC stats
		if c2_flag is True:
			plt.figure(4)
			c2_stddev1 = self.lc.t['c2_mean']+self.lc.t['c2_stddev']
			c2_stddev2 = self.lc.t['c2_mean']-self.lc.t['c2_stddev']

			if makecuts_apply is True:
				sp, plot2, dplot2 = dataPlot(lc_MJD_bad, lc_uJy_bad, dy=lc_duJy_bad)
				matlib.setp(plot2,mfc='white',ms=4,color='r')
				sp, plot, dplot = dataPlot(lc_MJD, lc_uJy, dy=lc_duJy)
				matlib.setp(plot,ms=4,color='r')
			else:
				sp, plot, dplot = dataPlot(lc_MJD, lc_uJy, dy=lc_duJy)
				matlib.setp(plot,ms=4,color='r')
			plt.fill_between(self.lc.t['MJD'],c2_stddev1,c2_stddev2)
			# plot data points in addition to fill_between
			'''
			sp, c2_stddev1, dplot_c2_stddev1 = dataPlot(lc_MJD,c2_stddev1)
			matlib.setp(c2_stddev1,ms=4,color='k')
			sp, c2_stddev2, dplot_c2_stddev2 = dataPlot(lc_MJD,c2_stddev2)
			matlib.setp(c2_stddev2,ms=4,color='k')
			'''

			plt.title('%s with control lc nans cut' % self.t.at[SNindex,'tnsname'])
			plt.xlabel('MJD')
			plt.ylabel('uJy')
			if makecuts_apply is True:
				plot_legend = self.t.at[SNindex,'tnsname']+' cleaned'
				plot2_legend = self.t.at[SNindex,'tnsname']+' cut data'
				#plt.legend((plot,plot2,c2_stddev1,c2_stddev2),(plot_legend,plot2_legend,'c2_stddev1','c2_stddev2'))
				plt.legend((plot,plot2),(plot_legend,plot2_legend))
			else:
				#plt.legend((plot,c2_stddev1,c2_stddev2),(self.t.at[SNindex,'tnsname'],'c2_stddev1','c2_stddev2'))
				plt.legend((plot),(self.t.at[SNindex,'tnsname']))
			plt.axhline(linewidth=1,color='k')
			#plt.xlim(58900,59250) # delete me
			plt.ylim(min(lc_uJy)*1.1,max(lc_uJy)*1.1) # delete me

			plotfilename = self.lcbasename(SNindex=SNindex)+'.mask_nan.png'
			print('Plot file name: ',plotfilename)
			plt.savefig(plotfilename)

	def plotlcloop(self,args,SNindex):
		print('###################################\nPlotting LCs...\n###################################')
		self.plot_lc(args,SNindex)
		
		# decide if plotting controlLC stats
		c1_flag = False
		c2_flag = False
		if self.cfg.params['plotlc']['plot_mask4mjd'] is True:
			c1_flag = True
		if self.cfg.params['plotlc']['plot_mask_nan'] is True:
			c2_flag = True
		if (c1_flag is True) and (c2_flag is True):
			print('Plotting mask4mjd and mask_nan control LC...')
			self.plot_lc_controlLC(args,SNindex,c1_flag=True,c2_flag=True)
		elif c1_flag is True:
			print('Plotting mask4mjd control LC...')
			self.plot_lc_controlLC(args,SNindex,c1_flag=True)
		elif c2_flag is True:
			print('Plotting mask_nan control LC...')
			self.plot_lc_controlLC(args,SNindex,c2_flag=True)

if __name__ == '__main__':

	plotlc = plotlcclass()
	parser = plotlc.define_options()
	args = parser.parse_args()

	SNindexlist = plotlc.initialize(args)

	if args.filt is None:
		print('Looping through c and o filters...')
		for filt in ['o','c']:
			print('### FILTER SET: %s' % filt)
			plotlc.filt = filt
			for SNindex in SNindexlist:
				print('Plotting lc for ',plotlc.t.at[SNindex,'tnsname'],', index %i/%i' % (SNindex,len(plotlc.t)))
				print(SNindex,plotlc.t.at[SNindex,'tnsname']) # delete me
				plotlc.plotlcloop(args,SNindex)
			print('Finished with filter %s!' % filt)
	else:
		print('### FILTER SET: %s' % args.filt)
		plotlc.filt = args.filt
		for SNindex in SNindexlist:
			print('Plotting lc for ',plotlc.t.at[SNindex,'tnsname'],', index %i/%i' % (SNindex,len(plotlc.t)))
			print(SNindex,plotlc.t.at[SNindex,'tnsname']) # delete me
			plotlc.plotlcloop(args,SNindex)
