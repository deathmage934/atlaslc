#!/usr/bin/env python

# detect_bumps.py 2010zzz --sim_gaussian 58000,58325 19,20,20.5,21 30

from SNloop import SNloopclass
from statistics import median
import sigmacut
import sys
import numpy as np
import pandas as pd
from astropy.time import Time
import matplotlib.pyplot as plt
import pylab as matlib
from asym_gaussian import gauss2lc
from plot_lc import dataPlot

# for ignoring matplotlib deprecation warnings
import warnings
warnings.filterwarnings("ignore")

class detectbumpsclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)

		self.filt_dict = None

	def define_options(self, **kwargs):
		parser = SNloopclass.define_options(self, **kwargs)
		parser.add_argument('--sim_gaussian', nargs=3, default=None, help=('Comma-separated peakMJD list, peak_appmag, gaussian_sigma: add a Gaussian at peakMJD with a peak apparent magnitude of peak_appmag and a sigma of gaussian_sigma in days.'))
		return(parser)

	def applyrolling_gaussian(self,SNindex,controlindex=0,MJDbinsize=1.0,gaussian_sigma_days=30.0,
							  simparams=None):
		
		self.load_lc(SNindex,controlindex=controlindex,MJDbinsize=MJDbinsize)

		indices = self.lc.ix_unmasked('Mask',self.flag_day_bad)
		#badindices = self.lc.ix_masked('Mask',self.flag_day_bad)
		
		# make sure there are no lingering simulations!
		dropcols=[]
		for col in ['uJysim','SNRsim','simLC','SNRsimsum']:
			if col in self.t.columns:
				dropcols.append(col)
		if len(dropcols)>0:
			self.lc.t.drop(columns=dropcols,inplace=True)
		
		self.lc.t['SNR']=0.0
		self.lc.t.loc[indices,'SNR']=(self.lc.t.loc[indices,self.flux_colname])/self.lc.t.loc[indices,self.dflux_colname]

		title = self.t.loc[SNindex,'tnsname']
		if not(simparams is None):
			peakMJDs = simparams['sim_peakMJD'].split(',')
			print('bbbb',peakMJDs,simparams['sim_peakMJD'])
			
			title = 'SIM Gaussian: mag: %.2f sigma-+:(%.1f,%.1f) peak MJDs:' % (simparams['sim_appmag'],simparams['sim_sigma_minus'],simparams['sim_sigma_plus'])


			# get teh simulated gaussian
			mjds = self.lc.t.loc[indices,'MJD']
			self.lc.t.loc[indices,'uJysim']=self.lc.t.loc[indices,'uJy']
			self.lc.t['simLC'] = 0.0
			for peakMJD in peakMJDs:
				peakMJD = float(peakMJD)
				print('ADDING SIM Gaussian: peak MJD: %.0f Apparent mag: %.2f sigma-: %.1f sigma+:%.1f' % (peakMJD,simparams['sim_appmag'],simparams['sim_sigma_minus'],simparams['sim_sigma_plus']))
 
				simflux = gauss2lc(mjds,peakMJD,simparams['sim_sigma_minus'],simparams['sim_sigma_plus'],app_mag=simparams['sim_appmag'])
				
				# add simflux to light curve
				self.lc.t.loc[indices,'uJysim']+=simflux

				# save the sim LC, for all MJDs
				simflux_all = gauss2lc(self.lc.t['MJDbin'],peakMJD,simparams['sim_sigma_minus'],simparams['sim_sigma_plus'],app_mag=simparams['sim_appmag'])
				self.lc.t['simLC'] += simflux_all

				title += ' %.1f' % peakMJD

			# Make sure all bad rows have SNRsim=0.0, so that they have no impact on the rolling SNRsum
			self.lc.t['SNRsim']=0.0
			# include simflux in the SNR
			self.lc.t.loc[indices,'SNRsim']=(self.lc.t.loc[indices,'uJysim'])/self.lc.t.loc[indices,self.dflux_colname]

		
		gaussian_sigma = round(gaussian_sigma_days/MJDbinsize)
		windowsize = int(3 * gaussian_sigma * 2)
		halfwindowsize = int(windowsize*0.5)+1
		print('sigma(days)=%.2f, MJD binsize=%.2f, sigma(bins)=%d, window size(bins)=%d' % (gaussian_sigma_days,MJDbinsize,gaussian_sigma,windowsize))

		dataindices = np.array(range(len(self.lc.t))+np.full(len(self.lc.t),halfwindowsize))

		# Calculate the rolling SNR sum
		temp = pd.Series(np.zeros(len(self.lc.t)+2*halfwindowsize),name='SNR', dtype=np.float64)
		temp[dataindices] = self.lc.t['SNR']
		SNRsum = temp.rolling(windowsize,center=True,win_type='gaussian').sum(std=gaussian_sigma)
		self.lc.t['SNRsum']=list(SNRsum.loc[dataindices])

		if not(simparams is None):
			# Calculate the rolling SNR sum for SNR including simflux
			temp = pd.Series(np.zeros(len(self.lc.t)+2*halfwindowsize),name='SNRsim', dtype=np.float64)
			temp[dataindices] = self.lc.t['SNRsim']
			SNRsimsum = temp.rolling(windowsize,center=True,win_type='gaussian').sum(std=gaussian_sigma)
			self.lc.t['SNRsimsum']=list(SNRsimsum.loc[dataindices])
			
		if not(simparams is None):
			# maginfo is to add to filenames!
			maginfo='.sim%.1fmag' % (simparams['sim_appmag'])
			if args.savelc is True:
				self.save_lc(SNindex=SNindex,controlindex=controlindex,filt=self.filt,MJDbinsize=MJDbinsize,addsuffix=maginfo,overwrite=True)
			outbasefilename = self.lcbasename(SNindex=SNindex,controlindex=controlindex,MJDbinsize=MJDbinsize,addsuffix=maginfo)
		else:
			outbasefilename = self.lcbasename(SNindex=SNindex,controlindex=controlindex,MJDbinsize=MJDbinsize)

		"""
		plt.close("all")
		if not(simparams is None):
			ax1 = self.lc.plot( x='MJDbin',y='uJysim',yerr='duJy',kind='scatter',color='cyan',legend=True,title=title)
		else:
			ax1 = self.lc.plot(x='MJDbin',y='uJy',yerr='duJy',kind='scatter',color='red')
		if not(simparams is None):
			self.lc.plot(x='MJDbin',y='simLC',color='cyan',ax=ax1)
		if not(simparams is None):
			outfile='%s.simLC.png' % outbasefilename
		else:
			outfile='%s.LC.png' % outbasefilename
		print('Saving ',outfile)
		plt.savefig(outfile)

		if not(simparams is None):
			ax2 =  self.lc.plot(indices=indices,x='MJDbin',y='SNRsim',kind='scatter',color='cyan',title=title)
		else:
			ax2 = self.lc.plot(indices=indices,x='MJDbin',y='SNR',kind='scatter',color='red')
		if not(simparams is None):
			self.lc.plot(x='MJDbin',y='SNRsimsum',ax=ax2,color='cyan')
		self.lc.plot(x='MJDbin',y='SNRsum',ax=ax2,color='red')
		if not(simparams is None):
			outfile='%s.simSNR.png' % outbasefilename
		else:
			outfile='%s.SNR.png' % outbasefilename
		print('Saving ',outfile)
		plt.savefig(outfile)
		"""

		# Updated version by Sofia
		
		# plot og lc with simulated bumps (if added)
		plt.figure()
		plt.axhline(linewidth=1,color='k')
		plt.rcParams['font.family'] = 'serif'
		plt.rcParams["font.serif"] = 'times'
		if not(simparams is None):
			sp, plot_uJysim, dplot_uJysim = dataPlot(x=self.lc.t.loc[indices,'MJDbin'],y=self.lc.t.loc[indices,'uJysim'],dy=self.lc.t.loc[indices,'duJy'],fmt='co',ecolor='c')
			matlib.setp(plot_uJysim,ms=3,alpha=1)
		sp, plot_uJy, dplot_uJy = dataPlot(x=self.lc.t.loc[indices,'MJDbin'],y=self.lc.t.loc[indices,'uJy'],dy=self.lc.t.loc[indices,'duJy'],fmt='ro',ecolor='r')
		matlib.setp(plot_uJy,ms=3,alpha=1)
		if not(simparams is None):
			sp, plot, dplot = dataPlot(x=self.lc.t['MJDbin'],y=self.lc.t['simLC'],fmt='c')
		# set title and legend
		if controlindex == 0:
			if not(simparams is None):
				plt.title('SN %s and Simulated Eruptions (2 Gaussians of width = 30 days)' % self.t.loc[SNindex,'tnsname'])
			else:
				plt.title('SN %s' % self.t.loc[SNindex,'tnsname'])
			if not(simparams is None):
				plt.legend((plot_uJysim,plot_uJy,plot),('SN %s + Simulated Gaussian'% self.t.loc[SNindex,'tnsname'],'SN %s' % self.t.loc[SNindex,'tnsname'],'Gaussian Models'))
		else:
			if not(simparams is None):
				plt.title('Control LC %d and Simulated Eruptions (2 Gaussians of width = 30 days)' % controlindex)
			else:
				plt.title('Control LC %d' % controlindex)
			if not(simparams is None):
				plt.legend((plot_uJysim,plot_uJy,plot),('Control LC %d + Simulated Gaussian'% controlindex,'Control LC %d' % controlindex,'Gaussian Models'))
		plt.xlabel('MJD')
		plt.ylabel('Flux ($\mu$Jy)')
		# get x and y limits from args; else, leave as is
		xlim_lower, xlim_upper = plt.xlim()
		if not(args.xlim_lower is None): 
			xlim_lower = args.xlim_lower
		if not(args.xlim_upper is None):
			xlim_upper = args.xlim_upper
		ylim_lower, ylim_upper = plt.xlim()
		if not(args.ylim_lower is None): 
			ylim_lower = args.ylim_lower
		if not(args.ylim_upper is None): 
			ylim_upper = args.ylim_upper
		plt.xlim(xlim_lower,xlim_upper)
		plt.ylim(ylim_lower,ylim_upper)
		if not(simparams is None):
			outfile = '%s.simLC.png' % outbasefilename
		else:
			outfile = '%s.png' % outbasefilename
		print('Saving ',outfile)
		plt.savefig(outfile,dpi=200)
		plt.close()
		
		# plot snr
		plt.figure()
		plt.axhline(linewidth=1,color='k')
		if not(simparams is None):
			sp, plot_SNRsim, dplot_SNRsim = dataPlot(x=self.lc.t.loc[indices,'MJDbin'],y=self.lc.t.loc[indices,'SNRsim'],fmt='co')
			matlib.setp(plot_SNRsim,ms=3,alpha=1)
		sp, plot_SNR, dplot_SNR = dataPlot(x=self.lc.t.loc[indices,'MJDbin'],y=self.lc.t.loc[indices,'SNR'],fmt='ro')
		matlib.setp(plot_SNR,ms=3,alpha=1)
		if not(simparams is None):
			sp, plot_simsum, dplot_simsum = dataPlot(x=self.lc.t['MJDbin'],y=self.lc.t['SNRsimsum'],fmt='c')
		sp, plot_sum, dplot_sum = dataPlot(x=self.lc.t['MJDbin'],y=self.lc.t['SNRsum'],fmt='r')
		# set title and legend
		if controlindex == 0:
			plt.title('SN %s S/N and Gaussian Weighted Rolling Sum of S/N' % self.t.loc[SNindex,'tnsname'])
			if not(simparams is None):
				plt.legend((plot_SNRsim,plot_simsum,plot_SNR,plot_sum),('SN %s + Eruption Model S/N' % self.t.loc[SNindex,'tnsname'],'Gaussian Weighted Rolling Sum','SN %s S/N' % self.t.loc[SNindex,'tnsname'],'Gaussian Weighted Rolling Sum'))
			else:
				plt.legend((plot_SNR,plot_sum),('SN %s S/N' % self.t.loc[SNindex,'tnsname'],'Gaussian Weighted Rolling Sum'))
		else:
			plt.title('Control LC %d S/N and Gaussian Weighted Rolling Sum of S/N' % controlindex)
			if not(simparams is None):
				plt.legend((plot_SNRsim,plot_simsum,plot_SNR,plot_sum),('Control LC %d + Eruption Model S/N' % controlindex,'Gaussian Weighted Rolling Sum','Control LC %d S/N' % controlindex,'Gaussian Weighted Rolling Sum'))
			else:
				plt.legend((plot_SNR,plot_sum),('Control LC %d S/N' % controlindex,'Gaussian Weighted Rolling Sum'))
		plt.xlabel('MJD')
		plt.ylabel('S/N')
		# get x and y limits from args; else, leave as is
		xlim_lower, xlim_upper = plt.xlim()
		if not(args.xlim_lower is None): 
			xlim_lower = args.xlim_lower
		if not(args.xlim_upper is None):
			xlim_upper = args.xlim_upper
		ylim_lower, ylim_upper = plt.xlim()
		if not(args.ylim_lower is None): 
			ylim_lower = args.ylim_lower
		if not(args.ylim_upper is None): 
			ylim_upper = args.ylim_upper
		plt.xlim(xlim_lower,xlim_upper)
		plt.ylim(ylim_lower,ylim_upper)
		if not(simparams is None):
			outfile = '%s.simSNR.png' % outbasefilename
		else:
			outfile = '%s.snr.png' % outbasefilename
		print('Saving ',outfile)
		plt.savefig(outfile,dpi=200)
		plt.close()

		# make plot of all SNR
		plt.figure(self.filt_dict[self.filt])
		plt.axhline(linewidth=1,color='k')
		if controlindex == 0:
			sp, plotn_sum, dplotn_sum = dataPlot(x=self.lc.t['MJDbin'],y=self.lc.t['SNRsum'],fmt='r')
		else:
			sp, plot_sum, dplot_sum = dataPlot(x=self.lc.t['MJDbin'],y=self.lc.t['SNRsum'],fmt='c')
		#plt.close(self.filt_dict[self.filt])
		
		# TEMPORARY - DELETE -----------------------------
		"""
		plt.figure(0)
		plt.rcParams['font.family'] = 'serif'
		plt.rcParams["font.serif"] = 'times'
		
		sp, plot_uJy, dplot_uJy = dataPlot(x=self.lc.t.loc[indices,'MJDbin'],y=self.lc.t.loc[indices,'uJy'],dy=self.lc.t.loc[indices,'duJy'],fmt='ro',ecolor='r')
		matlib.setp(plot_uJy,ms=3,alpha=1)
		plt.title('%s' % self.t.loc[SNindex,'tnsname'])
		plt.xlabel('MJD')
		plt.ylabel('Flux ($\mu$Jy)')
		outfile = '%s.lc.png' % outbasefilename
		print('Saving ',outfile)
		plt.savefig(outfile,dpi=200)
		
		plt.figure(1)
		#sp, plot_SNR, dplot_SNR = dataPlot(x=self.lc.t.loc[indices,'MJDbin'],y=self.lc.t.loc[indices,'SNR'],fmt='ro')
		#matlib.setp(plot_SNR,ms=3,alpha=1)
		if controlindex == 0:
			sp, plotn_sum, dplotn_sum = dataPlot(x=self.lc.t['MJDbin'],y=self.lc.t['SNRsum'],fmt='r')
		else:
			sp, plot_sum, dplot_sum = dataPlot(x=self.lc.t['MJDbin'],y=self.lc.t['SNRsum'],fmt='c')
			#plt.legend((plotn_sum,plot_sum),('%s S/N Gaussian Weighted Rolling Sum' % self.t.loc[SNindex,'tnsname'],'Control LCs Gaussian Weighted Rolling Sum'))
		plt.title('%s Gaussian Weighted Rolling Sum of S/N' % self.t.loc[SNindex,'tnsname'])
		plt.xlabel('MJD')
		plt.ylabel('S/N')
		outfile = '%s.snr.png' % outbasefilename
		print('Saving ',outfile)
		plt.savefig(outfile,dpi=200)
		"""
		# TEMPORARY - DELETE -----------------------------

	def detectbumpsloop(self,SNindex,MJDbinsize=1.0,simparams=None):
		print('###################################\nDetecting Bumps ...\n###################################')
		
		# loop through SN and control lcs
		#for controlindex in range(0,len(self.RADECtable.t)):
		for controlindex in range(len(self.RADECtable.t)-1,-1,-1):
			# stop loop if only SN should be done
			if (not self.cfg.params['detectBumps']['apply2offsets']) and (controlindex>0):
				break
			# average the light curve by MJDbinsize
			self.applyrolling_gaussian(SNindex,controlindex=controlindex,MJDbinsize=MJDbinsize,gaussian_sigma_days=self.cfg.params['detectBumps']['gaussian_sigma'],simparams=simparams)

		# save allsnr plot
		outbasefilename = self.lcbasename(SNindex=SNindex,MJDbinsize=MJDbinsize)
		plt.figure(self.filt_dict[self.filt])
		plt.title('%s Gaussian Weighted Rolling Sum of S/N' % self.t.loc[SNindex,'tnsname'])
		plt.xlabel('MJD')
		plt.ylabel('S/N')
		# get x and y limits from args; else, leave as is
		xlim_lower, xlim_upper = plt.xlim()
		if not(args.xlim_lower is None): 
			xlim_lower = args.xlim_lower
		if not(args.xlim_upper is None):
			xlim_upper = args.xlim_upper
		ylim_lower, ylim_upper = plt.xlim()
		if not(args.ylim_lower is None): 
			ylim_lower = args.ylim_lower
		if not(args.ylim_upper is None): 
			ylim_upper = args.ylim_upper
		plt.xlim(xlim_lower,xlim_upper)
		plt.ylim(ylim_lower,ylim_upper)
		outfile = '%s.allsnr.png' % outbasefilename
		print('Saving ',outfile)
		plt.savefig(outfile,dpi=200)

		# TEMPORARY - DELETE -----------------------------
		"""
		outbasefilename = self.lcbasename(SNindex=SNindex,MJDbinsize=MJDbinsize)

		plt.figure(0)
		plt.title('SN %s with mags 19, 20, 20.5, 21' % self.t.loc[SNindex,'tnsname'])
		plt.xlabel('MJD')
		plt.ylabel('Flux ($\mu$Jy)')
		plt.xlim(57250,58800)
		plt.ylim(-190,200)
		outfile='%s.simLC000000.png' % outbasefilename
		print('Saving ',outfile)
		plt.savefig(outfile,dpi=200)

		plt.figure(1)
		plt.title('SN %s S/N and Gaussian Weighted Rolling Sum of S/N, with mags 19, 20, 20.5, 21' % self.t.loc[SNindex,'tnsname'])
		plt.xlabel('MJD')
		plt.ylabel('S/N')
		plt.xlim(57250,58800)
		plt.ylim(-100,200)
		outfile='%s.simSNR000000.png' % outbasefilename
		print('Saving ',outfile)
		plt.savefig(outfile,dpi=200)
		"""
		# TEMPORARY - DELETE -----------------------------

		print('### Detecting bumps in LCs done')

if __name__ == '__main__':

	detectbumps = detectbumpsclass()
	parser = detectbumps.define_options()
	args = parser.parse_args()

	SNindexlist = detectbumps.initialize(args)

	if args.filt is None:
		print('Looping through c and o filters...')
		detectbumps.filt_dict = {'o':1,'c':2}
		for filt in detectbumps.filt_dict:
			detectbumps.filt = filt
			if args.sim_gaussian is None:
				simparams=None
				for SNindex in SNindexlist:
					print('Detecting bumps for ',detectbumps.t.at[SNindex,'tnsname'],', index %i/%i' % (SNindex,len(detectbumps.t)))
					detectbumps.loadRADEClist(SNindex)
					detectbumps.detectbumpsloop(SNindex,MJDbinsize=args.MJDbinsize,simparams=simparams)
			else:
				if ',' in args.sim_gaussian[1]:
					appmags = args.sim_gaussian[1].split(',')
					print('Multiple mags input: ',appmags)
				else:
					appmags = [args.sim_gaussian[1]]
					print('Only 1 mag input: ',appmags)
				for appmag in appmags:
					print('Mag set to: ',appmag)
					simparams = {'sim_peakMJD':args.sim_gaussian[0],'sim_appmag':float(appmag),'sim_sigma_minus':float(args.sim_gaussian[2]),'sim_sigma_plus':float(args.sim_gaussian[2])}
					for SNindex in SNindexlist:
						print('Detecting bumps for ',detectbumps.t.at[SNindex,'tnsname'],', index %i/%i' % (SNindex,len(detectbumps.t)))
						detectbumps.loadRADEClist(SNindex)
						detectbumps.detectbumpsloop(SNindex,MJDbinsize=args.MJDbinsize,simparams=simparams)
			print('Finished with filter %s!' % filt)
	else:
		print('### FILTER SET: %s' % args.filt)
		detectbumps.filt = args.filt
		if args.sim_gaussian is None:
			simparams=None
			for SNindex in SNindexlist:
				print('Detecting bumps for ',detectbumps.t.at[SNindex,'tnsname'],', index %i/%i' % (SNindex,len(detectbumps.t)))
				detectbumps.loadRADEClist(SNindex)
				detectbumps.detectbumpsloop(SNindex,MJDbinsize=args.MJDbinsize,simparams=simparams)
		else:
			if ',' in args.sim_gaussian[1]:
				appmags = args.sim_gaussian[1].split(',')
				print('Multiple mags input: ',appmags)
			else:
				appmags = [args.sim_gaussian[1]]
				print('Only 1 mag input: ',appmags)
			for appmag in appmags:
				print('Mag set to: ',appmag)
				simparams = {'sim_peakMJD':args.sim_gaussian[0],'sim_appmag':float(appmag),'sim_sigma_minus':float(args.sim_gaussian[2]),'sim_sigma_plus':float(args.sim_gaussian[2])}
				for SNindex in SNindexlist:
					print('Detecting bumps for ',detectbumps.t.at[SNindex,'tnsname'],', index %i/%i' % (SNindex,len(detectbumps.t)))
					detectbumps.loadRADEClist(SNindex)
					detectbumps.detectbumpsloop(SNindex,MJDbinsize=args.MJDbinsize,simparams=simparams)
