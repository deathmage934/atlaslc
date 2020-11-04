#!/usr/bin/env python

# S. Rest

from SNloop import SNloopclass
import numpy as np
import sigmacut
import sys
import pandas as pd

class averagelcclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)

		self.flag_uncertainty = 0x1
		self.flag_dynamic = 0x2
		self.flag_static = 0x4

	def calcaveragelc(self, args, SNindex, lc_uJy, lc_duJy, lc_MJD, MJDbinsize=None, offsetindex=None, cuts_indices=None):
		if not(cuts_indices is None):
			mjdindices = self.lc.ix_remove_null(indices=cuts_indices)
			#print(mjdindices)
		else:
			mjdindices = list(self.lc.ix_remove_null())

		MJD = np.amin(self.lc.t.loc[mjdindices,'MJD'])
		MJDmax = np.amax(self.lc.t.loc[mjdindices,'MJD'])

		while MJD <= MJDmax:
			windowindices = np.where(np.logical_and(lc_MJD>=MJD, lc_MJD<(MJD+MJDbinsize))==True)
			windowindices = list(windowindices[0])
			uJy_windowdata = lc_uJy.iloc[windowindices]
			duJy_windowdata = lc_duJy.iloc[windowindices]
			MJD_windowdata = lc_MJD.iloc[windowindices]

			calcaverage=sigmacut.calcaverageclass()
			calcaverage.calcaverage_sigmacutloop(uJy_windowdata,noise=duJy_windowdata,verbose=2,Nsigma=3.0,median_firstiteration=True,saveused=True)
			if calcaverage.Nused<=0:
				if self.verbose>1:
					print('No data values. Skipping...')
				MJD += MJDbinsize
				continue
			lcaverageindex = len(self.averagelctable.t)
			df = pd.DataFrame([[self.RADECtable.t.loc[offsetindex,'OffsetID'],0,calcaverage.mean,calcaverage.mean_err,calcaverage.stdev,calcaverage.X2norm,calcaverage.Nused,calcaverage.Nskipped,0,0]],columns=['OffsetID','MJD',self.flux_colname,self.dflux_colname,'stdev','X2norm','Nused','Nclipped','MJDNused','MJDNskipped'])
			self.averagelctable.t = self.averagelctable.t.append(df,ignore_index=True)
			
			mask=np.logical_not(calcaverage.use)
			calcaverage.calcaverage_sigmacutloop(MJD_windowdata,noise=duJy_windowdata,mask=mask,verbose=2,Nsigma=0.0,median_firstiteration=True,saveused=True)
			self.averagelctable.t.at[lcaverageindex,'MJD'] = calcaverage.mean # FIX TO INT
			self.averagelctable.t.at[lcaverageindex,'MJDNused'] = calcaverage.Nused
			self.averagelctable.t.at[lcaverageindex,'MJDNskipped'] = calcaverage.Nskipped

			if self.verbose>1:
				print('windowindices: ',windowindices)
				print("mean: %f, uncertainty: %f" % (calcaverage.mean,calcaverage.mean_err))
				print('lcaverageindex and calcaverage: ',lcaverageindex,calcaverage.use)

			MJD += MJDbinsize

		#print(self.averagelctable.write())
		return(0)
	
	def averagelcs(self, args, SNindex, offsetindex, lc_uJy, lc_duJy, lc_MJD, MJDbinsize=None,fileformat=None,overwrite=True,cuts_indices=None):
		if offsetindex>0:
			#self.averagelctable.t.remove_rows(slice(0,len(self.averagelctable.t)))
			self.averagelctable.t.drop(self.averagelctable.ix_inrange(lowlim=0,uplim=len(self.averagelctable.t)))
		result = self.calcaveragelc(args, SNindex, lc_uJy, lc_duJy, lc_MJD, offsetindex=offsetindex,MJDbinsize=MJDbinsize,cuts_indices=cuts_indices)

		filename = self.lcbasename(SNindex=SNindex,filt=self.filt,offsetindex=offsetindex,MJDbinsize=MJDbinsize)+'.txt'
		if fileformat is None: 
			fileformat = self.cfg.params['output']['fileformat']

		if result == 0: 
			self.averagelctable.write(filename,overwrite=True,verbose=True)
		else:
			print('Removing file because length of light curve is 0...')
			rmfile(filename)

	def averagelcloop(self,args,SNindex,offsetindex):
		MJDbinsize = self.cfg.params['output']['MJDbinsize']
		if not(args.MJDbinsize is None): MJDbinsize = args.MJDbinsize

		# load lc
		self.load_lc(SNindex, filt=self.filt, offsetindex=offsetindex, MJDbinsize=None)
		#print(self.lc.t)
		print('Length of self.lc.t: ',len(self.lc.t))
		if len(self.lc.t) == 0:
			return(1)

		makecuts_apply = self.cfg.params['averagelc']['makecuts']
		if not(args.avg_makecuts) is None:
			if args.avg_makecuts is True:
				makecuts_apply = True
			else:
				makecuts_apply = False

		# data set to lc table cuts_indices
		if makecuts_apply == True:
			if not('Mask' in self.lc.t.columns):
				raise RuntimeError('No "Mask" column exists! Please run "cleanup_lc.py %s" beforehand.' % self.t.at[SNindex,'tnsname'])
			lc_uJy, lc_duJy, lc_MJD, cuts_indices, bad_data = self.makecuts_indices(SNindex, offsetindex, 'averagelc')
			print('Calculating average_lc table for offset index %d' % offsetindex)
			self.averagelcs(args, SNindex, offsetindex, lc_uJy, lc_duJy, lc_MJD, MJDbinsize=MJDbinsize, cuts_indices=cuts_indices)
		# data set to lc table
		else:
			print('Skipping makecuts using mask column...')
			lc_uJy = self.lc.t[self.flux_colname]
			lc_duJy = self.lc.t[self.dflux_colname]
			lc_MJD = self.lc.t['MJD']
			print('Calculating average_lc table for offset index %d' % offsetindex) 
			self.averagelcs(args, SNindex, offsetindex, lc_uJy, lc_duJy, lc_MJD, MJDbinsize=MJDbinsize)
		
if __name__ == '__main__':

	averagelc = averagelcclass()
	parser = averagelc.define_options()
	args = parser.parse_args()

	SNindexlist = averagelc.initialize(args)

	MJDbinsize = averagelc.cfg.params['output']['MJDbinsize']
	if not(args.MJDbinsize is None):
		MJDbinsize = args.MJDbinsize

	for SNindex in SNindexlist:
		averagelc.loadRADEClist(SNindex, filt=averagelc.filt)
		for offsetindex in range(len(averagelc.RADECtable.t)):
			averagelc.averagelcloop(args,SNindex,offsetindex)

	print('\n')
