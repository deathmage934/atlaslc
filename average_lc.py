#!/usr/bin/env python

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

	def calcaveragelc(self, SNindex, lc_uJy, lc_duJy, lc_MJD, MJDbinsize=None, offsetindex=None):
		#self.averagelctable.formattable(formatMapping={'OffsetID':'3d','MJD':'.5f','uJy':'.2f','duJy':'.2f','stdev':'.2f','X2norm':'.3f'})
		MJD = int(np.amin(lc_MJD))
		MJDmax = np.amax(lc_MJD)

		while MJD <= MJDmax:
			windowindices = np.where(np.logical_and(lc_MJD>=MJD, lc_MJD<(MJD+MJDbinsize))==True)
			windowindices = list(windowindices[0])
			print('windowindices: ',windowindices)
			uJy_windowdata = lc_uJy.iloc[windowindices]
			duJy_windowdata = lc_duJy.iloc[windowindices]
			MJD_windowdata = lc_MJD.iloc[windowindices]

			calcaverage=sigmacut.calcaverageclass()
			calcaverage.calcaverage_sigmacutloop(uJy_windowdata,noise=duJy_windowdata,verbose=2,Nsigma=3.0,median_firstiteration=True,saveused=True)

			if calcaverage.Nused<=0:
				print('No data values. Skipping...')
				MJD += MJDbinsize
				continue

			print("mean: %f, uncertainty: %f" % (calcaverage.mean,calcaverage.mean_err))
			lcaverageindex = len(self.averagelctable.t)

			df = pd.DataFrame([[self.RADECtable.t.loc[offsetindex,'OffsetID'],calcaverage.mean,calcaverage.mean_err,calcaverage.stdev,calcaverage.X2norm,calcaverage.Nused,calcaverage.Nskipped]],columns=['OffsetID','uJy','duJy','stdev','X2norm','Nused','Nclipped'])
			self.averagelctable.t = self.averagelctable.t.append(df,ignore_index=True)
			# delete me
			#self.averagelctable.t.add_row({'OffsetID':self.RADECtable.t['OffsetID'][offsetindex],'uJy':calcaverage.mean,'duJy':calcaverage.mean_err, 'stdev':calcaverage.stdev, 'X2norm':calcaverage.X2norm, 'Nused':calcaverage.Nused, 'Nclipped':calcaverage.Nskipped})
			#print(self.averagelctable.write()) # delete me


			print('lcaverageindex and calcaverage: ',lcaverageindex,calcaverage.use)
			mask=np.logical_not(calcaverage.use)
			calcaverage.calcaverage_sigmacutloop(MJD_windowdata,noise=duJy_windowdata,mask=mask,verbose=2,Nsigma=0.0,median_firstiteration=True,saveused=True)
			self.averagelctable.t.at[lcaverageindex,'MJD'] = calcaverage.mean # FIX TO INT!!!
			self.averagelctable.t.at[lcaverageindex,'MJDNused'] = calcaverage.Nused
			self.averagelctable.t.at[lcaverageindex,'MJDNskipped'] = calcaverage.Nskipped

			MJD += MJDbinsize

		print(self.averagelctable.write())
		return(0)
	
	def averagelcs(self, SNindex, offsetindex, lc_uJy, lc_duJy, lc_MJD, MJDbinsize=None,fileformat=None,overwrite=True):
		if offsetindex>0:
			#self.averagelctable.t.remove_rows(slice(0,len(self.averagelctable.t)))
			self.averagelctable.t.drop(self.averagelctable.ix_inrange(lowlim=0,uplim=len(self.averagelctable.t)))
		result = self.calcaveragelc(SNindex, lc_uJy, lc_duJy, lc_MJD, offsetindex=offsetindex,MJDbinsize=MJDbinsize)

		filename = self.lcbasename(SNindex,filt=self.filt,offsetindex=offsetindex,MJDbinsize=MJDbinsize)+'.txt'
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
		print(self.lc.t)
		print('Length of self.lc.t: ',len(self.lc.t))
		if len(self.lc.t) == 0:
			return(1)

		makecuts_apply = self.cfg.params['averagelc']['makecuts']
		if args.skip_makecuts:
			makecuts_apply = False

		# data set to lc table cuts_indices
		if makecuts_apply == True:
			if not('Mask' in self.lc.t.columns):
				raise RuntimeError('No "Mask" column exists! Please run "cleanup_lc.py %s" beforehand.' % self.t[SNindex]['tnsname'])
			lc_uJy, lc_duJy, lc_MJD, totalcut = self.makecuts_indices(SNindex, offsetindex)
		# data set to lc table
		else:
			print('Skipping makecuts using mask column...')
			lc_uJy = self.lc.t['uJy']
			lc_duJy = self.lc.t['duJy']
			lc_MJD = self.lc.t['MJD'] 

		self.averagelcs(SNindex, offsetindex, lc_uJy, lc_duJy, lc_MJD, MJDbinsize=MJDbinsize)

if __name__ == '__main__':

	averagelc = averagelcclass()
	parser = averagelc.define_options()
	args = parser.parse_args()

	SNindexlist = averagelc.initialize(args)

	MJDbinsize = averagelc.cfg.params['output']['MJDbinsize']
	if not(args.MJDbinsize is None):
		MJDbinsize = args.MJDbinsize

	for SNindex in SNindexlist:
		averagelc.loadRADEClist(SNindex=SNindex, filt=averagelc.filt)
		for offsetindex in range(len(averagelc.RADECtable.t)):
			averagelc.averagelcloop(args,SNindex,offsetindex)

	print('\n')
