#!/usr/bin/env python

from SNloop import SNloopclass
import numpy as np
import sigmacut
import sys

class averagelcclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)

		self.flag_uncertainty = 0x1
		self.flag_dynamic = 0x2
		self.flag_static = 0x4

	def makecuts_indeces(self,SNindex,offsetindex):
		flags = self.cfg.params['averagelc']['flags']
		print('Setting indeces using flags: %x' % flags)
		
		mask2 = np.bitwise_and(self.lc.t['Mask'], flags, out=self.lc.t['Mask'])
		print(self.lc.t)

		print('Points to cut: ',np.where(self.lc.t['Mask']>0))
		cuts_indeces = np.where(self.lc.t['Mask']==0)

		lc_uJy = self.lc.t['uJy'][cuts_indeces]
		lc_duJy = self.lc.t['duJy'][cuts_indeces]
		lc_MJD = self.lc.t['MJD'][cuts_indeces]

		return(lc_uJy, lc_duJy, lc_MJD)

	def calcaveragelc(self, SNindex, lc_uJy, lc_duJy, lc_MJD, MJDbinsize=None, offsetindex=None):
		self.averagelctable.formattable(formatMapping={'OffsetID':'3d','MJD':'.5f','uJy':'.2f','duJy':'.2f','stdev':'.2f','X2norm':'.3f'})

		MJD = int(np.amin(lc_MJD))
		MJDmax = np.amax(lc_MJD)
		#print('MJD: ',MJD)

		while MJD <= MJDmax:
			windowindeces = np.logical_and(lc_MJD>=MJD, lc_MJD<(MJD+MJDbinsize))
			uJy_windowdata = lc_uJy[windowindeces]
			#print(uJy_windowdata)
			duJy_windowdata = lc_duJy[windowindeces]
			#print(duJy_windowdata)
			MJD_windowdata = lc_MJD[windowindeces]
			#print(MJD_windowdata)

			calcaverage=sigmacut.calcaverageclass()
			calcaverage.calcaverage_sigmacutloop(uJy_windowdata,noise=duJy_windowdata,verbose=2,Nsigma=3.0,median_firstiteration=True,saveused=True)

			if calcaverage.Nused<=0:
				print('No data values. Skipping...')
				MJD += MJDbinsize
				continue

			print("mean:%f (uncertainty:%f) " % (calcaverage.mean,calcaverage.mean_err))
			lcaverageindex = len(self.averagelctable.t)
			self.averagelctable.t.add_row({'OffsetID':self.RADECtable.t['OffsetID'][offsetindex],'uJy':calcaverage.mean,'duJy':calcaverage.mean_err, 'stdev':calcaverage.stdev, 'X2norm':calcaverage.X2norm, 'Nused':calcaverage.Nused, 'Nclipped':calcaverage.Nskipped})
			
			print('lcaverageindex and calcaverage: ',lcaverageindex,calcaverage.use)
			mask=np.logical_not(calcaverage.use)
			calcaverage.calcaverage_sigmacutloop(MJD_windowdata,noise=duJy_windowdata,mask=mask,verbose=2,Nsigma=0.0,median_firstiteration=True,saveused=True)
			self.averagelctable.t['MJD'][lcaverageindex] = calcaverage.mean
			self.averagelctable.t['MJDNused'][lcaverageindex] = calcaverage.Nused
			self.averagelctable.t['MJDNskipped'][lcaverageindex] = calcaverage.Nskipped

			MJD += MJDbinsize

		print(self.averagelctable.t)
		return(0)
	
	def averagelcs(self, SNindex, offsetindex, lc_uJy, lc_duJy, lc_MJD, MJDbinsize=None,fileformat=None,overwrite=True):
		if offsetindex>0:
			self.averagelctable.t.remove_rows(slice(0,len(self.averagelctable.t)))
		result = self.calcaveragelc(SNindex, lc_uJy, lc_duJy, lc_MJD, offsetindex=offsetindex,MJDbinsize=MJDbinsize)

		filename = self.lcbasename(SNindex,filt=self.filt,offsetindex=offsetindex,MJDbinsize=MJDbinsize)+'.txt'
		if fileformat is None: 
			fileformat = self.cfg.params['output']['fileformat']

		if result == 0: 
			self.averagelctable.write(filename,format=fileformat,overwrite=True,verbose=(self.verbose>0))
		else:
			print('Removing file because length of self.lc.t is 0...')
			rmfile(filename)

	def averagelcloop(self,args,SNindex,offsetindex):
		MJDbinsize = self.cfg.params['output']['MJDbinsize']
		if not(args.MJDbinsize is None): MJDbinsize = args.MJDbinsize
		#self.loadRADEClist(SNindex)

		# load lc
		self.load_lc(SNindex, filt=self.filt, offsetindex=offsetindex, MJDbinsize=None)
		print(self.lc.t)
		print('Length of self.lc.t: ',len(self.lc.t))
		if len(self.lc.t) == 0:
			return(1)

		makecuts_apply = self.cfg.params['averagelc']['makecuts']
		if args.skip_makecuts:
			makecuts_apply = False

		# data set to lc table cuts_indeces
		if makecuts_apply == True:
			lc_uJy, lc_duJy, lc_MJD = self.makecuts_indeces(SNindex, offsetindex)
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
