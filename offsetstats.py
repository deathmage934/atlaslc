#!/usr/bin/env python
'''
@author: S. Rest
'''

from SNloop import SNloopclass
from statistics import median
from copy import deepcopy
import sigmacut
import sys
import numpy as np
import pandas as pd

class offsetstatsclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)

	def cutandaveragelc(self,SNindex,offsetindex=0,MJDbinsize=1):
		daymax_Nclip = self.cfg.params['offsetstats']['daymax_Nclip']
		daymin_Ngood = self.cfg.params['offsetstats']['daymin_Ngood']
		daymax_X2norm = self.cfg.params['offsetstats']['daymax_X2norm']
		print('Max Nskipped: %d, min Nused: %d, max X2norm: %f' % (daymax_Nclip,daymin_Ngood,daymax_X2norm))

		MJD = np.amin(self.lc.t['MJD'])
		MJDmax = np.amax(self.lc.t['MJD'])

		# clear averagelc table and set dtypes
		self.averagelc.t = self.averagelc.t.iloc[0:0]
		self.averagelc.t['OffsetID'] = pd.Series([], dtype=np.int32)
		self.averagelc.t['MJD'] = pd.Series([], dtype=np.float64)
		self.averagelc.t['uJy'] = pd.Series([], dtype=np.float64)
		self.averagelc.t['duJy'] = pd.Series([], dtype=np.float64)
		self.averagelc.t['stdev'] = pd.Series([], dtype=np.float64)
		self.averagelc.t['X2norm'] = pd.Series([], dtype=np.float64)
		self.averagelc.t['Nclip'] = pd.Series([], dtype=np.int64)
		self.averagelc.t['Ngood'] = pd.Series([], dtype=np.int64)
		self.averagelc.t['Nexcluded'] = pd.Series([], dtype=np.int64)
		#self.averagelc.t['Mask'] = pd.Series([], dtype=np.int32)

		while MJD <= MJDmax:
			# get 4 measurements
			if self.verbose>1: 
				print('MJD range: ',MJD,' to ',MJD+1)
			ix1 = self.lc.ix_inrange(colnames=['MJD'],lowlim=MJD,uplim=MJD+MJDbinsize,exclude_uplim=True)
			if len(ix1)==0:
				if self.verbose>1: 
					print('Length of MJD range = 0, skipping MJD range...')
				MJD += MJDbinsize
				continue
			else:
				if self.verbose>2: 
					self.lc.write(indices=ix1) # delete me
			
			# get good and ok measurements (from offsetstats) out of 4
			ix2 = self.lc.ix_unmasked('Mask',maskval=self.flag_o2_bad|self.flag_o0_X2norm,indices=ix1)
			if len(ix2)==0:
				if self.verbose>1: 
					print('Length of good and ok measurements = 0, skipping MJD range...')
				MJD += MJDbinsize
				continue
			else:
				if self.verbose>2: 
					print('Good and ok indices (ix2): ',ix2)

			# sigmacut good and ok measurements
			self.lc.calcaverage_sigmacutloop('uJy',noisecol='duJy',indices=ix2,verbose=1,Nsigma=3.0,median_firstiteration=True)
			fluxstatparams = deepcopy(self.lc.statparams)
			if self.verbose>1: 
				print('Nclip: {}, Ngood: {}, X2norm: {}'.format(fluxstatparams['Nclip'],fluxstatparams['Ngood'],fluxstatparams['X2norm']))
			if fluxstatparams['mean'] is None or len(fluxstatparams['ix_good'])<1:
				if self.verbose>1: 
					print('Mean uJy is None OR length index good < 1, flagging bad day and skipping MJD range...')
				self.lc.t.loc[ix1,'Mask'] = np.bitwise_or(self.lc.t.loc[ix1,'Mask'],self.flag_daybad)
				MJD += MJDbinsize
				continue

			# get average mjd
			self.lc.calcaverage_sigmacutloop('MJD',noisecol='duJy',indices=fluxstatparams['ix_good'],verbose=1,Nsigma=0,median_firstiteration=False)
			averagemjd = self.lc.statparams['mean']
			if self.verbose>1:
				print('Average MJD: ',averagemjd,' (check MJD range: ',MJD,' to ',MJD+1,')')
			
			# add row to averagelc table
			df = {'OffsetID':offsetindex,'MJD':averagemjd,self.flux_colname:fluxstatparams['mean'],self.dflux_colname:fluxstatparams['mean_err'],'stdev':fluxstatparams['stdev'],'X2norm':fluxstatparams['X2norm'],'Nclip':fluxstatparams['Nclip'],'Ngood':fluxstatparams['Ngood'],'Nexcluded':len(ix1)-len(ix2)}
			lcaverageindex = self.averagelc.newrow(df)
			self.averagelc.t.at[lcaverageindex,'Mask'] = 0
			if self.verbose>2: 
				print(self.averagelc.t) 

			# flag clipped sigmacut measurements in original lc
			ix_bad = fluxstatparams['ix_clip']
			if self.verbose>1:
				print('Sigmacut clipped indices (ix_bad): ',ix_bad)
			if len(ix_bad) > 0:
				flag_array = np.full(len(ix_bad),self.flag_daysigma)
				self.lc.t.loc[ix_bad,'Mask'] = np.bitwise_or(self.lc.t.loc[ix_bad,'Mask'],flag_array)

			badflag = 0
			if len(ix2)<=2:
				# flag in original and averaged lc
				flag_array = np.full(len(ix2),self.flag_daysmallnumber)
				self.lc.t.loc[ix2,'Mask'] = np.bitwise_or(self.lc.t.loc[ix2,'Mask'],flag_array)
				self.averagelc.t.at[lcaverageindex,'Mask'] = np.bitwise_or(self.averagelc.t.at[lcaverageindex,'Mask'],self.flag_daysmallnumber)
			else:
				# check sigmacut stats and, if badflag, flag as bad day in original and averaged lc
				if fluxstatparams['Ngood'] < daymin_Ngood: 
					badflag = 1
				if fluxstatparams['Nclip'] > daymax_Nclip: 
					badflag = 1
				if not(fluxstatparams['X2norm'] is None) and (fluxstatparams['X2norm'] > daymax_X2norm): 
					badflag = 1
			if badflag == 1:
				if self.verbose>1:
					print('# Flagged as bad day!')
				flag_array = np.full(len(ix1),self.flag_daybad)
				self.lc.t.loc[ix1,'Mask'] = np.bitwise_or(self.lc.t.loc[ix1,'Mask'],flag_array)
				self.averagelc.t.at[lcaverageindex,'Mask'] = np.bitwise_or(self.averagelc.t.at[lcaverageindex,'Mask'],self.flag_daybad)
			else:
				if self.verbose>1:
					print('# Flagged as good day!')

			MJD += MJDbinsize

		avglcfilename = self.lcbasename(SNindex=SNindex,filt=self.filt,offsetindex=offsetindex,MJDbinsize=MJDbinsize)+'.txt'
		print("Saving averaged light curve: %s" % avglcfilename)
		#self.averagelc.write()
		self.averagelc.write(avglcfilename,overwrite=True,verbose=True)

	def offsetstatsloop(self,SNindex):
		# load main lc
		self.load_lc(SNindex,offsetindex=0,filt=self.filt)
		if self.verbose==1: print('Length of self.lc.t: ',len(self.lc.t))
		if len(self.lc.t)==0: return(1)

		# get sigmacut info and flag for 4-day bins
		print('Making cuts based on day measurement statistics...')
		self.cutandaveragelc(SNindex,offsetindex=0)

		# save lc
		self.save_lc(SNindex=SNindex,offsetindex=0,filt=self.filt,overwrite=True)

		if self.cfg.params['offsetstats']['apply2offsets'] is True:
			# get o1 and o2 masks from SN lc for copying to offset lc mask column
			flags = self.flag_o1_good | self.flag_o2_bad | self.flag_o2_ok | self.flag_o2_good
			flags_array = np.full(self.lc.t['Mask'].shape,flags)
			omask = self.lc.t['Mask'] & flags_array

			# loop through offset lcs
			for offsetindex in range(1,len(self.RADECtable.t)):
				# load offset lc
				self.load_lc(SNindex,offsetindex=offsetindex,filt=self.filt)
				if self.verbose==1: print('Length of self.lc.t: ',len(self.lc.t))
				if len(self.lc.t)==0: return(1)

				# copy over SN o1 and o2 masks to offset mask column
				self.lc.t['Mask'] = np.bitwise_or(self.lc.t['Mask'],omask)

				# get sigmacut info and flag for 4-day bins
				print('Making cuts based on day measurement statistics...')
				self.cutandaveragelc(SNindex,offsetindex=offsetindex)

				# save lc
				self.save_lc(SNindex=SNindex,offsetindex=offsetindex,filt=self.filt,overwrite=True)

if __name__ == '__main__':

	offsetstats = offsetstatsclass()
	parser = offsetstats.define_options()
	args = parser.parse_args()

	SNindexlist = offsetstats.initialize(args)

	for SNindex in SNindexlist:
		offsetstats.loadRADEClist(SNindex,filt=offsetstats.filt)
		offsetstats.offsetstatsloop(SNindex)
