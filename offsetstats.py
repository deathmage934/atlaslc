#!/usr/bin/env python

from SNloop import SNloopclass
from statistics import median
import sigmacut
import sys
import numpy as np
import pandas as pd

class offsetstatsclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)

	def offsetstatsloop(self,SNindex,filt):
		# load main lc
		self.load_lc(SNindex,offsetindex=0,filt=self.filt)
		if self.verbose:
			print('Length of self.lc.t: ',len(self.lc.t))
		if len(self.lc.t) == 0:
			return(1)

		N_lc = len(self.RADECtable.t)-1
		MJD_SN = self.lc.t['MJD']
		# debug
		#MJD_SN.append(pd.Series(59096.351849))	
		N_MJD = len(self.lc.t['MJD'])

		uJy = np.full((N_lc,N_MJD),np.nan,dtype=np.int64)
		duJy = np.full((N_lc,N_MJD),np.nan,dtype=np.int64)
		Mask = np.full((N_lc,N_MJD),np.nan,dtype=np.int32)

		for offsetindex in range(1,len(self.RADECtable.t)):
			self.load_lc(SNindex,offsetindex=offsetindex,filt=self.filt)
			if self.verbose:
				print('Length of self.lc.t: ',len(self.lc.t))
			if len(self.lc.t) == 0:
				return(1)

			# make sure MJD_SN is the same as self.lc.t['MJD']
			if (len(self.lc.t) != N_MJD) or (np.array_equal(MJD_SN, self.lc.t['MJD']) is False):
				if self.verbose:
					print(MJD_SN,'\n',self.lc.t['MJD']) # delete me
				print('WARNING: Offset lc not equal to SN lc')
				counter = 0
				for mjd_index in range(N_MJD):
					if counter >= len(self.lc.t['MJD']):
						break
					foundflag = False
					if MJD_SN[mjd_index] == self.lc.t['MJD'][counter]:
						foundflag = True
					else:
						if MJD_SN[mjd_index] > self.lc.t['MJD'][counter]:
							while self.lc.t['MJD'][counter]<MJD_SN[mjd_index] and counter<len(self.lc.t['MJD'])-1:
								counter += 1
							if MJD_SN[mjd_index] == self.lc.t['MJD'][counter]:
								foundflag = True
					if foundflag:
						uJy[offsetindex-1,mjd_index] = self.lc.t[self.flux_colname][counter]
						duJy[offsetindex-1,mjd_index] = self.lc.t[self.dflux_colname][counter]
						Mask[offsetindex-1,mjd_index] = self.lc.t['Mask'][counter]
						counter += 1
					else:
						uJy[offsetindex-1,mjd_index] = 0
						duJy[offsetindex-1,mjd_index] = 0
						Mask[offsetindex-1,mjd_index] = 0x8
			else:
				uJy[offsetindex-1,:] = self.lc.t[self.flux_colname]
				duJy[offsetindex-1,:] = self.lc.t[self.dflux_colname]
				Mask[offsetindex-1,:] = self.lc.t['Mask']
		if self.verbose>1:
			print('%s: ',uJy % self.flux_colname)
			print('%s: ',duJy % self.dflux_colname)
			print('Mask: ',Mask)

		self.load_lc(SNindex,offsetindex=0,filt=self.filt)
		if self.verbose:
			print('Length of self.lc.t: ',len(self.lc.t))
		if len(self.lc.t) == 0:
			return(1)

		if self.cfg.params['offsetstats']['procedure'] == 'mask1':
			print('Procedure mask1') # FIX
		else:
			print('Procedure mask2') # FIX
		for index in range(N_MJD):
			uJy4MJD = uJy[:,index]
			duJy4MJD = duJy[:,index]
			Mask4MJD = Mask[:,index]

			calcaverage=sigmacut.calcaverageclass()
			if self.cfg.params['offsetstats']['procedure'] == 'mask1':
				calcaverage.calcaverage_sigmacutloop(uJy4MJD,noise=duJy4MJD,mask=Mask4MJD,verbose=2,Nsigma=3.0,median_firstiteration=True,saveused=True)
				# add columns to self.lc.t
				self.lc.t.at[index,'o1_mean'] = calcaverage.mean
				self.lc.t.at[index,'o1_mean_err'] = calcaverage.mean_err
				self.lc.t.at[index,'o1_stddev'] = calcaverage.stdev
				self.lc.t.at[index,'o1_chi/N'] = calcaverage.X2norm
				self.lc.t.at[index,'o1_Nused'] = calcaverage.Nused
				self.lc.t.at[index,'o1_Nskipped'] = calcaverage.Nskipped
			elif self.cfg.params['offsetstats']['procedure'] == 'mask2':
				mask = np.bitwise_and(Mask4MJD, 0x8)
				calcaverage.calcaverage_sigmacutloop(uJy4MJD,noise=duJy4MJD,mask=mask,verbose=2,Nsigma=0.0,median_firstiteration=True,saveused=True)
				# add columns to self.lc.t
				self.lc.t.at[index,'o2_mean'] = calcaverage.mean
				self.lc.t.at[index,'o2_mean_err'] = calcaverage.mean_err
				self.lc.t.at[index,'o2_stddev'] = calcaverage.stdev
				self.lc.t.at[index,'o2_chi/N'] = calcaverage.X2norm
				self.lc.t.at[index,'o2_Nused'] = calcaverage.Nused
				self.lc.t.at[index,'o2_Nskipped'] = calcaverage.Nskipped
			else:
				raise RuntimeError("Mask procedure must be mask1 or mask2 in config file!")
			if calcaverage.Nused<=0:
				if self.verbose>2:
					print('No data values...')

		# round o1 or o2 data
		#self.lc.t = self.lc.formattable(roundingMapping={'o1_mean':3,'o1_mean_err':3,'o1_stddev':3,'o1_chi/N':4})#,dtypeMapping={'o1_mean':np.float64,'o1_mean_err':np.float64,'o1_stddev':np.float64,'o1_chi/N':np.float64,'o1_Nused':np.int64,'o1_Nskipped':np.int64})
		if 'o1_mean' in self.lc.t.columns:
			self.lc.t = self.lc.t.round({'o1_mean':3,'o1_mean_err':3,'o1_stddev':3,'o1_chi/N':4})
		if 'o2_mean' in self.lc.t.columns:
			self.lc.t = self.lc.t.round({'o2_mean':3,'o2_mean_err':3,'o2_stddev':3,'o2_chi/N':4})
		self.save_lc(SNindex,offsetindex=0,filt=self.filt,overwrite=True)

if __name__ == '__main__':

	offsetstats = offsetstatsclass()
	parser = offsetstats.define_options()
	args = parser.parse_args()

	SNindexlist = offsetstats.initialize(args)

	for SNindex in SNindexlist:
		offsetstats.loadRADEClist(SNindex, filt=offsetstats.filt)
		offsetstats.offsetstatsloop(SNindex,filt=offsetstats.filt)
