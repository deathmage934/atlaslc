#!/usr/bin/env python

# S. Rest

from SNloop import SNloopclass
from statistics import median
import sigmacut
import sys
import numpy as np
import pandas as pd

class offsetstatsclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)

		self.o1_nanindexlist = []
		self.o2_nanindexlist = []

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

		apply_mask4mjd = False
		apply_mask_nan = False
		if self.cfg.params['offsetstats']['procedure'] == 'mask4mjd':
			print('Procedure set to mask4mjd')
			apply_mask4mjd = True
		elif self.cfg.params['offsetstats']['procedure'] == 'mask_nan':
			print('Procedure set to mask_nan')
			apply_mask_nan = True
		elif self.cfg.params['offsetstats']['procedure'] == 'both':
			print('Procedures set to mask4mjd and mask_nan')
			apply_mask4mjd = True
			apply_mask_nan = True
		else:
			raise RuntimeError("Mask procedure must be 'mask4mjd', 'mask_nan', or 'both' in precursor.cfg!")
		
		for index in range(N_MJD):
			uJy4MJD = uJy[:,index]
			duJy4MJD = duJy[:,index]
			Mask4MJD = Mask[:,index]

			calcaverage=sigmacut.calcaverageclass()
			self.lc.t.at[index,'o0_Nmasked'] = np.count_nonzero(Mask4MJD)
			if apply_mask_nan:
				mask = np.bitwise_and(Mask4MJD, 0x8)
				calcaverage.calcaverage_sigmacutloop(uJy4MJD,noise=duJy4MJD,mask=mask,verbose=2,Nsigma=0.0,median_firstiteration=True,saveused=True)
				# add columns to self.lc.t
				self.lc.t.at[index,'o1_mean'] = calcaverage.mean
				self.lc.t.at[index,'o1_mean_err'] = calcaverage.mean_err
				self.lc.t.at[index,'o1_stddev'] = calcaverage.stdev
				self.lc.t.at[index,'o1_X2norm'] = calcaverage.X2norm
				self.lc.t.at[index,'o1_Nvalid'] = calcaverage.Nused
				self.lc.t.at[index,'o1_Nnan'] = calcaverage.Nskipped
				self.lc.t.at[index,'Noffsetlc'] = self.lc.t.at[index,'o1_Nvalid'] + self.lc.t.at[index,'o1_Nnan']
			if apply_mask4mjd is True:
				calcaverage.calcaverage_sigmacutloop(uJy4MJD,noise=duJy4MJD,mask=Mask4MJD,verbose=2,Nsigma=3.0,median_firstiteration=True,saveused=True)
				# add columns to self.lc.t
				self.lc.t.at[index,'o2_mean'] = calcaverage.mean
				self.lc.t.at[index,'o2_mean_err'] = calcaverage.mean_err
				self.lc.t.at[index,'o2_stddev'] = calcaverage.stdev
				self.lc.t.at[index,'o2_X2norm'] = calcaverage.X2norm
				self.lc.t.at[index,'o2_Nused'] = calcaverage.Nused
				self.lc.t.at[index,'o2_Nskipped'] = calcaverage.Nskipped
				self.lc.t.at[index,'o2_Nin'] = self.lc.t.at[index,'Noffsetlc'] - self.lc.t.at[index,'o0_Nmasked']
			
			use_o1 = False
			use_o2 = False
			# if just one measurement got cut, either poisson noise, cosmic ray hit, or some other issue that only affects 1 measurement
			# if more than 1 measurement are cut, then something is probably not right
			if self.lc.t.at[index,'o2_Nin'] - self.lc.t.at[index,'o2_Nused'] <= 1:
				use_o2 = True
			else:
				use_o1 = True

			Nsigma = 5 # should be from 3-5? check
			X2norm_max = 3 # check
			if use_o2 is True:
				if (np.isnan(self.lc.t.at[index,'o2_mean'])) or (np.isnan(self.lc.t.at[index,'o2_mean_err'])):
					#print('nans detected for ',index,' index! Skipping o2 masking for this measurement...')
					self.o2_nanindexlist = []
					self.o2_nanindexlist.append(index)
				else:
					if abs(self.lc.t.at[index,'o2_mean'] / self.lc.t.at[index,'o2_mean_err']) > Nsigma:
						self.lc.t.loc[index,'Mask'] = np.bitwise_or(self.lc.t.at[index,'Mask'],self.flag_o2_meannorm)
					if self.lc.t.at[index,'o2_X2norm'] > X2norm_max:
						self.lc.t.loc[index,'Mask'] = np.bitwise_or(self.lc.t.at[index,'Mask'],self.flag_o2_X2norm)
			else:
				if (np.isnan(self.lc.t.at[index,'o1_mean'])) or (np.isnan(self.lc.t.at[index,'o1_mean_err'])):
					#print('nans detected for ',index,' index! Skipping o1 masking for this measurement...')
					self.o1_nanindexlist.append(index)
				else:
					if abs(self.lc.t.at[index,'o1_mean'] / self.lc.t.at[index,'o1_mean_err']) > Nsigma:
						self.lc.t.loc[index,'Mask'] = np.bitwise_or(self.lc.t.at[index,'Mask'],self.flag_o1_meannorm)
					if self.lc.t.at[index,'o1_X2norm'] > X2norm_max:
						self.lc.t.loc[index,'Mask'] = np.bitwise_or(self.lc.t.at[index,'Mask'],self.flag_o1_X2norm)

			if calcaverage.Nused<=0:
				if self.verbose>2:
					print('No data values...')

		if len(self.o1_nanindexlist)==0:
			print('No o1 nans detected!')
		else:
			print('o1 nans detected! Index list: ',self.o1_nanindexlist)

		if len(self.o2_nanindexlist)==0:
			print('No o2 nans detected!')
		else:
			print('o2 nans detected! Index list: ',self.o2_nanindexlist)

		# round o1 or o2 data
		if 'o1_mean' in self.lc.t.columns:
			self.lc.t = self.lc.t.round({'o1_mean':3,'o1_mean_err':3,'o1_stddev':3,'o1_X2norm':4})
		if 'o2_mean' in self.lc.t.columns:
			self.lc.t = self.lc.t.round({'o2_mean':3,'o2_mean_err':3,'o2_stddev':3,'o2_X2norm':4})
		self.save_lc(SNindex=SNindex,offsetindex=0,filt=self.filt,overwrite=True)

if __name__ == '__main__':

	offsetstats = offsetstatsclass()
	parser = offsetstats.define_options()
	args = parser.parse_args()

	SNindexlist = offsetstats.initialize(args)

	for SNindex in SNindexlist:
		offsetstats.loadRADEClist(SNindex, filt=offsetstats.filt)
		offsetstats.offsetstatsloop(SNindex,filt=offsetstats.filt)
