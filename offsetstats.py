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

		self.apply_mask_nan = False
		self.apply_mask4mjd = False

	def cleanmask(self,indices=None):
		# cleans mask data if prior o1 and/or o2 columns detected
		if (self.apply_mask4mjd is True) and (self.apply_mask_nan is True):
			if ('o2_mean' in self.lc.t.columns) and ('o1_mean' in self.lc.t.columns):
				self.lc.t.loc[indices,'Mask'] = np.bitwise_xor(self.lc.t.loc[indices,'Mask'],self.flag_o1_good+self.flag_o2_good+self.flag_o2_ok+self.flag_o2_bad)
		elif self.apply_mask4mjd is True:
			if 'o2_mean' in self.lc.t.columns:
				self.lc.t.loc[indices,'Mask'] = np.bitwise_xor(self.lc.t.loc[indices,'Mask'],self.flag_o2_good+self.flag_o2_ok+self.flag_o2_bad)
		elif self.apply_mask_nan is True:
			if 'o1_mean' in self.lc.t.columns:
				self.lc.t.loc[indices,'Mask'] = np.bitwise_xor(self.lc.t.loc[indices,'Mask'],self.flag_o1_good)
		else:
			print('No prior mask4mjd or mask_nan data detected!')

	def calcstats(self,N_MJD=None,uJy=None,duJy=None,mask=None):
		if N_MJD is None:
			N_MJD = len(self.lc.t['MJD'])

		for index in range(N_MJD):
			uJy4MJD = uJy[:,index]
			duJy4MJD = duJy[:,index]
			Mask4MJD = mask[:,index]
			# sigmacut and get statistics
			calcaverage=sigmacut.calcaverageclass()
			self.lc.t.at[index,'o0_Nmasked'] = np.count_nonzero(Mask4MJD)
			if self.apply_mask_nan:
				mask1 = np.bitwise_and(Mask4MJD, 0x8)
				calcaverage.calcaverage_sigmacutloop(uJy4MJD,noise=duJy4MJD,mask=mask1,verbose=2,Nsigma=0.0,median_firstiteration=True,saveused=True)
				# add columns to self.lc.t
				self.lc.t.at[index,'o1_mean'] = calcaverage.mean
				self.lc.t.at[index,'o1_mean_err'] = calcaverage.mean_err
				self.lc.t.at[index,'o1_stddev'] = calcaverage.stdev
				self.lc.t.at[index,'o1_X2norm'] = calcaverage.X2norm
				self.lc.t.at[index,'o1_Nvalid'] = calcaverage.Nused
				self.lc.t.at[index,'o1_Nnan'] = calcaverage.Nskipped
				self.lc.t.at[index,'Noffsetlc'] = self.lc.t.at[index,'o1_Nvalid'] + self.lc.t.at[index,'o1_Nnan']
			if self.apply_mask4mjd is True:
				calcaverage.calcaverage_sigmacutloop(uJy4MJD,noise=duJy4MJD,mask=Mask4MJD,verbose=2,Nsigma=3.0,median_firstiteration=True,saveused=True)
				# add columns to self.lc.t
				self.lc.t.at[index,'o2_mean'] = calcaverage.mean
				self.lc.t.at[index,'o2_mean_err'] = calcaverage.mean_err
				self.lc.t.at[index,'o2_stddev'] = calcaverage.stdev
				self.lc.t.at[index,'o2_X2norm'] = calcaverage.X2norm
				self.lc.t.at[index,'o2_Nused'] = calcaverage.Nused
				self.lc.t.at[index,'o2_Nskipped'] = calcaverage.Nskipped
				self.lc.t.at[index,'o2_Nin'] = self.lc.t.at[index,'Noffsetlc'] - self.lc.t.at[index,'o0_Nmasked']

	def makecuts(self,N_MJD=None):
		if N_MJD is None:
			N_MJD = len(self.lc.t['MJD'])

		# check self.flag_cut0_uncertainty = 0x1 and self.flag_cut0_X2norm_static = 0x4. TO DO: CUT0 UNCERTAINTIES
		max_X2norm = self.cfg.params['cleanlc']['chi/N']['max_chi2norm']
		o1max_X2norm = self.cfg.params['offsetstats']['o1max_X2norm']
		o1max_meannorm = self.cfg.params['offsetstats']['o1max_meannorm']
		o2max_Nclipped = self.cfg.params['offsetstats']['o2max_Nclipped']
		o2max_Nused = self.cfg.params['offsetstats']['o2max_Nused']

		for index in range(N_MJD):
			if self.lc.t.at[index,'chi/N'] < max_X2norm:
				# if o1_x2norm < 2.5 and o1_mean_err < 3.0 : good. else : o2
				if (self.lc.t.at[index,'o1_X2norm']<o1max_X2norm) and (self.lc.t.at[index,'o1_mean_err']<o1max_meannorm):
					self.lc.t.at[index,'Mask'] = np.bitwise_or(self.lc.t.at[index,'Mask'],self.flag_o1_good)
				else:
					# if o2_Nskipped = 0 and o2_X2norm < 2.5 : ok1.
					if (self.lc.t.at[index,'o2_Nskipped']==0) and (self.lc.t.at[index,'o2_X2norm']<2.5):
						self.lc.t.at[index,'Mask'] = np.bitwise_or(self.lc.t.at[index,'Mask'],self.flag_o2_good)
					# if o2_X2norm < 2.5 and o2_Nskipped <= 1 and o2_Nused >= 3 : ok2.
					elif (self.lc.t.at[index,'o2_X2norm']<2.5) and (self.lc.t.at[index,'o2_Nskipped'].astype(int)<=o2max_Nclipped) and (self.lc.t.at[index,'o2_Nused'].astype(int)>=o2max_Nused):
						self.lc.t.at[index,'Mask'] = np.bitwise_or(self.lc.t.at[index,'Mask'],self.flag_o2_ok)
					# else: bad.
					else:
						self.lc.t.at[index,'Mask'] = np.bitwise_or(self.lc.t.at[index,'Mask'],self.flag_o2_bad)
			# else: bad. do nothing bc measurements already flagged with cut0

		'''
		# if chi/N < max_chi2norm : o1. else: bad
		indices = np.where(self.lc.t['chi/N']<max_chi2norm)
		indices = list(indices[0])
		# if o1_x2norm < 2.5 and o1_mean_err < 3.0 : good. else : o2
		if (self.lc.t.loc[indices,'o1_X2norm']<o1max_chi2norm) and (self.lc.t.loc[indices,'o1_mean_err']<o1max_meannorm):
			self.lc.t.loc[indices,'Mask'] = np.bitwise_or(self.lc.t.loc[indices,'Mask'],self.flag_o1_good)
		else:
			# if o2_Nskipped = 0 and o2_X2norm < 2.5 : ok1.
			if (self.lc.t.loc[indices,'o2_Nskipped']==0) and (self.lc.t.loc[indices,'o2_X2norm']<2.5):
				self.lc.t.loc[indices,'Mask'] = np.bitwise_or(self.lc.t.loc[indices,'Mask'],self.flag_o2_good)
			# if o2_X2norm < 2.5 and o2_Nskipped <= 1 and o2_Nused >= 3 : ok2.
			elif (self.lc.t.loc[indices,'o2_X2norm']<2.5) and (self.lc.t.loc[indices,'o2_Nskipped'].astype(int)<=1) and (self.lc.t.loc[indices,'o2_Nused'].astype(int)>=3):
				self.lc.t.loc[indices,'Mask'] = np.bitwise_or(self.lc.t.loc[indices,'Mask'],self.flag_o2_ok)
			# else: bad.
			else:
				self.lc.t.loc[indices,'Mask'] = np.bitwise_or(self.lc.t.loc[indices,'Mask'],self.flag_o2_bad)
		'''

	def daycuts(self,SNindex):
		daymax_Nskipped = self.cfg.params['offsetstats']['dayflagging']['daymax_Nskipped']
		daymin_Nused = self.cfg.params['offsetstats']['dayflagging']['daymin_Nused']
		daymax_X2norm = self.cfg.params['offsetstats']['dayflagging']['daymax_X2norm']

		MJD = np.amin(self.lc.t.loc['MJD'])
		MJDmax = self.t.at[SNindex,'MJDpreSN']
		
		while MJD <= MJDmax:
			# get 4 measurements
			ix1 = self.lc.ix_inrange(colnames=['MJD'],lowlim=MJD,uplim=MJD+1)
			
			# get good and ok measurements (from offsetstats) out of 4
			ix2 = self.lc.ix_unmasked('Mask',maskval=self.flag_o2_bad|self.flag_cut0_X2norm_static,indices=ix1)
			
			# sigmacut the 4 measurements
			self.lc.calcaverage_sigmacutloop('uJy',noisecol='duJy',indices=ix2,verbose=2,Nsigma=3.0,median_firstiteration=True)

			badflag = 0
			if calcaverage.Nskipped > daymax_Nskipped:
				badflag = 1
			if calcaverage.Nused < daymin_Nused:
				badflag = 1
			if calcaverage.X2norm > daymax_X2norm:
				badflag = 1

			if badflag == 1:
				ix_bad = statparams['ix_bad']
				print(ix_bad) # delete me
				sys.exit(0) # delete me
				self.lc.t.at[ix_bad,'Mask'] = np.bitwise_or(self.lc.t.at[ix_bad,'Mask'],self.flag_daysigma)

			MJD += 4

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

		if self.cfg.params['offsetstats']['procedure'] == 'mask4mjd':
			print('Procedure set to mask4mjd')
			self.apply_mask4mjd = True
		elif self.cfg.params['offsetstats']['procedure'] == 'mask_nan':
			print('Procedure set to mask_nan')
			self.apply_mask_nan = True
		elif self.cfg.params['offsetstats']['procedure'] == 'both':
			print('Procedures set to mask4mjd and mask_nan')
			apply_mask4mjd = True
			self.apply_mask_nan = True
		else:
			raise RuntimeError("Mask procedure must be 'mask4mjd', 'mask_nan', or 'both' in precursor.cfg!")

		indices = self.lc.getindices()
		# clear o1 and o2 masks
		self.cleanmask(indices=indices)
		# calculate offset stats
		self.calcstats(N_MJD=N_MJD,uJy=uJy,duJy=duJy,mask=Mask)
		# make cuts on good, bad, and ok measurements
		self.makecuts(N_MJD=N_MJD)
		# get sigmacut info and flag for 4-day bins
		self.daycuts(SNindex)

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
