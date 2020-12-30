#!/usr/bin/env python
'''
@author: S. Rest
'''

from SNloop import SNloopclass
from statistics import median
from astropy.table import Table, Column, MaskedColumn
from copy import deepcopy
import sigmacut
import sys
import numpy as np
import pandas as pd

class cleanuplcclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)

	def o0_PSF_uncertainty_cut(self):
		# define vars
		o0_Nmedian = self.cfg.params['cleanlc']['o0']['PSF_uncertainty']['o0_Nmedian']
		a = o0_Nmedian * median(self.lc.t[self.dflux_colname])
		print('Flagging all measurements with %s bigger than %i...' % (self.dflux_colname, a))

		# define indices
		a_indices = np.where(self.lc.t[self.dflux_colname]>a)
		a_indices = list(a_indices[0])
		if self.verbose>1: print('Indices: ',a_indices)
		if len(self.lc.t.loc[a_indices,self.dflux_colname])>0:
			if self.verbose: print('%s above %i: ' % (self.dflux_colname, a), len(self.lc.t.loc[a_indices,self.dflux_colname]))
		else:
			if self.verbose: print('No measurements flagged!')

		# update 'Mask' column
		flag_o0_uncertainty = np.full(self.lc.t.loc[a_indices,'Mask'].shape, self.flag_o0_uncertainty)
		self.lc.t.loc[a_indices,'Mask'] = np.bitwise_or(self.lc.t.loc[a_indices,'Mask'],flag_o0_uncertainty) 
	
	# ended up not working, will not use but keep just in case
	'''
	def flag_bigchi_dynamic(self):
		# define vars
		Nsigma = self.cfg.params['cleanlc']['chi/N']['Nsigma']
		chi_median = median(self.lc.t['chi/N'])
		chi_stddev = np.std(self.lc.t['chi/N'])
		a = int(chi_median+(Nsigma*chi_stddev)) # !! CURRENTLY ROUNDS DOWN
		print('Flagging all measurements with chi/N bigger than %i...' % a)

		# define indices
		a_indices = np.where(self.lc.t['chi/N']>a)
		a_indices = list(a_indices[0])
		if self.verbose>1: print('Indices: ',a_indices) 
		if len(self.lc.t.loc[a_indices,'chi/N'])>0:
			if self.verbose: print('chi/N above %i: ' % a,len(self.lc.t.loc[a_indices,'chi/N']))
		else:
			if self.verbose: print('No measurements flagged!')

		# update 'Mask' column
		flag_cut0_X2norm_dynamic = np.full(self.lc.t.loc[a_indices,'Mask'].shape, self.flag_cut0_X2norm_dynamic)
		self.lc.t.loc[a_indices,'Mask'] = np.bitwise_or(self.lc.t.loc[a_indices,'Mask'], flag_cut0_X2norm_dynamic)
		print('Nsigma: %.1f, chi_median: %f, chi_stddev: %f' % (Nsigma, chi_median, chi_stddev))
	'''
	
	def o0_PSF_X2norm_cut(self):
		# define vars
		o0max_X2norm = self.cfg.params['cleanlc']['o0']['PSF_X2norm']['o0max_X2norm']
		print('Flagging all measurements with chi/N bigger than %i...' % o0max_X2norm)
		
		# define indices
		a_indices = np.where(self.lc.t['chi/N']>o0max_X2norm)
		a_indices = list(a_indices[0])
		if self.verbose>1: print('Indices: ',a_indices) 
		if len(self.lc.t.loc[a_indices,'chi/N'])>0:
			if self.verbose: print('chi/N above %i: ' % o0max_X2norm,len(self.lc.t.loc[a_indices,'chi/N']))
		else:
			if self.verbose: print('No measurements flagged!')		

		# update 'Mask' column
		flag_o0_X2norm = np.full(self.lc.t.loc[a_indices,'Mask'].shape, self.flag_o0_X2norm)
		self.lc.t.loc[a_indices,'Mask'] = np.bitwise_or(self.lc.t.loc[a_indices,'Mask'], flag_o0_X2norm)

	# no use for this now, but keep just in case
	'''
	def cleanmask(self,indices=None):
		# cleans mask data if prior o1 and/or o2 columns detected
		if ('o2_mean' in self.lc.t.columns) and ('o1_mean' in self.lc.t.columns):
			self.lc.t.loc[indices,'Mask'] = np.bitwise_xor(self.lc.t.loc[indices,'Mask'],self.flag_o1_good+self.flag_o2_good+self.flag_o2_ok+self.flag_o2_bad)
		else:
			print('No prior mask4mjd or mask_nan data detected!')
	'''
	# to initialize:
	#print('Cleaning mask column...')
	#indices = self.lc.getindices()
	#self.cleanmask(indices=indices)

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

			# mask_nan (o1) sigmacut
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
			
			# mask4mjd (o2) sigmacut
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

		o0max_X2norm = self.cfg.params['cleanlc']['o0']['PSF_X2norm']['o0max_X2norm']
		o1max_X2norm = self.cfg.params['cleanlc']['o1']['o1max_X2norm']
		o1max_meannorm = self.cfg.params['cleanlc']['o1']['o1max_meannorm']
		o2max_Nclipped = self.cfg.params['cleanlc']['o2']['o2max_Nclipped']
		o2max_Nused = self.cfg.params['cleanlc']['o2']['o2max_Nused']

		for index in range(N_MJD):
			# check self.flag_o0_uncertainty = 0x1 and self.flag_o0_X2norm = 0x2. TO DO: CUT0 UNCERTAINTIES
			if self.lc.t.at[index,'chi/N'] < o0max_X2norm:
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
			# else: bad. do nothing because measurements already flagged with o0		

	def cleanuplcloop(self,args,SNindex):
		# o0 - mask lcs based on PSF X2norm and uncertainty

		for offsetindex in range(len(cleanuplc.RADECtable.t)):
			# load lc
			self.load_lc(SNindex, filt=self.filt, offsetindex=offsetindex, MJDbinsize=None)
			print('Length of self.lc.t: ',len(self.lc.t))
			if len(self.lc.t) == 0:
				return(1)

			# clear mask column
			self.lc.t['Mask'] = 0

			if self.cfg.params['cleanlc']['o0']['PSF_uncertainty']['apply'] is True:
				print('Applying uncertainty cleanup...')
				self.o0_PSF_uncertainty_cut()
			else:
				print('Skipping uncertainty cleanup...')
			if self.cfg.params['cleanlc']['o0']['PSF_X2norm']['apply'] is True:
				print('Applying chi/N static cleanlc')
				self.o0_PSF_X2norm_cut()
			else:
				print('Skipping chi/N cleanup...')

			self.save_lc(SNindex=SNindex,filt=self.filt,overwrite=True,offsetindex=offsetindex)

		# o1 and o2 - mask lcs based on offset sigmacut

		if self.cfg.params['cleanlc']['apply_o1o2'] is True:
			print('o1 and o2 masking in progress...')
			# load main SN lc
			self.load_lc(SNindex,offsetindex=0,filt=self.filt)
			if self.verbose:
				print('Length of self.lc.t: ',len(self.lc.t))
			if len(self.lc.t) == 0:
				return(1)

			N_lc = len(self.RADECtable.t)
			MJD_SN = self.lc.t['MJD']
			N_MJD = len(self.lc.t['MJD'])

			# construct arrays for offset data
			uJy = np.full((N_lc,N_MJD),np.nan,dtype=np.int64)
			duJy = np.full((N_lc,N_MJD),np.nan,dtype=np.int64)
			Mask = np.full((N_lc,N_MJD),np.nan,dtype=np.int32)

			for offsetindex in range(1,len(self.RADECtable.t)):
				# load offset lc
				self.load_lc(SNindex,offsetindex=offsetindex,filt=self.filt)
				if self.verbose>=1:
					print('Length of offset light curve: ',len(self.lc.t))
				if len(self.lc.t) == 0:
					return(1)

				# make sure MJD_SN is the same as self.lc.t['MJD'], then fill array of offset uJy, duJy, Mask
				if (len(self.lc.t) != N_MJD) or (np.array_equal(MJD_SN, self.lc.t['MJD']) is False):
					if self.verbose>1: 
						print(MJD_SN,'\n',self.lc.t['MJD'])
					print('WARNING: Offset lc not equal to SN lc')
					counter = 0
					for mjd_index in range(N_MJD):
						if counter >= len(self.lc.t['MJD']):
							break
						foundflag = False
						if MJD_SN[mjd_index] == self.lc.t.at[counter,'MJD']:
							foundflag = True
						else:
							if MJD_SN[mjd_index] > self.lc.t.at[counter,'MJD']:
								while self.lc.t.at[counter,'MJD']<MJD_SN[mjd_index] and counter<len(self.lc.t['MJD'])-1:
									counter += 1
								if MJD_SN[mjd_index] == self.lc.t.at[counter,'MJD']:
									foundflag = True
						if foundflag:
							uJy[offsetindex,mjd_index] = self.lc.t.at[counter,self.flux_colname]
							duJy[offsetindex,mjd_index] = self.lc.t.at[counter,self.dflux_colname]
							Mask[offsetindex,mjd_index] = self.lc.t.at[counter,'Mask']
							counter += 1
						else:
							uJy[offsetindex,mjd_index] = 0
							duJy[offsetindex,mjd_index] = 0
							Mask[offsetindex,mjd_index] = 0x8
				else:
					uJy[offsetindex,:] = self.lc.t[self.flux_colname]
					duJy[offsetindex,:] = self.lc.t[self.dflux_colname]
					Mask[offsetindex,:] = self.lc.t['Mask']
			
			if self.verbose>1:
				print(self.flux_colname,': ',uJy)
				print(self.dflux_colname,': ',duJy)
				print('Mask: ',Mask)

			# load main lc
			self.load_lc(SNindex,offsetindex=0,filt=self.filt)
			if self.verbose==1: print('Length of self.lc.t: ',len(self.lc.t))
			if len(self.lc.t)==0: return(1)
			
			# calculate offset stats
			print('Calculating offset statistics...')
			self.calcstats(N_MJD=N_MJD,uJy=uJy,duJy=duJy,mask=Mask)	
			
			# make cuts on good, bad, and ok measurements
			print('Making cuts based on offset statistics...')
			self.makecuts(N_MJD=N_MJD)

			# save lc
			self.save_lc(SNindex=SNindex,offsetindex=0,filt=self.filt,overwrite=True)

		else:
			print('Skipping o1 and o2 masking...')

if __name__ == '__main__':

	cleanuplc = cleanuplcclass()
	parser = cleanuplc.define_options()
	args = parser.parse_args()

	SNindexlist = cleanuplc.initialize(args)

	for SNindex in SNindexlist:
		cleanuplc.loadRADEClist(SNindex=SNindex,filt=cleanuplc.filt)
		cleanuplc.cleanuplcloop(args,SNindex)

	print('\n')
