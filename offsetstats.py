#!/usr/bin/env python

# S. Rest

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

	def cleanmask(self,indices=None):
		# cleans mask data if prior o1 and/or o2 columns detected
		if ('o2_mean' in self.lc.t.columns) and ('o1_mean' in self.lc.t.columns):
			self.lc.t.loc[indices,'Mask'] = np.bitwise_xor(self.lc.t.loc[indices,'Mask'],self.flag_o1_good+self.flag_o2_good+self.flag_o2_ok+self.flag_o2_bad)
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

	def cutandaveragelc(self,SNindex,offsetindex=0,MJDbinsize=1):
		daymax_Nclip = self.cfg.params['offsetstats']['dayflagging']['daymax_Nclip']
		daymin_Ngood = self.cfg.params['offsetstats']['dayflagging']['daymin_Ngood']
		daymax_X2norm = self.cfg.params['offsetstats']['dayflagging']['daymax_X2norm']
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
			ix2 = self.lc.ix_unmasked('Mask',maskval=self.flag_o2_bad|self.flag_cut0_X2norm_static,indices=ix1)
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
				self.lc.t.loc[ix1,'Mask'] = np.bitwise_or(self.lc.t.loc[ix1,'Mask'],self.flag_badday)
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
				if fluxstatparams['Nclip'] > daymax_Nclip: 
					badflag = 1
				if fluxstatparams['Ngood'] < daymin_Ngood: 
					badflag = 1
				if (fluxstatparams['X2norm'] > daymax_X2norm): 
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

	def offsetstatsloop(self,SNindex,filt):
		# load main SN lc
		self.load_lc(SNindex,offsetindex=0,filt=self.filt)
		if self.verbose:
			print('Length of self.lc.t: ',len(self.lc.t))
		if len(self.lc.t) == 0:
			return(1)

		N_lc = len(self.RADECtable.t)
		MJD_SN = self.lc.t['MJD']
		N_MJD = len(self.lc.t['MJD'])

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

			# make sure MJD_SN is the same as self.lc.t['MJD']
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
		
		# clear o1 and o2 masks
		print('Cleaning mask column...')
		indices = self.lc.getindices()
		self.cleanmask(indices=indices)
		
		# calculate offset stats
		print('Calculating offset statistics...')
		self.calcstats(N_MJD=N_MJD,uJy=uJy,duJy=duJy,mask=Mask)	
		
		# make cuts on good, bad, and ok measurements
		print('Making cuts based on offset statistics...')
		self.makecuts(N_MJD=N_MJD)

		# get sigmacut info and flag for 4-day bins
		print('Making cuts based on day measurement statistics...')
		self.cutandaveragelc(SNindex,offsetindex=0)

		# round o1 or o2 data and save lc
		#if 'o1_mean' in self.lc.t.columns:
			#self.lc.t = self.lc.t.round({'o1_mean':3,'o1_mean_err':3,'o1_stddev':3,'o1_X2norm':4})
		#if 'o2_mean' in self.lc.t.columns:
			#self.lc.t = self.lc.t.round({'o2_mean':3,'o2_mean_err':3,'o2_stddev':3,'o2_X2norm':4})
		self.save_lc(SNindex=SNindex,offsetindex=0,filt=self.filt,overwrite=True)

		if self.cfg.params['offsetstats']['dayflagging']['apply2offsets'] is True:
			# get o1 and o2 masks from SN lc for copying to offset lc mask column
			flags = self.flag_o1_good | self.flag_o2_bad | self.flag_o2_ok | self.flag_o2_good
			flags_array = np.full(self.lc.t['Mask'].shape,flags)
			omask = self.lc.t['Mask'] & flags_array

			# loop through offset lcs
			for offsetindex in range(1,len(self.RADECtable.t)):
				#if offsetindex>=1:
					#sys.exit(0)
				# load offset lc
				self.load_lc(SNindex,offsetindex=offsetindex,filt=self.filt)
				if self.verbose==1: 
					print('Length of self.lc.t: ',len(self.lc.t))
				if len(self.lc.t)==0: return(1)

				# copy over SN o1 and o2 masks to offset mask column
				self.lc.t['Mask'] = np.bitwise_or(self.lc.t['Mask'],omask)

				# get sigmacut info and flag for 4-day bins
				print('Making cuts based on day measurement statistics...')
				self.cutandaveragelc(SNindex,offsetindex=offsetindex)

				if len(self.o1_nanindexlist)==0:
					print('No o1 nans detected!')
				else:
					print('o1 nans detected! Index list: ',self.o1_nanindexlist)

				if len(self.o2_nanindexlist)==0:
					print('No o2 nans detected!')
				else:
					print('o2 nans detected! Index list: ',self.o2_nanindexlist)

				# round o1 or o2 data and save lc
				if 'o1_mean' in self.lc.t.columns:
					self.lc.t = self.lc.t.round({'o1_mean':3,'o1_mean_err':3,'o1_stddev':3,'o1_X2norm':4})
				if 'o2_mean' in self.lc.t.columns:
					self.lc.t = self.lc.t.round({'o2_mean':3,'o2_mean_err':3,'o2_stddev':3,'o2_X2norm':4})
				self.save_lc(SNindex=SNindex,offsetindex=offsetindex,filt=self.filt,overwrite=True)

if __name__ == '__main__':

	offsetstats = offsetstatsclass()
	parser = offsetstats.define_options()
	args = parser.parse_args()

	SNindexlist = offsetstats.initialize(args)

	for SNindex in SNindexlist:
		offsetstats.loadRADEClist(SNindex, filt=offsetstats.filt)
		offsetstats.offsetstatsloop(SNindex,filt=offsetstats.filt)
