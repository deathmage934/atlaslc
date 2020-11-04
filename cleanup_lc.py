#!/usr/bin/env python

from SNloop import SNloopclass
from statistics import median
import numpy as np
from astropy.table import Table, Column, MaskedColumn

class cleanuplcclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)

	def flag_biguncertainty_duJy(self,SNindex,offsetindex):
		# define vars
		Nmedian = self.cfg.params['cleanlc']['uncertainty']['Nmedian']
		a = Nmedian * median(self.lc.t[self.dflux_colname])
		print('Flagging all measurements with %s bigger than %i...' % (self.dflux_colname, a))

		# define indices
		a_indices = np.where(self.lc.t[self.dflux_colname]>a)
		a_indices = list(a_indices[0])
		print('Indices: ',a_indices)
		if len(self.lc.t.loc[a_indices,self.dflux_colname])>0:
			if self.verbose:
				print('%s above %i: ' % (self.dflux_colname, a), len(self.lc.t.loc[a_indices,self.dflux_colname]))
		else:
			if self.verbose:
				print('No measurements flagged!')

		# update 'Mask' column
		flag_cut0_uncertainty = np.full(self.lc.t.loc[a_indices,'Mask'].shape, self.flag_cut0_uncertainty)
		self.lc.t.loc[a_indices,'Mask'] = np.bitwise_or(self.lc.t.loc[a_indices,'Mask'],flag_cut0_uncertainty) 

	def flag_bigchi_dynamic(self,SNindex,offsetindex):
		# define vars
		Nsigma = self.cfg.params['cleanlc']['chi/N']['Nsigma']
		chi_median = median(self.lc.t['chi/N'])
		chi_stddev = np.std(self.lc.t['chi/N'])
		a = int(chi_median+(Nsigma*chi_stddev)) # !! CURRENTLY ROUNDS DOWN
		print('Flagging all measurements with chi/N bigger than %i...' % a)

		# define indices
		a_indices = np.where(self.lc.t['chi/N']>a)
		a_indices = list(a_indices[0])
		print('Indices: ',a_indices) 
		if len(self.lc.t.loc[a_indices,'chi/N'])>0:
			if self.verbose:
				print('chi/N above %i: ' % a,len(self.lc.t.loc[a_indices,'chi/N']))
		else:
			if self.verbose:
				print('No measurements flagged!')

		# update 'Mask' column
		flag_cut0_X2norm_dynamic = np.full(self.lc.t.loc[a_indices,'Mask'].shape, self.flag_cut0_X2norm_dynamic)
		self.lc.t.loc[a_indices,'Mask'] = np.bitwise_or(self.lc.t.loc[a_indices,'Mask'], flag_cut0_X2norm_dynamic)

		if self.verbose:
			print('Nsigma: %.1f, chi_median: %f, chi_stddev: %f' % (Nsigma, chi_median, chi_stddev))

	def flag_bigchi_static(self,SNindex,offsetindex):
		# define vars
		chi_max = self.cfg.params['cleanlc']['chi/N']['max_chi2norm']
		print('Flagging all measurements with chi/N bigger than %i...' % chi_max)
		
		# define indices
		a_indices = np.where(self.lc.t['chi/N']>chi_max)
		a_indices = list(a_indices[0])
		if len(self.lc.t.loc[a_indices,'chi/N'])>0:
			if self.verbose:
				print('chi/N above %i: ' % chi_max,len(self.lc.t.loc[a_indices,'chi/N']))
		else:
			if self.verbose:
				print('No measurements flagged!')		

		# update 'Mask' column
		flag_cut0_X2norm_static = np.full(self.lc.t.loc[a_indices,'Mask'].shape, self.flag_cut0_X2norm_static)
		self.lc.t.loc[a_indices,'Mask'] = np.bitwise_or(self.lc.t.loc[a_indices,'Mask'],flag_cut0_X2norm_static)		

		if self.verbose:
			print('chi_max: ',chi_max)

	def cleanuplcloop(self,args,SNindex,offsetindex,filt):
		# load lc
		self.load_lc(SNindex, filt=self.filt, offsetindex=offsetindex, MJDbinsize=None)
		print('Length of self.lc.t: ',len(self.lc.t))
		if len(self.lc.t) == 0:
			return(1)

		# add or replace mask column
		if 'Mask' in self.lc.t.columns:
			if self.verbose:
				print('Replacing existing Mask column...')
			for i in range(len(self.lc.t)):
				self.lc.t.loc[i,'Mask'] = 0
		else: 
			mask = np.zeros((len(self.lc.t)), dtype=int)
			self.lc.t = self.lc.t.assign(Mask=mask)

		# determine if using uncertainty cleanup
		uncert_apply = self.cfg.params['cleanlc']['uncertainty']['apply']
		if args.skip_uncert: uncert_apply = False
		if uncert_apply == True:
			print('Applying uncertainty cleanup...')
			self.flag_biguncertainty_duJy(SNindex,offsetindex)
		else:
			print('Skipping uncertainty cleanup...')

		# determine if using chi/N cleanup, and if type is dynamic or static
		chi_apply = self.cfg.params['cleanlc']['chi/N']['apply']
		if args.skip_chi: chi_apply = False
		if chi_apply == True:
			chi_type = self.cfg.params['cleanlc']['chi/N']['type']
			if chi_type == 'dynamic':
				print('Applying chi/N cleanlc type: %s ' % chi_type)
				self.flag_bigchi_dynamic(SNindex,offsetindex=offsetindex)
			elif chi_type == 'static':
				print('Applying chi/N cleanlc type: %s ' % chi_type)
				self.flag_bigchi_static(SNindex,offsetindex=offsetindex)
			else:
				print('Cleanup type error--type: ',chi_type)
				raise RuntimeError('chi/N cleanup type must be dynamic or static!')
		else:
			print('Skipping chi/N cleanup...')

		# save lc with additional mask column
		self.save_lc(SNindex=SNindex,filt=self.filt,overwrite=True,offsetindex=offsetindex)

if __name__ == '__main__':

	cleanuplc = cleanuplcclass()
	parser = cleanuplc.define_options()
	args = parser.parse_args()

	SNindexlist = cleanuplc.initialize(args)

	for SNindex in SNindexlist:
		cleanuplc.loadRADEClist(SNindex=SNindex, filt=cleanuplc.filt)
		for offsetindex in range(len(cleanuplc.RADECtable.t)):
			cleanuplc.cleanuplcloop(args,SNindex,offsetindex=offsetindex,filt=cleanuplc.filt)

	print('\n')
