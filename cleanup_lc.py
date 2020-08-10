#!/usr/bin/env python

from SNloop import SNloopclass
from statistics import median
import numpy as np

class cleanuplcclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)

	def uncertainty(self,SNindex,offsetindex):
		# load lc
		self.load_lc(SNindex, filt=self.filt, offsetindex=offsetindex, MJDbinsize=None)
		print('Length of self.lc.t: ',len(self.lc.t))
		if len(self.lc.t) == 0:
			return(1)
		# define vars
		# define indeces

	def dynamic(self,SNindex,offsetindex):
		# load lc
		self.load_lc(SNindex, filt=self.filt, offsetindex=offsetindex, MJDbinsize=None)
		print('Length of self.lc.t: ',len(self.lc.t))
		if len(self.lc.t) == 0:
			return(1)
		# define vars
		Nsigma = self.cfg.params['cleanlc']['chi/N']['Nsigma']
		chi_median = median(self.lc.t['chi/N'])
		chi_stddev = np.std(self.lc.t['chi/N'])
		print('Nsigma: %.1f, chi_median: %f, chi_stddev: %f' % (Nsigma, chi_median, chi_stddev))
		# define indeces
		a = int(chi_median+Nsigma*chi_stddev) # !! ROUND DOWN WITH INT OR ROUND UP?
		print('Removing all measurements with chi/N bigger than %i...' % a)
		a_indeces = np.where(self.lc.t['chi/N']<a)
		print('chi/N below %i: ' % a) # delete me
		print(self.lc.t['chi/N'][a_indeces]) # delete me

		Nremoved = len(self.lc.t['chi/N']) - len(self.lc.t['chi/N'][a_indeces])
		print('Length of cleaned up chi/N: ',len(self.lc.t['chi/N'][a_indeces]),', removed %i measurements' % Nremoved)

	def static(self,SNindex,offsetindex):
		# load lc
		self.load_lc(SNindex, filt=self.filt, offsetindex=offsetindex, MJDbinsize=None)
		print('Length of self.lc.t: ',len(self.lc.t))
		if len(self.lc.t) == 0:
			return(1)
		# define vars
		chi_max = self.cfg.params['cleanlc']['chi/N']['max_chi2norm']
		print('chi_max: ',chi_max)
		print('Removing all measurements with chi/N bigger than %i...' % chi_max)
		# define indeces
		a_indeces = np.where(self.lc.t['chi/N']<chi_max)
		print('chi/N below %i: ' % chi_max) # delete me
		print(self.lc.t['chi/N'][a_indeces]) # delete me

		Nremoved = len(self.lc.t['chi/N']) - len(self.lc.t['chi/N'][a_indeces])
		print('Length of cleaned up chi/N: %i, removed %i measurements' % (len(self.lc.t['chi/N'][a_indeces]), Nremoved))

	def cleanuplcloop(self,args,SNindex,offsetindex):
		uncert_apply = self.cfg.params['cleanlc']['uncertainty']['apply']
		if args.skip_uncert: uncert_apply = False
		if uncert_apply == True:
			print('Applying uncertainty cleanup...')
			self.uncertainty(SNindex,offsetindex)
		else:
			print('Skipping uncertainty cleanup...')

		chi_apply = self.cfg.params['cleanlc']['chi/N']['apply']
		if args.skip_chi: chi_apply = False
		if chi_apply == True:
			chi_type = self.cfg.params['cleanlc']['chi/N']['type']
			if chi_type == 'dynamic':
				print('Applying chi/N cleanlc type: %s ' % chi_type)
				self.dynamic(SNindex,offsetindex=offsetindex)
			elif chi_type == 'static':
				print('Applying chi/N cleanlc type: %s ' % chi_type)
				self.static(SNindex,offsetindex=offsetindex)
			else:
				print('Cleanup type error--type: ',chi_type)
				raise RuntimeError('chi/N cleanup type must be dynamic or static!')
		else:
			print('Skipping chi/N cleanup...')

if __name__ == '__main__':

	cleanuplc = cleanuplcclass()
	parser = cleanuplc.define_options()
	args = parser.parse_args()

	SNindexlist = cleanuplc.initialize(args)

	for SNindex in SNindexlist:
		cleanuplc.loadRADEClist(SNindex=SNindex, filt=cleanuplc.filt)
		for offsetindex in range(len(cleanuplc.RADECtable.t)):
			cleanuplc.cleanuplcloop(args,SNindex,offsetindex=offsetindex)
