#!/usr/bin/env python

from SNloop import SNloopclass
import sys
import numpy as np
import pandas as pd

class getcleanlcclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)

	def getcleanlcfile(self,args,SNindex):
		# for original lcs
		self.load_lc(SNindex,filt=self.filt,controlindex=0)
		self.lctype = 'og'
		indices = self.getusableindices()
		#indices = self.getgoodindices()
		self.save_lc(SNindex=SNindex,indices=indices,filt=self.filt,overwrite=True,controlindex=0,addsuffix='.clean')

		# for average lcs
		self.load_lc(SNindex,filt=self.filt,controlindex=0,MJDbinsize=args.MJDbinsize)
		self.lctype = 'avg'
		indices = self.getusableindices()
		#indices = self.getgoodindices()
		self.save_lc(SNindex=SNindex,indices=indices,filt=self.filt,overwrite=True,controlindex=0,MJDbinsize=args.MJDbinsize,addsuffix='.clean')

if __name__ == '__main__':

	getcleanlc = getcleanlcclass()
	parser = getcleanlc.define_options()
	args = parser.parse_args()

	SNindexlist = getcleanlc.initialize(args)

	if args.filt is None:
		print('Looping through c and o filters...')
		filtlist = ['c','o']
	else:
		filtlist = [args.filt]
	
	for filt in filtlist:
			print('### FILTER SET: %s' % filt)
			getcleanlc.filt = filt
			for SNindex in SNindexlist:
				print('Getting and saving clean (usable) lc for ',getcleanlc.t.at[SNindex,'tnsname'],', index %i/%i' % (SNindex,len(getcleanlc.t)))
				getcleanlc.loadRADEClist(SNindex=SNindex,filt=getcleanlc.filt)
				getcleanlc.getcleanlcfile(args,SNindex)
			print('Finished with filter %s!' % filt)
	"""
	else:
		print('### FILTER SET: %s' % args.filt)
		getcleanlc.filt = args.filt
		for SNindex in SNindexlist:
			print('Getting and saving clean lc for ',getcleanlc.t.at[SNindex,'tnsname'],', index %i/%i' % (SNindex,len(getcleanlc.t)))
			getcleanlc.loadRADEClist(SNindex=SNindex,filt=getcleanlc.filt)
			getcleanlc.getcleanlcfile(args,SNindex)
	"""