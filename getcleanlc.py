#!/usr/bin/env python

from SNloop import SNloopclass
import sys
import numpy as np
import pandas as pd

class getcleanlcclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)

	def getcleanlcfile(self,args,SNindex):
		self.load_lc(SNindex,filt=self.filt,controlindex=0)
		indices = self.getusableindices()
		if self.lctype == 'avg':
			self.save_lc(SNindex=SNindex,indices=indices,filt=self.filt,overwrite=True,controlindex=0,MJDbinsize=args.MJDbinsize,addsuffix='.clean')
		else:
			self.save_lc(SNindex=SNindex,indices=indices,filt=self.filt,overwrite=True,controlindex=0,addsuffix='.clean')

		self.load_lc(SNindex,filt=self.filt,controlindex=0,MJDbinsize=args.MJDbinsize)
		indices = self.getusableindices()
		if self.lctype == 'avg':
			self.save_lc(SNindex=SNindex,indices=indices,filt=self.filt,overwrite=True,controlindex=0,MJDbinsize=args.MJDbinsize,addsuffix='.clean')
		else:
			self.save_lc(SNindex=SNindex,indices=indices,filt=self.filt,overwrite=True,controlindex=0,addsuffix='.clean')

if __name__ == '__main__':

	getcleanlc = getcleanlcclass()
	parser = getcleanlc.define_options()
	args = parser.parse_args()

	SNindexlist = getcleanlc.initialize(args)

	if args.filt is None:
		print('Looping through c and o filters...')
		for filt in ['o','c']:
			print('### FILTER SET: %s' % filt)
			getcleanlc.filt = filt
			for SNindex in SNindexlist:
				print('Getting and saving clean lc for ',getcleanlc.t.at[SNindex,'tnsname'],', index %i/%i' % (SNindex,len(getcleanlc.t)))
				getcleanlc.loadRADEClist(SNindex=SNindex,filt=getcleanlc.filt)
				getcleanlc.getcleanlcfile(args,SNindex)
			print('Finished with filter %s!' % filt)
	else:
		print('### FILTER SET: %s' % args.filt)
		getcleanlc.filt = args.filt
		for SNindex in SNindexlist:
			print('Getting and saving clean lc for ',getcleanlc.t.at[SNindex,'tnsname'],', index %i/%i' % (SNindex,len(getcleanlc.t)))
			getcleanlc.loadRADEClist(SNindex=SNindex,filt=getcleanlc.filt)
			getcleanlc.getcleanlcfile(args,SNindex)