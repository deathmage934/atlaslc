#!/usr/bin/env python

# S. Rest

from SNloop import SNloopclass
from statistics import median
import sigmacut
import sys
import numpy as np
import pandas as pd

class detectbumpsclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)

	def detectbumpsloop(self,SNindex):
		self.load_lc(SNindex,controlindex=0,filt=self.filt)
		if self.verbose==1:
			print('Length of self.lc.t: ',len(self.lc.t))
		if len(self.lc.t) == 0:
			return(1)

if __name__ == '__main__':

	detectbumps = detectbumpsclass()
	parser = detectbumps.define_options()
	args = parser.parse_args()

	SNindexlist = detectbumps.initialize(args)

	for SNindex in SNindexlist:
		detectbumps.detectbumpsloop(SNindex)