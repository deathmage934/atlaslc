#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 12:32:14 2020

@author: arest
"""

from SNloop import SNloopclass

class verifyMJDclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)

if __name__ == '__main__':

	verifyMJD = verifyMJDclass()
	parser = verifyMJD.define_options()
	args = parser.parse_args()

	SNindexlist = verifyMJD.initialize(args)

	for SNindex in SNindexlist:
		verifyMJD.loadRADEClist(SNindex, filt=verifyMJD.filt)
		verifyMJD.verifyMJD(SNindex,filt=verifyMJD.filt)
 