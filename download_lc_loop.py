#!/usr/bin/env python

import numpy as np
import math
import sys,socket,os,re
#import astropy.table as at
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import Angle
#import matplotlib.pyplot as plt
#import pylab as matlib
#from datetime import datetime as dt
import argparse
#from tools import yamlcfgclass
from tools import rmfile
from tools import RaInDeg
from tools import DecInDeg
#import astrotable
from astrotable import astrotableclass
from download_atlas_lc import download_atlas_lc_class
from SNloop import SNloopclass
from cleanup_lc import cleanuplcclass
from plot_lc import plotlcclass, dataPlot
from average_lc import averagelcclass
import sigmacut

class downloadlcloopclass(cleanuplcclass,plotlcclass,averagelcclass):
	def __init__(self):
		cleanuplcclass.__init__(self)
		plotlcclass.__init__(self)
		averagelcclass.__init__(self)
		self.download_atlas_lc = download_atlas_lc_class()

	def downloadlc(self, SNindex, lookbacktime_days=None, savelc=False, overwrite=False, fileformat=None, offsetindex=None):
		if offsetindex is None:
		   RA = self.t[SNindex]['ra']
		   Dec = self.t[SNindex]['dec']
		else:
			RA = self.RADECtable.t['RaNew'][offsetindex]
			Dec = self.RADECtable.t['DecNew'][offsetindex]

		self.download_atlas_lc.verbose = self.verbose
		self.download_atlas_lc.debug = self.debug
		self.download_atlas_lc.get_lc(RA,Dec,lookbacktime_days=lookbacktime_days)

		# read the lc into an ascii table
		self.lc.load_spacesep(self.download_atlas_lc.lcraw, formatMapping={'MJD':'%.6f','chi/N':'%.2f','x':'%.2f','y':'%.2f','m':'%.3f','dm':'%.3f'})

		# save the lc file with the automatic output filename
		if savelc:
			self.save_lc(SNindex,overwrite=overwrite,fileformat=fileformat, offsetindex=offsetindex)
			for filt in ['c','o']:
				filename = self.lcbasename(SNindex,filt=filt, offsetindex=offsetindex)+'.txt'
				if fileformat is None: fileformat = self.cfg.params['output']['fileformat']
				detections4filt=np.where(self.lc.t['F']==filt)
				if len(detections4filt[0]) is 0:
					print('Removing %s because nothing to save...' % filename)
					rmfile(filename)
				else: 
					self.lc.write(filename, indeces=detections4filt, format=fileformat, overwrite=overwrite, verbose=(self.verbose>0))

	def defineRADEClist(self,RA,Dec,pattern):
		if pattern=='circular':
			n = self.cfg.params['forcedphotpatterns']['circular']['n']
			radii = self.cfg.params['forcedphotpatterns']['circular']['radii']
			self.RADECtable.t.add_row({'OffsetID':0,'Ra':RaInDeg(RA),'Dec':DecInDeg(Dec),'RaNew':RaInDeg(RA),'DecNew':DecInDeg(Dec)})
			OffsetID=1
			for radius in radii:
				R = Angle(radius,u.arcsec)
				print('R = %f arcsec, %f deg, %f rad' % (R.arcsec,R.degree,R.radian))
				for i in range(n):
					angle = Angle(i*360.0/n, u.degree)
					if self.verbose>1:
						print('\nAngle: %f deg' % angle.degree)

					Dec_angle = Angle(Dec, u.degree)
					cosdec = math.cos(Dec_angle.radian)
					if self.verbose>1:
						print('Dec: %f deg' % Dec_angle.degree)

					RAdistance = Angle(R.degree*math.cos(angle.radian),u.degree)
					RAoffset = RAdistance*(1.0/cosdec)
					if self.verbose>1:
						print('RA distance: %f arcsec' % RAdistance.arcsec)
					if self.verbose>1:
						print('RA offset to be added to RA: %f arcsec' % RAoffset.arcsec)

					RAnew = RaInDeg(RA) + RAoffset.degree
					if self.verbose>1:
						print('RA new: %f deg' % RAnew)

					DECoffset = Angle(R.degree*math.sin(angle.radian),u.degree)
					if self.verbose>1:
						print('Dec offset: %f arcsec' % DECoffset.arcsec)
					DECnew = DecInDeg(Dec) + DECoffset.degree
					if self.verbose>1:
						print('Dec new: %f deg' % DECnew)

					if self.verbose:
						print('Angle: %.1f, new RA and Dec: %f, %f' % (angle.degree, RAnew, DECnew))

					index = self.RADECtable.t.add_row({'OffsetID':OffsetID,'Ra':RaInDeg(RA),'Dec':DecInDeg(Dec),'RaNew':RAnew,'DecNew':DECnew,'RaDistance':RAdistance.arcsec,'RaOffset':RAoffset.arcsec,'DecOffset':DECoffset.arcsec,'Radius':radius})
					OffsetID += 1
		elif pattern=='box':
			pass
		else:
			print("Pattern %s is not defined" % pattern)
			self.RADECtable.t.add_row({'OffsetID':0,'Ra':RaInDeg(RA),'Dec':DecInDeg(Dec),'RaNew':RaInDeg(RA),'DecNew':DecInDeg(Dec)})

		self.RADECtable.formattable(formatMapping={'OffsetID':'3d','Ra':'.8f','Dec':'.8f','RaNew':'.8f','DecNew':'.8f','RaDistance':'.2f','RaOffset':'.2f','DecOffset':'.2f','Ndet':'4d','Ndet_o':'4d','Ndet_c':'4d'})

		if self.verbose>1:
			print(self.RADECtable.t)

	def downloadoffsetlc(self, SNindex, forcedphot_offset=False, lookbacktime_days=None, savelc=False, overwrite=False, fileformat=None,pattern=None):
		print('Offset status: ',forcedphot_offset)
		if forcedphot_offset=='True':
			# add SN data in new row
			RA = self.t[SNindex]['ra']
			Dec = self.t[SNindex]['dec']
			if pattern is None: pattern=self.cfg.params['forcedphotpatterns'][pattern]
			self.defineRADEClist(RA,Dec,pattern)
			# print(self.RADECtable.t)
			self.RADECtable.formattable(formatMapping={'OffsetID':'3d','Ra':'.8f','Dec':'.8f','RaNew':'.8f','DecNew':'.8f','RaDistance':'.2f','DecOffset':'.2f','Ndet':'4d','Ndet_o':'4d','Ndet_c':'4d'})

			# add new row for each offset using data from RADECtable
			for i in range(len(self.RADECtable.t)):
				print(self.RADECtable.t['OffsetID', 'Ra', 'Dec'][i])
				self.downloadlc(SNindex,
							 lookbacktime_days=lookbacktime_days,
							 savelc=savelc,
							 overwrite=overwrite,
							 fileformat=fileformat,
							 offsetindex=i)
				
				print('Length of lc: ',len(self.lc.t))
				self.RADECtable.t['Ndet'][i]=len(self.lc.t)

				ofilt = np.where(self.lc.t['F']=='o')
				self.RADECtable.t['Ndet_o'][i]=len(ofilt[0])

				cfilt = np.where(self.lc.t['F']=='c')
				self.RADECtable.t['Ndet_c'][i]=len(cfilt[0])

			if self.verbose>1: 
				print(self.RADECtable.t)
			if savelc:
				self.saveRADEClist(SNindex)
				self.saveRADEClist(SNindex,filt='c')
				self.saveRADEClist(SNindex,filt='o')
		else:
			print('Skipping forcedphot offsets lc...')
			# only add SN data in new row without offsets
			RA = self.t[SNindex]['ra']
			Dec = self.t[SNindex]['dec']
			self.defineRADEClist(RA,Dec,pattern=None)
			self.RADECtable.formattable(formatMapping={'OffsetID':'3d','Ra':'.8f','Dec':'.8f','RaNew':'.8f','DecNew':'.8f','RaDistance':'.2f','DecOffset':'.2f','Ndet':'4d','Ndet_o':'4d','Ndet_c':'4d'})

			for i in range(len(self.RADECtable.t)):
				print(self.RADECtable.t['OffsetID', 'Ra', 'Dec'][i])
				self.downloadlc(SNindex,
							 lookbacktime_days=lookbacktime_days,
							 savelc=savelc,
							 overwrite=overwrite,
							 fileformat=fileformat,
							 offsetindex=i)
				
				print('Length of lc: ',len(self.lc.t))
				self.RADECtable.t['Ndet'][i]=len(self.lc.t)

				ofilt = np.where(self.lc.t['F']=='o')
				self.RADECtable.t['Ndet_o'][i]=len(ofilt[0])

				cfilt = np.where(self.lc.t['F']=='c')
				self.RADECtable.t['Ndet_c'][i]=len(cfilt[0])

			if self.verbose>1: 
				print(self.RADECtable.t)
			if savelc:
				self.saveRADEClist(SNindex)
				self.saveRADEClist(SNindex,filt='c')
				self.saveRADEClist(SNindex,filt='o')

if __name__ == '__main__':

	downloadlc = downloadlcloopclass()
	parser = downloadlc.download_atlas_lc.define_optional_args()
	parser = downloadlc.define_options(parser=parser)
	args = parser.parse_args()

	SNindexlist = downloadlc.initialize(args)

	username = downloadlc.cfg.params['username']
	if username is False:
		username = os.environ['USER']
	if not(args.user is None):
		username = args.user
	print('ATLAS username: ', username)
	password=args.passwd
	print('ATLAS password length: ',len(password))

	downloadlc.download_atlas_lc.connect(args.atlasmachine,username,password)

	for SNindex in SNindexlist:
		downloadlc.downloadoffsetlc(SNindex,
				lookbacktime_days=args.lookbacktime_days,
				savelc=args.savelc,
				overwrite=args.overwrite,
				fileformat=args.fileformat,
				pattern=args.pattern,
				forcedphot_offset=args.forcedphot_offset)
		for offsetindex in range(len(downloadlc.RADECtable.t)):
			downloadlc.cleanuplcloop(args,SNindex,offsetindex=offsetindex)
		if args.plot: 
			downloadlc.plotlcloop(args,SNindex)
		if args.averagelc: 
			downloadlc.averagelcloop(args,SNindex,offsetindex=offsetindex)
