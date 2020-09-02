#!/usr/bin/env python

import numpy as np
import math
import sys,socket,os,re
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import Angle
from astropy.table import Table
import argparse
from tools import rmfile
from tools import RaInDeg
from tools import DecInDeg
import pandas as pd
from astrotable import astrotableclass
from download_atlas_lc import download_atlas_lc_class
from SNloop import SNloopclass
from cleanup_lc import cleanuplcclass
from plot_lc import plotlcclass, dataPlot
from average_lc import averagelcclass
import sigmacut
import io

class downloadlcloopclass(cleanuplcclass,plotlcclass,averagelcclass):
	def __init__(self):
		cleanuplcclass.__init__(self)
		plotlcclass.__init__(self)
		averagelcclass.__init__(self)
		self.download_atlas_lc = download_atlas_lc_class()

	def downloadlc(self, SNindex, lookbacktime_days=None, savelc=False, overwrite=False, fileformat=None, offsetindex=None):
		if offsetindex is None:
			RA = self.t.loc[SNindex,'ra']
			Dec = self.t.loc[SNindex,'dec']
		else:
			RA = self.RADECtable.t.loc[offsetindex,'RaNew']
			Dec = self.RADECtable.t.loc[offsetindex,'DecNew']

		self.download_atlas_lc.verbose = self.verbose
		self.download_atlas_lc.debug = self.debug
		self.download_atlas_lc.get_lc(RA,Dec,lookbacktime_days=lookbacktime_days)

		# read the lc into a pandas table
		self.lc.t = pd.read_csv(io.StringIO('\n'.join(self.download_atlas_lc.lcraw)),delim_whitespace=True,skipinitialspace=True)

		# save the lc file with the output filename
		if savelc:
			self.save_lc(SNindex,overwrite=overwrite,fileformat=fileformat,offsetindex=offsetindex)
			for filt in ['c','o']:
				filename = self.lcbasename(SNindex,filt=filt, offsetindex=offsetindex)+'.txt'
				if fileformat is None: 
					fileformat = self.cfg.params['output']['fileformat']
				detections4filt=np.where(self.lc.t['F']==filt)
				if len(detections4filt[0]) is 0:
					print('Saving blank light curve: %s' % filename)
					self.lc.write(filename,index=False,indices=detections4filt,overwrite=True,verbose=False,columns=['MJD','m','dm','uJy','duJy','F','err','chi/N','RA','Dec','x','y','maj','min','phi','apfit','Sky','ZP','Obs','Mask'])
				else: 
					print('Saving light curve: %s' % filename)
					self.lc.write(filename, index=False, indices=detections4filt, overwrite=True, verbose=False)

	def defineRADEClist(self,RA,Dec,pattern=None):
		if not(pattern is None):
			# number of offsets and radii depends on pattern specified in precursor.cfg or in args
			if pattern=='circular':
				n = self.cfg.params['forcedphotpatterns']['circular']['n']
				radii = self.cfg.params['forcedphotpatterns']['circular']['radii']
			elif pattern=='box':
				n = 4
				radii = [self.cfg.params['forcedphotpatterns']['box']['sidelength']]
			elif pattern=='galaxy':
				n = self.cfg.params['forcedphotpatterns']['galaxy']['n']
				galRA = self.t.at[SNindex,'galRA']
				galDec = self.t.at[SNindex,'galDec']
				RA = self.t.at[SNindex,'ra']
				Dec = self.t.at[SNindex,'dec']
				rdist = 1 # !! FIX ASAP
				radii = [rdist,rdist+5]
			else:
				raise RuntimeError("Pattern %s is not defined" % pattern)

			# sets up RADECtable, fills in OffsetID, Ra, Dec, RaNew, DecNew for (n*len(radii)) offsets
			df = pd.DataFrame([[0,RaInDeg(RA),DecInDeg(Dec),RaInDeg(RA),DecInDeg(Dec),0,0,0,0,0,0,0]], columns=['OffsetID','Ra','Dec','RaNew','DecNew','RaDistance','RaOffset','DecOffset','Radius','Ndet','Ndet_c','Ndet_o'])
			self.RADECtable.t = self.RADECtable.t.append(df, ignore_index=True)
			OffsetID=1
			for radius in radii:
				R = Angle(radius,u.arcsec)
				print('R = %f arcsec, %f deg, %f rad' % (R.arcsec,R.degree,R.radian))
				for i in range(n):
					angle = Angle(i*360.0/n, u.degree)
					Dec_angle = Angle(Dec, u.degree)
					cosdec = math.cos(Dec_angle.radian)
					RAdistance = Angle(R.degree*math.cos(angle.radian),u.degree)
					RAoffset = RAdistance*(1.0/cosdec)
					RAnew = RaInDeg(RA) + RAoffset.degree
					DECoffset = Angle(R.degree*math.sin(angle.radian),u.degree)
					DECnew = DecInDeg(Dec) + DECoffset.degree
					if self.verbose>1:
						print('\nAngle: %f deg' % angle.degree)
						print('Dec: %f deg' % Dec_angle.degree)
						print('RA distance: %f arcsec' % RAdistance.arcsec)
						print('RA offset to be added to RA: %f arcsec' % RAoffset.arcsec)
						print('RA new: %f deg' % RAnew)
						print('Dec offset: %f arcsec' % DECoffset.arcsec)
						print('Dec new: %f deg' % DECnew)
					if self.verbose:
						print('Angle: %.1f, new RA and Dec: %f, %f' % (angle.degree, RAnew, DECnew))
					df = pd.DataFrame([[OffsetID,RaInDeg(RA),DecInDeg(Dec),RAnew, DECnew,RAdistance.arcsec,RAoffset.arcsec,DECoffset.arcsec,radius,0,0,0]],columns=['OffsetID','Ra','Dec','RaNew','DecNew','RaDistance','RaOffset','DecOffset','Radius','Ndet','Ndet_c','Ndet_o'])
					self.RADECtable.t = self.RADECtable.t.append(df,ignore_index=True)
					OffsetID += 1
		else:
			df = pd.DataFrame([[0,RaInDeg(RA),DecInDeg(Dec),RaInDeg(RA),DecInDeg(Dec),0,0,0,0,0,0,0]], columns=['OffsetID','Ra','Dec','RaNew','DecNew','RaDistance','RaOffset','DecOffset','Radius','Ndet','Ndet_c','Ndet_o'])
			self.RADECtable.t = self.RADECtable.t.append(df, ignore_index=True)
		#print(self.RADECtable.write(index=True,overwrite=False))

	def downloadoffsetlc(self, SNindex, forcedphot_offset=False, lookbacktime_days=None, savelc=False, overwrite=False, fileformat=None,pattern=None):
		print('Offset status: ',forcedphot_offset)
		if forcedphot_offset =='True':
			# add SN data in new row
			RA = self.t.at[SNindex,'ra']
			Dec = self.t.at[SNindex,'dec']
			if pattern is None: pattern=self.cfg.params['forcedphotpatterns'][pattern]
			self.defineRADEClist(RA,Dec,pattern=pattern)
			
			# add new row for each offset using data from RADECtable
			for i in range(len(self.RADECtable.t)):
				if self.verbose>1:
					print(self.RADECtable.write(indices=i, columns=['OffsetID','Ra','Dec']))
				self.downloadlc(SNindex,
							 lookbacktime_days=lookbacktime_days,
							 savelc=savelc,
							 overwrite=overwrite,
							 fileformat=fileformat,
							 offsetindex=i)
				
				print('Length of lc: ',len(self.lc.t))
				self.RADECtable.t.loc[i,'Ndet']=len(self.lc.t)

				ofilt = np.where(self.lc.t['F']=='o')
				self.RADECtable.t.loc[i,'Ndet_o']=len(ofilt[0])

				cfilt = np.where(self.lc.t['F']=='c')
				self.RADECtable.t.loc[i,'Ndet_c']=len(cfilt[0])
			print(self.RADECtable.write(index=True,overwrite=False))
			if savelc:
				self.saveRADEClist(SNindex,filt='c')
				self.saveRADEClist(SNindex,filt='o')
		else:
			print('Skipping forcedphot offsets lc...')
			# only add SN data in new row without offsets
			RA = self.t.at[SNindex, 'ra']
			Dec = self.t.at[SNindex, 'dec']
			self.defineRADEClist(RA,Dec,pattern=None)
			#self.RADECtable.formattable(formatMapping={'OffsetID':'3d','Ra':'.8f','Dec':'.8f','RaNew':'.8f','DecNew':'.8f','RaDistance':'.2f','DecOffset':'.2f','Ndet':'4d','Ndet_o':'4d','Ndet_c':'4d'})

			for i in range(len(self.RADECtable.t)):
				if self.verbose:
					print(self.RADECtable.t.loc[i,['OffsetID', 'RaNew', 'DecNew']])
				self.downloadlc(SNindex,
							 lookbacktime_days=lookbacktime_days,
							 savelc=savelc,
							 overwrite=overwrite,
							 fileformat=fileformat,
							 offsetindex=i)
				
				print('Length of lc: ',len(self.lc.t))
				self.RADECtable.t.loc[i,'Ndet']=len(self.lc.t)

				ofilt = np.where(self.lc.t['F']=='o')
				self.RADECtable.t.loc[i,'Ndet_o']=len(ofilt[0])

				cfilt = np.where(self.lc.t['F']=='c')
				self.RADECtable.t.loc[i,'Ndet_c']=len(cfilt[0])

			if self.verbose>1: 
				print(self.RADECtable)
			if savelc:
				#self.saveRADEClist(SNindex)
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
			downloadlc.cleanuplcloop(args,SNindex,offsetindex=offsetindex,filt=downloadlc.filt)
			if args.averagelc: 
				downloadlc.averagelcloop(args,SNindex,offsetindex=offsetindex)
		if args.plot: 
			downloadlc.plotlcloop(args,SNindex)
