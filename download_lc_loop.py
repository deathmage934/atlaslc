#!/usr/bin/env python

import numpy as np
import math
import sys,socket,os,re
import astropy.table as at
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import Angle
import matplotlib.pyplot as plt
import pylab as matlib
from datetime import datetime as dt
import argparse
from tools import yamlcfgclass
from tools import rmfile
from tools import RaInDeg
from tools import DecInDeg
import astrotable
from astrotable import astrotableclass
from download_atlas_lc import download_atlas_lc_class
from SNloop import SNloopclass
from download_atlas_lc_loop import downloadloopclass
import sigmacut

def dataPlot(x, y, dx=None, dy=None, sp=None, fmt='bo', ecolor='k', elinewidth=None, barsabove = False, capsize=1, logx=False, logy=False):
	if sp == None:
		sp = matlib.subplot(111)
	if dx is None and dy is None:
		if logy:
			if logx:
				plot, = sp.loglog(x, y, fmt)
			else:
				 plot, = sp.semilogy(x, y, fmt)
		elif logx:
			plot, = sp.semilogx(x, y, fmt)
		else:
			if barsabove:
				plot, dplot,dummy = sp.errorbar(x, y, fmt=fmt, capsize=capsize, barsabove=barsabove)
			else:
				plot, = sp.plot(x, y, fmt)
		return sp, plot, None
	else:
		if logy:
			sp.set_yscale("log", nonposx='clip')
		if logx:
			sp.set_xscale("log", nonposx='clip')
		plot, dplot,dummy = sp.errorbar(x, y, xerr=dx, yerr=dy, fmt=fmt, ecolor=ecolor, elinewidth=elinewidth, capsize=capsize, barsabove=barsabove)
		return sp, plot, dplot

class downloadlcloopclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)
		self.download_atlas_lc = download_atlas_lc_class()

	def downloadlc4SN(self, SNindex, lookbacktime_days=None, savelc=False, overwrite=False, fileformat=None, offsetindex=None):
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

		if self.verbose:
			print(self.RADECtable.t)

	def downloadoffsetlc4SN(self, SNindex,forcedphot_offset=False,lookbacktime_days=None, savelc=False, overwrite=False, fileformat=None,pattern=None):
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
				self.downloadlc4SN(SNindex,
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
			#print(self.RADECtable.t)
			self.RADECtable.formattable(formatMapping={'OffsetID':'3d','Ra':'.8f','Dec':'.8f','RaNew':'.8f','DecNew':'.8f','RaDistance':'.2f','DecOffset':'.2f','Ndet':'4d','Ndet_o':'4d','Ndet_c':'4d'})

			for i in range(len(self.RADECtable.t)):
				print(self.RADECtable.t['OffsetID', 'Ra', 'Dec'][i])
				self.downloadlc4SN(SNindex,
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

			if self.verbose: 
				print(self.RADECtable.t)
			if savelc:
				self.saveRADEClist(SNindex)
				self.saveRADEClist(SNindex,filt='c')
				self.saveRADEClist(SNindex,filt='o')

	def plot_lc(self,SNindex,sp=None,plotoffsetlcflag=False):
		print('Plotting SN lc and offsets...')

		MJD_offsetlc = []
		uJy_offsetlc = []
		duJy_offsetlc = []

		if sp is None:
			sp = matlib.subplot(111)

		self.loadRADEClist(SNindex)

		for i in range(len(self.RADECtable.t)-1,-1,-1):
			self.load_lc(SNindex, filt=self.filt, offsetindex=self.RADECtable.t['OffsetID'][i])
			print('Offset index: ',i)

			if i==0:
				(sp, plotSN, dplotSN)=dataPlot(self.lc.t['MJD'], self.lc.t['uJy'], dy=self.lc.t['duJy'],sp=sp)
				matlib.setp(plotSN,ms=5,color='r')
			else:
				MJD_offsetlc.extend(self.lc.t['MJD'])
				uJy_offsetlc.extend(self.lc.t['uJy'])
				duJy_offsetlc.extend(self.lc.t['duJy'])
				(sp, plotOffset, dplotOffset)=dataPlot(self.lc.t['MJD'],self.lc.t['uJy'],dy=self.lc.t['duJy'],sp=sp)
				matlib.setp(plotOffset,ms=5,color='b')

		if len(self.cfg.params['forcedphotpatterns']['circular']['radii'])==1:
			plt.legend(('SN', '%s %s" Offset' % (self.cfg.params['forcedphotpatterns']['circular']['n'], self.cfg.params['forcedphotpatterns']['circular']['radii'][0])))
		else:
			plt.legend(('SN', '%s %s" Offset and %s %s" Offset' % (self.cfg.params['forcedphotpatterns']['circular']['n'], self.cfg.params['forcedphotpatterns']['circular']['radii'][0],self.cfg.params['forcedphotpatterns']['circular']['n'], self.cfg.params['forcedphotpatterns']['circular']['radii'][1])))
		plt.axhline(linewidth=1,color='k')
		#plt.xlim(59000,59050)
		#plt.ylim(-5000,5000)
		plt.xlabel('MJD')
		plt.ylabel('uJy')

	def calcaveragelc(self, SNindex, MJDbinsize=None, offsetindex=None):
		self.load_lc(SNindex, filt=self.filt, offsetindex=offsetindex, MJDbinsize=None)
		self.averagelctable.formattable(formatMapping={'OffsetID':'3d','MJD':'.5f','uJy':'.2f','duJy':'.2f','stdev':'.2f','X2norm':'.3f'})
		print('Length of self.lc.t: ',len(self.lc.t))
		if len(self.lc.t) == 0:
			return(1)

		MJD = int(np.amin(self.lc.t['MJD']))
		MJDmax = np.amax(self.lc.t['MJD'])
		print('MJD: ',MJD)

		while MJD <= MJDmax:
			windowindeces = np.logical_and(self.lc.t['MJD']>=MJD, self.lc.t['MJD']<(MJD+MJDbinsize))
			uJy_windowdata = self.lc.t['uJy'][windowindeces]
			print(uJy_windowdata)
			duJy_windowdata = self.lc.t['duJy'][windowindeces]
			print(duJy_windowdata)
			MJD_windowdata = self.lc.t['MJD'][windowindeces]
			print(MJD_windowdata)

			calcaverage=sigmacut.calcaverageclass()
			calcaverage.calcaverage_sigmacutloop(uJy_windowdata,noise=duJy_windowdata,verbose=2,Nsigma=3.0,median_firstiteration=True,saveused=True)

			if calcaverage.Nused<=0:
				print('No data values. Skipping...')
				MJD += MJDbinsize
				continue

			print("mean:%f (uncertainty:%f) " % (calcaverage.mean,calcaverage.mean_err))
			lcaverageindex = len(self.averagelctable.t)
			self.averagelctable.t.add_row({'OffsetID':self.RADECtable.t['OffsetID'][offsetindex],'uJy':calcaverage.mean,'duJy':calcaverage.mean_err, 'stdev':calcaverage.stdev, 'X2norm':calcaverage.X2norm, 'Nused':calcaverage.Nused, 'Nclipped':calcaverage.Nskipped})
			
			print('lcaverageindex and calcaverage: ',lcaverageindex,calcaverage.use)
			mask=np.logical_not(calcaverage.use)
			calcaverage.calcaverage_sigmacutloop(MJD_windowdata,noise=duJy_windowdata,mask=mask,verbose=2,Nsigma=0.0,median_firstiteration=True,saveused=True)
			self.averagelctable.t['MJD'][lcaverageindex] = calcaverage.mean
			self.averagelctable.t['MJDNused'][lcaverageindex] = calcaverage.Nused
			self.averagelctable.t['MJDNskipped'][lcaverageindex] = calcaverage.Nskipped

			MJD += MJDbinsize

		print(self.averagelctable.t)
		return(0)
	
	def averagelcs(self, SNindex, MJDbinsize=None,fileformat=None,overwrite=True):
		for offsetindex in range(len(self.RADECtable.t)):
			if offsetindex>0:
				self.averagelctable.t.remove_rows(slice(0,len(self.averagelctable.t)))
			result = self.calcaveragelc(SNindex,offsetindex=offsetindex,MJDbinsize=MJDbinsize)

			filename = self.lcbasename(SNindex,filt=self.filt,offsetindex=offsetindex,MJDbinsize=MJDbinsize)+'.txt'
			if fileformat is None: 
				fileformat = self.cfg.params['output']['fileformat']

			if result == 0: 
				self.averagelctable.write(filename,format=fileformat,overwrite=True,verbose=(self.verbose>0))
			else:
				print('Removing file because length of self.lc.t is 0...')
				rmfile(filename)


	def define_options(self, parser=None, usage=None, conflict_handler='resolve'):
		parser = downloadloopclass.define_options(self,parser=parser, usage=usage, conflict_handler=conflict_handler)
		parser.add_argument('--forcedphot_offset', default=False)
		parser.add_argument('--pattern', default='circular',
							help=('offset pattern, defined in the config file, (default=%(default)s)'))
		parser.add_argument('--plot', default=False, help=('plot lcs'))
		parser.add_argument('--averagelc',default=False,help=('average lcs'))
		return(parser)

	def indexlistloop(self,args,indexlist):
		indexlist = self.initialize(args)
		for i in indexlist:
			self.downloadoffsetlc4SN(i,
									lookbacktime_days=args.lookbacktime_days,
									savelc=args.savelc,
									overwrite=args.overwrite,
									fileformat=args.fileformat,
									pattern=args.pattern,
									forcedphot_offset=args.forcedphot_offset)

			if args.plot:
				self.plot_lc(i)
				plotfilename = self.lcbasename(i)+'.png'
				print('Plot file name: ',plotfilename)
				plt.savefig(plotfilename)

			if args.averagelc:
				MJDbinsize = self.cfg.params['output']['MJDbinsize']
				if not(args.MJDbinsize is None): MJDbinsize = args.MJDbinsize
				self.loadRADEClist(i)
				self.averagelcs(i, MJDbinsize=MJDbinsize)

if __name__ == '__main__':

	downloadlc = downloadlcloopclass()
	parser = downloadlc.download_atlas_lc.define_optional_args()
	parser = downloadlc.define_options(parser=parser)
	args = parser.parse_args()

	indexlist = downloadlc.initialize(args)

	# get username from config file
	username = downloadlc.cfg.params['username']
	print('params user: ',username)
	# if config file username is false, overwrite with environment username
	if username is False:
		username = os.environ['USER']
		print('environ user: ',username)
	# if args.user in command, overwrite
	if not(args.user is None):
		username = args.user
		print('args user: ',username)
	print('ATLAS username: ', username)

	password=args.passwd
	print('ATLAS password length: ',len(password))

	downloadlc.download_atlas_lc.connect(args.atlasmachine,username,password)
	downloadlc.indexlistloop(args,indexlist)
