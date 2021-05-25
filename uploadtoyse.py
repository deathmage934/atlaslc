#!/usr/bin/env python

# Code adapted from Q. Wang by S. Rest

from SNloop import SNloopclass
from download_lc_loop import downloadlcloopclass
from uploadTransientData import upload, runDBcommand, DBOps
from autoadd import autoaddclass
from pdastro import pdastroclass, AandB
from download_atlas_lc import download_atlas_lc_class
from tools import RaInDeg
from tools import DecInDeg

import pandas as pd
import numpy as np
import sigmacut
import sys,socket,os,re,io
import optparse
import configparser
import argparse
import math
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.table import Table

"""
for command input: uploadtoyse.py -t 2020lse --user USERNAME --passwd 'PASSWORD'
for TNSlistfile input: uploadtoyse.py -n tnslist.txt --user USERNAME --passwd 'PASSWORD'
for YSE list input: uploadtoyse.py --user USERNAME --passwd 'PASSWORD'

you must specify source directory, output root directory (has tnslist.txt if you want to use it), and output subdirectory (has all the downloaded data).
you can specify these in the commandline using --sourcedir, --outrootdir, and --outsubdir
OR
1. you can set the source directory and the output root directory in atlaslc.sourceme
2. then set the output subdirectory in precursor.cfg next to 'yse_outsubdir' (currently set to default 'ysetest')
"""

class uploadtoyseclass(downloadlcloopclass,autoaddclass):
	def __init__(self):
		downloadlcloopclass.__init__(self)
		autoaddclass.__init__(self)
		#upload.__init__(self)
		self.download_atlas_lc = download_atlas_lc_class()

		# important vars
		self.sourcedir = None
		self.outrootdir = None
		self.outsubdir = None
		self.flux_colname = None
		self.dflux_colname = None
		#self.existflag = False

		# tables and table/list names
		self.YSEtable = pdastroclass()
		self.RADECtable = pdastroclass()
		self.averagelctable = pdastroclass()
		self.TNSnamelist = None
		self.TNSlistfilename = None
		self.TNSlistfile = pdastroclass(columns=['TNSname','ra','dec'])
		self.lc = pdastroclass()
		self.api = False

	def YSE_list(self):
		all_cand = pd.read_csv('https://ziggy.ucolick.org/yse/explorer/147/download?format=csv')
		all_cand = all_cand.drop_duplicates(subset='name')
		df = pd.DataFrame()
		df['Name'] = all_cand['name'] 
		df['RA'] = all_cand['transient_RA']
		df['Dec'] = all_cand['transient_Dec']
		df['Disc date'] = all_cand['disc_date']
		return df

	def checkTNSlistfile(self,TNSname):
		# search for TNSname in TNSlistfile
		onTNSlistfile = False
		for name in self.TNSlistfile.t['TNSname']:
			if name == TNSname:
				onTNSlistfile = True
		return(onTNSlistfile)

	def yselcfilename(self,TNSname,controlindex,filt,MJDbinsize=None):
		oindex = '%03d' % controlindex # why not just do this on the actual line?? like %03d.%s.lc.txt?? too afraid to fix lol
		SNID = TNSname
		if not(MJDbinsize is None):
			filename = '%s/%s/%s/%s_i%s.%s.%.2fdays.lc.txt' % (self.outrootdir,self.outsubdir,SNID,SNID,oindex,filt,MJDbinsize)
		else:
			filename = '%s/%s/%s/%s_i%s.%s.lc.txt' % (self.outrootdir,self.outsubdir,SNID,SNID,oindex,filt)
		return(filename)

	def saveyselc(self,TNSname,controlindex,filt=None,indices=None,overwrite=False,MJDbinsize=None):
		#oindex = '%03d' % controlindex
		if filt is None:
			filt = self.filt
		filename = self.yselcfilename(TNSname,controlindex,filt,MJDbinsize)
		self.lc.write(filename,indices=indices,overwrite=overwrite,verbose=True)
		return(0)

	def loadyselc(self,TNSname,controlindex,filt=None,MJDbinsize=None):
		if filt is None:
			filt = self.filt
		filename = self.yselcfilename(TNSname,controlindex,filt,MJDbinsize)
		self.lc.load_spacesep(filename,delim_whitespace=True)
		return(0)

	def atlas2yse(self,TNSname,outname,ra,dec,atlas_data_file,filt):
	    t = ascii.read(atlas_data_file)
	    outname = atlas_data_file[:-9]+'%s.yse.csv' % filt
	    filter_fict = {'o':'orange-ATLAS', 'c':'cyan-ATLAS'}
	   
	    with open(outname, 'w+') as f:
	        f.write('SNID: '+TNSname+' \nRA: '+str(ra)+'     \nDECL: '+str(dec)+' \n \nVARLIST:  MJD  FLT  FLUXCAL   FLUXCALERR MAG     MAGERR DQ \n')
	        for k in t:
	        	#if args.cleanlc:
	        		#if not(k['Mask']>0): #if k['Mask']==0:
	        			# copy and paste
	        	#else:
	        		# copy and paste
	            flt = filter_fict[k['F']]
	            if k['m']>0:
	                flux = 10**(-0.4*(k['m']-27.5))
	                fluxerr= k['dm']*10**(-0.4*(k['m']-27.5))
	            else:
	                flux = -10**(-0.4*(-k['m']-27.5))
	                fluxerr= k['dm']*10**(-0.4*(-k['m']-27.5))
	            mag = k['m']
	            magerr = k['dm']
	            f.write('OBS: ' + str(k['MJD']) +' '+ flt+' '+ str(flux)+ ' '+str(fluxerr)+' '+ str(mag)+' '+ str(magerr)+' 0 \n')
	    print("Converted ATLAS lc to YSE format: %s" % outname)
	    return(outname)

	def uploadtoyse(self,filename):
		os.system('python %s/uploadTransientData.py -e -s %s/settings.ini -i %s --instrument ACAM1 --fluxzpt 27.5' % (self.sourcedir,self.sourcedir,filename))

	def saveRADECtable(self,TNSname,filt):
		RADECtablefilename = '%s/%s/%s/%s.RADECtable.txt' % (self.outrootdir,self.outsubdir,TNSname,TNSname)
		print('Saving RADECtable: %s' % RADECtablefilename)
		self.RADECtable.write(RADECtablefilename,overwrite=True,verbose=True)
		return(0)

	def loadRADECtable(self,TNSname):
		# get RADECtable from already existing file
		RADECtablefilename = '%s/%s/%s/%s.RADECtable.txt' % (self.outrootdir,self.outsubdir,TNSname,TNSname)
		print('Loading RADECtable: %s' % RADECtablefilename)
		self.RADECtable.load_spacesep(RADECtablefilename, delim_whitespace=True)
		#print(self.RADECtable.write())
		return(0)

	def defineRADECtable(self,RA,Dec,pattern=None):
		self.RADECtable.t = self.RADECtable.t[0:0]
		if not(pattern is None):
			pattern_list = pattern
			print('Pattern(s) set to ',pattern_list)
			ControlID = 1
			foundflag = False
			RA = Angle(RaInDeg(RA),u.degree)
			Dec = Angle(DecInDeg(Dec),u.degree)
			
			for pattern in pattern_list: 
				# number of offsets and radii depends on pattern specified in precursor.cfg or in args
				if pattern=='circle':
					PatternID = 1
					n = self.cfg.params['forcedphotpatterns']['circle']['n']
					radii = self.cfg.params['forcedphotpatterns']['circle']['radii']
				elif pattern=='box':
					PatternID = 2
					n = 4
					radii = [self.cfg.params['forcedphotpatterns']['box']['sidelength']]
				elif pattern=='closebright':
					PatternID = 3
					mindist = self.cfg.params['forcedphotpatterns']['closebright']['mindist']
					n = self.cfg.params['forcedphotpatterns']['closebright']['n']

					# query panstarrs for closest bright object and add to snlist.txt
					if self.cfg.params['forcedphotpatterns']['closebright']['autosearch'] is True:
						results = self.autosearch(RA.degree, Dec.degree, 20)
						print('Close bright objects found: \n',results)
						sys.exit(0)
						cbRA = '_' # FIX
						cbDec = '_' # FIX
						self.t.at[SNindex,'closebrightRA'] = Angle(RaInDeg(cbRA),u.degree)
						self.t.at[SNindex,'closebrightDec'] = Angle(RaInDeg(cbDec),u.degree)
					# use coordinates listed in snlist.txt
					else:
						cbRA = Angle(RaInDeg(self.t.at[SNindex,'closebrightRA']),u.degree)
						cbDec = Angle(DecInDeg(self.t.at[SNindex,'closebrightDec']),u.degree)

					# radius is distance between SN and bright object
					c1 = SkyCoord(RA, Dec, frame='fk5')
					c2 = SkyCoord(cbRA, cbDec, frame='fk5')
					sep = c1.separation(c2)
					r1 = sep.arcsecond
					radii = [r1]

					if self.verbose:
						print('Minimum distance from SN to control LC: ',mindist)
				else:
					raise RuntimeError("Pattern %s is not defined" % pattern)

				# sets up RADECtable, fills in ControlID, Ra, Dec, RaNew, DecNew for (n*len(radii)) offsets
				if foundflag is False:
					foundflag = True
					df = pd.DataFrame([[0,0,RA.degree,Dec.degree,0,0,0,0,0,0]], columns=['ControlID','PatternID','Ra','Dec','RaOffset','DecOffset','Radius','Ndet','Ndet_c','Ndet_o'])
					self.RADECtable.t = self.RADECtable.t.append(df, ignore_index=True)
				for radius in radii:
					R = Angle(radius,u.arcsec)

					# define center coordinates
					if pattern=='closebright':
						RAcenter = Angle(cbRA.degree,u.degree)
						Deccenter = Angle(cbDec.degree,u.degree)
					elif (pattern=='circle') or (pattern=='box'):	
						RAcenter = Angle(RA.degree,u.degree)
						Deccenter = Angle(Dec.degree,u.degree)
					else:
						raise RuntimeError('Pattern must be circle, box, or closebright!')

					#print('CENTER: ',RAcenter.degree,Deccenter.degree) # delete me
					#print('R = %f arcsec or %f deg or %f rad' % (R.arcsec,R.degree,R.radian)) # delete me

					for i in range(n):
						angle = Angle(i*360.0/n, u.degree)
						RAdistance = Angle(R.degree*math.cos(angle.radian),u.degree)
						RAoffset = Angle(RAdistance.degree*(1.0/math.cos(Deccenter.radian)),u.degree)
						RAnew = Angle(RAcenter.degree+RAoffset.degree,u.degree)
						DECoffset = Angle(R.degree*math.sin(angle.radian),u.degree)
						DECnew = Angle(Deccenter.degree+DECoffset.degree,u.degree)

						# check to see if offset location is within mindist arcsec from SN
						if pattern=='closebright':
							c1 = SkyCoord(RA, Dec, frame='fk5')
							#print(RAnew,DECnew) # delete me
							c2 = SkyCoord(RAnew, DECnew, frame='fk5')
							offset_sep = c1.separation(c2)
							offset_sep = offset_sep.arcsecond
							if offset_sep < mindist:
								print('Control LC with ControlID %d too close to SN with distance of %f arcsec, skipping...' % (ControlID, offset_sep))
								continue
						df = pd.DataFrame([[ControlID,PatternID,RAnew.degree,DECnew.degree,RAdistance.arcsec,DECoffset.arcsec,radius,0,0,0]],columns=['ControlID','PatternID','Ra','Dec','RaOffset','DecOffset','Radius','Ndet','Ndet_c','Ndet_o'])
						self.RADECtable.t = self.RADECtable.t.append(df,ignore_index=True)

						if self.verbose>1:
							print('#\nAngle: %f deg' % angle.degree)
							print('#RA center: %f deg' % RAcenter.degree)
							print('#Dec center: %f deg' % Deccenter.degree)
							print('#RA distance: %f arcsec' % RAdistance.arcsec)
							print('#RA offset to be added to RA: %f arcsec' % RAoffset.arcsec)
							print('#Dec offset: %f arcsec' % DECoffset.arcsec)
						if self.verbose:
							print('#Angle: %.1f, new RA and Dec: %f, %f' % (angle.degree, RAnew.degree, DECnew.degree))
						ControlID += 1
			print(self.RADECtable.write(index=True,overwrite=False))
		else:
			RA = Angle(RaInDeg(RA),u.degree)
			Dec = Angle(DecInDeg(Dec),u.degree)
			df = pd.DataFrame([[0,0,RA.degree,Dec.degree,0,0,0,0,0,0]], columns=['ControlID','PatternID','Ra','Dec','RaOffset','DecOffset','Radius','Ndet','Ndet_c','Ndet_o'])
			self.RADECtable.t = self.RADECtable.t.append(df, ignore_index=True)

	def downloadyselc(self,args,ra,dec,controlindex,lookbacktime_days=None):
		self.download_atlas_lc.verbose = 1
		if self.api is True:
			print('Connecting to API...')
			token_header = self.download_atlas_lc.connect_atlas(args.user,args.passwd)
			print('TOKEN HEADER: ',token_header)
			self.lc.t = self.download_atlas_lc.get_result(RaInDeg(ra), DecInDeg(dec), token_header, lookbacktime_days=lookbacktime_days)
		else:
			print('Connecting to SSH...')
			self.download_atlas_lc.connect(args.atlasmachine,args.user,args.passwd)
			self.download_atlas_lc.get_lc(ra,dec,lookbacktime_days=lookbacktime_days)
			self.lc.t = pd.read_csv(io.StringIO('\n'.join(self.download_atlas_lc.lcraw)),delim_whitespace=True,skipinitialspace=True)
		
		mask = np.zeros((len(self.lc.t)), dtype=int)
		self.lc.t = self.lc.t.assign(Mask=mask)

		self.lc.t = self.lc.t.sort_values(by=['MJD'],ignore_index=True)
		indices = self.lc.ix_remove_null(colnames='uJy')

		# save the lc file with the output filename
		self.saveyselc(TNSname,controlindex,indices=indices)

		# split the lc file into 2 separate files by filter
		for filt in ['c','o']:
			filename = self.yselcfilename(TNSname,controlindex,filt)
			fileformat = self.cfg.params['output']['fileformat']
			detections4filt=np.where(self.lc.t['F']==filt)
			newindices = AandB(indices,detections4filt)
			if len(detections4filt[0]) is 0:
				print('Saving blank light curve: %s' % filename)
				self.lc.write(filename,index=False,indices=newindices,overwrite=True,verbose=False,columns=['MJD','m','dm',self.flux_colname,self.dflux_colname,'F','err','chi/N','RA','Dec','x','y','maj','min','phi','apfit','Sky','ZP','Obs','Mask'])
			else: 
				print('Saving light curve: %s' % filename)
				self.lc.write(filename, index=False, indices=newindices, overwrite=True, verbose=False)

	def downloadYSEcontrollc(self,args,TNSname,ra,dec,pattern=None,lookbacktime_days=None):
		print('Offset status: ',args.forcedphot_offset)
		if args.forcedphot_offset == 'True':
			self.defineRADECtable(ra,dec,pattern=pattern)

			for i in range(len(self.RADECtable.t)):
				print(self.RADECtable.write(indices=i, columns=['ControlID', 'Ra', 'Dec']))
				self.downloadyselc(args,ra,dec,i,lookbacktime_days=lookbacktime_days)
				print('Length of lc: ',len(self.lc.t))

				self.RADECtable.t.loc[i,'Ndet']=len(self.lc.t)
				ofilt = np.where(self.lc.t['F']=='o')
				self.RADECtable.t.loc[i,'Ndet_o']=len(ofilt[0])
				cfilt = np.where(self.lc.t['F']=='c')
				self.RADECtable.t.loc[i,'Ndet_c']=len(cfilt[0])
			
			self.saveRADECtable(TNSname,'c')
			self.saveRADECtable(TNSname,'o')
		else:
			print('Skipping forcedphot offsets lc...')
			self.defineRADECtable(ra,dec,pattern=None)

			print(self.RADECtable.write(index=True,overwrite=False))
			for i in range(len(self.RADECtable.t)):
				self.downloadyselc(args,ra,dec,i,lookbacktime_days=lookbacktime_days)
				print('Length of lc: ',len(self.lc.t))

				self.RADECtable.t.loc[i,'Ndet']=len(self.lc.t)
				ofilt = np.where(self.lc.t['F']=='o')
				self.RADECtable.t.loc[i,'Ndet_o']=len(ofilt[0])
				cfilt = np.where(self.lc.t['F']=='c')
				self.RADECtable.t.loc[i,'Ndet_c']=len(cfilt[0])
			
			self.saveRADECtable(TNSname,'c')
			self.saveRADECtable(TNSname,'o')

	def averageyselc(self,args,TNSname):
		self.loadRADECtable(TNSname)
		for filt in ['c','o']:
			self.loadyselc(TNSname,0,filt)
			
			# add masks to mask column

			# average light curve and add points to new averaged df

	"""
	def cleanupYSEcontrollc(self,args,TNSname):
		self.loadRADECtable(TNSname)
		for controlindex in range(len(self.RADECtable.t)):
			for filt in ['c','o']:
				# load lc
				#oindex = '%03d' % controlindex
				self.loadyselc(TNSname,controlindex,filt)
				print('Loaded %s light curve! Length of self.lc.t: ' % filt,len(self.lc.t))
				if len(self.lc.t) == 0: return(1)

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
					self.flag_biguncertainty_duJy(TNSname,controlindex)
				else:
					print('Skipping uncertainty cleanup...')

				# determine if using chi/N cleanup, and if type is dynamic or static
				chi_apply = self.cfg.params['cleanlc']['chi/N']['apply']
				if args.skip_chi: chi_apply = False
				if chi_apply == True:
					chi_type = self.cfg.params['cleanlc']['chi/N']['type']
					if chi_type == 'dynamic':
						print('Applying chi/N cleanlc type: %s ' % chi_type)
						self.flag_bigchi_dynamic(TNSname,controlindex)
					elif chi_type == 'static':
						print('Applying chi/N cleanlc type: %s ' % chi_type)
						self.flag_bigchi_static(TNSname,controlindex)
					else:
						print('Cleanup type error--type: ',chi_type)
						raise RuntimeError('chi/N cleanup type must be dynamic or static!')
				else:
					print('Skipping chi/N cleanup...')

				# save lc with additional mask column
				self.saveyselc(TNSname,controlindex,filt,overwrite=True)
	"""
	"""
	def getcleanlc(self,args):
		flags = self.cfg.params['upltoyse']['flags']
		print('Setting indices using flags: %x' % flags)
		
		mask=np.bitwise_and(self.lc.t['Mask'], flags)
		cuts_indices = np.where(mask==0)
		cuts_indices = cuts_indices[0]
		bad_data = np.where(mask!=0)
		bad_data = bad_data[0]

		datacut = len(self.lc.t)-len(self.lc.t.loc[cuts_indices])
		print('Length original lc: ',len(self.lc.t),', length cleaned lc: ',len(self.lc.t.loc[cuts_indices]),', data points cut: ',datacut)

		lc_m = self.lc.t.loc[cuts_indices, 'm']
		lc_dm = self.lc.t.loc[cuts_indices,'dm']
		lc_mjd = self.lc.t.loc[cuts_indices,'MJD']
		#return(lc_m, lc_dm, lc_mjd, cuts_indices)
	"""
	"""
	def averageYSElc(self,args,TNSname):
		self.loadRADECtable(TNSname)

		# decide whether or not to make cuts on data to average
		makecuts_apply = self.cfg.params['averagelc']['makecuts']
		if not(args.avg_makecuts) is None:
			if args.avg_makecuts is True:
				makecuts_apply = True
			else:
				makecuts_apply = False

		for controlindex in range(len(self.RADECtable.t)):
			MJDbinsize = self.cfg.params['output']['MJDbinsize']
			if not(args.MJDbinsize is None): MJDbinsize = args.MJDbinsize
			
			# load lc
			self.loadyselc(TNSname,controlindex,filt)
			print('Loaded %s light curve! Length of self.lc.t: ' % filt,len(self.lc.t))
			if len(self.lc.t) == 0: return(1)

			# AHHHHHHHHHHHHH 
			if makecuts_apply == True:
				if not('Mask' in self.lc.t.columns):
					raise RuntimeError('No "Mask" column exists! Please run "cleanup_lc.py %s" beforehand.' % self.t.at[SNindex,'tnsname'])
				lc_uJy, lc_duJy, lc_MJD, cuts_indices, bad_data = self.makecuts_indices(SNindex, controlindex, 'averagelc')
				print('Calculating average_lc table for offset index %d' % controlindex)
				self.averagelcs(args, SNindex, controlindex, lc_uJy, lc_duJy, lc_MJD, MJDbinsize=MJDbinsize, cuts_indices=cuts_indices)
			# data set to lc table
			else:
				print('Skipping makecuts using mask column...')
				lc_uJy = self.lc.t[self.flux_colname]
				lc_duJy = self.lc.t[self.dflux_colname]
				lc_MJD = self.lc.t['MJD']
				print('Calculating average_lc table for offset index %d' % controlindex) 
				self.averagelcs(args, SNindex, controlindex, lc_uJy, lc_duJy, lc_MJD, MJDbinsize=MJDbinsize)
			'''
	"""
	def uploadloop(self,args,TNSname,overwrite=False):
		# GET RA AND DEC
		if args.tnsnamelist:
			# get ra and dec automatically
			ra, dec = self.getradec(TNSname)
		elif args.tnslistfilename:
			onTNSlistfile = self.checkTNSlistfile(TNSname)
			if onTNSlistfile is False:
				# get ra and dec automatically, then append to TNSlistfile
				ra, dec = self.getradec(TNSname)
				df = pd.DataFrame([[TNSname,ra,dec]], columns=['TNSname','RA','Dec'])
				self.TNSlistfile.t = self.TNSlistfile.t.append(df, ignore_index=True)
			else:
				# get ra and dec from TNSlistfile
				index = self.TNSlistfile.ix_equal('TNSname',val=TNSname)
				if len(index)>0:
					ra = self.TNSlistfile.t.at[index[0],'RA']
					dec = self.TNSlistfile.t.at[index[0],'Dec']
				else:
					raise RuntimeError('Something went wrong: TNSname does not exist!')
		else:
			# get ra and dec from yse table
			index = self.YSEtable.ix_equal('Name',val=TNSname)
			if len(index)>0:
				ra = self.YSEtable.t.at[index[0],'RA']
				dec = self.YSEtable.t.at[index[0],'Dec']
			else:
				raise RuntimeError('Something went wrong: TNSname does not exist!')

		#self.existflag = False
		if args.overwrite is False:
			for filt in ['c','o']:
				filename = self.yselcfilename(TNSname,0,filt)
				if os.path.exists(filename):
					print("Data for %s with filter %s already exists" % (TNSname,filt))
					self.existflag = True
				else:
					print("Found no data for %s with filter %s, downloading full lc..." % (TNSname,filt))
					self.existflag = False

		# set offset pattern
		if not(args.pattern is None):
			pattern = args.pattern
		else:
			pattern = self.cfg.params['forcedphotpatterns']['patterns_to_use']

		# download lc and, if forcedphot_offset, offset lcs
		if args.lookbacktime_days:
			lookbacktime_days = args.lookbacktime_days
		else:
			lookbacktime_days = 60
		self.downloadYSEcontrollc(args,TNSname,ra,dec,pattern=pattern,lookbacktime_days=lookbacktime_days)
		
		#if args.averagelc is True:
			#self.averageYSElc(args,TNSname)
		
		# upload to YSE-PZ
		for filt in ['c','o']:
			filename = self.yselcfilename(TNSname,0,filt)
			outname = self.atlas2yse(TNSname,filename,ra,dec,filename,filt)
			self.uploadtoyse(outname)

			if args.averagelc:
				filename = self.yselcfilename(TNSname,0,filt,args.MJDbinsize)
				outname = self.atlas2yse(TNSname,outname,ra,dec,atlas_data_file,filt)
				self.uploadtoyse(outname)

if __name__ == '__main__':

	upltoyse = uploadtoyseclass()

	# add arguments
	parser = argparse.ArgumentParser(conflict_handler='resolve')
	parser.add_argument('-t','--tnsnamelist', default=None, nargs='+', help='name of transients to download and upload')
	parser.add_argument('-n','--tnslistfilename', default=None, help='address of file containing TNS names, ra, and dec')
	parser.add_argument('-o','--overwrite',default=True,help='overwrite existing files')
	parser.add_argument('--sourcedir', default=None, help='source code directory')
	parser.add_argument('--outrootdir', default=None, help='output root directory')
	parser.add_argument('--outsubdir', default=None, help='output subdirectory')
	parser.add_argument('--api', default=False, help=('use API instead of SSH to get light curves from ATLAS'))
	parser.add_argument('--forcedphot_offset', default=False, help=("download offsets (settings in config file)"))
	parser.add_argument('--pattern', choices=['circle','box','closebright'], help=('offset pattern, defined in the config file; options are circle, box, or closebright'))
	#parser.add_argument('--plotlc', default=False, help=('plot lcs')) # to add
	#parser.add_argument('--cleanlc', default=False, help=('upload only clean data to yse')) # in progress
	parser.add_argument('--averagelc', default=False, help=('average lcs')) # to add
	#parser.add_argument('--skip_uncert', default=False, help=('skip cleanup lcs using uncertainties'))
	#parser.add_argument('--skip_chi', default=False, help=('skip cleanup lcs using chi/N'))
	parser.add_argument('-m','--MJDbinsize', default=None, help=('specify MJD bin size'),type=float) # to add

	# add config file and atlaslc arguments
	cfgfile = upltoyse.defineoptions()
	parser = upltoyse.download_atlas_lc.define_optional_args(parser=parser)
	args = parser.parse_args()

	# set up source code directory, output root directory, and output subdirectory
	# get either from arguments, atlaslc.sourceme (sets up directories automatically), or precursor.cfg (config file)
	upltoyse.loadcfgfile(cfgfile)
	if args.sourcedir: 
		upltoyse.sourcedir = args.sourcedir
	else: 
		upltoyse.sourcedir = os.environ['ATLASLC_SOURCEDIR']
	if args.outrootdir: 
		upltoyse.outrootdir = args.outrootdir
	else: 
		upltoyse.outrootdir = upltoyse.cfg.params['output']['outrootdir']
	if args.outsubdir: 
		upltoyse.outsubdir = args.outsubdir
	else: 
		upltoyse.outsubdir = upltoyse.cfg.params['output']['yse_outsubdir']

	# GET TNSNAMELIST
	if not(args.tnsnamelist is None):
		# TNSnamelist set to objects put in command line
		upltoyse.TNSnamelist = args.tnsnamelist
		print("TNSnamelist from command: ",upltoyse.TNSnamelist)
	elif not(args.tnslistfilename is None):
		upltoyse.TNSlistfilename = '%s/%s' % (upltoyse.outrootdir,args.tnslistfilename)
		# check if TNSlistfilename exists; if not, make TNSlistfile
		if os.path.exists(upltoyse.TNSlistfilename):
			upltoyse.TNSlistfile.load_spacesep(upltoyse.TNSlistfilename, delim_whitespace=True)
		else: 
			raise RuntimeError("%s does not exist! Please create a file with the headers 'TNSname', 'RA', and 'Dec', then try again.")
		# TNSnamelist set to objects in TNSlistfilename
		upltoyse.TNSnamelist = upltoyse.TNSlistfile.t['TNSname'].values
		print("TNSnamelist from TNSlistfile: ",upltoyse.TNSnamelist)
		print("TNSlistfilename: ",upltoyse.TNSlistfilename)
	else:
		upltoyse.YSEtable.t = upltoyse.YSE_list()
		# TNSnamelist set to ojects in YSE list
		upltoyse.TNSnamelist = upltoyse.YSEtable.t['Name'].values
		print("TNSnamelist from YSE: ",upltoyse.TNSnamelist) # change me

	# SET UP TABLES
	upltoyse.flux_colname = upltoyse.cfg.params['flux_colname']
	upltoyse.dflux_colname = upltoyse.cfg.params['dflux_colname']
	upltoyse.RADECtable = pdastroclass(columns=['ControlID','PatternID','Ra','Dec','RaOffset','DecOffset','Radius','Ndet','Ndet_c','Ndet_o'])
	upltoyse.averagelctable = pdastroclass(columns=['ControlID','MJD',upltoyse.flux_colname,upltoyse.dflux_colname,'stdev','X2norm','Nused','Nclipped','MJDNused','MJDNskipped'])
	upltoyse.RADECtable.default_formatters = {'ControlID':'{:3d}'.format,'PatternID':'{:2d}'.format,'Ra':'{:.8f}'.format,'Dec':'{:.8f}'.format,'RaOffset':'{:.2f}'.format,'DecOffset':'{:.2f}'.format,'Radius':'{:.2f}'.format,'Ndet':'{:4d}'.format,'Ndet_c':'{:4d}'.format,'Ndet_o':'{:4d}'.format}
	upltoyse.averagelctable.default_formatters = {'ControlID':'{:3d}'.format,'MJD':'{:.5f}'.format,upltoyse.flux_colname:'{:.2f}'.format,upltoyse.dflux_colname:'{:.2f}'.format,'stdev':'{:.2f}'.format,'X2norm':'{:.3f}'.format,'Nused':'{:4d}'.format,'Nclipped':'{:4d}'.format,'MJDNused':'{:4d}'.format,'MJDNskipped':'{:4d}'.format}

	# api
	upltoyse.api = upltoyse.cfg.params['api']
	if args.api: 
		upltoyse.api = True

	#for TNSname in upltoyse.TNSnamelist:
	for index in range(0,len(upltoyse.TNSnamelist)):
		TNSname = upltoyse.TNSnamelist[index]
		print("\nUploading and/or downloading data for %s, TNSnamelist index %d/%d" % (TNSname,index,len(upltoyse.TNSnamelist)))
		upltoyse.uploadloop(args,TNSname,overwrite=True)

	if not(args.tnslistfilename is None):
		upltoyse.TNSlistfile.write(upltoyse.TNSlistfilename,overwrite=True)

