#!/usr/bin/env python

# Code adapted from Q. Wang by S. Rest

from SNloop import SNloopclass
from download_lc_loop import downloadlcloopclass
from uploadTransientData import upload, runDBcommand, DBOps
from autoadd import autoaddclass
from pdastro import pdastroclass, pdastrostatsclass, AnotB, AandB
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
from copy import deepcopy
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

		# tables and table/list names
		self.YSEtable = pdastrostatsclass()
		self.RADECtable = pdastrostatsclass()
		self.averagelctable = pdastrostatsclass()
		self.TNSnamelist = None
		self.TNSlistfilename = None
		self.TNSlistfile = pdastrostatsclass(columns=['TNSname','ra','dec'])
		self.lc = pdastrostatsclass(hexcols=['Mask'])
		
		self.api = False
		self.verbose = 0

	def YSE_list(self):
		all_cand = pd.read_csv(self.cfg.params['upltoyse']['yse_list'])
		all_cand = all_cand.drop_duplicates(subset='name')
		df = pd.DataFrame()
		df['Name'] = all_cand['name'] 
		df['RA'] = all_cand['transient_RA']
		df['Dec'] = all_cand['transient_Dec']
		df['Disc date'] = all_cand['disc_date']
		print(df)
		sys.exit(0)
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

	def saveyselc(self,TNSname,controlindex,filt=None,indices=None,overwrite=True,MJDbinsize=None):
		#oindex = '%03d' % controlindex
		if filt is None:
			filt = self.filt
		filename = self.yselcfilename(TNSname,controlindex,filt,MJDbinsize)
		if not(MJDbinsize is None):
			self.averagelctable.write(filename,indices=indices,overwrite=overwrite,verbose=True)
		else:
			self.lc.write(filename,indices=indices,overwrite=overwrite,verbose=True)
		return(0)

	def loadyselc(self,TNSname,controlindex,filt=None,MJDbinsize=None):
		if filt is None:
			filt = self.filt
		filename = self.yselcfilename(TNSname,controlindex,filt,MJDbinsize)
		if not(MJDbinsize is None):
			self.averagelctable.load_spacesep(filename,delim_whitespace=True)
		else:
			self.lc.load_spacesep(filename,delim_whitespace=True)
		return(0)

	def atlas2yse(self,TNSname,outname,ra,dec,atlas_data_file,filt):
		t = ascii.read(atlas_data_file)
		
		# different output names for regular lcs vs. averaged lcs
		if 'days' in atlas_data_file: # if averaged
			outname = atlas_data_file[:-7]+'.yse.csv'
		else:
			outname = atlas_data_file[:-9]+'.%s.yse.csv' % filt
		
		filter_fict = {'o':'orange-ATLAS', 'c':'cyan-ATLAS'}
		with open(outname, 'w+') as f:
			f.write('SNID: '+TNSname+' \nRA: '+str(ra)+'     \nDECL: '+str(dec)+' \n \nVARLIST:  MJD  FLT  FLUXCAL   FLUXCALERR MAG     MAGERR DQ \n')
			for row in t:
				"""
				flt = filter_fict[row['F']]
				if row['m']>0:
					flux = 10**(-0.4*(row['m']-27.5)) # mag = (-1/0.4)*log(row[self.flux_colname])+27.5
					fluxerr= row['dm']*10**(-0.4*(row['m']-27.5))
				else:
					flux = -10**(-0.4*(-row['m']-27.5)) # mag = (1/0.4)*log(-row[self.flux_colname])+27.5
					fluxerr= row['dm']*10**(-0.4*(-row['m']-27.5))
				mag = row['m']
				magerr = row['dm']
				"""
				f.write('OBS: '+str(row['MJD'])+' '+filter_fict[row['F']]+' '+str(row[self.flux_colname])+' '+str(row[self.dflux_colname])+' '+str(row['m'])+' '+str(row['dm'])+' 0 \n')
		print("Converted ATLAS lc to YSE format: %s" % outname)
		return(outname)

	def uploadtoyse(self,filename):
		os.system('python %s/uploadTransientData.py -e -s %s/settings.ini -i %s --instrument ACAM1 --fluxzpt 23.9' % (self.sourcedir,self.sourcedir,filename))

	def saveRADECtable(self,args,TNSname,filt):
		RADECtablefilename = '%s/%s/%s/%s.RADECtable.txt' % (self.outrootdir,self.outsubdir,TNSname,TNSname)
		print('Saving RADECtable: %s' % RADECtablefilename)
		self.RADECtable.write(RADECtablefilename,overwrite=args.overwrite,verbose=True)
		return(0)

	def loadRADECtable(self,args,TNSname):
		# get RADECtable from already existing file
		RADECtablefilename = '%s/%s/%s/%s.RADECtable.txt' % (self.outrootdir,self.outsubdir,TNSname,TNSname)
		print('Loading RADECtable: %s' % RADECtablefilename)
		self.RADECtable.load_spacesep(RADECtablefilename, delim_whitespace=True)
		#print(self.RADECtable.write())
		return(0)

	def defineRADECtable(self,args,RA,Dec,pattern=None):
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
						print('Autosearch not functional yet!!')
						sys.exit(0)
						"""
						results = self.autosearch(RA.degree, Dec.degree, 20)
						print('Close bright objects found: \n',results)
						cbRA = 
						cbDec = 
						self.t.at[SNindex,'closebrightRA'] = Angle(RaInDeg(cbRA),u.degree)
						self.t.at[SNindex,'closebrightDec'] = Angle(RaInDeg(cbDec),u.degree)
						"""
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
					#print('Minimum distance from SN to control LC: ',mindist)
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
			print(self.RADECtable.write(index=True,overwrite=args.overwrite))
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

		# split the lc file into 2 separate files by filter
		for filt in ['c','o']:
			filename = self.yselcfilename(TNSname,controlindex,filt)
			fileformat = self.cfg.params['output']['fileformat']
			detections4filt=np.where(self.lc.t['F']==filt)
			newindices = AandB(indices,detections4filt)
			if len(detections4filt[0]) is 0:
				print('Saving blank light curve: %s' % filename)
				self.lc.write(filename,index=False,indices=newindices,overwrite=args.overwrite,verbose=False,columns=['MJD','m','dm',self.flux_colname,self.dflux_colname,'F','err','chi/N','RA','Dec','x','y','maj','min','phi','apfit','Sky','ZP','Obs','Mask'])
			else: 
				print('Saving light curve: %s' % filename)
				self.lc.write(filename, index=False,indices=newindices,overwrite=args.overwrite,verbose=False)

	def downloadYSEcontrollc(self,args,TNSname,ra,dec,pattern=None,lookbacktime_days=None):
		print('Offset status: ',args.forcedphot_offset)
		if args.forcedphot_offset == 'True':
			self.defineRADECtable(args,ra,dec,pattern=pattern)

			for i in range(len(self.RADECtable.t)):
				print(self.RADECtable.write(indices=i, columns=['ControlID', 'Ra', 'Dec']))
				self.downloadyselc(args,ra,dec,i,lookbacktime_days=lookbacktime_days)
				print('Length of lc: ',len(self.lc.t))

				self.RADECtable.t.loc[i,'Ndet']=len(self.lc.t)
				ofilt = np.where(self.lc.t['F']=='o')
				self.RADECtable.t.loc[i,'Ndet_o']=len(ofilt[0])
				cfilt = np.where(self.lc.t['F']=='c')
				self.RADECtable.t.loc[i,'Ndet_c']=len(cfilt[0])
			
			self.saveRADECtable(args,TNSname,'c')
			self.saveRADECtable(args,TNSname,'o')
		else:
			print('Skipping forcedphot offsets lc...')
			self.defineRADECtable(args,ra,dec,pattern=None)

			print(self.RADECtable.write(index=True,overwrite=False))
			for i in range(len(self.RADECtable.t)):
				self.downloadyselc(args,ra,dec,i,lookbacktime_days=lookbacktime_days)
				print('Length of lc: ',len(self.lc.t))

				self.RADECtable.t.loc[i,'Ndet']=len(self.lc.t)
				ofilt = np.where(self.lc.t['F']=='o')
				self.RADECtable.t.loc[i,'Ndet_o']=len(ofilt[0])
				cfilt = np.where(self.lc.t['F']=='c')
				self.RADECtable.t.loc[i,'Ndet_c']=len(cfilt[0])
			
			self.saveRADECtable(args,TNSname,'c')
			self.saveRADECtable(args,TNSname,'o')

	def averageyselc(self,args,TNSname,filt,MJDbinsize=1.0):
		self.loadRADECtable(args,TNSname)
		self.loadyselc(TNSname,0,filt)

		# clear averagelctable and set columns
		self.averagelctable.t = self.averagelctable.t.iloc[0:0]
		self.averagelctable.t['MJD'] = pd.Series([], dtype=np.float64)
		self.averagelctable.t['MJDbin'] = pd.Series([], dtype=np.float64)
		self.averagelctable.t[self.flux_colname] = pd.Series([], dtype=np.float64)
		self.averagelctable.t[self.dflux_colname] = pd.Series([], dtype=np.float64)
		self.averagelctable.t['m'] = pd.Series([], dtype=np.float64)
		self.averagelctable.t['dm'] = pd.Series([], dtype=np.float64)
		self.averagelctable.t['stdev'] = pd.Series([], dtype=np.float64)
		self.averagelctable.t['X2norm'] = pd.Series([], dtype=np.float64)
		self.averagelctable.t['Nclipped'] = pd.Series([], dtype=np.int64)
		self.averagelctable.t['Nused'] = pd.Series([], dtype=np.int64)

		# add masks to mask column by cleaning lc - x2norm and uncertainty masks
		self.lc.t['Mask'] = 0
		# drop old columns
		dropcols=[]
		if 'Noffsetlc' in self.lc.t.columns: dropcols.append('Noffsetlc')
		for col in self.lc.t.columns:
			if re.search('^c\d_',col): dropcols.append(col)
			if re.search('^o\d_',col): dropcols.append(col)
		if len(dropcols)>0: self.lc.t.drop(columns=dropcols,inplace=True)
		# masks
		self.c0_PSF_uncertainty_cut(self.cfg.params['cleanlc']['cut0']['N_dflux_max'])
		self.c0_PSF_X2norm_cut(self.cfg.params['cleanlc']['cut0']['PSF_X2norm_max'])
		self.saveyselc(TNSname,0,filt=filt,overwrite=args.overwrite)

		# average light curve add points to new averaged df
		MJD = int(np.amin(self.lc.t['MJD']))
		MJDmax = int(np.amax(self.lc.t['MJD']))+1
		while MJD <= MJDmax:
			# get measurements for MJD range
			mjd_ix = self.lc.ix_inrange(colnames=['MJD'],lowlim=MJD,uplim=MJD+MJDbinsize,exclude_uplim=True)
			if self.cfg.params['upltoyse']['use_cut0']:
				# get measurements without x2norm and uncertainty masks
				mjd_ix = self.lc.ix_unmasked('Mask',maskval=self.flag_c0_X2norm|self.flag_c0_uncertainty,indices=mjd_ix)

			if len(mjd_ix)>0:
				# add row to averagelc table
				df = {'MJDbin':MJD+0.5*MJDbinsize,'F':filt}
				lcaverageindex = self.averagelctable.newrow(df)

			if len(mjd_ix)==0: 
				if self.verbose>1: print('No data in MJD range = 0, skipping MJD range...')
				MJD += MJDbinsize
				continue

			# sigmacut indices
			self.lc.calcaverage_sigmacutloop(self.flux_colname,noisecol=self.dflux_colname,indices=mjd_ix,verbose=1,Nsigma=3.0,median_firstiteration=True)
			fluxstatparams = deepcopy(self.lc.statparams)
			if self.verbose>1: print('Nclip: {}, Ngood: {}, X2norm: {}'.format(fluxstatparams['Nclip'],fluxstatparams['Ngood'],fluxstatparams['X2norm']))

			# get average mjd
			self.lc.calcaverage_sigmacutloop('MJD',noisecol=self.dflux_colname,indices=fluxstatparams['ix_good'],verbose=1,Nsigma=0,median_firstiteration=False)
				
			# add row to averagelc table
			df = {'MJD':self.lc.statparams['mean'],self.flux_colname:fluxstatparams['mean'],self.dflux_colname:fluxstatparams['mean_err'],
				  'stdev':fluxstatparams['stdev'],'X2norm':fluxstatparams['X2norm'],'Nclipped':fluxstatparams['Nclip'],'Nused':fluxstatparams['Ngood']}
			self.averagelctable.add2row(lcaverageindex,df)

			MJD += MJDbinsize

		# convert flux to magnitude
		self.averagelctable.flux2mag(self.flux_colname,self.dflux_colname,'m','dm',zpt=23.9,upperlim_Nsigma=3)

		self.averagelctable.write()

		self.saveyselc(TNSname,0,filt=filt,overwrite=args.overwrite,MJDbinsize=MJDbinsize)

	def uploadloop(self,args,TNSname,overwrite=True,skipdownload=False):
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

		"""
		print('Overwrite status: ',args.overwrite)
		if args.overwrite is False:
			for filt in ['c','o']:
				filename = self.yselcfilename(TNSname,0,filt)
				if os.path.exists(filename):
					print("Data for %s with filter %s already exists" % (TNSname,filt))
					self.existflag = True
				else:
					print("Found no data for %s with filter %s, downloading full lc..." % (TNSname,filt))
					self.existflag = False
		"""

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
		#self.downloadYSEcontrollc(args,TNSname,ra,dec,pattern=pattern,lookbacktime_days=lookbacktime_days)
		
        if not skipdownload:
            self.downloadYSEcontrollc(args,TNSname,ra,dec,pattern=pattern,lookbacktime_days=lookbacktime_days)

		# upload to YSE-PZ
		for filt in ['c','o']:
			print('### FILTER SET: ',filt)
			filename = self.yselcfilename(TNSname,0,filt)
			outname = self.atlas2yse(TNSname,filename,ra,dec,filename,filt)
			self.uploadtoyse(outname)

			if args.averagelc:
				filename = self.yselcfilename(TNSname,0,filt,args.MJDbinsize)
				self.averageyselc(args,TNSname,filt,MJDbinsize=args.MJDbinsize)
				outname = self.atlas2yse(TNSname,filename,ra,dec,filename,filt)
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
    parser.add_argument('--averagelc', default=False, help=('average lcs'))
    parser.add_argument('--skipdownload', default=False, help=('askip downloading'))
    parser.add_argument('-m','--MJDbinsize', default=1.0, help=('specify MJD bin size'),type=float)

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
    upltoyse.verbose = args.verbose

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

    # set up tables
    upltoyse.flux_colname = upltoyse.cfg.params['flux_colname']
    upltoyse.dflux_colname = upltoyse.cfg.params['dflux_colname']
    upltoyse.RADECtable = pdastroclass(columns=['ControlID','PatternID','Ra','Dec','RaOffset','DecOffset','Radius','Ndet','Ndet_c','Ndet_o'])
    upltoyse.RADECtable.default_formatters = {'ControlID':'{:3d}'.format,'PatternID':'{:2d}'.format,'Ra':'{:.8f}'.format,'Dec':'{:.8f}'.format,'RaOffset':'{:.2f}'.format,'DecOffset':'{:.2f}'.format,'Radius':'{:.2f}'.format,'Ndet':'{:4d}'.format,'Ndet_c':'{:4d}'.format,'Ndet_o':'{:4d}'.format}
    upltoyse.averagelctable = pdastroclass(columns=['MJDbin','MJD','F',upltoyse.flux_colname,upltoyse.dflux_colname,'m','dm','stdev','X2norm','Nused','Nclipped'])
    # the following line caused a formatting error and I didn't want to figure out why... this is a problem for future sofia:
    #upltoyse.averagelctable.default_formatters = {'MJDbin':'{:.1f}'.format,'MJD':'{:.5f}'.format,upltoyse.flux_colname:'{:.2f}'.format,upltoyse.dflux_colname:'{:.2f}'.format,'m':'{:.3f}'.format,'dm':'{:.3f}'.format,'stdev':'{:.2f}'.format,'X2norm':'{:.3f}'.format,'Nused':'{:4d}'.format,'Nclipped':'{:4d}'.format,'MJDNused':'{:4d}'.format,'MJDNskipped':'{:4d}'.format}

    # api
    upltoyse.api = upltoyse.cfg.params['api']
    if args.api: 
        upltoyse.api = True

    #for TNSname in upltoyse.TNSnamelist:
    for index in range(0,len(upltoyse.TNSnamelist)):
        TNSname = upltoyse.TNSnamelist[index]
        print("\nUploading and/or downloading data for %s, TNSnamelist index %d/%d" % (TNSname,index+1,len(upltoyse.TNSnamelist)))
        upltoyse.uploadloop(args,TNSname,overwrite=True,skipdownload=args.skipdownload)

    if not(args.tnslistfilename is None):
        upltoyse.TNSlistfile.write(upltoyse.TNSlistfilename,overwrite=True)

