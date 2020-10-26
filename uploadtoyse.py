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
from astropy.io import ascii

# for command input: uploadtoyse.py -t 2020lse --user USERNAME --passwd 'PASSWORD'
# for TNSlistfile input: uploadtoyse.py -n tnslist.txt --user USERNAME --passwd 'PASSWORD'
# for YSE list input: uploadtoyse.py --user USERNAME --passwd 'PASSWORD'

# you must specify source directory, output root directory (has tnslist.txt if you want to use it), and output subdirectory (has all the downloaded data)
# you can specify these in the commandline using --sourcedir, --outrootdir, and --outsubdir
# OR
# 1. you can set the source directory and the output root directory in atlaslc.sourceme
# 2. you can set the output subdirectory in precursor.cfg next to 'yse_outsubdir' (currently set to default 'ysetest')

class uploadtoyseclass(downloadlcloopclass,autoaddclass):
	def __init__(self):
		downloadlcloopclass.__init__(self)
		autoaddclass.__init__(self)
		#upload.__init__(self)
		self.download_atlas_lc = download_atlas_lc_class()

		self.sourcedir = None
		self.outrootdir = None
		self.outsubdir = None
		self.flux_colname = None
		self.dflux_colname = None

		self.YSEtable = pdastroclass()
		self.TNSnamelist = None
		self.TNSlistfilename = None
		self.TNSlistfile = pdastroclass(columns=['TNSname','ra','dec'])
		self.yselc = pdastroclass()

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

	def yselcfilename(self,TNSname,filt):
		SNID = TNSname
		filt = filt
		if not(filt is None):
			filename = '%s/%s/%s/%s.%s.lc.txt' % (self.outrootdir,self.outsubdir,SNID,SNID,filt)
		else:
			filename = '%s/%s/%s/%s.lc.txt' % (self.outrootdir,self.outsubdir,SNID,SNID)
		return(filename)

	def downloadyselc(self,args,ra,dec,lookbacktime_days=None,fileformat=None):
		self.download_atlas_lc.verbose = 1
		self.download_atlas_lc.connect(args.atlasmachine,args.user,args.passwd)
		self.download_atlas_lc.get_lc(ra,dec,lookbacktime_days=lookbacktime_days)

		# read the lc into a pandas table, sort by MJD, and remove nans
		self.yselc.t = pd.read_csv(io.StringIO('\n'.join(self.download_atlas_lc.lcraw)),delim_whitespace=True,skipinitialspace=True)
		mask = np.zeros((len(self.yselc.t)), dtype=int)
		self.yselc.t = self.yselc.t.assign(Mask=mask)
		self.yselc.t = self.yselc.t.sort_values(by=['MJD'],ignore_index=True)
		indices = self.yselc.ix_remove_null(colnames='uJy')

		# save the lc file with the output filename
		self.saveyselc(TNSname,indices=indices)

		# split the lc file into 2 separate files by filter
		for filt in ['c','o']:
			filename = self.yselcfilename(TNSname=TNSname,filt=filt)
			if fileformat is None: 
				fileformat = self.cfg.params['output']['fileformat']
			detections4filt=np.where(self.yselc.t['F']==filt)
			newindices = AandB(indices,detections4filt)
			if len(detections4filt[0]) is 0:
				print('Saving blank light curve: %s' % filename)
				self.yselc.write(filename,index=False,indices=newindices,overwrite=True,verbose=False,columns=['MJD','m','dm',self.flux_colname,self.dflux_colname,'F','err','chi/N','RA','Dec','x','y','maj','min','phi','apfit','Sky','ZP','Obs','Mask'])
			else: 
				print('Saving light curve: %s' % filename)
				self.yselc.write(filename, index=False, indices=newindices, overwrite=True, verbose=False)

	def saveyselc(self,TNSname,filt=None,indices=None,overwrite=False):
		# write table and save lc as file
		filename = self.yselcfilename(TNSname,filt)
		self.yselc.write(filename,indices=indices,overwrite=overwrite,verbose=True)
		return(0)

	def atlas2yse(self, TNSname, outname, ra, dec, atlas_data_file, filt):
	    t = ascii.read(atlas_data_file)
	    outname = atlas_data_file[:-9]+'%s.yse.csv' % filt
	    filter_fict = {'o':'orange-ATLAS', 'c':'cyan-ATLAS'}
	   
	    with open(outname, 'w+') as f:
	        f.write('SNID: '+TNSname+' \nRA: '+str(ra)+'     \nDECL: '+str(dec)+' \n \nVARLIST:  MJD  FLT  FLUXCAL   FLUXCALERR MAG     MAGERR DQ \n')
	        for k in t:
	            flt = filter_fict[k['F']]
	            if k[self.flux_colname]>0:
	                flux = 10**(-0.4*(k[self.flux_colname]-27.5))
	                fluxerr= k[self.dflux_colname]*10**(-0.4*(k[self.flux_colname]-27.5))
	            else:
	                flux = -10**(-0.4*(-k[self.flux_colname]-27.5))
	                fluxerr= k[self.dflux_colname]*10**(-0.4*(-k[self.flux_colname]-27.5))
	            mag = k[self.flux_colname]
	            magerr = k[self.dflux_colname]
	            f.write('OBS: ' + str(k['MJD']) +' '+ flt+' '+ str(flux)+ ' '+str(fluxerr)+' '+ str(mag)+' '+ str(magerr)+' 0 \n')
	    print("Converted ATLAS lc to YSE format: %s" % outname)
	    return(outname)

	def uploadtoyse(self,filename):
		os.system('python %s/uploadTransientData.py -e -s %s/settings.ini -i %s --instrument ACAM1 --fluxzpt 27.5' % (self.sourcedir,self.sourcedir,filename))

	def uploadloop(self,args,TNSname,overwrite=False):
		self.flux_colname = self.cfg.params['flux_colname']
		self.dflux_colname = self.cfg.params['dflux_colname']

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
				#ra = self.TNSlistfile.t.at['RA',TNSname]
				#dec = self.TNSlistfile.t.at['Dec',TNSname]
		else:
			# get ra and dec from yse table
			index = self.YSEtable.ix_equal('Name',val=TNSname)
			if len(index)>0:
				ra = self.YSEtable.t.at[index[0],'RA']
				dec = self.YSEtable.t.at[index[0],'Dec']
			else:
				raise RuntimeError('Something went wrong: TNSname does not exist!')

		ra = RaInDeg(ra)
		dec = DecInDeg(dec)

		# check if sn already has lc downloaded; if not, download
		existflag = True
		for filt in ['c','o']:
			filename = self.yselcfilename(TNSname,filt)
			if os.path.exists(filename):
				print("Data for %s with filter %s already exists" % (TNSname,filt))
			else:
				print("Found no data for %s with filter %s, downloading full lc..." % (TNSname,filt))
				existflag = False
		
		if existflag == False:
			if args.lookbacktime_days:
				lookbacktime_days = args.lookbacktime_days
			else:
				lookbacktime_days = 60
			self.downloadyselc(args,ra,dec,lookbacktime_days=lookbacktime_days)

		for filt in ['c','o']:
			filename = self.yselcfilename(TNSname,filt)
			outname = self.atlas2yse(TNSname,filename,ra,dec,filename,filt)
			self.uploadtoyse(outname)

if __name__ == '__main__':

	upltoyse = uploadtoyseclass()

	# add arguments
	parser = argparse.ArgumentParser(conflict_handler='resolve')
	parser.add_argument('-t','--tnsnamelist', default=None, nargs='+', help='name of transients to download and upload')
	parser.add_argument('-n','--tnslistfilename', default=None, help='address of file containing TNS names, ra, and dec')
	parser.add_argument('--sourcedir', default=None, help='source code directory')
	parser.add_argument('--outrootdir', default=None, help='output root directory')
	parser.add_argument('--outsubdir', default=None, help='output subdirectory')

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

	for TNSname in upltoyse.TNSnamelist:
		print("Uploading and/or downloading data for ",TNSname)
		upltoyse.uploadloop(args,TNSname,overwrite=True)

	if not(args.tnslistfilename is None):
		upltoyse.TNSlistfile.write(upltoyse.TNSlistfilename,overwrite=True)

