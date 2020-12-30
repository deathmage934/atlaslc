#!/usr/bin/env python
'''
@author: S. Rest
'''

import argparse
import numpy as np
from astropy.time import Time
from pdastro import pdastroclass
from SNloop import SNloopclass
from tools import yamlcfgclass
from tools import RaInDeg
from tools import DecInDeg
import pandas as pd
from lxml import html
import requests
import sys,os

class autoaddclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)

		self.cfg = yamlcfgclass()
		self.debug = False
		self.outrootdir = None
		self.snlist = pdastroclass(columns=['atlasdesignation','otherdesignation','ra','dec','spectraltype','earliestmjd','tnsname','closebrightRA','closebrightDec','MJDpreSN'])
		self.snlistsurvey = pdastroclass(columns=['atlasdesignation','otherdesignation','ra','dec',	'spectraltype','earliestmjd','tnsname','closebrightRA','closebrightDec','MJDpreSN'])

	def defineoptions(self):
		if 'ATLASLC_SOURCEDIR' in os.environ and os.environ['ATLASLC_SOURCEDIR'] != '':
			outrootdir = os.environ['ATLASLC_DATA']
			cfgfile = '%s/precursor.cfg' % os.environ['ATLASLC_SOURCEDIR']
		else:
			cfgfile = None
			outrootdir = None
		return(cfgfile)
	
	def setoutputdir(self,outrootdir=None):
		if outrootdir is None:
			basedir = self.cfg.params['output']['outrootdir']
		else:
			basedir = outrootdir
		self.outrootdir = basedir
	
	def loadcfgfile(self, *pargs, **kwargs):
		result=self.cfg.loadcfgfiles(*pargs, **kwargs)
		return(result)

	def getradec(self, tnsname):
		obj_name=tnsname
		print('Getting data from https://wis-tns.weizmann.ac.il/object/'+obj_name)
		page = requests.get('https://wis-tns.weizmann.ac.il/object/'+obj_name)
		tree = html.fromstring(page.content)
		
		radec = tree.xpath('//div[@class="field field-radec"]//div[@class="value"]/text()')
		print(radec)

		def my_split(seq):
			for item in seq:
				a, b, c = item.partition(' ')
				yield a
				yield b+c

		datasplit = list(my_split(radec))
		ra = datasplit[0]
		dec = datasplit[1]
		print('RA: %s, Dec: %s' % (ra, dec))
		print('RA in degrees: %s, Dec in degrees: %s' % (RaInDeg(ra),DecInDeg(dec)))
		return(ra, dec)

	def getdisc_date(self,tnsname):
		obj_name=tnsname
		print('Getting data from https://wis-tns.weizmann.ac.il/object/'+obj_name)
		page = requests.get('https://wis-tns.weizmann.ac.il/object/'+obj_name)
		tree = html.fromstring(page.content)
		
		discoverydate = tree.xpath('//div[@class="field field-discoverydate"]//div[@class="value"]/b/text()')
		print(discoverydate)

		def my_split(seq):
			for item in seq:
				a, b, c = item.partition(' ')
				yield a
				yield b+c

		datasplit2 = list(my_split(discoverydate))
		date = datasplit2[0]
		time = datasplit2[1][1:]
		#time = time[1:]
		print('Date: %s, Time: %s' % (date, time))

		disc_date_format = date+"T"+time
		dateobjects = Time(disc_date_format, format='isot', scale='utc')
		disc_date = dateobjects.mjd
		print("MJD: %.4f" % disc_date)
		return(disc_date)

	def addrow2snlist(self, tnsname, ra, dec, MJDpreSN, closebrightRA=None, closebrightDec=None):
		if os.path.exists(snlistfilename):
			print('snlist.txt OLD:')
			print(self.snlist.write())

		if closebrightRA is None:
			df = pd.DataFrame([['<NA>','<NA>',ra,dec,'<NA>','NaN',tnsname,'<NA>','<NA>',disc_date]], columns=['atlasdesignation','otherdesignation','ra','dec','spectraltype','earliestmjd','tnsname','closebrightRA','closebrightDec','MJDpreSN'])
		else: 
			df = pd.DataFrame([['<NA>','<NA>',ra,dec,'<NA>','NaN',tnsname,closebrightRA,closebrightDec,disc_date]], columns=['atlasdesignation','otherdesignation','ra','dec','spectraltype','earliestmjd','tnsname','closebrightRA','closebrightDec','MJDpreSN'])
		self.snlist.t = self.snlist.t.append(df, ignore_index=True)
		
		print('snlist.txt NEW:')
		print(self.snlist.write())

	def getcmd(self, tnsname):
		print('Your command is: download_lc_loop.py %s -v -o -s --forcedphot_offset True --plot True --averagelc True --MJDbinsize 1 -l 70 --user %s --passwd "X"' % (tnsname, self.cfg.params['username']))

if __name__ == '__main__':
	autoadd = autoaddclass()

	parser = argparse.ArgumentParser()
	parser.add_argument('tnsname', help="TNS name; type 'survey' for survey MJDpreSN", type=str)
	parser.add_argument('--ra', help="RA position", default=None, type=str)
	parser.add_argument('--dec', help="Dec position", default=None, type=str)
	parser.add_argument('--disc_date', help="Discovery date of SN", default=None, type=float)
	parser.add_argument('--autosearch',help="Search automatically for the closest bright object",default=False)
	args = parser.parse_args()

	cfgfile=autoadd.defineoptions()
	autoadd.loadcfgfile(cfgfile)
	autoadd.setoutputdir()

	# special survey procedure
	if args.tnsname == 'survey':
		snlistfilename = autoadd.cfg.params['snlistsurveyfilename']
		autoadd.snlistsurvey.default_formatters = {'MJDpreSN':'{:.4f}'.format}
		autoadd.snlistsurvey.load_spacesep(snlistfilename, delim_whitespace=True)

		indexlist = autoadd.snlistsurvey.getindices()
		#indexlist = range(len(autoadd.snlistsurvey.t))
		for index in indexlist:
			try: 
				tnsname = autoadd.snlistsurvey.t.at[index,'tnsname']
				print('TNSname: %s, index: %d/%d' % (tnsname, index, max(indexlist)))
				if not(isinstance(tnsname,str)):
					print('Found nan, skipping...')
					continue
				else:
					disc_date = autoadd.getdisc_date(tnsname)
					autoadd.snlistsurvey.t.at[index,'MJDpreSN'] = disc_date.astype(float)
					print('ADDED: %.4f' % autoadd.snlistsurvey.t.at[index,'MJDpreSN']) # delete me
			except:
				autoadd.snlistsurvey.write(snlistfilename,overwrite=True,verbose=True)
				raise RuntimeError('ERROR at index %d! Saving table...' % index)
		autoadd.snlistsurvey.write(snlistfilename,overwrite=True,verbose=True)
	
	# regular procedure
	else:
		if not args.ra:
			print(args.tnsname)
			ra, dec = autoadd.getradec(args.tnsname)
			if not args.disc_date:
				disc_date = autoadd.getdisc_date(args.tnsname)
			else:
				disc_date = args.disc_date
		else:
			ra=args.ra
			dec=args.dec
			if not args.disc_date:
				disc_date = autoadd.getdisc_date(args.tnsname)
			else:
				disc_date = args.disc_date
		
		deltat = autoadd.cfg.params['deltat']
		disc_date = disc_date.astype(float) - deltat
		print(disc_date, type(disc_date)) # delete me

		snlistfilename = autoadd.cfg.params['snlistfilename']
		if os.path.exists(snlistfilename):
			autoadd.snlist.load_spacesep(snlistfilename, delim_whitespace=True)
		else: 
			autoadd.snlist.t = pd.DataFrame({'atlasdesignation':[],'otherdesignation':[],'ra':[],'dec':[],'spectraltype':[],'earliestmjd':[],'tnsname':[],'closebrightRA':[],'closebrightDec':[],'MJDpreSN':[]})

		if args.autosearch:
			print('Running autosearch for nearest close bright object\nObject will be within 20 arcsec, brighter than 18 mag, and could be a star or galaxy')
			results = autoadd.autosearch(RaInDeg(ra), DecInDeg(dec), 20)
			print('Close bright objects found: \n',results)
			print(RaInDeg('19:58:26.6356'),DecInDeg('+62:08:05.497')) # delete me
			closebrightRA = results['182562996106752735']['raMean']
			closebrightDec = results['182562996106752735']['decMean']
			print(closebrightRA,closebrightDec)
			autoadd.addrow2snlist(tnsname=args.tnsname,ra=ra,dec=dec,MJDpreSN=disc_date,closebrightRA=closebrightRA,closebrightDec=closebrightDec)
		else:
			autoadd.addrow2snlist(tnsname=args.tnsname,ra=ra,dec=dec,MJDpreSN=disc_date)
		
		autoadd.snlist.write(snlistfilename,overwrite=True,verbose=True)
		autoadd.getcmd(tnsname=args.tnsname)
