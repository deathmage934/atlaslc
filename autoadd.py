#!/usr/bin/env python
'''
@author: S. Rest
'''

import argparse,requests,sys,os,time,json
import numpy as np
from astropy.time import Time
from pdastro import pdastroclass
from SNloop import SNloopclass
from tools import yamlcfgclass
from tools import RaInDeg
from tools import DecInDeg
import pandas as pd
from lxml import html
from collections import OrderedDict

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

	def my_split(self,seq):
		for item in seq:
			a, b, c = item.partition(' ')
			yield a
			yield b+c

	def getdata(self,tnsname):
		try:
			get_obj = [("objname",tnsname), ("objid",""), ("photometry","1"), ("spectra","1")]
			get_url = 'https://www.wis-tns.org/api/get/object'
			json_file = OrderedDict(get_obj)
			get_data = {'api_key':'2eca323a16b17d78fbc99cd6f1f801699a81a91c','data':json.dumps(json_file)}
			response = requests.post(get_url, data=get_data, headers={'User-Agent':'tns_marker{"tns_id":104739,"type": "bot", "name":"Name and Redshift Retriever"}'})
			json_data = json.loads(response.text,object_pairs_hook=OrderedDict) #self.format(response.text)
			return json_data
		except Exception as e:
			return [None,'Error message : \n'+str(e)]

	def getradec(self,tnsname):
		json_data = self.getdata(tnsname)
		ra = json_data['data']['reply']['ra']
		dec = json_data['data']['reply']['dec']
		print('In sexagesimal: RA: %s, Dec: %s' % (ra,dec))
		print('In decimal: RA %0.8f Dec: %0.8f' % (RaInDeg(ra),DecInDeg(dec)))
		return ra, dec
		"""
		try:
			ra = json_data['data']['reply']['ra']
			dec = json_data['data']['reply']['dec']
			print('In sexagesimal: RA: %s, Dec: %s' % (ra,dec))
			print('In decimal: RA %0.8f Dec: %0.8f' % (RaInDeg(ra),DecInDeg(dec)))
			return ra, dec
		except KeyError:
			print('No results found, skipping...')
			return np.nan, np.nan
		"""

	def getdisc_date(self,tnsname):
		json_data = self.getdata(tnsname)
		discoverydate = json_data['data']['reply']['discoverydate']

		date = list(discoverydate.partition(' '))[0]
		time = list(discoverydate.partition(' '))[2]

		disc_date_format = date+"T"+time
		dateobjects = Time(disc_date_format, format='isot', scale='utc')
		disc_date = dateobjects.mjd
		print("Discovery Date MJD: %.4f" % disc_date)
		return disc_date

	def addrow2snlist(self, tnsname, ra, dec, MJDpreSN, closebrightRA=None, closebrightDec=None):
		if closebrightRA is None:
			df = pd.DataFrame([['<NA>','<NA>','%0.8f'%RaInDeg(ra),'%0.8f'%DecInDeg(dec),'<NA>','NaN',tnsname,'<NA>','<NA>',disc_date]], columns=['atlasdesignation','otherdesignation','ra','dec','spectraltype','earliestmjd','tnsname','closebrightRA','closebrightDec','MJDpreSN'])
		else: 
			df = pd.DataFrame([['<NA>','<NA>','%0.8f'%RaInDeg(ra),'%0.8f'%DecInDeg(dec),'<NA>','NaN',tnsname,closebrightRA,closebrightDec,disc_date]], columns=['atlasdesignation','otherdesignation','ra','dec','spectraltype','earliestmjd','tnsname','closebrightRA','closebrightDec','MJDpreSN'])
		print('Adding row: \n',df)
		self.snlist.t = self.snlist.t.append(df, ignore_index=True)
		#print(self.snlist.t)

	def getcmd(self, tnsname):
		print('Your command is: download_lc_loop.py %s -v -o -s --api --forcedphot_offset --plot --user %s --passwd "X"' % (tnsname, self.cfg.params['username']))

if __name__ == '__main__':
	autoadd = autoaddclass()

	parser = argparse.ArgumentParser()
	parser.add_argument('tnsname', nargs='+', help="TNS name; type 'survey' for survey MJDpreSN", type=str)
	parser.add_argument('--ra', help="RA position", default=None, type=str)
	parser.add_argument('--dec', help="Dec position", default=None, type=str)
	parser.add_argument('--disc_date', help="Discovery date of SN", default=None, type=float)
	parser.add_argument('--autosearch', action="store_true", help="Search automatically for the closest bright object",default=False)
	args = parser.parse_args()

	cfgfile = autoadd.defineoptions()
	autoadd.loadcfgfile(cfgfile)
	autoadd.setoutputdir()

	# special survey procedure
	if args.tnsname == 'survey':
		snlistfilename = autoadd.cfg.params['snlistsurveyfilename']
		autoadd.snlistsurvey.default_formatters = {'MJDpreSN':'{:.4f}'.format}
		autoadd.snlistsurvey.load_spacesep(snlistfilename, delim_whitespace=True)

		indexlist = autoadd.snlistsurvey.getindices()
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
		print(len(args.tnsname))
		for i in range(0,len(args.tnsname)):
			tnsname = args.tnsname[i]
			if not args.ra:
				print('Searching for '+tnsname+', index %d/%d' % (i,len(args.tnsname)-1))
				try:
					ra, dec = autoadd.getradec(tnsname)
				except KeyError:
					print('No results found, skipping...')
					ra = np.nan 
					dec = np.nan 
				if isinstance(ra,str): #not(np.isnan(ra)):
					if not args.disc_date:
						disc_date = autoadd.getdisc_date(tnsname)
					else:
						disc_date = args.disc_date
				else:
					disc_date = np.nan
			else:
				ra=args.ra
				dec=args.dec
				if not args.disc_date:
					disc_date = autoadd.getdisc_date(tnsname)
				else:
					disc_date = args.disc_date
			
			if isinstance(disc_date,str):
				deltat = autoadd.cfg.params['deltat']
				disc_date = disc_date.astype(float) - deltat

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
				if isinstance(ra,str):
					autoadd.addrow2snlist(tnsname=tnsname,ra=ra,dec=dec,MJDpreSN=disc_date,closebrightRA=closebrightRA,closebrightDec=closebrightDec)
			else:
				if isinstance(ra,str):
					autoadd.addrow2snlist(tnsname=tnsname,ra=ra,dec=dec,MJDpreSN=disc_date)
			
			autoadd.snlist.write(snlistfilename,overwrite=True,verbose=True)
			
			if isinstance(ra,str):
				autoadd.getcmd(tnsname=tnsname)

			print('sleeping...')
			time.sleep(20)
			print('done')

