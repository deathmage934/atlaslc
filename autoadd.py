#!/usr/bin/env python

import argparse
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
		self.snlist = pdastroclass(columns=['atlasdesignation','otherdesignation','ra','dec','spectraltype','earliestmjd','tnsname','closebrightRA','closebrightDec'])

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
		data = tree.xpath('//div[@class="value"]/text()')

		def my_split(seq):
			for item in seq:
				a, b, c = item.partition(' ')
				yield a
				yield b+c

		datasplit = list(my_split(data))
		ra = datasplit[0]
		dec = datasplit[1]
		print('RA: %s, Dec: %s' % (ra, dec))
		print('RA in degrees: %s, Dec in degrees: %s' % (RaInDeg(ra),DecInDeg(dec)))
		return(ra, dec)

	def addrow2snlist(self, tnsname, ra, dec, closebrightRA=None, closebrightDec=None):
		if os.path.exists(snlistfilename):
			print('snlist.txt OLD:')
			print(self.snlist.write())

		if closebrightRA is None:
			df = pd.DataFrame([['x','x',ra,dec,'x','x',tnsname,'x','x']], columns=['atlasdesignation','otherdesignation','ra','dec','spectraltype','earliestmjd','tnsname','closebrightRA','closebrightDec'])
		else: 
			df = pd.DataFrame([['x','x',ra,dec,'x','x',tnsname,closebrightRA,closebrightDec]], columns=['atlasdesignation','otherdesignation','ra','dec','spectraltype','earliestmjd','tnsname','closebrightRA','closebrightDec'])
		self.snlist.t = self.snlist.t.append(df, ignore_index=True)
		
		print('snlist.txt NEW:')
		print(self.snlist.write())

	def getcmd(self, tnsname):
		print('Your command is: download_lc_loop.py %s -v -o -s --forcedphot_offset True --plot True --averagelc True --MJDbinsize 1 -l 70 --user %s --passwd "X"' % (tnsname, self.cfg.params['username']))

if __name__ == '__main__':
	autoadd = autoaddclass()

	parser = argparse.ArgumentParser()
	parser.add_argument('tnsname', help="TNS name", type=str)
	parser.add_argument('--ra', help="RA position", default=None, type=str)
	parser.add_argument('--dec', help="Dec position", default=None, type=str)
	parser.add_argument('--autosearch',help="Search automatically for the closest bright object",default=False)
	args = parser.parse_args()

	cfgfile=autoadd.defineoptions()
	autoadd.loadcfgfile(cfgfile)
	autoadd.setoutputdir()

	if not args.ra:
		print(args.tnsname)
		ra, dec = autoadd.getradec(args.tnsname)
	else:
		ra=args.ra
		dec=args.dec

	snlistfilename = autoadd.cfg.params['snlistfilename']
	if os.path.exists(snlistfilename):
		autoadd.snlist.load_spacesep(snlistfilename, delim_whitespace=True)
	else: 
		autoadd.snlist.t = pd.DataFrame({'atlasdesignation':[],'otherdesignation':[],'ra':[],'dec':[],'spectraltype':[],'earliestmjd':[],'tnsname':[],'closebrightRA':[],'closebrightDec':[]})

	if args.autosearch:
		print('Running autosearch for nearest close bright object\nObject will be within 20 arcsec, brighter than 18 mag, and could be a star or galaxy')
		results = autoadd.autosearch(RaInDeg(ra), DecInDeg(dec), 20)
		print('Close bright objects found: \n',results)
		print(RaInDeg('19:58:26.6356'),DecInDeg('+62:08:05.497')) # delete me
		closebrightRA = results['182562996106752735']['raMean']
		closebrightDec = results['182562996106752735']['decMean']
		print(closebrightRA,closebrightDec)
		sys.exit(0)
		autoadd.addrow2snlist(tnsname=args.tnsname,ra=ra,dec=dec, closebrightRA=closebrightRA, closebrightDec=closebrightDec)
	else:
		autoadd.addrow2snlist(tnsname=args.tnsname,ra=ra,dec=dec)
	
	autoadd.snlist.write(snlistfilename,overwrite=True,verbose=True)
	autoadd.getcmd(tnsname=args.tnsname)
