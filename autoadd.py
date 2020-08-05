#!/usr/bin/env python

from astropy.io import ascii
import argparse
from astrotable import astrotableclass
from tools import yamlcfgclass
import pandas as pd
from lxml import html
import requests
import sys,os

class autoaddclass(astrotableclass):
	def __init__(self):
		astrotableclass.__init__(self)

		self.cfg = yamlcfgclass()
		self.debug = False
		self.outrootdir = None
		self.snlist = None

	def define_options(self):
		if 'ATLASLC_SOURCEDIR' in os.environ and os.environ['ATLASLC_SOURCEDIR'] != '':
			outrootdir = os.environ['ATLASLC_DATA']
			cfgfile = '%s/precursor.cfg' % os.environ['ATLASLC_SOURCEDIR']
		else:
			cfgfile = None
			outrootdir = None
		return(cfgfile)

	def setoutdir(self,outrootdir=None):
		if outrootdir is None:
			basedir = self.cfg.params['output']['outrootdir']
		else:
			basedir = outrootdir
		self.outrootdir = basedir

	def loadcfgfiles(self, *pargs, **kwargs):
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
		return(ra, dec)

	def addrow2snlist(self, tnsname, ra, dec):
		print('snlist.txt OLD: \n',self.snlist)
		self.snlist.add_row({'atlasdesignation':'x','otherdesignation':'x','ra':ra,'dec':dec,'spectraltype':'x','earliestmjd':'x','tnsname':tnsname})
		print('snlist.txt NEW: \n',self.snlist)

	def getcmd(self, tnsname):
		print('Your command is: download_lc_loop_sofia.py %s -v -o -s --forcedphot_offset True --plot True --averagelc True --MJDbinsize 1 -l 70' % tnsname)

if __name__ == '__main__':
	autoadd = autoaddclass()

	parser = argparse.ArgumentParser()
	parser.add_argument('tnsname', help="TNS name", type=str)
	parser.add_argument('--ra', help="RA position", default=None, type=str)
	parser.add_argument('--dec', help="Dec position", default=None, type=str)
	args = parser.parse_args()

	cfgfile=autoadd.define_options()
	autoadd.loadcfgfiles(cfgfile)
	autoadd.setoutdir()
	snlistfilename = autoadd.cfg.params['snlistfilename']
	autoadd.snlist = ascii.read(snlistfilename)

	if not args.ra:
		print(args.tnsname)
		ra, dec = autoadd.getradec(args.tnsname)
	else:
		ra=args.ra
		dec=args.dec

	autoadd.addrow2snlist(tnsname=args.tnsname,ra=ra,dec=dec)

	ascii.write([autoadd.snlist['atlasdesignation'], autoadd.snlist['otherdesignation'], autoadd.snlist['ra'], autoadd.snlist['dec'], autoadd.snlist['spectraltype'], autoadd.snlist['earliestmjd'], autoadd.snlist['tnsname']], 
				snlistfilename, names=['atlasdesignation', 'otherdesignation', 'ra', 'dec', 'spectraltype', 'earliestmjd', 'tnsname'], 
				format='fixed_width_two_line', overwrite=True)

	autoadd.getcmd(tnsname=args.tnsname)
