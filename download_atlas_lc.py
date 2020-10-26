#!/usr/bin/env python

import astropy.table as at
import sys,socket
from astropy.io import ascii
from datetime import datetime as dt
import requests, re
from jumpssh import SSHSession
import argparse
from astropy.time import Time
from astrotable import astrotableclass
from tools import DecInDeg,RaInDeg

class download_atlas_lc_class:
	def __init__(self):
		self.verbose = 0
		self.outputbasedir = '.'
		self.remote_session = None
		self.lcraw=None
		
	def define_args(self, parser=None, usage=None,conflict_handler='resolve'):
		if parser is None:
			parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)
		parser.add_argument('RA',help="RA position")
		parser.add_argument('Dec',help="Dec position")
		
		return(parser)
		
	def define_optional_args(self, parser=None, usage=None,conflict_handler='resolve'):
		if parser is None:
			parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)
		
		parser.add_argument('--verbose', '-v', default=0, action='count')
		parser.add_argument('-d', '--debug', help="debug", action='count')
		parser.add_argument('--outputbasedir', default=None,
							help=('Define the output basedirectory (default=%(default)s)'))
		parser.add_argument('-l','--lookbacktime_days', type=int, default=None, help=("lookback time in days"))
		parser.add_argument('--atlasmachine', default='atlas-base-sc01.ifa.hawaii.edu',
							help=('address for atlas machine (default=%(default)s)'))
		parser.add_argument('--user', default=None,
							help=('user name for atlas machine (default=%(default)s)'))
		parser.add_argument('--passwd', default=None,
							help=('password for atlas machine'))
		parser.add_argument('-s','--savelc', help=("save lc with this name"))
		parser.add_argument('-f','--fileformat', default='fixed_width_two_line', choices=['basic','csv','rdb','tab','fixed_width','fixed_width_two_line'],
							help=("specify the file format (https://docs.astropy.org/en/stable/io/ascii/index.html#supported-formats) (default=%(default)s)"))
		parser.add_argument('-o','--overwrite', help="overwrite existing file if saved",
							action="store_true", default=False)
		
		return(parser)
	   
	def connect(self, address, username=None, passwd=None):
		if self.verbose:
			print('Connecting session...')

		self.gateway_session = SSHSession(address, username, password=passwd).open()

		if self.verbose:
			print('session %s@%s' % ( username, address ))
		
		self.remote_session = self.gateway_session.get_remote_session(address, password=passwd)
		
		return(self.check_connection())

	def check_connection(self):
		if not self.remote_session.is_active():
			print("Something went wrong with connection!!!")
			return(False)
		
		if self.verbose:
			print ('Connection established!!')

		return(True)
		
	def build_cmd(self, ra, dec, dodb=1, parallel=20, lookbacktime_days=None):
		print(RaInDeg(ra),DecInDeg(dec),ra,dec)
		cmd = 'force.sh %f %f' % (RaInDeg(ra), DecInDeg(dec))
		
		if dodb>0:
			cmd+=' dodb=%d' % dodb 
		
		if parallel>0:
			cmd+=' parallel=%d' % parallel
		
		if lookbacktime_days!=None:
			mjd = Time.now().mjd
			if self.verbose>1:
				print('Today\'s MJD: %.0f' % mjd)
			cmd += ' m0=%d' % int(mjd - lookbacktime_days)
		
		return(cmd)

	def exec_cmd(self,cmd):
		if self.remote_session is None:
			raise(RuntimeError,"Remote session not initiated!!")
		
		tmp=self.remote_session.get_cmd_output(cmd)
		self.lcraw = tmp.split('\r\n')
		Nlines = len(self.lcraw)-1

		if self.verbose>1:
			print('LC has %d lines of data' % Nlines)
		del(tmp)

		# remove leading #
		if Nlines>=0:
			self.lcraw[0]= re.sub('^\#+','',self.lcraw[0])
		
		return(0)

	def get_lc(self, ra, dec, lookbacktime_days=None):
		# check connection
		if not self.check_connection():
			raise(RuntimeError,"No Connection!!!")
		
		# build the cmd
		cmd = self.build_cmd(ra, dec, lookbacktime_days=lookbacktime_days)
		print('Command for photometry: %s' % cmd)

		# run the cmd
		self.exec_cmd(cmd)
		if self.verbose>2:
			print('Raw ATLAS lc:')
			for s in self.lcraw:print(s)

	def save_lc(self, filename, overwrite=False, fileformat='fixed_width_two_line'):
		self.lc.write(filename, format=fileformat, overwrite=overwrite, verbose=(self.verbose>0))
			
if __name__ == "__main__":
	download_atlas_lc = download_atlas_lc_class()
	parser = download_atlas_lc.define_args()
	parser = download_atlas_lc.define_optional_args(parser=parser)
	args = parser.parse_args()

	download_atlas_lc.verbose = args.verbose
	download_atlas_lc.debug = args.debug
	
	download_atlas_lc.connect(args.atlasmachine,args.user,args.passwd)
	download_atlas_lc.get_lc(args.RA,args.Dec,
							lookbacktime_days=args.lookbacktime_days,
							overwrite=args.overwrite,
							fileformat=args.fileformat)
	
	if not(args.savefile is None):
		download_atlas_lc.save_lc(args.savelc, overwrite=args.overwrite, fileformat=args.fileformat)
