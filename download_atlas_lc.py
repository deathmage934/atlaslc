#!/usr/bin/env python

# Adapted from Q. Wang by S. Rest

import astropy.table as at
import sys,socket
from astropy.io import ascii
from datetime import datetime as dt
import requests,re,io,sys
import sqlite3
from jumpssh import SSHSession
import argparse
import time
from astropy.time import Time
#from astrotable import astrotableclass
from tools import DecInDeg,RaInDeg
import pandas as pd

class download_atlas_lc_class:
	def __init__(self):
		self.verbose = 0
		self.outputbasedir = '.'
		self.remote_session = None
		self.lcraw = None
		self.baseurl = 'https://fallingstar-data.com/forcedphot'
		
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

	"""
	def connect_api(self,username,password):
		#ra, dec, name = sys.argv[1:]
		token_header = connect_atlas(username,password)
		data = get_result(ra, dec, token_header)
		#ascii.write(data, name + '.csv', format = 'csv', overwrite = True)
	"""

	def connect_atlas(self,username,password):
		resp = requests.post(url=f"{self.baseurl}/api-token-auth/",data={'username':username,'password':password})
		if resp.status_code == 200:
			token = resp.json()['token']
			print(f'Your token is {token}')
			headers = {'Authorization':f'Token {token}','Accept':'application/json'}
		else:
			print(f'ERROR {resp.status_code}')
			print(resp.json())
		return headers

	# API GUIDE: https://fallingstar-data.com/forcedphot/apiguide/
	def get_result(self,ra,dec,headers,lookbacktime_days=None,mjd_max=None):
		today = dt.today()
		con = sqlite3.connect(":memory:")

		if not(lookbacktime_days is None):
			lookbacktime_days = int(Time.now().mjd - lookbacktime_days)
			#lookbacktime_days = '  '+str(list(con.execute("select julianday('"+today.strftime("%Y-%m-%d")+"')"))[0][0]-lookbacktime_days-2400000)
		else:
			lookbacktime_days = int(Time.now().mjd - 1900)

		task_url = None
		while not task_url:
			with requests.Session() as s:
				resp = s.post(f"{self.baseurl}/queue/",headers=headers,data={'ra':ra,'dec':dec,'send_email':False,"mjd_min":lookbacktime_days,"mjd_max":mjd_max})
				if resp.status_code == 201:  # success
					task_url = resp.json()['url']
					print(f'The task URL is {task_url}')
				elif resp.status_code == 429:  # throttled
					message = resp.json()["detail"]
					print(f'{resp.status_code} {message}')
					t_sec = re.findall(r'available in (\d+) seconds', message)
					t_min = re.findall(r'available in (\d+) minutes', message)
					if t_sec:
						waittime = int(t_sec[0])
					elif t_min:
						waittime = int(t_min[0]) * 60
					else:
						waittime = 10
					print(f'Waiting {waittime} seconds')
					time.sleep(waittime)
				else:
					print(f'ERROR {resp.status_code}')
					print(resp.json())
					sys.exit()
		result_url = None
		
		while not result_url:
			with requests.Session() as s:
				resp = s.get(task_url, headers=headers)
				if resp.status_code == 200:  # HTTP OK
					if resp.json()['finishtimestamp']:
						result_url = resp.json()['result_url']
						print(f"Task is complete with results available at {result_url}")
						break
					elif resp.json()['starttimestamp']:
						print(f"Task is running (started at {resp.json()['starttimestamp']})")
					else:
						print("Waiting for job to start. Checking again in 10 seconds...")
					time.sleep(10)
				else:
					print(f'ERROR {resp.status_code}')
					print(resp.json())
					sys.exit()
		
		with requests.Session() as s:
			result = s.get(result_url, headers=headers).text
		
		dfresult = pd.read_csv(io.StringIO(result.replace("###", "")), delim_whitespace=True)
		return dfresult

# don't need the following main:
"""
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
"""
