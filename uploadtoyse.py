#!/usr/bin/env python

# Code adapted from Q. Wang and D. Jones by S. Rest

import sigmacut
from SNloop import SNloopclass
from download_lc_loop import downloadlcloopclass
from uploadTransientData import upload, runDBcommand, DBOps
from autoadd import autoaddclass
from pdastro import pdastroclass, pdastrostatsclass, AnotB, AandB
from download_atlas_lc import download_atlas_lc_class
from tools import RaInDeg,DecInDeg

import pandas as pd
import numpy as np
import sys,socket,os,re,io,math
import optparse,configparser,argparse
from copy import deepcopy
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
import requests,json,urllib.request,urllib,ast,datetime,time,coreapi
from requests.exceptions import MissingSchema
from requests.auth import HTTPBasicAuth

"""
get & upload original ATLAS light curves from daily YSE list: 
    uploadtoyse.py --user USERNAME --passwd 'PASSWORD' --api
get & upload specific original ATLAS light curves: 
    uploadtoyse.py -t 2020lse 2019vxm --user USERNAME --passwd 'PASSWORD' --api
get & upload original ATLAS light curves from user-made table of TNS names, RA, and Dec: 
    uploadtoyse.py -n tnslist.txt --user USERNAME --passwd 'PASSWORD' --api
get & upload original ATLAS light curves AND averaged light curves: 
    add --averagelc to the command
        uploadtoyse.py --user USERNAME --passwd 'PASSWORD' --api --averagelc
    change the MJD bin size to 0.04 days instead of the default 1.00 day when averaging:
        add --MJDbinsize 0.04 to the command (must be averaging light curve with the command --averagelc)
            uploadtoyse.py --user USERNAME --passwd 'PASSWORD' --api --averagelc --MJDbinsize 0.04
set the lookback time in days to 500 instead of the default 60:
    add -l 500 to the command
        uploadtoyse.py --user USERNAME --passwd 'PASSWORD' --api -l 500

you must specify source directory, output root directory (has tnslist.txt if you want to use it), and output subdirectory (has all the downloaded data).
you can specify these in the commandline using --sourcedir, --outrootdir, and --outsubdir
OR
1. you can set the source directory and the output root directory in atlaslc.sourceme
2. then set the output subdirectory in precursor.cfg through the variable 'yse_outsubdir' (currently set to default 'ysetest')
"""

class uploadtoyseclass(downloadlcloopclass,autoaddclass):
    def __init__(self):
        downloadlcloopclass.__init__(self)
        autoaddclass.__init__(self)

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

        #flags
        self.flag_day_bad = 0x800000

    def add_options(self,parser=None,usage=None,config=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage,conflict_handler="resolve")
        parser.add_argument('-s','--settingsfile', default='%s/settings.ini'%self.sourcedir, type=str,help='settings file (login/password info)')
        parser.add_argument('-v', '--verbose', action="count", dest="verbose",default=1)
        parser.add_argument('--clobber', default=False, action="store_true",help='clobber output file')
        parser.add_argument('--status', default='Following', type=str,help='transient status (new, follow, etc')
        parser.add_argument('--obsgroup', default='Foundation', type=str,help='group who observed this transient')
        parser.add_argument('--permissionsgroup', default='', type=str,help='group that has permission to view this photometry on YSE_PZ')
        parser.add_argument('--inputformat', default='basic', type=str,help="input file format, can be 'basic' or 'snana' (photometry only) ")
        parser.add_argument('--instrument', default='ACAM1', type=str,help="instrument name")
        parser.add_argument('--forcedphot', default=1, type=int,help="set to 1 if forced photometry")
        parser.add_argument('--diffim', default=1, type=int,help="set to 1 if difference imaged")
        parser.add_argument('--fluxzpt', default=23.9, type=float,help="flux zero point")
        parser.add_argument('-e','--onlyexisting', default=True, action="store_true",help="only add light curves for existing objects")
        parser.add_argument('-m','--mjdmatchmin', default=0.01, type=float,help="""if clobber flag not set, photometric observation with MJD separation less than this, in the same filter/instrument are treated as the same data. Allows updates to the photometry""")
        parser.add_argument('--spectrum', default=False, action="store_true",help='input file is a spectrum') 
        return(parser)

    def db_init_params(self,args):
        self.dblogin = args.dblogin
        self.dbpassword = args.dbpassword
        self.dburl = args.dburl
        self.baseposturl = "http --ignore-stdin -a %s:%s POST %s"%(self.dblogin,self.dbpassword,self.dburl)
        self.basegeturl = "http --ignore-stdin -a %s:%s GET %s"%(self.dblogin,self.dbpassword,self.dburl)
        self.baseputurl = "http --ignore-stdin -a %s:%s PUT %s"%(self.dblogin,self.dbpassword,self.dburl)
        self.basegetobjurl = "http --ignore-stdin -a %s:%s GET "%(self.dblogin,self.dbpassword)

    def add_db_options(self,parser=None,usage=None,config=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage,conflict_handler="resolve")
        parser.add_argument('-p','--photometry',default=False,action="store_true",help='input file is photometry')
        parser.add_argument('--login', default=config.get('main','login'), type=str,help='gmail login (default=%default)')
        parser.add_argument('--password', default=config.get('main','password'), type=str,help='gmail password (default=%default)')
        parser.add_argument('--dblogin', default=config.get('main','dblogin'), type=str,help='gmail login (default=%default)')
        parser.add_argument('--dbpassword', default=config.get('main','dbpassword'), type=str,help='gmail password (default=%default)')
        parser.add_argument('--dburl', default=config.get('main','dburl'), type=str,help='base URL to POST/GET,PUT to/from a database (default=%default)')
        parser.add_argument('--transientapi', default=config.get('main','transientapi'), type=str,help='URL to POST transients to a database (default=%default)')
        parser.add_argument('--internalsurveyapi', default=config.get('main','internalsurveyapi'), type=str,help='URL to POST transients to a database (default=%default)')
        parser.add_argument('--transientclassesapi', default=config.get('main','transientclassesapi'), type=str,help='URL to POST transients classes to a database (default=%default)')
        parser.add_argument('--hostapi', default=config.get('main','hostapi'), type=str,help='URL to POST transients to a database (default=%default)')
        parser.add_argument('--photometryapi', default=config.get('main','photometryapi'), type=str,help='URL to POST transients to a database (default=%default)')
        parser.add_argument('--photdataapi', default=config.get('main','photdataapi'), type=str,help='URL to POST transients to a database (default=%default)')
        parser.add_argument('--hostphotometryapi', default=config.get('main','hostphotometryapi'), type=str,help='URL to POST transients to a database (default=%default)')
        parser.add_argument('--hostphotdataapi', default=config.get('main','hostphotdataapi'), type=str,help='URL to POST transients to a database (default=%default)')
        parser.add_argument('--obs_groupapi', default=config.get('main','obs_groupapi'), type=str,help='URL to POST group to a database (default=%default)')
        parser.add_argument('--statusapi', default=config.get('main','statusapi'), type=str,help='URL to POST status to a database (default=%default)')
        parser.add_argument('--instrumentapi', default=config.get('main','instrumentapi'), type=str,help='URL to POST instrument to a database (default=%default)')
        parser.add_argument('--bandapi', default=config.get('main','bandapi'), type=str,help='URL to POST band to a database (default=%default)')
        parser.add_argument('--observatoryapi', default=config.get('main','observatoryapi'), type=str,help='URL to POST observatory to a database (default=%default)')
        parser.add_argument('--telescopeapi', default=config.get('main','telescopeapi'), type=str,help='URL to POST telescope to a database (default=%default)')
        return(parser)

    def get_transient_from_DB(self,fieldname,debug=False):
        if debug: tstart = time.time()
        tablename = 'transients'
        auth = coreapi.auth.BasicAuthentication(username=self.dblogin,password=self.dbpassword)
        client = coreapi.Client(auth=auth)
        try:
            schema = client.get('%s%s'%(self.dburl.replace('/api','/get_transient'),fieldname))
        except:
            raise RuntimeError('Error : couldn\'t get schema!')
        if not schema['transient']:
            return None
        return(schema['transient']['url'])

    def parsePhotHeaderData(self,args,tnsname,ra,dec):
        transientdict = {'obs_group':args.obsgroup,'status':args.status,'name':tnsname,'ra':ra,'dec':dec,'groups':args.permissionsgroup}
        photdict = {'instrument':args.instrument,'obs_group':args.obsgroup,'transient':tnsname,'groups':args.permissionsgroup}
        return(transientdict,photdict)

    def uploadatlasphotometry(self,args,tnsname,ra,dec,filt,MJDbinsize=None):
        filename = self.yselcfilename(TNSname,0,filt,MJDbinsize=MJDbinsize)
        self.lc.load_spacesep(filename,delim_whitespace=True,raiseError=True)
        
        # SET ANY NANS IN DM AND MJD TO NONE BECAUSE NAN BREAKS IT
        ix_dm = self.lc.ix_null(colnames=['dm'])
        self.lc.t.loc[ix_dm,'dm'] = 5.0
        # remove any nans in the MJD column
        ix_notmjd = AnotB(self.lc.getindices(),self.lc.ix_null(colnames=['MJD']))
        self.lc.t = self.lc.t.loc[ix_notmjd]

        transid = self.get_transient_from_DB(tnsname)
        if args.onlyexisting and not transid:
            raise RuntimeError('Object %s not found! Exiting...' % tnsname)
        print('Uploading object %s...' % tnsname)
        transientdict,photdict = self.parsePhotHeaderData(args,tnsname,ra,dec)
        PhotUploadAll = {'transient':transientdict,'photheader':photdict}

        for mjd,flux,fluxerr,mag,magerr,flt,i,mask in zip(self.lc.t['MJD'],self.lc.t[self.flux_colname],self.lc.t[self.dflux_colname],self.lc.t['m'],self.lc.t['dm'],self.lc.t['F'],range(len(self.lc.t['F'])),self.lc.t['Mask']):
            print(mjd)
            obsdate = Time(mjd,format='mjd').isot
            if flt == 'o': flt = 'orange-ATLAS'
            elif flt == 'c': flt = 'cyan-ATLAS'
            else: 
                raise RuntimeError('Error converting filter %s' % flt)
            
            PhotUploadDict = {'obs_date':obsdate,'flux':flux,'flux_err':fluxerr,'forced':args.forcedphot,'diffim':args.diffim,'band':flt,'groups':[],'flux_zero_point':args.fluxzpt,'discovery_point':0}
            if flux > 0:
                PhotUploadDict['mag'] = mag
                PhotUploadDict['mag_err'] = magerr
            else:
                PhotUploadDict['mag'] = None
                PhotUploadDict['mag_err'] = None
            if not(mask==0):
                PhotUploadDict['data_quality'] = 1
            else:
                PhotUploadDict['data_quality'] = 0
            PhotUploadAll['%s_%i'%(obsdate,i)] = PhotUploadDict
            PhotUploadAll['header'] = {'clobber':args.clobber,'mjdmatchmin':args.mjdmatchmin}

        url = '%s' % args.dburl.replace('/api','/add_transient_phot')
        def myconverter(obj):
            if obj is np.nan:
                return None
            elif isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, datetime.datetime):
                return obj.__str__()

        r = requests.post(url=url,data=json.dumps(PhotUploadAll,default=myconverter),auth=HTTPBasicAuth(args.dblogin,args.dbpassword))
        print('YSE_PZ says: %s'%json.loads(r.text)['message'])

    def yseupload(self,tnsname,ra,dec,filt,MJDbinsize=None,parser=None):
        # load in args, parse config file, then load in db args with config file passed
        parser = self.add_options(parser=parser)
        args = parser.parse_args()
        config = configparser.ConfigParser()
        config.read(args.settingsfile)
        parser = self.add_db_options(parser=parser,config=config)
        args = parser.parse_args()

        self.db_init_params(args)
        self.uploadatlasphotometry(args,tnsname,ra,dec,filt,MJDbinsize=MJDbinsize)

    def YSE_list(self,args):
        yse_list = self.cfg.params['upltoyse']['yse_list']
        if not(args.ysequery is None):
            yse_list = 'https://ziggy.ucolick.org/yse/explorer/'+str(args.ysequery)+'/download?format=csv'
        print('Obtaining YSE list from '+yse_list+'...')
        all_cand = pd.read_csv(yse_list)
        all_cand = all_cand.drop_duplicates(subset='name')
        df = pd.DataFrame()
        df['Name'] = all_cand['name']
        df['RA'] = all_cand['ra']
        df['Dec'] = all_cand['dec']
        df['Disc_date'] = all_cand['disc_date']
        print(df)
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
        filename = self.yselcfilename(TNSname,controlindex,filt,MJDbinsize=MJDbinsize)
        if not(MJDbinsize is None):
            self.averagelctable.load_spacesep(filename,delim_whitespace=True,raiseError=True)
        else:
            self.lc.load_spacesep(filename,delim_whitespace=True,raiseError=True)
        return(0)

    """
    def atlas2yseold(self,TNSname,outname,ra,dec,atlas_data_file,filt):
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
                f.write('OBS: '+str(row['MJD'])+' '+filter_fict[row['F']]+' '+str(row[self.flux_colname])+' '+str(row[self.dflux_colname])+' '+str(row['m'])+' '+str(row['dm'])+' 0 \n')
        print("Converted ATLAS lc to YSE format: %s" % outname)
        return(outname)
    """

    def atlas2yse(self,TNSname,outname,ra,dec,lc,filt):
        lc.t['dummy']='OBS: '
        ix_o =  lc.ix_equal('F','o')
        ix_c =  lc.ix_equal('F','c')
        lc.t.loc[ix_o,'F']='orange-ATLAS'
        lc.t.loc[ix_c,'F']='cyan-ATLAS'
        
        ix = lc.ix_null('dm')
        lc.t.loc[ix,'dm'] = np.nan
        
        ix = lc.ix_null('Mask')
        lc.t.loc[ix,'Mask']=int(0)
        lc.t['Mask']=lc.t['Mask'].astype(int)
        
        with open(outname, 'w+') as f:
            f.write('SNID: '+TNSname+' \nRA: '+str(ra)+'     \nDECL: '+str(dec)+' \n \nVARLIST:  MJD  FLT  FLUXCAL   FLUXCALERR MAG     MAGERR DQ \n')
            lines = lc.t.to_string(header=False, index=False,columns=['dummy','MJD','F',self.flux_colname,self.dflux_colname,'m','dm','Mask'])
            f.writelines(lines)
        f.close()
        print("Converted ATLAS lc to YSE format: %s" % outname)
        return(outname)

    def uploadtoyse(self,filename):
        os.system('python %s/uploadTransientData.py -e -s %s/settings.ini -i %s --instrument ACAM1 --fluxzpt 23.9 --forcedphot 1 --diffim 1' % (self.sourcedir,self.sourcedir,filename))

    def saveRADECtable(self,args,TNSname):
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
                        raise RuntimeError('Autosearch for closest bright object is not functional yet!')
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
                        cbRA = Angle(RaInDeg(self.t.loc[SNindex,'closebrightRA']),u.degree)
                        cbDec = Angle(DecInDeg(self.t.loc[SNindex,'closebrightDec']),u.degree)

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
            try:
                print('Lookbacktime in days: ',lookbacktime_days,'. Maximum lookbacktime in days: ',args.lookbacktime_days_max)
                self.lc.t = self.download_atlas_lc.get_result(RaInDeg(ra), DecInDeg(dec), token_header, lookbacktime_days=lookbacktime_days, mjd_max=args.lookbacktime_days_max) # HERE
            except MissingSchema:
                print('WARNING: NO NEW DATA AVAILABLE. Skipping cleaning, averaging, and uploading of this light curve...')
                return(1)
        else:
            print('Connecting to SSH...')
            self.download_atlas_lc.connect(args.atlasmachine,args.user,args.passwd)
            self.download_atlas_lc.get_lc(ra,dec,lookbacktime_days=lookbacktime_days)
            self.lc.t = pd.read_csv(io.StringIO('\n'.join(self.download_atlas_lc.lcraw)),delim_whitespace=True,skipinitialspace=True)
        
        mask = np.zeros((len(self.lc.t)), dtype=int)
        self.lc.t = self.lc.t.assign(Mask=mask)

        # delete all magnitudes and dmagnitudes by setting them to nan
        self.lc.t['m'] = np.nan
        self.lc.t['dm'] = np.nan
        # convert flux to magnitude
        self.lc.flux2mag(self.flux_colname,self.dflux_colname,'m','dm',zpt=23.9,upperlim_Nsigma=3)

        self.lc.t = self.lc.t.sort_values(by=['MJD'],ignore_index=True)
        indices = self.lc.ix_remove_null(colnames='uJy')

        # split the lc file into 2 separate files by filter
        for filt in ['c','o']:
            filename = self.yselcfilename(TNSname,controlindex,filt)
            fileformat = self.cfg.params['output']['fileformat']
            detections4filt=np.where(self.lc.t['F']==filt)
            newindices = AandB(indices,detections4filt)
            if len(detections4filt[0]) == 0:
                print('Saving blank light curve: %s' % filename)
                self.lc.write(filename,index=False,indices=newindices,overwrite=args.overwrite,verbose=False,columns=['MJD','m','dm',self.flux_colname,self.dflux_colname,'F','err','chi/N','RA','Dec','x','y','maj','min','phi','apfit','Sky','ZP','Obs','Mask'])
            else: 
                print('Saving light curve: %s' % filename)
                self.lc.write(filename, index=False,indices=newindices,overwrite=args.overwrite,verbose=False)
        return(0)

    def downloadYSEcontrollc(self,args,TNSname,ra,dec,pattern=None,lookbacktime_days=None):
        print('Offset status: ',args.forcedphot_offset)
        if args.forcedphot_offset == 'True':
            self.defineRADECtable(args,ra,dec,pattern=pattern)

            for i in range(len(self.RADECtable.t)):
                print(self.RADECtable.write(indices=i, columns=['ControlID', 'Ra', 'Dec']))
                result = self.downloadyselc(args,ra,dec,i,lookbacktime_days=lookbacktime_days)
                if result == 0:
                    print('Length of lc: ',len(self.lc.t))
                    self.RADECtable.t.loc[i,'Ndet']=len(self.lc.t)
                    ofilt = np.where(self.lc.t['F']=='o')
                    self.RADECtable.t.loc[i,'Ndet_o']=len(ofilt[0])
                    cfilt = np.where(self.lc.t['F']=='c')
                    self.RADECtable.t.loc[i,'Ndet_c']=len(cfilt[0])
                else:
                    return(1)
            self.saveRADECtable(args,TNSname)
        else:
            print('Skipping forcedphot offsets lc...')
            self.defineRADECtable(args,ra,dec,pattern=None)

            print(self.RADECtable.write(index=True,overwrite=False))
            for i in range(len(self.RADECtable.t)):
                result = self.downloadyselc(args,ra,dec,i,lookbacktime_days=lookbacktime_days)
                if result == 0:
                    print('Length of lc: ',len(self.lc.t))
                    self.RADECtable.t.loc[i,'Ndet']=len(self.lc.t)
                    ofilt = np.where(self.lc.t['F']=='o')
                    self.RADECtable.t.loc[i,'Ndet_o']=len(ofilt[0])
                    cfilt = np.where(self.lc.t['F']=='c')
                    self.RADECtable.t.loc[i,'Ndet_c']=len(cfilt[0])
                else:
                    return(1)
            self.saveRADECtable(args,TNSname)
        return(0)

    def averageyselc(self,args,TNSname,filt,MJDbinsize=1.0):
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
        self.averagelctable.t['Mask'] = pd.Series([], dtype=np.int64)

        # drop old columns
        dropcols=[]
        if 'Noffsetlc' in self.lc.t.columns: dropcols.append('Noffsetlc')
        for col in self.lc.t.columns:
            if re.search('^c\d_',col): dropcols.append(col)
            if re.search('^o\d_',col): dropcols.append(col)
        if len(dropcols)>0: self.lc.t.drop(columns=dropcols,inplace=True)
        
        # first flag measurements using uncertainty and x2norm cuts (cut0)
        self.lc.t['Mask'] = 0
        self.c0_PSF_uncertainty_cut(self.cfg.params['cleanlc']['cut0']['N_dflux_max'])
        self.c0_PSF_X2norm_cut(self.cfg.params['cleanlc']['cut0']['PSF_X2norm_max'])

        self.saveyselc(TNSname,0,filt=filt,overwrite=args.overwrite)

        # set maximums for flagging
        Nclip_max = self.cfg.params['upltoyse']['Nclip_max']
        Ngood_min = self.cfg.params['upltoyse']['Ngood_min']
        X2norm_max = self.cfg.params['upltoyse']['X2norm_max']

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
                df = {'MJDbin':MJD+0.5*MJDbinsize,'F':filt,'Mask':0}
                lcaverageindex = self.averagelctable.newrow(df)

            if len(mjd_ix)==0: 
                if self.verbose>1: print('No data in MJD range = 0, skipping MJD range...')
                MJD += MJDbinsize
                continue

            # sigmacut indices
            self.lc.calcaverage_sigmacutloop(self.flux_colname,noisecol=self.dflux_colname,indices=mjd_ix,verbose=1,Nsigma=3.0,median_firstiteration=True)
            fluxstatparams = deepcopy(self.lc.statparams)
            if self.verbose>1: print('Nclip: {}, Ngood: {}, X2norm: {}'.format(fluxstatparams['Nclip'],fluxstatparams['Ngood'],fluxstatparams['X2norm']))

            if fluxstatparams['mean'] is None or len(fluxstatparams['ix_good'])<1:
                if self.verbose>1: print('No measurements used, skipping MJD bin...')
                self.averagelctable.t = self.averagelctable.t.drop(index=lcaverageindex)
            else:
                # get average mjd
                self.lc.calcaverage_sigmacutloop('MJD',noisecol=self.dflux_colname,indices=fluxstatparams['ix_good'],verbose=1,Nsigma=0,median_firstiteration=False)
                    
                # add row to averagelc table
                df = {'MJD':self.lc.statparams['mean'],self.flux_colname:fluxstatparams['mean'],self.dflux_colname:fluxstatparams['mean_err'],'stdev':fluxstatparams['stdev'],'X2norm':fluxstatparams['X2norm'],'Nclipped':fluxstatparams['Nclip'],'Nused':fluxstatparams['Ngood']}
                self.averagelctable.add2row(lcaverageindex,df)

                # if any of the badflags are true, flag as bad
                badflag1 = fluxstatparams['Ngood'] < Ngood_min
                badflag2 = fluxstatparams['Nclip'] > Nclip_max
                badflag3 = not(fluxstatparams['X2norm'] is None) and (fluxstatparams['X2norm'] > X2norm_max)
                if badflag1 or badflag2 or badflag3:
                    self.averagelctable.t.loc[lcaverageindex,'Mask'] = int(self.averagelctable.t.loc[lcaverageindex,'Mask'])|self.flag_day_bad

            MJD += MJDbinsize

        # convert flux to magnitude
        self.averagelctable.flux2mag(self.flux_colname,self.dflux_colname,'m','dm',zpt=23.9,upperlim_Nsigma=3)
        if '__tmp_SN' in self.averagelctable.t.columns:
            self.averagelctable.t = self.averagelctable.t.drop(columns=['__tmp_SN'])

        self.saveyselc(TNSname,0,filt=filt,overwrite=args.overwrite,MJDbinsize=MJDbinsize)

    def uploadloop(self,args,TNSname,overwrite=True,skipdownload=False,parser=None):
        # GET RA AND DEC
        if args.ra and args.dec:
            ra = args.ra
            dec = args.dec
        elif args.tnsnamelist:
            try:
                # look in TNS
                print('Obtaining RA and Dec from TNS...')
                ra, dec = self.getradec(TNSname)
            except KeyError:
                # look in default YSE SQL query list
                print(TNSname+' not found in TNS, loading YSE list and searching for matching entry...')
                self.YSEtable.t = upltoyse.YSE_list(args)
                indices = self.YSEtable.ix_equal('Name',val=TNSname)
                if len(indices) == 0:
                    raise RuntimeError('Something went wrong: SN name '+TNSname+' not found in TNS or default YSE list!')
                else:
                    print('Entry or entries for SN name found in YSE list, obtaining RA and Dec from first index found...')
                    ra = self.YSEtable.t.loc[indices[0],'RA']
                    dec = self.YSEtable.t.loc[indices[0],'Dec']
        elif args.tnslistfilename:
            onTNSlistfile = self.checkTNSlistfile(TNSname)
            if onTNSlistfile:
                print('Obtaining RA and Dec from given TNS list file...')
                # get ra and dec from TNSlistfile
                index = self.TNSlistfile.ix_equal('TNSname',val=TNSname)
                if len(index)>0:
                    ra = self.TNSlistfile.t.loc[index[0],'RA']
                    dec = self.TNSlistfile.t.loc[index[0],'Dec']
                else:
                    raise RuntimeError('Something went wrong: TNSname does not exist!')
            else: 
                try:
                    # get ra and dec automatically from TNS, then append to TNSlistfile
                    print('Obtaining RA and Dec from TNS, then appending to given TNS list file...')
                    ra, dec = self.getradec(TNSname)
                    df = pd.DataFrame([[TNSname,ra,dec]], columns=['TNSname','RA','Dec'])
                    self.TNSlistfile.t = self.TNSlistfile.t.append(df, ignore_index=True)
                except KeyError:
                    # look in default YSE SQL query list
                    print(TNSname+' not found in TNS, loading YSE list and searching for matching entry...')
                    self.YSEtable.t = upltoyse.YSE_list(args)
                    indices = self.YSEtable.ix_equal('Name',val=TNSname)
                    if len(indices) == 0:
                        raise RuntimeError('Something went wrong: SN name '+TNSname+' not found in TNS or default YSE list!')
                    else:
                        print('Entry or entries for SN name found in YSE list, obtaining RA and Dec from first index found...')
                        ra = self.YSEtable.t.loc[indices[0],'RA']
                        dec = self.YSEtable.t.loc[indices[0],'Dec']
        else:
            # get ra and dec from yse table
            index = self.YSEtable.ix_equal('Name',val=TNSname)
            if len(index)>0:
                ra = self.YSEtable.t.loc[index[0],'RA']
                dec = self.YSEtable.t.loc[index[0],'Dec']
            else:
                raise RuntimeError('Something went wrong: TNSname does not exist!')

        # set offset pattern and lookback time in days
        if not(args.pattern is None):
            pattern = args.pattern
        else:
            pattern = self.cfg.params['forcedphotpatterns']['patterns_to_use']
        if args.lookbacktime_days:
            lookbacktime_days = args.lookbacktime_days
        else:
            lookbacktime_days = 60
        
        # download lc and, if forcedphot_offset, offset lcs
        result = 0   
        if not skipdownload:
            result = self.downloadYSEcontrollc(args,TNSname,ra,dec,pattern=pattern,lookbacktime_days=lookbacktime_days)

        if result == 0:
            # upload to YSE-PZ
            for filt in ['c','o']:
                print('### FILTER SET: ',filt)
                self.loadRADECtable(args,TNSname)
                self.loadyselc(TNSname,0,filt)
                if len(self.lc.t)>0:
                    # load single measurement light curve and calculate the average lc; also applies cut0 to single measurement light curve
                    self.averageyselc(args,TNSname,filt,MJDbinsize=args.MJDbinsize)
                    if args.averagelc:
                        if len(self.averagelctable.t)==0:
                            print('### No good averaged measurements, skipping uploading...')
                        else:
                            print('### Uploading averaged measurements for filter ',filt)
                            self.yseupload(TNSname,ra,dec,filt,MJDbinsize=args.MJDbinsize,parser=parser)
                    else:
                        print('### Uploading single measurements for filter ',filt)
                        self.yseupload(TNSname,ra,dec,filt,parser=parser)
                else:
                    print('### WARNING: empty light curve, skipping averaging and uploading...')

if __name__ == '__main__':
    upltoyse = uploadtoyseclass()

    # add arguments
    parser = argparse.ArgumentParser(conflict_handler='resolve')
    parser.add_argument('-t','--tnsnamelist', default=None, nargs='+', help='name of transients to download and upload')
    parser.add_argument('-n','--tnslistfilename', default=None, help='address of file containing TNS names, ra, and dec')
    parser.add_argument('-o','--overwrite',default=True,help='overwrite existing files')
    parser.add_argument('-i','--index', default=0, type=int, help='start at a certain index of the TNSnamelist')
    parser.add_argument('--sourcedir', default=None, help='source code directory')
    parser.add_argument('--outrootdir', default=None, help='output root directory')
    parser.add_argument('--outsubdir', default=None, help='output subdirectory')
    parser.add_argument('--api', default=False, action="store_true", help=('use API instead of SSH to get light curves from ATLAS'))
    parser.add_argument('--forcedphot_offset', default=False, action="store_true", help=("download offsets (settings in config file)"))
    parser.add_argument('--pattern', choices=['circle','box','closebright'], help=('offset pattern, defined in the config file; options are circle, box, or closebright'))
    parser.add_argument('--averagelc', default=False, action="store_true", help=('average lcs'))
    parser.add_argument('-m','--MJDbinsize', default=1.0, help=('specify MJD bin size'),type=float)
    #parser.add_argument('--skipsingleupload', default=False, action="store_true", help=('skip uploading single-measurement light curves'))
    parser.add_argument('--skipdownload', default=False, action="store_true", help=('skip downloading'))
    parser.add_argument('--ysequery', default=147, help=('enter the query number for the desired YSE list'))
    parser.add_argument('--lookbacktime_days_max',default=None, type=int, help='enter maximum lookbacktime in days for forced photometry')
    parser.add_argument('--ra', help="RA position", default=None, type=str)
    parser.add_argument('--dec', help="Dec position", default=None, type=str)

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
        print("TNSnamelist from command: \n",upltoyse.TNSnamelist)
    elif not(args.tnslistfilename is None):
        upltoyse.TNSlistfilename = '%s/%s' % (upltoyse.outrootdir,args.tnslistfilename)
        # check if TNSlistfilename exists; if not, make TNSlistfile
        if os.path.exists(upltoyse.TNSlistfilename):
            upltoyse.TNSlistfile.load_spacesep(upltoyse.TNSlistfilename, delim_whitespace=True)
        else: 
            raise RuntimeError("%s does not exist! Please create a file with the headers 'TNSname', 'RA', and 'Dec', then try again.")
        # TNSnamelist set to objects in TNSlistfilename
        upltoyse.TNSnamelist = upltoyse.TNSlistfile.t['TNSname'].values
        print("TNSnamelist from TNSlistfile: \n",upltoyse.TNSnamelist)
        print("TNSlistfilename: \n",upltoyse.TNSlistfilename)
    else:
        upltoyse.YSEtable.t = upltoyse.YSE_list(args)
        # TNSnamelist set to ojects in YSE list
        upltoyse.TNSnamelist = upltoyse.YSEtable.t['Name'].values
        print("TNSnamelist from YSE: \n",upltoyse.TNSnamelist) # change me

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
    for index in range(args.index,len(upltoyse.TNSnamelist)):
        TNSname = upltoyse.TNSnamelist[index]
        print("\nUploading and/or downloading data for %s, index %d out of %d in TNSnamelist" % (TNSname,index,len(upltoyse.TNSnamelist)-1))
        upltoyse.uploadloop(args,TNSname,overwrite=True,skipdownload=args.skipdownload,parser=parser)

    if not(args.tnslistfilename is None):
        upltoyse.TNSlistfile.write(upltoyse.TNSlistfilename,overwrite=True)

