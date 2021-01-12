#!/usr/bin/env python

# S. Rest

from SNloop import SNloopclass
from statistics import median
import sigmacut
import sys
import numpy as np
import pandas as pd
from astropy.time import Time
import matplotlib.pyplot as plt
from asym_gaussian import gauss2lc

class detectbumpsclass(SNloopclass):
    def __init__(self):
        SNloopclass.__init__(self)

    def applyrollingDATETIME_OLDNOTUSED(self,SNindex,controlindex=0,MJDbinsize=1,gaussian_sigma_days=30.0,rolling_win_type='gaussian'):
        self.load_lc(SNindex,controlindex=controlindex,MJDbinsize=MJDbinsize)

        indices = self.lc.ix_unmasked('Mask',self.flag_day_bad)

        # DELME!!!!!!!!!
        self.lc.t.loc[indices,self.flux_colname]=0.0
        print()
        

        self.lc.t.loc[indices,'SNR']=self.lc.t.loc[indices,self.flux_colname]/self.lc.t.loc[indices,self.dflux_colname]

        # DELME!!!!!!!!!
        self.lc.t.loc[indices[20:30],'SNR']=1.0
        self.lc.t.loc[indices[50:90],'SNR']=2.0

        mjds = self.lc.t.loc[indices,'MJD']
        datelist = Time(mjds,format='mjd', scale='utc')
        
        temp = pd.DataFrame({'SNR':list(self.lc.t.loc[indices,'SNR'])},index=pd.to_datetime(datelist.iso))
        #temp = pd.DataFrame(self.lc.t['SNR'])

        windowsize = gaussian_sigma_days*24*3600
        windowsize_sec = 5*24*3600
        windowsize_days = gaussian_sigma_days
        print(temp)
        SNRsum = temp.rolling('%dd' % windowsize_days).sum()
        self.lc.t.loc[indices,'SNRsum']=list(SNRsum['SNR'])
        print('bbb',list(SNRsum['SNR']))
        #sys.exit(0)
        self.lc.write(indices = indices)
        
        plt.close("all")
        #self.lc.plot(indices = indices, x='MJD',y=['SNR','SNRsum'], kind="scatter")
        ax = self.lc.plot(indices = indices, x='MJD',y='SNR',kind='scatter',color='blue')
        self.lc.plot(indices = indices, x='MJD',y='SNRsum',kind='scatter',color='red',ax=ax)
        plt.show()
        sys.exit(0)
        if gaussian_sigma_days is None:
            SNRsum = temp.rolling('%ds' % windowsize, win_type=rolling_win_type).sum()
        else:
            SNRsum = temp.rolling('%ds' % windowsize, win_type=rolling_win_type).sum(std=gaussian_sigma_days)
            #SNRsum = temp.rolling('%ds' % windowsize, win_type=rolling_win_type).sum()
        #SNRsum = temp.rolling('2s').sum()
        print('cccc')
        print(SNRsum)
        temp['SNRsum']=SNRsum
        sys.exit(0)

    def define_options(self, **kwargs):
        parser = SNloopclass.define_options(self, **kwargs)
        parser.add_argument('--sim_gaussian', nargs=3, default=None, help=(' comma-separated peakMJD list, peak_appmag, gaussian_sigma): add a Gaussian at peakMJD with a peak apparent magnitude of peak_appmag and a sigma of gaussian_sigma in days.'))
        return(parser)

    def applyrolling_gaussian(self,SNindex,controlindex=0,MJDbinsize=1.0,gaussian_sigma_days=30.0,
                              simparams=None):
        
        self.load_lc(SNindex,controlindex=controlindex,MJDbinsize=MJDbinsize)

        indices = self.lc.ix_unmasked('Mask',self.flag_day_bad)
        #badindices = self.lc.ix_masked('Mask',self.flag_day_bad)
        
        # make sure there are no lingering simulations!
        dropcols=[]
        for col in ['uJysim','SNRsim','simLC','SNRsimsum']:
            if col in self.t.columns:
                dropcols.append(col)
        if len(dropcols)>0:
            self.lc.t.drop(columns=dropcols,inplace=True)
        
        self.lc.t['SNR']=0.0
        self.lc.t.loc[indices,'SNR']=(self.lc.t.loc[indices,self.flux_colname])/self.lc.t.loc[indices,self.dflux_colname]

        if not(simparams is None):
            peakMJDs = simparams['sim_peakMJD'].split(',')
            print('bbbb',peakMJDs,simparams['sim_peakMJD'])
            
            title = 'SIM Gaussian: mag: %.2f sigma-+:(%.1f,%.1f) peak MJDs:' % (simparams['sim_appmag'],simparams['sim_sigma_minus'],simparams['sim_sigma_plus'])


            # get teh simulated gaussian
            mjds = self.lc.t.loc[indices,'MJD']
            self.lc.t.loc[indices,'uJysim']=self.lc.t.loc[indices,'uJy']
            self.lc.t['simLC'] = 0.0
            for peakMJD in peakMJDs:
                peakMJD = float(peakMJD)
                print('ADDING SIM Gaussian: peak MJD: %.0f Apparent mag: %.2f sigma-: %.1f sigma+:%.1f' % (peakMJD,simparams['sim_appmag'],simparams['sim_sigma_minus'],simparams['sim_sigma_plus']))
 
                simflux = gauss2lc(mjds,peakMJD,simparams['sim_sigma_minus'],simparams['sim_sigma_plus'],app_mag=simparams['sim_appmag'])
                
                # add simflux to light curve
                self.lc.t.loc[indices,'uJysim']+=simflux

                # save the sim LC, for all MJDs
                simflux_all = gauss2lc(self.lc.t['MJDbin'],peakMJD,simparams['sim_sigma_minus'],simparams['sim_sigma_plus'],app_mag=simparams['sim_appmag'])
                self.lc.t['simLC'] += simflux_all

                title += ' %.0f' % peakMJD

            # Make sure all bad rows have SNRsim=0.0, so that they have no impact on the rolling SNRsum
            self.lc.t['SNRsim']=0.0
            # include simflux in the SNR
            self.lc.t.loc[indices,'SNRsim']=(self.lc.t.loc[indices,'uJysim'])/self.lc.t.loc[indices,self.dflux_colname]

        
        gaussian_sigma = round(gaussian_sigma_days/MJDbinsize)
        windowsize = int(3 * gaussian_sigma * 2)
        halfwindowsize = int(windowsize*0.5)+1
        print('sigma(days)=%.2f, MJD binsize=%.2f, sigma(bins)=%d, window size(bins)=%d' % (gaussian_sigma_days,MJDbinsize,gaussian_sigma,windowsize))

        dataindices = np.array(range(len(self.lc.t))+np.full(len(self.lc.t),halfwindowsize))

        # Calculate the rolling SNR sum
        temp = pd.Series(np.zeros(len(self.lc.t)+2*halfwindowsize),name='SNR', dtype=np.float64)
        temp[dataindices] = self.lc.t['SNR']
        SNRsum = temp.rolling(windowsize,center=True,win_type='gaussian').sum(std=gaussian_sigma)
        self.lc.t['SNRsum']=list(SNRsum.loc[dataindices])

        if not(simparams is None):
            # Calculate the rolling SNR sum for SNR including simflux
            temp = pd.Series(np.zeros(len(self.lc.t)+2*halfwindowsize),name='SNRsim', dtype=np.float64)
            temp[dataindices] = self.lc.t['SNRsim']
            SNRsimsum = temp.rolling(windowsize,center=True,win_type='gaussian').sum(std=gaussian_sigma)
            self.lc.t['SNRsimsum']=list(SNRsimsum.loc[dataindices])
            
        # maginfo is to add to filenames!
        maginfo='.sim%.0fmag' % (simparams['sim_appmag'])

        # save lc
        self.save_lc(SNindex=SNindex,controlindex=controlindex,filt=self.filt,MJDbinsize=MJDbinsize,addsuffix=maginfo,overwrite=True)

        
        outbasefilename = self.lcbasename(SNindex=SNindex,controlindex=controlindex,MJDbinsize=MJDbinsize,addsuffix=maginfo)

        plt.close("all")
        ax1 = self.lc.plot( x='MJDbin',y='uJysim',yerr='duJy',kind='scatter',color='cyan',legend=True,title=title)
        self.lc.plot(x='MJDbin',y='uJy',yerr='duJy',kind='scatter',color='red',ax=ax1)
        self.lc.plot(x='MJDbin',y='simLC',color='cyan',ax=ax1)
        outfile='%s.simLC.png' % outbasefilename
        print('Saving ',outfile)
        plt.savefig(outfile)

         
        ax2 =  self.lc.plot(indices=indices,x='MJDbin',y='SNRsim',kind='scatter',color='cyan',title=title)
        self.lc.plot(indices=indices,x='MJDbin',y='SNR',kind='scatter',color='red',ax=ax2)
        self.lc.plot(x='MJDbin',y='SNRsimsum',ax=ax2,color='cyan')
        self.lc.plot(x='MJDbin',y='SNRsum',ax=ax2,color='red')
        outfile='%s.simSNR.png' % outbasefilename
        print('Saving ',outfile)
        plt.savefig(outfile)

    def detectbumpsloop(self,SNindex,MJDbinsize=1.0,simparams=None):
        print('###################################\Detecting Bumps ...\n###################################')
        
        # loop through SN and control lcs
        for controlindex in range(0,len(self.RADECtable.t)):    
            # stop loop if only SN should be done
            if (not self.cfg.params['detectBumps']['apply2offsets']) and (controlindex>0):
                break

            # average the light curve by MJDbinszie
            self.applyrolling_gaussian(SNindex,controlindex=controlindex,MJDbinsize=MJDbinsize,
                              gaussian_sigma_days=self.cfg.params['detectBumps']['gaussian_sigma'],
                              simparams=simparams)

        print('### Averaging LCs done')


if __name__ == '__main__':

    detectbumps = detectbumpsclass()
    parser = detectbumps.define_options()
    args = parser.parse_args()
    #print(pd.__version__)
    #sys.exit(0)

    SNindexlist = detectbumps.initialize(args)
    
    
    if not(args.sim_gaussian is None):
        simparams = {'sim_peakMJD':args.sim_gaussian[0],
                     'sim_appmag':float(args.sim_gaussian[1]),
                     'sim_sigma_minus':float(args.sim_gaussian[2]),
                     'sim_sigma_plus':float(args.sim_gaussian[2])}
    else:
        simparams=None

    for SNindex in SNindexlist:
        print('Averaging lc for ',detectbumps.t.at[SNindex,'tnsname'],', index %i/%i' % (SNindex,len(detectbumps.t)))
        detectbumps.loadRADEClist(SNindex)
        detectbumps.detectbumpsloop(SNindex,MJDbinsize=args.MJDbinsize,simparams=simparams)