#!/usr/bin/env python

from SNloop import SNloopclass
import pylab as matlib
import matplotlib.pyplot as plt

def dataPlot(x, y, dx=None, dy=None, sp=None, label=None, fmt='bo', ecolor='k', elinewidth=None, barsabove = False, capsize=1, logx=False, logy=False):
	if sp == None:
		sp = matlib.subplot(111)
	if dx is None and dy is None:
		if logy:
			if logx:
				plot, = sp.loglog(x, y, fmt)
			else:
				 plot, = sp.semilogy(x, y, fmt)
		elif logx:
			plot, = sp.semilogx(x, y, fmt)
		else:
			if barsabove:
				plot, dplot,dummy = sp.errorbar(x, y, label=label, fmt=fmt, capsize=capsize, barsabove=barsabove)
			else:
				plot, = sp.plot(x, y, fmt)
		return sp, plot, None
	else:
		if logy:
			sp.set_yscale("log", nonposx='clip')
		if logx:
			sp.set_xscale("log", nonposx='clip')
		plot, dplot,dummy = sp.errorbar(x, y, xerr=dx, yerr=dy, label=label, fmt=fmt, ecolor=ecolor, elinewidth=elinewidth, capsize=capsize, barsabove=barsabove)
		return sp, plot, dplot

class plotlcclass(SNloopclass):
	def __init__(self):
		SNloopclass.__init__(self)

		self.flag_uncertainty = 0x1
		self.flag_dynamic = 0x2
		self.flag_static = 0x4

	def plot_lc(self,args,SNindex,sp=None,plotoffsetlcflag=False):
		print('Plotting SN lc and offsets...')
		if sp is None:
			sp = matlib.subplot(111)
		self.loadRADEClist(SNindex)
		totalcuts = 0

		# plot SN in red first, then loop through offsets to plot in blue
		for offsetindex in range(len(self.RADECtable.t)-1,-1,-1):
			self.load_lc(SNindex, filt=self.filt, offsetindex=self.RADECtable.t.at[offsetindex,'OffsetID'])
			print('Offset index: ',offsetindex)

			makecuts_apply = self.cfg.params['plotlc']['makecuts']
			if args.skip_makecuts:
				makecuts_apply = False
			if makecuts_apply == True:
				if not('Mask' in self.lc.t.columns):
					raise RuntimeError('No "Mask" column exists! Please run "cleanup_lc.py %s" beforehand.' % self.t.at[SNindex,'tnsname'])
				lc_uJy, lc_duJy, lc_MJD, datacut = self.makecuts_indices(SNindex, offsetindex=offsetindex)
				totalcuts += datacut
			else:
				print('Skipping makecuts using mask column...')
				lc_uJy = self.lc.t['uJy']
				lc_duJy = self.lc.t['duJy']
				lc_MJD = self.lc.t['MJD']

			if offsetindex==0:
				#sp, plotSNog, dplotSNog = dataPlot(self.lc.t['MJD'],self.lc.t['uJy'],dy=self.lc.t['duJy'],sp=sp)
				#matlib.setp(plotSN,ms=5,color='r')
				sp, plotSN, dplotSN = dataPlot(lc_MJD,lc_uJy,dy=lc_duJy,sp=sp)
				matlib.setp(plotSN,ms=5,color='r')
			elif len(self.RADECtable.t)==1:
				print('No offsets, skipping plotting...')
			else: 
				sp, plotOffset, dplotOffset = dataPlot(lc_MJD,lc_uJy,dy=lc_duJy,sp=sp)
				matlib.setp(plotOffset,ms=5,color='b')
		print('Total data points cut: ',totalcuts)

		# determine legend
		if len(self.RADECtable.t)>1:
			if len(self.cfg.params['forcedphotpatterns']['circular']['radii'])==1:
				offsetlabel = '%s %s" Offset' % (self.cfg.params['forcedphotpatterns']['circular']['n'], self.cfg.params['forcedphotpatterns']['circular']['radii'][0])
			else:
				offsetlabel = '%s %s" Offset and %s %s" Offset' % (self.cfg.params['forcedphotpatterns']['circular']['n'], self.cfg.params['forcedphotpatterns']['circular']['radii'][0],self.cfg.params['forcedphotpatterns']['circular']['n'], self.cfg.params['forcedphotpatterns']['circular']['radii'][1])
			plt.legend((plotSN,plotOffset),('SN %s' % self.t.at[SNindex,'tnsname'],offsetlabel))
			plt.title('SN %s' % self.t.at[SNindex,'tnsname'])
		else:
			plt.title('SN %s' % self.t.at[SNindex,'tnsname'])

		plt.axhline(linewidth=1,color='k')
		plt.xlabel('MJD')
		plt.ylabel('uJy')
		if not(len(lc_MJD)==0):
			plt.ylim(min(lc_uJy)*1.1,max(lc_uJy)*1.1)

	def plotlcloop(self,args,SNindex):
		self.plot_lc(args,SNindex)
		plotfilename = self.lcbasename(SNindex)+'.png'
		print('Plot file name: ',plotfilename)
		plt.savefig(plotfilename)

if __name__ == '__main__':

	plotlc = plotlcclass()
	parser = plotlc.define_options()
	args = parser.parse_args()

	SNindexlist = plotlc.initialize(args)

	for SNindex in SNindexlist:
		plotlc.plotlcloop(args,SNindex)

	print('\n')
