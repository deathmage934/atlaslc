#!/usr/bin/env python

import sys
from download_lc_loop import downloadlcloopclass

if __name__ == '__main__':
	downloadlc = downloadlcloopclass()
	parser = downloadlc.download_atlas_lc.define_optional_args()
	parser = downloadlc.define_options(parser=parser)
	args = parser.parse_args()
	
	SNindexlist = downloadlc.initialize(args)
	if len(SNindexlist)<1:
		print('No matching SN found, exiting!!')
		sys.exit(0)

	downloadlc.download_atlas_lc.connect(args.atlasmachine,'sofia','Starswirl1410!@')
	for SNindex in SNindexlist:
		downloadlc.downloadoffsetlc(SNindex,
				lookbacktime_days=args.lookbacktime_days,
				savelc=args.savelc,
				overwrite=args.overwrite,
				fileformat=args.fileformat,
				pattern=args.pattern,
				forcedphot_offset=args.forcedphot_offset)
		for offsetindex in range(len(downloadlc.RADECtable.t)):
			downloadlc.cleanuplcloop(args,SNindex,offsetindex=offsetindex)
		if args.plot: 
			downloadlc.plotlcloop(args,SNindex)
		if args.averagelc: 
			downloadlc.averagelcloop(args,SNindex)