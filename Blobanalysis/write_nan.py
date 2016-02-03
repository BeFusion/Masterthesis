#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from datetime import datetime			# for date and time printing in output file

def write_nan(fu,shift_Block, Measurement,maxlen, NewData, Radial, WriteHeader, DenCase, EmCase, BlockCase, Counterbad, blub, t_start,t_end, x,y,Refdec_ind, shift, BlobCount, Delx, tau_B, Blobf,vr, rval, vdmax, vimax, B0, Lp, q95, Te0, Ti0, ne0, omegaCi, rhoS, endtime):

	Counterbad = Counterbad +1		# Number of bad events
	blub=0
	BlobCount = np.NaN
	Delx = np.NaN
	tau_B = np.NaN
	Blobf = np.NaN
	vr = np.NaN
	vdmax = np.NaN
	vimax = np.NaN
	rval = np.NaN
	if NewData:					# Write only, if needed
	#					print('memory used:', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
		if not Radial:
			if not WriteHeader:
				print('no header is written for outputfile!')
			if not DenCase or BlockCase:
				fu = open('Emdata.txt', 'a')		# use only writing (option 'w') if needed (ATTENTION: discards the whole file, appending uses 'a')
			if DenCase:
				fu = open('Dendata.txt', 'a')		# use only writing (option 'w') if needed (ATTENTION: discards the whole file, appending uses 'a')
			if BlockCase:
				fu = open('Blockdata.txt', 'a')		# use only writing (option 'w') if needed (ATTENTION: discards the whole file, appending uses 'a')


		# data is list consisting of tubles (!) with data and format specifier
		if not BlockCase and BlobCount>=0:
			data = [("00"+str(Measurement), ""), (t_start,".2f"), (t_end-t_start, ".2f"), (x[Refdec_ind], ".2f"), (y[shift*2],".2f"), (BlobCount, "d"), (Delx, ".3f"), (tau_B,".4f"), (Blobf,".1f"), (vr, ".1f"), (rval,".1f"), (vdmax,".1f"), (vimax,".1f"), (B0,".2f"), (Lp,".2f"), (q95,".2f"), (Te0,".2f"), (Ti0,".2f"), (ne0,".2e"), (omegaCi,".2e"), (rhoS,".6f"), (endtime,".2f")]
		if not BlockCase and np.isnan(BlobCount):
			data = [("00"+str(Measurement), ""), (t_start,".2f"), (t_end-t_start, ".2f"), (x[Refdec_ind], ".2f"), (y[shift*2],".2f"), (BlobCount, ".1f"), (Delx, ".1f"), (tau_B,".1f"), (Blobf,".1f"), (vr, ".1f"), (rval,".1f"), (vdmax,".1f"), (vimax,".1f"), (B0,".2f"), (Lp,".2f"), (q95,".2f"), (Te0,".2f"), (Ti0,".2f"), (ne0,".2e"), (omegaCi,".2e"), (rhoS,".6f"), (endtime,".2f")]
		if BlockCase and BlobCount>=0:
			data = [("00"+str(Measurement), ""), (t_start,".2f"), (t_end-t_start, ".2f"), (x[Refdec_ind], ".2f"), (y[shift_Block*2],".2f"), (BlobCount, "d"), (Delx, ".3f"), (tau_B,".4f"), (Blobf,".1f"), (vr, ".1f"), (rval,".1f"), (vdmax,".1f"), (vimax,".1f"), (B0,".2f"), (Lp,".2f"), (q95,".2f"), (Te0,".2f"), (Ti0,".2f"), (ne0,".2e"), (omegaCi,".2e"), (rhoS,".6f"), (endtime,".2f")]
		if np.isnan(BlobCount) and BlockCase:
			data = [("00"+str(Measurement), ""), (t_start,".2f"), (t_end-t_start, ".2f"), (x[Refdec_ind], ".2f"), (y[shift_Block*2],".2f"), (BlobCount, ".1f"), (Delx, ".1f"), (tau_B,".1f"), (Blobf,".1f"), (vr, ".1f"), (rval,".1f"), (vdmax,".1f"), (vimax,".1f"), (B0,".2f"), (Lp,".2f"), (q95,".2f"), (Te0,".2f"), (Ti0,".2f"), (ne0,".2e"), (omegaCi,".2e"), (rhoS,".6f"), (endtime,".2f")]

		myrow="{0:%a_%d/%m/%Y_%H:%M:%S}".format(datetime.now())		# initialize format for data output (first entry is going to be the date and time)
		i=0
		for elem in data:						# go through all tuples in data and format every element
			myrow+=("{"+str(0)+":^"+ str(maxlen) +  elem[1] + "}").format(elem[0])	# format all element of data with specified format in tupel elem[1]
		myrow+="\n"

		fu.write(myrow)							# write to file
		if not Radial:
			fu.close()


		#			Set all values to zero for multiple runs:
		Delx	= 0
		tau_B 	= 0
		Blobf	= 0
		vr	= 0
		rval	= 0
		vdmax	= 0
		vimax	= 0


	return(Counterbad, blub, Delx, tau_B, Blobf,vr, rval, vdmax, vimax)