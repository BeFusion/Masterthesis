#!/usr/bin/env python
# -*- coding: utf-8 -*-


########### Importing modules ###########################################################################

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#from sympy import *
from pylab import plot
import matplotlib.gridspec as gridspec		# define layout of grid
from mpl_toolkits.axes_grid1 import make_axes_locatable # for colorbars
from scipy import ndimage			# calculate com
from scipy.interpolate import spline		# for interpolate COM position for good speed calculation of blob
from scipy import stats, polyfit		# for linear interpolation data
from hdf5_read import Li_read_hdf5		# read hdf5 file in
from block_func import block_beam_func		# for calculating the block data
from write_nan import write_nan			# for writing in bad cases
from dect_func import dect_func
from Blob_params import *			# set a few parameters 
import timeit					# for timekeeping
from datetime import datetime			# for date and time printing in output file
import sys					# for exiting the program if errors occur
import pdb					# for error tracking (set_trace at right position)
import resource
#set for error tracking:
#pdb.set_trace()

#set for time keeping:
#tic1 = timeit.default_timer()
#tic2= timeit.default_timer()
#print('time needed for ...:', tic2-tic1)


###############################################################################################################################
########### COMMENTS ##########################################################################################################
###############################################################################################################################
#
# some parameters are set in Blob_params
#
# exit for bad results via emergency exit blub = 999
#
#
#
#
#
#

###############################################################################################################################
########### READING DATA AND PRELIMINARY settings #############################################################################
###############################################################################################################################


# Select cases and switches ###################################################################################################

# Evaluates the 'infinite' resolved beam emission, the density evaluation or the Block evaluation case

Measurement=4				# select number of measurment to treat (is used for hdf-selection etc)
# Case Loop
Counterbad = 0				# Number of bad events
for Case in range (3):
	if Case == 0:
		EmCase = True				# infinite emission case --> actually always true, since every different setting during analysis is made by if-statements with the other two cases --> can stay false
		DenCase = False				# density analysis case
		BlockCase = False			# block-averaging emission case
	if Case == 1:
		EmCase = False				# infinite emission case --> actually always true, since every different setting during analysis is made by if-statements with the other two cases --> can stay false
		DenCase = True				# density analysis case
		BlockCase = False			# block-averaging emission case
	if Case == 2:
		EmCase = False				# infinite emission case --> actually always true, since every different setting during analysis is made by if-statements with the other two cases --> can stay false
		DenCase = False				# density analysis case
		BlockCase = True			# block-averaging emission case


	if (EmCase and DenCase or EmCase and BlockCase or DenCase and BlockCase):		#if you do something wrong here ;)
		sys.exit('Accidentally set two cases for evaluation!')


	yvariation = True		# Calculate more values for y-variation
	
	Standardblob= True		# switch for calculations

	NewData = True			# Switch to write output

	Radial = True
	
	blub = 0			# emergency exit for bad results

	# Writing head of file for radial case ########################################################################################
	if Radial:	

		WriteHeader = True 

		if not WriteHeader:
			print('no header is written for outputfile!')
		if EmCase and NewData:
			fu = open('{0:03d}Emresults.txt'.format(Measurement), 'wb')		# use only writing (option 'w') if needed (ATTENTION: discards the whole file, appending uses 'a')
			print('results are stored in {0:03d}Emresults.txt'.format(Measurement))
		if DenCase and NewData:
			fu = open('{0:03d}Denresults.txt'.format(Measurement), 'wb')		# use only writing (option 'w') if needed (ATTENTION: discards the whole file, appending uses 'a')
			print('results are stored in {0:03d}Denresults.txt'.format(Measurement))
		if BlockCase and NewData:
			fu = open('{0:03d}Blockresults.txt'.format(Measurement), 'wb')		# use only writing (option 'w') if needed (ATTENTION: discards the whole file, appending uses 'a')	
			print('results are stored in {0:03d}Blockresults.txt'.format(Measurement))
	#	Putting everything into lists and pre-processing for output ###########################################################
		# Collumn index print for better human readability
		if WriteHeader and NewData:
			pre=[None]*23
			countblub=0
			for i in range(len(pre)):
				pre[i]= countblub
				countblub=countblub+1

			pre[0]="Col.-Index:  0"

			# header is list consisting of header titles:
			header = ['Date & Time', 'Measurement #', 'Starting time [ms]', 'Observed time [ms]', 'Refdec-pos[cm]', 'y-pos of beam [cm]', '# of Blobs', 'HWHM of Blob [cm]', 'tau_B [ms]','Blobfrequency [1/s]','Average speed vr [m/s]','Std. error vr [m/s]', 'max speed vdmax [m/s]', 'max i-pol. vimax [m/s]', 'B0 [T]', 'Lp [m]', 'q95 ', 'Te0 [eV]', 'Ti0 [eV]', 'ne0 [m-3]', 'omegaCi [1/s]', 'rhoS [m]', 'endtime [ms]']

			# determine maximum length of strings in header to determine space needed for each collumn
			maxlen = 0
			for elem in header:
				strlen = len(elem)
				if (strlen>maxlen):
					maxlen=strlen
			maxlen=maxlen+3							# safty distance from next collumn +x 
				
		#	Formating everything ########################################################################################################
		
			row_format=""							# initialize string
			for elem in range(len(header)):					# go through all elements in data-list
				row_format+="{"+ str(elem) + ":^"+ str(maxlen) +"}"	# format every element for header and pre (str(elem)) accoringly (centering:^) with space of maxlen

			row_format+= "\n"						# add new line
				
			fu.write(row_format.format(*pre))				# writes collumn index
			fu.write(row_format.format(*header))				# writes header



	# data reading  ################################################################################################################


									# will be used for file search and in output
	print('Data is going to be read from file: /hesel/data/ASDEX.{0:03d}.h5'.format(Measurement))

	# read in data from function Li_read_hdf5: Specify correct hdf5-file
	time, x, y, ne, LiTemp, Li2p1, B0, Lp, q95, Te0, Ti0, ne0, omegaCi, rhoS, endtime = Li_read_hdf5('/pfs/home/bschiess/public/hesel/data/ASDEX.{0:03d}.h5'.format(Measurement))

	# other settings and preparations  #############################################################################################

	# selected times, which are supposed to be investigated [us]
	t_start=0.2
	t_end = time[-1]
	print('time to be investigated betweem %.3gms and %.3gms:' % (t_start, t_end))

	# calculate time in [us] for snapshots:
	point=position*(time[11]-time[10])

	# some lengthes used during the program:
	mn=len(x)
	x_axis = x					# just useful sometimes to avoid confusion
	stepsize=len(y)					# misleading name, sorry
	timestep=len(time)				# misleading name, sorry
	avy=y[stepsize-1]/stepsize			
	avt=time[timestep-1]/timestep

	# some matrices we need copies from during the run:

	if DenCase:
		Li2p1Z = ne.copy()				 # copy matrix for use in calculation of Lip1Z (the zero-mean-fluctuation)
		LiBlobEm = ne.copy()				 # copy Li-Matrix for BlobCount use
	if not DenCase:
		Li2p1Z = Li2p1.copy()				 # copy matrix for use in calculation of Lip1Z (the zero-mean-fluctuation)
		LiBlobEm = Li2p1.copy()				 # copy Li-Matrix for BlobCount use

	if BlockCase:		# stuff which is always overwritten inside of Block-Case-loop and which has to be copied in...
		xcop = x.copy()
		Li2p1cop = Li2p1.copy()
		mncop = mn
	#define center of beam and shift for correct beam position:
	if yvariation:
		beamrange = np.arange (y_starting,y_ending,yResolution)
	else: 
		beamrange = [bcenter]
	for bcenter in beamrange:
		shift=int(bcenter/y[-1]*stepsize/2)+1		# shift of beam center in order to place the center at correct position
		shift2 = int(b/2/y[-1]*stepsize)		# shift*2-shift2 is now left index of beam edge and shift*2+shift2 the right one
		print('Position of Beam: {0:.2f}cm'.format(y[shift*2]))
		print('Beamwidth: {0:.2f}cm'.format(y[shift2+shift2]))	
			

		# use the block-emission data for the block case:
		if BlockCase:
			x_ori = xcop
		#	Li2p1 = block_func(Li2p1,timestep,mn,b,avy,shift,shift2)			# if original block-function should be used

			Li2p1_beam, Li2p1, x, blur = dect_func(Li2p1cop,timestep,mncop,b,avy,shift,shift2,xcop)	# if detector reduced block function should be used

			x_axis = xcop
			mn = len(x)
			shift_Block = shift
			shift = 0									# we are in the 1D-Case, so every time the entry at shift is called in matrix, this is set to the zero entry


		#Create zero-mean fluctuation matrix

		if not DenCase:
#			Li2p1Z = Li2p1.copy()				 # copy matrix for use in calculation of Lip1Z (the zero-mean-fluctuation)

			if BlockCase:					 # we are only treating the 1D-case here, since it has been averaged already in y-direction
				Li2p1Z = Li2p1.copy()
				mean = [np.NaN]*mn			 # create vector for mean in correct length
				mean = np.array(mean)			 # convert to array
				mean = mean.reshape(1,mn)		 # reshape to matrix for x- and y-dimension (y-dim is 1 in this case!) 
				for x_ind in range (0,mn):		 
					mean[0,x_ind] = np.mean(Li2p1[:,0,x_ind])
			if not BlockCase:
				mean = [np.NaN]*stepsize*mn		# create vector for mean in correct length
				mean = np.array(mean)			# convert to array
				mean = mean.reshape(stepsize,mn)	# reshape to matrix for x- and y-dimension 
				for y_ind in range (0, stepsize):	# do for loop over x-and y-positionand calculate the zero-mean-matrix (average over all timesteps)
					for x_ind in range (0,mn):		 
						mean[y_ind,x_ind] = np.mean(Li2p1[:,y_ind,x_ind])
					
			for t_ind in range (int(t_start/avt),int(t_end/avt)):
				Li2p1Z[t_ind,:,:]=Li2p1[t_ind,:,:]-mean  # subtract mean value for every timestep (takes most of the time for the program!)

		if DenCase:
#			Li2p1Z = ne.copy()				 # copy matrix for use in calculation of Lip1Z (the zero-mean-fluctuation)
			mean = [None]*stepsize*mn			 # create vector for mean in correct length
			mean = np.array(mean)				 # convert to array
			mean = mean.reshape(stepsize,mn)		 # reshape to matrix for x- and y-dimension 
			for y_ind in range (0, stepsize):		 # do for loop over x-and y-positionand calculate the zero-mean-matrix (average over all timesteps)
				for x_ind in range (0,mn):		 
					mean[y_ind,x_ind] = np.mean(ne[:,y_ind,x_ind])
					
			for t_ind in range (int(t_start/avt),int(t_end/avt)):
				Li2p1Z[t_ind,:,:]=ne[t_ind,:,:]-mean  # subtract mean value for every timestep (takes most of the time for the program!)


		# Calculate index and position of Reference detector on xaxis (indices are counted backwards x[len(x)-1]<0!)


		if Radial:
			if BlockCase:
				Startrange = int(len(x)/5)
				Resolution = 1
				
			if not BlockCase:
				Startrange = int(1.5/0.02)

			for Refdec_ind in range(Startrange, len(x),Resolution):			# Wall-region can be neglected (approximately 1.5-2cm) --> go through Refdec-Positions

				if not BlockCase:					# otherwise this position is determined later
					print('Position along beam-axis for Reference Detector: {0:.2f}cm'.format(x[Refdec_ind])) # check, if Reference detector is at correct position
					if (Refdec>x[0]):				# exit, if error with reference detector position occurs
						sys.exit('Reference Detector outside of observable region')

				if BlockCase:
					print('Reference detector has to be shifted to a position available for beam evaluation!')
					print('Position along beam-axis for Reference Detector in Block-Case: {0:.2f}cm'.format(x[Refdec_ind])) # check, if Reference detector is at correct position
					print('Position of Reference Detector set in Parameter-file: {0:.2f}cm'.format(Refdec)) # check, if Reference detector is at correct position
					if (Refdec>x[0]):				# exit, if error with reference detector position occurs
						sys.exit('Reference Detector outside of observable region')
				#Meshgrids
				X,Y=np.meshgrid(x,y)						# transform x,y-axis in matrix for contourf-plot of 2D-plots
				if BlockCase:
					X_Ori,Y = np.meshgrid(x_ori,y)				# only needed, if BlockCase changes x-axis
				X2,T2=np.meshgrid(x,time[int(t_start/avt):int(t_end/avt)]) 	# Meshgrid for timeplots:



				###################################################################################################################################
				########### FIGURE 1 ##############################################################################################################
				###################################################################################################################################

				if not Radial:
					# create fiugre with certain ratio and facecolor:
					
		
					f1=plt.figure(figsize=(8,8), facecolor = 'white')
					gs = gridspec.GridSpec(4, 1,				# ratio of grid space (2 plots per collumn, 3 per row)
							width_ratios=[1],		# width ratios of the 3 plots per row
							height_ratios=[1,1,1,1]		# height ratios of the 2 polots per collumn
							)
					

				
					#define limits for plots
					xmin=0.0
					xmax=2.2
					ymin=0
					ymax=8


					# create axes with the ratios specified in gs above. ax5, ax6 are the axes for the colorbars:
					ax1 = plt.subplot(gs[0])
					ax2 = plt.subplot(gs[1])
					ax3 = plt.subplot(gs[2])
					ax4 = plt.subplot(gs[3])

					# create a list in for simple modifications on all plots at the same time:
					allax=[ax1,ax2,ax3,ax4]


					# AX1: Density PLOT ############################################################################################################

					# first subplot: density
					if BlockCase:
						density = ax1.contourf(Y,X_Ori,ne[position,:,:],300)
					if not BlockCase:
						density = ax1.contourf(Y,X,ne[position,:,:],300)
					#ax1.invert_yaxis()
					ax1.axvline(bcenter,color='w', linestyle='-')
					ax1.set_xlabel('axis $y$ (cm)', labelpad= -10)			# switched of, since axis in row below
					ax1.set_ylabel('beam axis $x$ (cm)')
					ax1.set_title('input density $n_e$ (cm$^{-3}$) at %.3g ms' % (point))
					#ax1.set_xlim(xmin,xmax)
					#ax1.set_ylim(ymin,ymax)
					div1 =  make_axes_locatable(ax1)
					cax1=div1.append_axes("right",size="5%",pad=0.05)
					chbar_den = plt.colorbar(density,cax=cax1)
					chbar_den.set_label('density $n_e$ (cm$^{-3}$)')

					# AX2:emission profile #########################################################################################################

					if not DenCase:
						if BlockCase:
							emission2 = ax2.contourf(Y,X_Ori,Li2p1_beam[position,:,:],300)
						if not BlockCase:
							emission2 = ax2.contourf(Y,X,Li2p1[position,:,:],300)
						#ax2.invert_yaxis()
						ax2.axvline(bcenter,color='w', linestyle='-')
						ax2.set_ylabel('beam axis $x$ (cm)')
						ax2.set_xlabel('axis $y$ (cm)', labelpad=-10)
						ax2.set_title('beam emission spectrum at %.3g ms'% (point))
						#ax2.set_xlim(xmin,xmax)
						#ax2.set_ylim(ymin,ymax)
						div2 =  make_axes_locatable(ax2)
						cax2=div2.append_axes("right",size="5%",pad=0.05)
						cbar_emission2 = plt.colorbar(emission2,cax=cax2)
						cbar_emission2.set_label('Li$_{2p}$ emission (a.u.)')

					if DenCase:
						emission2 = ax2.contourf(Y,X,Li2p1Z[position,:,:],300)
						#ax2.invert_yaxis()
						ax2.axvline(bcenter,color='w', linestyle='-')
						ax2.set_ylabel('beam axis $x$ (cm)')
						ax2.set_xlabel('axis $y$ (cm)', labelpad=-10)
						ax2.set_title('Density fluctuation at %.3g ms'% (point))
						#ax2.set_xlim(xmin,xmax)
						#ax2.set_ylim(ymin,ymax)
						div2 =  make_axes_locatable(ax2)
						cax2=div2.append_axes("right",size="5%",pad=0.05)
						cbar_emission2 = plt.colorbar(emission2,cax=cax2)
						cbar_emission2.set_label('Density fluctuation $\delta n_e$ (cm$^{-3}$)')


					#plot of emission along beam over time
					#print(np.shape(T2), np.shape(X2), np.shape(Li2p1Z))


					# AX3: Plot of emission along beam over time ####################################################################################
					timeemission = ax3.contourf(T2,X2,Li2p1Z[int(t_start/avt):int(t_end/avt),shift*2,:],300)
					#ax3.invert_yaxis()
					ax3.axhline(x[Refdec_ind],color='w',linestyle='-')
					ax3.axvline(point,color='w', linestyle='-')
					ax3.set_ylabel('beam axis $x$ (cm)')			
					ax3.set_xlabel('time $t$ (ms)', labelpad=-10)
					if not DenCase:
						ax3.set_title('timeseries of emission fluctuation for beam position %.3g cm' % (bcenter))
					if DenCase:
						ax3.set_title('timeseries of density fluctuation for beam position %.3g cm' % (bcenter))
					ax3.set_xlim(t_start,t_end)
					#ax3.set_ylim(ymin,ymax)
					div3 =  make_axes_locatable(ax3)
					cax3=div3.append_axes("right",size="5%",pad=0.05)
					cbar_timeemission = plt.colorbar(timeemission,cax=cax3)
					if not DenCase:
						cbar_timeemission.set_label('Li$_{2p}$ emission fluctuation (a.u.)')
					if DenCase:
						cbar_timeemission.set_label('density $n_e$ (cm$^{-3}$)')



				# Calculations for conditional average ###########################################################################################

				# COMMENT: Settings for threshold etc are done in the top
				if not DenCase:
					LiAv = np.mean(Li2p1[int(t_start/avt):int(t_end/avt),shift*2,Refdec_ind])		# Average of Li value
					LiSig = np.std(Li2p1[int(t_start/avt):int(t_end/avt),shift*2,Refdec_ind])		# Standard deviation of Li value
				if DenCase:
					LiAv = np.mean(ne[int(t_start/avt):int(t_end/avt),shift*2,Refdec_ind])			# Average of Li value
					LiSig = np.std(ne[int(t_start/avt):int(t_end/avt),shift*2,Refdec_ind])			# Standard deviation of Li value
				LiCon = LiAv + ThresCon*LiSig									# Set condition
				TWind = TWindMult*avt							# Time window for conditional averaging [us] (Gregor's Publication: 500 us)

				if not DenCase:
					print('Average Emission Value at Reference Detector: {0:.4f} (a.u.) '.format(LiAv))
					print('Corresonding Standard Deviation: {0:4f} (a.u.)'.format(LiSig))
				if DenCase:
					print('Average Density at Reference Detector: {0:.2e}cm-3'.format(LiAv))
					print('Corresonding Standard Deviation: {0:.2e}cm-3'.format(LiSig))
				print('time window for conditional averaging is {0:.3f}ms consisting of {1:d} timesteps'.format(TWind, int(TWind/avt)))
				print('timewindow, in which too events are counted as one: {0:.3f}'.format(WinBar*avt))


				# Calculate number of blobs in time period:
				if BlockCase:
					LiBlobEm = Li2p1.copy()						# copy Li-Matrix for BlobCount use
				#if DenCase:
					#LiBlobEm = ne.copy()						# copy density-Matrix for BlobCount use

				# Change matrix to 0s and 1s (0 below, 1 above threshold) for finding the position of blob occurance
				for m in range (int(t_start/avt),int(t_end/avt)):					# go through every time point
					if LiBlobEm[m,shift*2,Refdec_ind]<LiCon:		# go in, if emission is below threshold
						LiBlobEm[m,shift*2,Refdec_ind]=0		# set all values to zero, if it is below the threshold
					else:							# go in, if emission is above threshold
						LiBlobEm[m,shift*2,Refdec_ind]=1		# set all values to 1, if it is above threshold
					m=m+1

				# Count Blobs at first time position, when threshold is passed, with condtion, that in LiBlobEm there should occur the event 0 1 1 or 0 1 0:
				BlobCount=0					# Counter counting the number of blobs
				remind=[np.NaN]*200				# index used for time window in order to not double count an index
				for m in range (int(t_start/avt),int(t_end/avt)-1):					# go through every time point
					if LiBlobEm[m,shift*2,Refdec_ind]>0 and LiBlobEm[m+1,shift*2,Refdec_ind]>0 and LiBlobEm[m-1,shift*2,Refdec_ind]==0: # condition for counting an event as a blob (0 1 1)
						BlobCount=BlobCount+1							# raise blob count
						if (remind[BlobCount-1]>=m-WinBar):					# safty time intervall to not count one blob as two
							BlobCount=BlobCount-1
						remind[BlobCount]=m
					if LiBlobEm[m,shift*2,Refdec_ind]>0 and LiBlobEm[m+1,shift*2,Refdec_ind]==0 and LiBlobEm[m-1,shift*2,Refdec_ind]==0:
						BlobCount=BlobCount+1							# raise blob count
						if (remind[BlobCount-1]>=m-WinBar):					# safty time intervall to not count one blob as two
							BlobCount=BlobCount-1
						remind[BlobCount]=m

				print('Number of Blobs in observed time period: ', BlobCount)

				if not Radial:
					# AX4: Plot of emission at Reference detector position #############################################################################

					if not DenCase:
						emission3 = ax4.plot(T2,Li2p1[int(t_start/avt):int(t_end/avt),shift*2,Refdec_ind])
					if DenCase:
						emission3 = ax4.plot(T2,ne[int(t_start/avt):int(t_end/avt),shift*2,Refdec_ind])
					ax4.axhline(LiAv,color='r', linestyle='-')
					ax4.axhline(LiCon,color='g', linestyle='-')
					ax4.axvline(point,color='k', linestyle='-')
					if not DenCase:
						ax4.set_ylabel('Li$_{2p}$ emission (a.u.)')
					if DenCase:
						ax4.set_ylabel('density $n_e$ (cm$^{-3}$)')
					ax4.set_xlabel('time $t$ (ms)',labelpad=-10)
					ax4.set_title('Series for Detector Position %.2g cm' % (Refdec))
					ax4.set_xlim(t_start,t_end)
					#ax4.set_ylim(ymin,ymax)
					div4 =  make_axes_locatable(ax4)
					cax4=div4.append_axes("right",size="5%",pad=0.05)

					for ylabel_i in cax4.get_yticklabels():
						ylabel_i.set_visible(False)
						ylabel_i.set_fontsize(0.0)
					for xlabel_i in cax4.get_xticklabels():
						xlabel_i.set_fontsize(0.0)
						xlabel_i.set_visible(False)


					##################################################################################################################################


				print('Some Dimensions to check shapes etc:')
				print('Li2p1, time, Y,X', np.shape(Li2p1), np.shape(time), np.shape(Y), np.shape(X))


				###################################################################################################################################
				########### Standard-Blob analysis (part 1) #######################################################################################
				###################################################################################################################################

				Delx = 0		# some stuff set to zero for case, that there are no blobs detected for first refdec-position
				tau_B = 0
				Blobf = 0
				vr = 0
				rval = 0
				vdmax = 0
				vimax = 0

				if (Standardblob==True and BlobCount>0):				# only occurs if switch is on and if there is certain number of blobs


				# 	Find emission values for standard blobs ###################################################################################

					Count2=0							# set blobcount to zero (now count again and find window around)
					LiConBlo = np.empty((BlobCount,int(TWind/avt),len(x)))		# Creating an empty, corectly shaped matrix for all blobs and their intensity values (later) [left time shift, right time shift, blobindex] and fill the seccond dim of the all-blob-matrix with zeros
					LiConBlo.fill(0)
					LiConAv = np.empty((int(TWind/avt), len(x)))			# array to hold values for standard blob later
					LiConAv.fill(0)
					print('LiConAv', np.shape(LiConAv))

					# find the time-windows of every blob event and put them into a matrix ######################################################
					remind2 = [np.NaN]*200			# reminding vector holding the indizes of blobs (maximum 200 blobs)
					intermediatestore = np.NaN
					for m in range (int(t_start/avt),int(t_end/avt)-1-int(TWind/avt)):					# go through every time point until one TWind before the edge --> a possible blob at the right edge-blob is not used
						if LiBlobEm[m,shift*2,Refdec_ind]>0 and LiBlobEm[m+1,shift*2,Refdec_ind]>0 and LiBlobEm[m-1,shift*2,Refdec_ind]==0 or \
						LiBlobEm[m,shift*2,Refdec_ind]>0 and LiBlobEm[m+1,shift*2,Refdec_ind]==0 and LiBlobEm[m-1,shift*2,Refdec_ind]==0: # condition for counting an event as a blob (0 1 1) or (0 1 0)
							Count2=Count2+1
							if (remind2[Count2-1]>=m-WinBar):								# safty time intervall to not count one blob as two
								intermediatestore = remind2[Count2-1]							# just for the entrance-condition to calculate emission time intervall
								Count2=Count2-1
							if not intermediatestore>=m-WinBar:								# enter only if that blob is not part of the first blob
								for l in range (0, int(TWind/avt)):							# go through time window
									for k in range (0, len(x)):							# go through all x-positions within time window
										LiConBlo[Count2-1,l,k]=Li2p1Z[(m-int(TWind/avt/2)+l),int(shift*2),k]	# store emission values within time window in LiConBlo
								print('Blob # {0} can be found at t={1:.2f}ms'.format(Count2,time[m]))
							remind2[Count2]=m

					# calculate average over all blob events --> Standard Blob in LiConAv #######################################################
					for m in range (0,len(x)):
						for l in range (0, int(TWind/avt)):
							LiConAv[l,m] = np.mean(LiConBlo[:,l,m])	 # calculate the values for the standard blob and put in array
							l=l+1
						m=m+1

					Delt=[None]*int(TWind/avt)				# create emty list for standard blob time window

					NullPos_ind = 0						# initialize to be zero (just for checking afterwards, that it is not any more
					for m in range (0,int(TWind/avt)):			# create a vector containing the times of the time window for standard blob
						Delt[m]=-TWind/2+avt*m
						if (Delt[m]==0.0):
							NullPos_ind = m				# NullPos_ind is index, used for calculating blobwidth at HWHM#
						m=m+1
					if not Radial:
						BloX, BloT = np.meshgrid(x,Delt)			# new meshgrid for standard blob
					if (NullPos_ind==0):
						sys.exit('The parameter for time window was not chosen to be even, NullPos_ind is not exact')

				# 	Define Time Window of Standard Blob Analysis ##############################################################################
					
					LookWinMin = Delt[0]
					LookWinMin_ind = 0
					LookWinMax = Delt[-1]
					LookWinMax_ind = len(Delt)-1



				# 	interpolate data for block-case evaluation for calculation of FWHM (correlation time)
					if BlockCase:
						x_inter = np.linspace(x[-1],x[0],300)
						x_here=[np.NaN]*len(x)							# introduce new x-axis going from negative to positive values for spline interpolation
						LiConAv_here_a = [np.NaN]*len(x)
						LiConAv_here_b = [np.NaN]*len(x)
						LiConAv_here_c = [np.NaN]*len(x)
						LiConAv_here_Null = [np.NaN]*len(x)
						for i in range(len(x)):
							x_here[len(x)-i-1]=x[i]
							LiConAv_here_a[len(x)-i-1] = LiConAv[a_an,i]
							LiConAv_here_b[len(x)-i-1] = LiConAv[b_an,i]
							LiConAv_here_c[len(x)-i-1] = LiConAv[c_an,i]
							LiConAv_here_Null[len(x)-i-1] = LiConAv[NullPos_ind,i]
						LiConAv_smooth_a = spline(x_here,LiConAv_here_a,x_inter)
						LiConAv_smooth_b = spline(x_here,LiConAv_here_b,x_inter)
						LiConAv_smooth_c = spline(x_here,LiConAv_here_c,x_inter)
						LiConAv_smooth_Null = spline(x_here,LiConAv_here_Null,x_inter)

				# 	Calculate HWHM of standard blob (more work, since curve not necessarily gausian!) ##########################################

					if BlockCase:
						LiConAv_max = max(LiConAv_smooth_Null[:])		# Maximum emission value
						for c in range (0,len(x_inter)):
							if (LiConAv_max==LiConAv_smooth_Null[c]):
								max_ind = c				# index on x-axis of maximum emission value

						# find closest index to half-max-position (lower index)	
						PeakC = [None]*10					# peak counter to determine for sure the position of half-max
						RandCount = 1
						for d in range (max_ind):
							if (LiConAv_smooth_Null[d]>LiConAv_max/2 and LiConAv_smooth_Null[d-1]<LiConAv_max/2):
								PeakC[RandCount]=d			# store index in peak counter array (in best case only one entry peakC[1]!)
								RandCount = RandCount+1
						low_index = max(PeakC)					# lower index is maximum for last time it surpasses HM

						# find closest index to half-max-position (upper index)
						PeakC = [None]*10					# peak counter array to determine for sure the position of half-max
						PeakC = np.array(PeakC)
						PeakC.fill(10000)					# fill with high values (higher than len(x), since we need the lowest one
						RandCount = 1
						for e in range (max_ind+1,len(x_inter)):
							if (LiConAv_smooth_Null[e-1]>LiConAv_max/2 and LiConAv_smooth_Null[e]<LiConAv_max/2):
								PeakC[RandCount]=e			# store index in peak counter array (in best case only one entry peakC[1]!)
								RandCount = RandCount+1	
						up_index = min(PeakC)					# up index is minimum index for first time it surpasses HM
						
						if (up_index>5000 or low_index<1):			# exit loop for this analysis if a proper calculation is not possible here
							blub=999
							print(blub, 'in up/low_index')
							Counterbad, blub, Delx, tau_B, Blobf,vr, rval, vdmax, vimax = write_nan(fu, shift_Block, Measurement,maxlen, NewData, Radial, WriteHeader, DenCase, EmCase, BlockCase, Counterbad, blub, t_start,t_end, x,y,Refdec_ind, shift, BlobCount, Delx, tau_B, Blobf,vr, rval, vdmax, vimax, B0, Lp, q95, Te0, Ti0, ne0, omegaCi, rhoS, endtime)
							continue
						
						Delx = abs(x_inter[low_index]-x_inter[up_index])/2.	# HWHM of intensity at Delt = 0.
						print('HWHM of standard blob (interpolated data): {0:.3f}cm'.format(Delx))

#					if blub==999:							# exit analysis because of bad result
#						print(blub, 'one level deeper')
#						break

					if not BlockCase:
						LiConAv_max = max(LiConAv[NullPos_ind,:])		# Maximum emission value
						for c in range (0,len(x)):
							if (LiConAv_max==LiConAv[NullPos_ind,c]):
								max_ind = c				# index on x-axis of maximum emission value

						# find closest index to half-max-position (lower index)	
						PeakC = [np.NaN]*10					# peak counter to determine for sure the position of half-max
						RandCount = 1
						for d in range (max_ind):
							if (LiConAv[NullPos_ind,d]>LiConAv_max/2 and LiConAv[NullPos_ind,d-1]<LiConAv_max/2):
								PeakC[RandCount]=d			# store index in peak counter array (in best case only one entry peakC[1]!)
								RandCount = RandCount+1
						low_index = np.nanmax(PeakC)					# lower index is maximum for last time it surpasses HM

						# find closest index to half-max-position (upper index)
						PeakC = [None]*10					# peak counter array to determine for sure the position of half-max
						PeakC = np.array(PeakC)
						PeakC.fill(10000)					# fill with high values (higher than len(x), since we need the lowest one
						RandCount = 1
						for e in range (max_ind+1,len(x)):
							if (LiConAv[NullPos_ind,e-1]>LiConAv_max/2 and LiConAv[NullPos_ind,e]<LiConAv_max/2):
								PeakC[RandCount]=e			# store index in peak counter array (in best case only one entry peakC[1]!)
								RandCount = RandCount+1	
						up_index = min(PeakC)					# up index is minimum index for first time it surpasses HM
						print('low index', low_index, 'up_index', up_index)
						
						if (up_index>5000 or low_index<1):			# exit loop for this analysis if a proper calculation is not possible here
							blub=999
							print(blub, 'in up/low_index')
							if not BlockCase:
								shift_Block = 0
							Counterbad, blub, Delx, tau_B, Blobf,vr, rval, vdmax, vimax = write_nan(fu, shift_Block, Measurement,maxlen, NewData, Radial, WriteHeader, DenCase, EmCase, BlockCase, Counterbad, blub, t_start,t_end, x,y,Refdec_ind, shift, BlobCount, Delx, tau_B, Blobf,vr, rval, vdmax, vimax, B0, Lp, q95, Te0, Ti0, ne0, omegaCi, rhoS, endtime)
							continue
						Delx = (x_axis[low_index]-x_axis[up_index])/2.		# HWHM of intensity at Delt = 0.
						print('HWHM of standard blob: {0:.3f}cm'.format(Delx))

#					if blub==999:							# exit analysis because of bad result
#						print(blub, 'one level deeper')
#						break

				#	Calculate self correlation time of standard blob a FWHM of I(x,Delt) at reference channel ####################################
					
					# interpolate data for better calculation of slope
					TWind_inter = np.linspace(LookWinMin,LookWinMax, 900)		
					LiConAv_self_smooth = spline(Delt,LiConAv[:,Refdec_ind],TWind_inter)

					LiConAv_self = max(LiConAv_self_smooth[:])
					RandCount=1 
					for c in range (0,len(TWind_inter)):
						if (LiConAv_self==LiConAv_self_smooth[c]):
							self_ind = c				# index on time-window-axis of maximum emission value



					# find closest index to half-max-position (lower index)	
					PeakC = [None]*10					# peak counter to determine for sure the position of half-max
					RandCount = 1
					for d in range (self_ind):
						if (LiConAv_self_smooth[d]>LiConAv_self/2 and LiConAv_self_smooth[d-1]<LiConAv_self/2):
							PeakC[RandCount]=d			# store index in peak counter array (in best case only one entry peakC[1]!)
							RandCount = RandCount+1
					left_index = max(PeakC)					# lower index is maximum for last time it surpasses HM

					# find closest index to half-max-position (upper index)
					PeakC = [None]*10					# peak counter array to determine for sure the position of half-max
					PeakC = np.array(PeakC)
					PeakC.fill(10000)					# fill with high values (higher than len(x), since we need the lowest one
					RandCount = 1
					for e in range (self_ind+1,len(TWind_inter)):
						if (LiConAv_self_smooth[e-1]>LiConAv_self/2 and LiConAv_self_smooth[e]<LiConAv_self/2):
							PeakC[RandCount]=e			# store index in peak counter array (in best case only one entry peakC[1]!)
							RandCount = RandCount+1	
					right_index = min(PeakC)				# up index is minimum index for first time it surpasses HM
					if (right_index>5000 or left_index==0):
						blub==999	
						print(blub)			# exit analysis because of bad result
						if not BlockCase:
							shift_Block = 0
						Counterbad, blub, Delx, tau_B, Blobf,vr, rval, vdmax, vimax = write_nan(fu, shift_Block, Measurement,maxlen, NewData, Radial, WriteHeader, DenCase, EmCase, BlockCase, Counterbad, blub, t_start,t_end, x,y,Refdec_ind, shift, BlobCount, Delx, tau_B, Blobf,vr, rval, vdmax, vimax, B0, Lp, q95, Te0, Ti0, ne0, omegaCi, rhoS, endtime)
						continue
					if (right_index==10000 or left_index==0):
						sys.exit('Please change time-window, because HWHM-calculation could not take place')
					tau_B = (TWind_inter[right_index]-TWind_inter[left_index])/2.	# HWHM of intensity at Delt = 0.
					print('Self correlation time of standard blob: {0:.3f}ms'.format(tau_B))
					self_TWind_inter = TWind_inter



				#	Calculate Blobfrequency #########################################################################################################
					Blobf = BlobCount/time[-1]*1000
					print('The Blobfrequency is: {0:d}'.format(int(Blobf)))


				# 	Calculate center of mass (COM) for blobs at different t in TWind ######################################################################

					# set all negative values of LiConAv to zero -- it should not count for calculation of COM
					for l in range (0,len(LiConAv[:,1])):
						for k in range (0,len(LiConAv[1,:])):
							if (LiConAv[l,k]<0):
								LiConAv[l,k]=0
					COM = [None]*len(LiConAv[:,1])		# create COM-vector with appropriate shape
					for m in range (0,int(TWind/avt)):
						Blubern = ndimage.measurements.center_of_mass(LiConAv[m,:])		# Provides tuple (x and y) --> Blubern is intermediate to hold tuple and not used, while COM[m] holds the index of the COM on the x-axis
						if (Blubern<0):
							print('There a negative index values for COM location in x[m] --> change time Window or position of detector and look into it -- something is weird!')
					#	print(Blubern[0])
						COM[m] = x[int(Blubern[0])]						# use only first part of tuple
					COM = np.array(COM)								# convert list to array




				###################################################################################################################################
				########### FIGURE 2 (PART1) ######################################################################################################
				###################################################################################################################################
					if not Radial:
						f2=plt.figure(figsize=(8,8), facecolor = 'white')
						gs = gridspec.GridSpec(2, 2,			# ratio of grid space (x plots per collumn [hight], y per row [width])
							width_ratios=[1,1],		# width ratios of the 3 plots per row
							height_ratios=[1,1]		# height ratios of the 2 polots per collumn
							)


						# create axes with the ratios specified in gs above. ax5, ax6 are the axes for the colorbars:
						ax21 = plt.subplot(gs[0])
						ax22 = plt.subplot(gs[1])
						ax23 = plt.subplot(gs[2])
						ax24 = plt.subplot(gs[3])



					#	AX21: Standard-Blob: representative image of radial-temporal data I(x,Deta t) ###############################################

						I2 = ax21.contourf(BloT, BloX, LiConAv,300)
						if BlockCase:
							ax21.plot([Delt[NullPos_ind],Delt[NullPos_ind]],[x_inter[up_index],x_inter[low_index]], linewidth=3, color='w')	# FWHM-plot
						if not BlockCase:
							ax21.plot([Delt[NullPos_ind],Delt[NullPos_ind]],[x[up_index],x[low_index]], linewidth=3, color='w')	# FWHM-plot
						ax21.plot([self_TWind_inter[left_index],self_TWind_inter[right_index]],[x[Refdec_ind],x[Refdec_ind]],linewidth=3, color = 'w')
						ax21.set_xlabel('time $\Delta$t', labelpad= -10)			
						ax21.set_ylabel('beam axis $x$ (cm)')
						ax21.set_title('Spectrum for standard blob')
						#time points for representatives blobs are drawn into plot
						ax21.axvline(Delt[a_an],color='r', linestyle='--')
						ax21.axvline(Delt[b_an],color='g', linestyle='--')
						ax21.axvline(Delt[c_an],color='b', linestyle='--')
						ax21.axvline(Delt[NullPos_ind], color='w', linestyle='-.')						# Null pos v-oline
						# ax1.set_xlim(xmin,xmax)
						# ax1.set_ylim(ymin,ymax)
						div21 =  make_axes_locatable(ax21)
						cax21=div21.append_axes("right",size="5%",pad=0.05)
						chbar_I = plt.colorbar(I2,cax=cax21)
						chbar_I.set_label('I(x,$\Delta$t) (a.u.)')




					#	AX22: Radial profiles of intensity response 
						I31 = ax22.plot(x,LiConAv[a_an,:], color='r')
						I32 = ax22.plot(x,LiConAv[b_an,:], color='g')
						I33 = ax22.plot(x,LiConAv[c_an,:], color='b')
						I34 = ax22.plot(x,LiConAv[NullPos_ind,:], color='k', linestyle ='-.')					# corresponding plot of NullPos v-line
						if BlockCase:
							I35 = ax22.plot(x_inter,LiConAv_smooth_a, color='r')
							I36 = ax22.plot(x_inter,LiConAv_smooth_b, color='g')
							I37 = ax22.plot(x_inter,LiConAv_smooth_c, color='b')
							I38 = ax22.plot(x_inter,LiConAv_smooth_Null, color='k', linestyle = '-.')
							ax22.plot([x_inter[up_index],x_inter[low_index]], [LiConAv_max/2,LiConAv_max/2], linewidth=3, color='k')		# FWHM
						if not BlockCase:
							ax22.plot([x[up_index],x[low_index]], [LiConAv_max/2,LiConAv_max/2], linewidth=3, color='k')		# FWHM
						# Center of mass positions drawn into plot:
						ax22.axvline(COM[a_an],color='r', linestyle='--')
						ax22.axvline(COM[b_an],color='g', linestyle='--')
						ax22.axvline(COM[c_an],color='b', linestyle='--')
						ax22.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			
						ax22.set_ylabel('I(x,$\Delta$t) (a.u.)')
						ax22.set_title('Blob representatives for dashed lined time positions in left plot')
						# ax1.set_xlim(xmin,xmax)
						# ax1.set_ylim(ymin,ymax)
						# div21 =  make_axes_locatable(ax21)
						# cax21=div21.append_axes("right",size="5%",pad=0.05)
						# chbar_I = plt.colorbar(LiConAv,cax=cax21)
						# chbar_I.set_label('I(x,$\Delta$t', size=10)



				###################################################################################################################################
				########### Standard Blob analysis (part 2) #######################################################################################
				###################################################################################################################################

				#	Position of Blob and speed ################################################################################################

					# calculate lienar regression statistics:
					v_r, v_intercept, r_value, p_value, std_err = stats.linregress(Delt[LookWinMin_ind:LookWinMax_ind],COM[LookWinMin_ind:LookWinMax_ind])

					#create plot information of linear regression
					(a_lgr,b_lgr) =polyfit(Delt[LookWinMin_ind:LookWinMax_ind],COM[LookWinMin_ind:LookWinMax_ind],1) 
					COM_lgr = np.polyval([a_lgr,b_lgr],Delt[LookWinMin_ind:LookWinMax_ind])

					# interpolate data for better calculation of slope
					TWind_inter = np.linspace(LookWinMin,LookWinMax,300)		
					XC_COM_smooth = spline(Delt,COM,TWind_inter)
					
					#calculate maximum slopes of real data and interpolated data:
					slope_data=[None]*(LookWinMax_ind-LookWinMin_ind)
					slope_inter=[None]*(len(TWind_inter)-1)
					slope_data=np.array(slope_data)						# convert list to array
					slope_inter=np.array(slope_inter)					# convert list to array
					slope_data[0]=0
					slope_inter[0]=0
					for q in range (LookWinMin_ind+1,LookWinMax_ind):			# Calculate slopes between every point
						slope_data[q] = (COM[q]-COM[q-1])/abs(Delt[q]-Delt[q-1])
						if (slope_data[q]<0):
							slope_data[q]=0

					v_dmax=np.amax(slope_data)						# find maximum slope of data 
					for p in range (0,len(TWind_inter)-1):					# Calculate slopes between every point
						slope_inter[p] = (XC_COM_smooth[p]-XC_COM_smooth[p-1])/abs(TWind_inter[p]-TWind_inter[p-1])
						if (slope_inter[p]<0):
							slope_inter[p]=0
					v_imax=np.amax(slope_inter)						# find maximum slope interpolation


					#convert into [m/s]
					vr = v_r*0.01*pow(10,3)
					rval = r_value*0.01*pow(10,3)
					vdmax = v_dmax*0.01*pow(10,3)
					vimax = v_imax*0.01*pow(10,3)

					print('Average speed of blob:{0:.1f}m/s with standard deviation of {1:.1f}m/s'.format(vr,rval))	# Average speed [m/s]
					print('Maximum speed of blob (real data): {0:.1f}m/s'.format(vdmax))	# Maximum speed from data [m/s]
					print('Maximum speed of blob (smoothed): {0:.1f}m/s:'.format(vimax))	# Maximum speed from smoothing [m/s]



				###################################################################################################################################
				########### FIGURE 2 (PART2) ######################################################################################################
				###################################################################################################################################
					if not Radial:
					# 	AX23: plot different data for blob position and speed

						XC = ax23.plot(Delt,COM, color = 'k')				# plot original data
						XC_inter = ax23.plot(TWind_inter, XC_COM_smooth, color = 'm')	# plot interpolated data
						XC_lgr = ax23.plot(Delt[LookWinMin_ind:LookWinMax_ind],COM_lgr, 'r--') # plot regression
						ax23.set_xlabel('time $\Delta$t', labelpad= -10)			
						ax23.set_ylabel('COM (cm)')
						ax23.set_title('COM of Blob')
						ax23.set_xlim(LookWinMin,LookWinMax)
						ax23.set_ylim(x[-1],x[0])
						div23 =  make_axes_locatable(ax23)
						cax23 = div23.append_axes("right",size="5%",pad=0.05)
					
						for ylabel_i in cax23.get_yticklabels():
							ylabel_i.set_visible(False)
							ylabel_i.set_fontsize(0.0)
						for xlabel_i in cax23.get_xticklabels():
							xlabel_i.set_fontsize(0.0)
							xlabel_i.set_visible(False)

					# 	AX24: plot emission data for standard blob at reference position

						Em_Stan = ax24.plot(Delt,LiConAv[:,Refdec_ind], color = 'k')				# plot original data
						Em_Stan_inter = ax24.plot(self_TWind_inter, LiConAv_self_smooth, color = 'm')	# plot interpolated data
						ax24.plot([self_TWind_inter[left_index-1],self_TWind_inter[right_index]],[LiConAv_self_smooth[left_index-1],LiConAv_self_smooth[right_index]],linewidth=3, color = 'm')
						ax24.set_xlabel('time $\Delta$t', labelpad= -10)			
						ax24.set_ylabel('I(x,$\Delta$t) (a.u.)')
						ax24.set_title('Standard blob at reference detector at %.2g cm' % (Refdec))
						ax24.set_xlim(LookWinMin,LookWinMax)
						# ax24.set_ylim(x[-1],x[0])
						# div24 =  make_axes_locatable(ax23)
						# cax24 = div23.append_axes("right",size="5%",pad=0.05)

				###################################################################################################################################
				########### PRINT TO FILE - OUTPUT STATEMENST #####################################################################################
				###################################################################################################################################
				if not NewData:
					print('output file is not produced!')
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
					

	if NewData:	
		fu.close()
print('Number of events left out because of bad time windows:', Counterbad)

