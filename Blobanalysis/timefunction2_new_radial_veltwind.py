#!/usr/bin/env python
# -*- coding: utf-8 -*-


########### Importing modules ###########################################################################

import numpy as np
import matplotlib
#matplotlib.use('Agg')
#default backend, can be changed by matplotlib = reload(matplotlib) and then plt.switch_backend('GTK') in e.g. pdb.
import matplotlib.pyplot as plt
#from sympy import *
from pylab import plot
import matplotlib.gridspec as gridspec		# define layout of grid
from mpl_toolkits.axes_grid1 import make_axes_locatable # for colorbars
from scipy import ndimage			# calculate com
from scipy.interpolate import spline		# for interpolate COM position for good speed calculation of blob
from scipy import stats, polyfit		# for linear interpolation data
from scipy_cookbook_smooth import smooth	# for interpolating noisy data
from scipy import signal			# for smoothing time signal
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

# Measurements to be analyzed:
#files = [004,108,109,110,111,112,113]

files = [113]
shots = [29315]
shots = [29302,29303, 29306, 29307, 29308, 29309, 29310,29311,29312,29315]

RealCase = False		# for reading in real emission data from ASDEX

thressweep = True


if RealCase:
	files = shots

storepref = '_veltwind'				# specify fileextension to be read in and to be stored (e.g '_dens_test', if file is called '004Emresults_dens_test')
if RealCase:
	storepref = 'Real_data'
if RealCase and thressweep:
	storepref = 'Real_data_thressweep'
if thressweep:
	storepref = '_thressweep'	

for flup in range (0,len(files)):
	if not RealCase:
		Measurement=files[flup]						# select number of measurment to treat (is used for hdf-selection etc)
	if RealCase:
		shot = files[flup]


	# Case Loop
	Counterbad = 0				# Number of bad events
	if not RealCase:
		ad = 0
		bd = 4

	if RealCase:
		ad = 5
		bd = 6
	for Case in range (ad,bd):

		EmCase = False			# infinite emission case --> actually always true, since every different setting during analysis is made by if-statements with the other two cases --> can stay false
		DenCase = False				# density analysis case
		BlockCase = False			# block-averaging emission case

		if Case == 0:
			EmCase = True				# infinite emission case --> actually always true, since every different setting during analysis is made by if-statements with the other two cases --> can stay false
			DenCase = False				# density analysis case
			BlockCase = False			# block-averaging emission case
		if Case == 1:
			EmCase = False				# infinite emission case --> actually always true, since every different setting during analysis is made by if-statements with the other two cases --> can stay false
			DenCase = True				# density analysis case
			BlockCase = False			# block-averaging emission case
		if Case == 2: # BlockCase with Noise!
			EmCase = False				# infinite emission case --> actually always true, since every different setting during analysis is made by if-statements with the other two cases --> can stay false
			DenCase = False				# density analysis case
			BlockCase = True			# block-averaging emission case
		if Case == 3: # BlockCase without Noise!
			EmCase = False				# infinite emission case --> actually always true, since every different setting during analysis is made by if-statements with the other two cases --> can stay false
			DenCase = False				# density analysis case
			BlockCase = True			# block-averaging emission case


		if (EmCase and DenCase or EmCase and BlockCase or DenCase and BlockCase):		#if you do something wrong here ;)
			sys.exit('Accidentally set two cases for evaluation!')

		Noise = False				# switch on additive white Gaussian noise --> Very time consuming!
		Smooth = False				# smoothing noisy data afterwards

		NoiseTime = False		# puts  noise on time signal
		SmoothTime = False		# smoothes time-signal
		realSNR = False		# reads in SNR profile for the Block Case

		if Case == 2:
			NoiseTime = True		# puts  noise on time signal
			SmoothTime = True		# smoothes time-signal
			realSNR = True	# reads in SNR profile for the Block Case

		yvariation = True		# Calculate more values for y-variation
			
		Standardblob= True		# switch for calculations

		NewData = True			# Switch to write output

		Radial = True
		
		blub = 0			# emergency exit for bad results

		# Writing head of file for radial case ########################################################################################
		'''
		if Noise:
			storename = storepref+'_SNR{0:03d}'.format(int(SNR*100))			# change file extension as needed!
		if Smooth:
			storename = storepref+'_SNR{0:03d}'.format(int(SNR*100))+'_Smooth'			# change file extension as needed!
		if NoiseTime:
			storename = storepref+'_TIME_SNR{0:03d}'.format(int(SNR*100))			# change file extension as needed!
		if SmoothTime:
			storename = storepref+'_TIME_SNR{0:03d}'.format(int(SNR*100))+'_Smooth'			# change file extension as needed!
		if not NoiseTime and not Noise:
			storename = storepref
		'''

		if realSNR and NoiseTime:
			storename = storepref+'_TIME_SNR_Comparison'
				
		if not NoiseTime and not Noise:
			storename = storepref


		if Radial:	

			WriteHeader = True 

			if not WriteHeader:
				print('no header is written for outputfile!')
			if EmCase and NewData:
				fu = open('{0:03d}Emresults{1:}.txt'.format(Measurement,storename), 'wb')		# use only writing (option 'w') if needed (ATTENTION: discards the whole file, appending uses 'a')
				print('results are stored in {0:03d}Emresults{1:}.txt'.format(Measurement,storename))
			if DenCase and NewData:
				fu = open('{0:03d}Denresults{1:}.txt'.format(Measurement,storename), 'wb')		# use only writing (option 'w') if needed (ATTENTION: discards the whole file, appending uses 'a')
				print('results are stored in {0:03d}Denresults{1:}.txt'.format(Measurement,storename))
			if BlockCase and NewData:
				fu = open('{0:03d}Blockresults{1:}.txt'.format(Measurement,storename), 'wb')		# use only writing (option 'w') if needed (ATTENTION: discards the whole file, appending uses 'a')	
				print('results are stored in {0:03d}Blockresults{1:}.txt'.format(Measurement,storename))
			if RealCase and NewData:
				fu = open('{0:}Realresults{1:}.txt'.format(shot,storename), 'wb')		# use only writing (option 'w') if needed (ATTENTION: discards the whole file, appending uses 'a')	
				print('results are stored in {0:}Realresults{1:}.txt'.format(shot,storename))

				
			# we do not need noise in DenCase:
			if DenCase:
				NoiseTime = False
				SmoothTime = False

		#	Putting everything into lists and pre-processing for output ###########################################################
			# Collumn index print for better human readability
			if WriteHeader and NewData:
				pre=[None]*24
				countblub=0
				for i in range(len(pre)):
					pre[i]= countblub
					countblub=countblub+1

				pre[0]="Col.-Index:  0"

				# header is list consisting of header titles:
				if not RealCase:
					header = ['Date & Time', 'Measurement #', 'Starting time [ms]', 'Observed time [ms]', 'Refdec-pos[cm]', 'y-pos of beam [cm]', '# of Blobs', 'HWHM of Blob [cm]', 'tau_B [ms]','Blobfrequency [1/s]','Average speed vr [m/s]','Std. error vr [m/s]', 'max speed vdmax [m/s]', 'max i-pol. vimax [m/s]', 'rel. flucutation', 'B0 [T]', 'Lp [m]', 'q95 ', 'Te0 [eV]', 'Ti0 [eV]', 'ne0 [m-3]', 'omegaCi [1/s]', 'rhoS [m]', 'endtime [ms]']
				if RealCase:
					header = ['Date & Time', 'Shot #', 'Starting time [ms]', 'Observed time [ms]', 'Refdec-pos[cm]', 'End time [ms]', '# of Blobs', 'HWHM of Blob [cm]', 'tau_B [ms]','Blobfrequency [1/s]','Average speed vr [m/s]','Std. error vr [m/s]', 'max speed vdmax [m/s]', 'max i-pol. vimax [m/s]', 'rel. flucutation', 'B0 [T]', 'Lp [m]', 'q95 ', 'Te0 [eV]', 'Ti0 [eV]', 'ne0 [m-3]', 'omegaCi [1/s]', 'rhoS [m]', 'endtime [ms]']

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


		if not RealCase:
			# will be used for file search and in output
			print('Data is going to be read from file: /hesel/data/ASDEX.{0:03d}.h5'.format(Measurement))
			# read in data from function Li_read_hdf5: Specify correct hdf5-file
		#	time, x, y, ne, LiTemp, Li2p1, B0, Lp, q95, Te0, Ti0, ne0, omegaCi, rhoS, endtime = Li_read_hdf5('/pfs/home/bschiess/public/hesel/data/ASDEX.{0:03d}.h5'.format(Measurement))
			time, x, y, ne, LiTemp, Li2p1, B0, Lp, q95, Te0, Ti0, ne0, omegaCi, rhoS, endtime = Li_read_hdf5('/pfs/home/bschiess/itmwork/hesel_data/ASDEX.{0:03d}.h5'.format(Measurement))

		if RealCase:

			print('Data is going to be read from file: {0:}_LiBes_signal_data.txt'.format(shot))
			Li2p1 = np.loadtxt('/afs/ipp-garching.mpg.de/u/bschiess/Analyzetools/Li-BES GUI/{0:}_LiBes_signal_data.txt'.format(shot))
			time = np.loadtxt('/afs/ipp-garching.mpg.de/u/bschiess/Analyzetools/Li-BES GUI/{0:}_LiBes_time_data.txt'.format(shot))
			indizes_beam_on = np.loadtxt('/afs/ipp-garching.mpg.de/u/bschiess/Analyzetools/Li-BES GUI/{0:}_LiBes_beam_on_indizes.txt'.format(shot))
			indizes_beam_on = indizes_beam_on.astype(int)				# convert to int
			indizes_beam_on = np.delete(indizes_beam_on,(-1),axis=0)		# last value is wrong
			shot_params = np.loadtxt('/afs/ipp-garching.mpg.de/u/bschiess/Analyzetools/Li-BES GUI/{0:}_LiBes_parameter.txt'.format(shot))
			#x1 = np.loadtxt('/afs/ipp-garching.mpg.de/u/bschiess/Analyzetools/Li-BES GUI/29302_LiBes_rhop_at2900ms.txt')
			
			#x = [None]*16
			#for x_ind in range(len(x)):
			#	x[x_ind] = x1[x_ind]

			B0 = shot_params[2]
			q95 = shot_params[0]
			Lp = shot_params[1]
			Ip = shot_params[3]
		

			Te0 = 0
			Ti0 = 0
			ne0 = 0
			endtime = 0
			omegaCi = 0
			rhoS = 0
			Measurement = shot
			y = 0
			shift = 0
			shift_Block = 0

			Li2p1 = Li2p1.transpose() 	# for right shape, as it was used before

		# other settings and preparations  #############################################################################################

		if not RealCase:
			# selected times, which are supposed to be investigated [ms]
			t_start=0.2
			t_end = time[-1]
			print('time to be investigated betweem %.3gms and %.3gms:' % (t_start, t_end))

			# calculate time in [us] for snapshots:
			point=position*(time[11]-time[10])

		if RealCase:
			# define total times to analyze
			t_start = 2.9		# starttime for the whole analysis (including standard-blob-analysis)
			t_end = 3.3		# endtime
			if shot == 29315:
				t_start = 2.6		# starttime for the whole analysis (including standard-blob-analysis)
				t_end = 3.0		# endtime
										
			# start and end indizes for time array, for the phase in between 2.9 and 3.3s 
			for i in range(len(time)):
				if time[i]>t_start:
					start_ind = i
					break
			for i in range(len(time)):
				if time[i]>t_end:
					end_ind = i-1
					break

			# start and end indizes for beam-on-phase, for the phase in between 2.9 and 3.3s 
			for j in range(len(indizes_beam_on)):
				if start_ind<indizes_beam_on[j]:
					ind_b_start = j
					break

			for j in range(len(indizes_beam_on)):
				if end_ind<indizes_beam_on[j] or j == len(indizes_beam_on)-1:
					ind_b_end = j-1
					break

			print('time to be investigated betweem %.3gs and %.3gs:' % (t_start, t_end))
			
			#approximated beam position:
		#       ch = np.arange(1,27)				# assumed number of channels: 26
		#        x = 0.652*(ch-1) +1.01				# equation from calibrate_functions.py
			
			#make a short version of x-axis:
			ch = np.arange(1,15)				# assumed number of channels: 26
			x_1 = 0.652*(ch-1) +1.01	

			x = np.array([0, 1.2,1.7,2.4,3,3.6,4.22,5.45,6.1,6.75,7.37,8.02,8.65,9.3,9.95,10.58]) 

		# some lengthes used during the program:
		mn=len(x)
		x_axis = x					# just useful sometimes to avoid confusion
		timestep=len(time)				# misleading name, sorry
		avt=(time[timestep-1]-time[0])/timestep
		if RealCase:
			avt = avt*10**3
			for i in range(len(time)):
				time[i]=time[i]*10**3
		if not RealCase:
			stepsize=len(y)					# misleading name, sorry
			avy=y[stepsize-1]/stepsize			
			

		
		if Noise:
			for x_ind in range(mn):
				for y_ind in range(stepsize):				
					noiLi = np.random.normal(0,np.mean(Li2p1[:,y_ind,x_ind])/SNR,timestep)
					for t_ind in range(timestep):
						Li2p1[t_ind,y_ind,x_ind] = Li2p1[t_ind,y_ind,x_ind]+noiLi[t_ind]
			Li2p1noise=Li2p1.copy()

			# smoothing noisy data
			if Smooth:
				Li2p1hn = Li2p1.copy()
			#	Li2p1hm = Li2p1.copy()
			#	Li2p1ba = Li2p1.copy()
			#	Li2p1bl = Li2p1.copy()

				for t_ind in range(timestep):
					for y_ind in range(stepsize):
						Li2p1hn[t_ind,y_ind,:] = smooth(Li2p1[t_ind,y_ind,:],window_len=smoothlen,window='hanning')
					#	Li2p1hm[t_ind,shift*2,:] = smooth(Li2p1[t_ind,shift*2,:],window_len=smoothlen,window='hamming')
					#	Li2p1ba[t_ind,shift*2,:] = smooth(Li2p1[t_ind,shift*2,:],window_len=smoothlen,window='bartlett')
					#	Li2p1bl[t_ind,shift*2,:] = smooth(Li2p1[t_ind,shift*2,:],window_len=smoothlen,window='blackman')
				Li2p1 = Li2p1hn.copy()
				Li2p1cop = Li2p1hn.copy()	# used for passing to detector function


		
		if NoiseTime and not BlockCase:
			print('The SNR is set to: {0:.2f}'.format(SNR))
			''' if only the quick smoothing analysis of the reference detector position is required 
			noiLi = np.random.normal(0,np.mean(Li2p1[:,shift*2,Refdec_ind])/SNR,timestep)
			for t_ind in range(int(t_start/avt),int(t_end/avt)):
				Li2p1[t_ind,shift*2,Refdec_ind] = Li2p1[t_ind,shift2,Refdec_ind]+noiLi[t_ind]
			Li2p1noise=Li2p1.copy()
			'''
			Li2p1_orio = Li2p1.copy()
			Li2p1noise=Li2p1.copy()
			for x_ind in range(mn):
				for y_ind in range(stepsize):			
					noiLi = np.random.normal(0,np.mean(Li2p1[:,y_ind,x_ind])/SNR,timestep)
					for t_ind in range(int(t_start/avt),int(t_end/avt)):
						Li2p1noise[t_ind,y_ind,x_ind] = Li2p1[t_ind,y_ind,x_ind]+noiLi[t_ind]

			# smoothing possible noisy data
			
			def butter_lowpass(cutoff,fs,order=5):
				nyq = 0.5*fs
				normal_cutoff = cutoff/nyq
				b, a = signal.butter(order, normal_cutoff, btype = 'low', analog = False)
				return b, a
			def butter_lowpass_filtfilt(data, cutoff, fs, order=5):
				b, a = butter_lowpass(cutoff, fs, order = order)
				y = signal.filtfilt(b, a, data)
				return y
			cutoff = 6000
			fs = 200000

			Li2p1hn = Li2p1.copy()
			Li2p1hn2 = Li2p1.copy()
			for x_ind in range(mn):
				for y_ind in range(stepsize):		
					Li2p1hn2[int(t_start/avt):int(t_end/avt),y_ind,x_ind] = butter_lowpass_filtfilt(Li2p1noise[int(t_start/avt):int(t_end/avt),y_ind,x_ind],cutoff,fs)
					Li2p1hn[int(t_start/avt):int(t_end/avt),y_ind,x_ind] = smooth(Li2p1noise[int(t_start/avt):int(t_end/avt),y_ind,x_ind],window_len=smoothlen,window='hanning')

	#		smooth x-direction afterwards for Em-Case (high resolution)
			Li2p1xs = Li2p1.copy()
			for t_ind in range(timestep):
				for y_ind in range(stepsize):
					Li2p1xs[t_ind,y_ind,:] = smooth(Li2p1hn[t_ind,y_ind,:],window_len=smoothlenx,window='hanning')
			

			# only use the new noise and smoothed data if selected
			Li2p1=Li2p1noise.copy()
			if SmoothTime:
				Li2p1 = Li2p1xs.copy()


		if RealCase:
			# smoothing noisy real data
			def butter_lowpass(cutoff,fs,order=5):
				nyq = 0.5*fs
				normal_cutoff = cutoff/nyq
				b, a = signal.butter(order, normal_cutoff, btype = 'low', analog = False)
				return b, a
			def butter_lowpass_filtfilt(data, cutoff, fs, order=5):
				b, a = butter_lowpass(cutoff, fs, order = order)
				y = signal.filtfilt(b, a, data)
				return y
			cutoff = 15000
			fs = 200000
			
			#Li2p1hn = Li2p1.copy()
			Li2p1hn2 = Li2p1.copy()
			for x_ind in range(len(x)):
				Li2p1hn2[indizes_beam_on[ind_b_start:ind_b_end],x_ind] = butter_lowpass_filtfilt(Li2p1[indizes_beam_on[ind_b_start:ind_b_end],x_ind],cutoff,fs)
				#Li2p1hn[indizes_beam_on[ind_b_start:ind_b_end],x_ind] = smooth(Li2p1[indizes_beam_on[ind_b_start:ind_b_end],x_ind],window_len=smoothlen,window='hanning')


		# some stuff we need copies from during the run:

		if BlockCase:						 # stuff which is always overwritten inside of Block-Case-loop and which has to be copied in...
			xcop = x.copy()
			Li2p1cop = Li2p1.copy()			         # used for passing in detector function
			mncop = mn
		if DenCase:
			Li2p1Z = ne.copy()				 # copy matrix for use in calculation of Lip1Z (the zero-mean-fluctuation)
			LiBlobEm = ne.copy()				 # copy Li-Matrix for BlobCount use
		if BlockCase or EmCase:
			Li2p1Z = Li2p1.copy()				 # copy matrix for use in calculation of Lip1Z (the zero-mean-fluctuation)
			LiBlobEm = Li2p1.copy()				 # copy Li-Matrix for BlobCount use
		if RealCase:
			Li2p1 = Li2p1hn2.copy()
			#calculate the matrix only, if it has not been calculated before or changes occurred:
			Li2p1Z = Li2p1.copy()			 # copy matrix for use in calculation of Lip1Z (the zero-mean-fluctuation)

		if RealCase:
			#find the different On-phases for the beam (usually around 5)
			num = 0
			for i in range(len(indizes_beam_on)):
				if indizes_beam_on[i] != indizes_beam_on[i-1]+1 and indizes_beam_on[i+1] == indizes_beam_on[i]+1 or indizes_beam_on[i] == 0:
					num = num+1
			Onswitch = [None]*num		
			Offswitch = [None]*num
			num = 0
			for i in range(len(indizes_beam_on)-1):
				if indizes_beam_on[i] != indizes_beam_on[i-1]+1 and indizes_beam_on[i+1] == indizes_beam_on[i]+1 or indizes_beam_on[i] == 0:
					Onswitch[num] = indizes_beam_on[i]		# index giving the start of the block
				if indizes_beam_on[i] != indizes_beam_on[i+1]-1 and indizes_beam_on[i] == indizes_beam_on[i-1]+1 or i == len(indizes_beam_on)-2:
					Offswitch[num] = indizes_beam_on[i]		# index giving the end of a block
					num = num+1

		# use realistic profile for noise for BlockCase
		if BlockCase:
			if realSNR:
				xSNR, rSNR = np.loadtxt('SNR_29302.txt', usecols = (0,1), unpack = True, skiprows=1)
				maxSNR = max(rSNR)
				for i in range(len(rSNR)):
					if maxSNR == rSNR[i]:
						xmaxSNR = xSNR[maxSNR]
						maxSNR_ind = i 

		if not RealCase:
			if yvariation:
				beamrange = np.arange (y_starting,y_ending,yResolution)
			else: 
				beamrange = [bcenter]

		if RealCase:
			beamrange = np.arange(0,num,1)		# do the whole measurement for one block alone

		for bcenter in beamrange:

			if RealCase:
				timepoint = bcenter	# rename for better understanding
				t_start = time[Onswitch[timepoint]]
				t_end = time[Offswitch[timepoint]]
											
				# start and end indizes for time array, for the phase in between 2.9 and 3.3s 
				for i in range(len(time)):
					if time[i]>=t_start:
						start_ind = i
						break
				for i in range(len(time)):
					if time[i]>t_end:
						end_ind = i-1
						break

				# start and end indizes for beam-on-phase, for the phase in between 2.9 and 3.3s 
				for j in range(len(indizes_beam_on)):
					if start_ind<=indizes_beam_on[j]:
						ind_b_start = j
						break

				for j in range(len(indizes_beam_on)):
					if end_ind<indizes_beam_on[j] or j == len(indizes_beam_on)-1:
						ind_b_end = j-1
						break

			if not RealCase:
				shift=int(bcenter/y[-1]*stepsize/2)+1		# shift of beam center in order to place the center at correct position
				shift2 = int(b/2/y[-1]*stepsize)		# shift*2-shift2 is now left index of beam edge and shift*2+shift2 the right one
				print('Position of Beam: {0:.2f}cm'.format(y[shift*2]))
			#	print('Beamwidth: {0:.2f}cm'.format(y[shift2+shift2]))	
				

			# use the block-emission data for the block case:
			if BlockCase:
				x_ori = xcop
			#	Li2p1 = block_func(Li2p1,timestep,mn,b,avy,shift,shift2)			# if original block-function should be used

				Li2p1_beam, Li2p1, x, blur = dect_func(Li2p1cop,timestep,mncop,b,avy,shift,shift2,xcop)	# if detector reduced block function should be used

				x_axis = xcop
				mn = len(x)
				shift_Block = shift
				shift = 0									# we are in the 1D-Case, so every time the entry at shift is called in matrix, this is set to the zero entry
			
				if realSNR:
					SNR = [None]*len(x)
					#find point with maximum emission, where also the SNR is highest
					for t_ind in range(int(t_start/avt),int(t_end/avt)):
						maxLi = max(Li2p1[t_ind,0,:])
						for x_ind in range(0,len(x)):
							if maxLi == Li2p1[t_ind,0,x_ind]:
								xmaxLi = x[x_ind]
								xmaxLi_ind = x_ind
						for x_ind in range(0,len(x)):
							SNR[x_ind] = rSNR[maxSNR_ind-xmaxLi_ind+x_ind]

				if NoiseTime:
					Li2p1noise=Li2p1.copy()
					for x_ind in range(len(x)):
						if realSNR:			# use only arry of SNR if real-SNR is read in
							noiLi = np.random.normal(0,np.mean(Li2p1[:,0,x_ind])/SNR[x_ind],timestep)
						if not realSNR:			# use no array
							noiLi = np.random.normal(0,np.mean(Li2p1[:,0,x_ind])/SNR,timestep)
						for t_ind in range(int(t_start/avt),int(t_end/avt)):
							Li2p1noise[t_ind,0,x_ind] = Li2p1[t_ind,0,x_ind]+noiLi[t_ind]

					# smoothing possible noisy data
					
					def butter_lowpass(cutoff,fs,order=5):
						nyq = 0.5*fs
						normal_cutoff = cutoff/nyq
						b, a = signal.butter(order, normal_cutoff, btype = 'low', analog = False)
						return b, a
					def butter_lowpass_filtfilt(data, cutoff, fs, order=5):
						b, a = butter_lowpass(cutoff, fs, order = order)
						y = signal.filtfilt(b, a, data)
						return y
					cutoff = 15000
					fs = 200000
					
					Li2p1hn = Li2p1.copy()
					Li2p1hn2 = Li2p1.copy()
					for x_ind in range(len(x)):		
						Li2p1hn2[int(t_start/avt):int(t_end/avt),0,x_ind] = butter_lowpass_filtfilt(Li2p1noise[int(t_start/avt):int(t_end/avt),0,x_ind],cutoff,fs)
						Li2p1hn[int(t_start/avt):int(t_end/avt),0,x_ind] = smooth(Li2p1noise[int(t_start/avt):int(t_end/avt),0,x_ind],window_len=smoothlen,window='hanning')
		
				# smoothing in x-direction for low resolution not necessary!
				#	Li2p1xs = Li2p1.copy()
				#	for t_ind in range(timestep):
				#		Li2p1xs[t_ind,0,:] = smooth(Li2p1hn[t_ind,0,:],window_len=smoothlenx,window='hanning')
			
					# only use the new noise and smoothed data if selected
					Li2p1=Li2p1noise.copy()
					if SmoothTime:
						Li2p1 = Li2p1hn2.copy()


			#Create zero-mean fluctuation matrix

			if EmCase or BlockCase:
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

			if RealCase:
				mean = [np.NaN]*mn			 # create vector for mean in correct length
				mean = np.array(mean)			 # convert to array
				for x_ind in range (0,mn):		 
					mean[x_ind] = np.mean(Li2p1[indizes_beam_on[ind_b_start:ind_b_end],x_ind])

				for x_ind in range(len(x)):
					for t_ind in range (start_ind,end_ind):
						Li2p1Z[t_ind,x_ind]=Li2p1[t_ind,x_ind]-mean[x_ind]  # subtract mean value for every timestep (takes most of the time for the program!)


			# Calculate index and position of Reference detector on xaxis (indices are counted backwards x[len(x)-1]<0!)


			if Radial:
				if RealCase:
					Startrange = 0
					Resolution = 1
				if BlockCase:
					Startrange = int(len(x)/5)
					Resolution = 1
					
				if DenCase or EmCase:
					Startrange = SetStartrange
					Resolution = SetResolution

				for Refdec_ind in range(Startrange, len(x),Resolution):			# Wall-region can be neglected (approximately 1.5-2cm) --> go through Refdec-Positions
					if RealCase:
						print('Position along beam-axis for Reference Detector: {0:.2f}cm'.format(x[Refdec_ind])) # check, if Reference detector is at correct position

	
					if EmCase or DenCase:					# otherwise this position is determined later
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
					if not RealCase:
						X,Y=np.meshgrid(x,y)						# transform x,y-axis in matrix for contourf-plot of 2D-plots
						if BlockCase:
							X_Ori,Y = np.meshgrid(x_ori,y)				# only needed, if BlockCase changes x-axis
						X2,T2=np.meshgrid(x,time[int(t_start/avt):int(t_end/avt)]) 	# Meshgrid for timeplots:
					#if RealCase:
					#	X2,T2=np.meshgrid(x,timelimview) 


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
						ax1.set_xlabel('axis $y$ (cm)')			# switched of, since axis in row below
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
							ax2.set_xlabel('axis $y$ (cm)')
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
							ax2.set_xlabel('axis $y$ (cm)')
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
						ax3.set_xlabel('time $t$ (ms)')
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
					if EmCase or BlockCase:
						LiAv = np.mean(Li2p1[int(t_start/avt):int(t_end/avt),shift*2,Refdec_ind])		# Average of Li value
						LiSig = np.std(Li2p1[int(t_start/avt):int(t_end/avt),shift*2,Refdec_ind])		# Standard deviation of Li value
					if DenCase:
						LiAv = np.mean(ne[int(t_start/avt):int(t_end/avt),shift*2,Refdec_ind])			# Average of Li value
						LiSig = np.std(ne[int(t_start/avt):int(t_end/avt),shift*2,Refdec_ind])			# Standard deviation of Li value
					
					if RealCase:
						LiAv = np.mean(Li2p1[indizes_beam_on[ind_b_start:ind_b_end],Refdec_ind])		# Average of Li value
						LiSig = np.std(Li2p1[indizes_beam_on[ind_b_start:ind_b_end],Refdec_ind])		# Standard deviation of Li value

					# choose the according thres-hold for the sweep-case:
					if thressweep:
						if x[0]<x[-1]:
							#sys.exit('The x-axis is inverted! Please change the thressweep-parameters accordingly!')
							for i in range(len(x)):
								if x[i]>=7.3 and x[i-1]<7.3:
									LCFS_ind = i
									break
							for i in range(len(x)):
								if x[i]>=2.9 and x[i-1]<2.9:
									Wall_ind = i
									break
							T1 = [2]*(Wall_ind)
							T1 = np.array(T1)
							T2 = np.linspace(2,1,abs(LCFS_ind-Wall_ind))
							T3 = [1]*(len(x)-LCFS_ind)
							T3 = np.array(T3)

							ThresCon_arr = np.append(T1,T2)
							ThresCon_arr = np.append(ThresCon_arr,T3)
						if x[0]>x[-1]:
							for i in range(len(x)):
								if x[i]<=4 and x[i-1]>4:
									Wall_ind = i
									break
							for i in range(len(x)):
								if x[i]<=0 and x[i-1]>0:
									LCFS_ind = i
									break
							T1 = [2]*(Wall_ind)
							T1 = np.array(T1)
							T2 = np.linspace(2,1,abs(LCFS_ind-Wall_ind))
							T3 = [1]*(len(x)-LCFS_ind)
							T3 = np.array(T3)

							ThresCon_arr = np.append(T1,T2)
							ThresCon_arr = np.append(ThresCon_arr,T3)
							#plt.plot(x,ThresCon_arr, linewidth = 2, color = 'b')
							#plt.xlabel('radial axis (cm)')
							#plt.ylabel('Threshold')
							#plt.ylim(0.7,2.3)
							#matplotlib = reload(matplotlib)
							#plt.switch_backend('GTK')
							#plt.show()
						ThresCon = ThresCon_arr[Refdec_ind]
						print('ThresCon for Refdec = {0:.2f} at {1:.2f}.'.format(x[Refdec_ind],ThresCon))
					LiCon = LiAv + ThresCon*LiSig									# Set condition
					TWind = TWindMult*avt							# Time window for conditional averaging [us] (Gregor's Publication: 500 us)
	# 				longer time window calculations ONLY for tau_B-calculation
					TWind_long = TWindMult_long*avt
					if RealCase: 
						TWind = 0.3
						if int(TWind/avt) % 2 == 0:
							TWind = TWind+avt
	


					if not DenCase:
						print('Average Emission Value at Reference Detector: {0:.4f} (a.u.) '.format(LiAv))
						print('Corresonding Standard Deviation: {0:4f} (a.u.)'.format(LiSig))
					if DenCase:
						print('Average Density at Reference Detector: {0:.2e}cm-3'.format(LiAv))
						print('Corresonding Standard Deviation: {0:.2e}cm-3'.format(LiSig))
					print('time window for conditional averaging is {0:.3f}ms consisting of {1:d} timesteps'.format(TWind, int(TWind/avt)))
					print('timewindow, in which too events are counted as one: {0:.3f}'.format(WinBar*avt))


					# Calculate number of blobs in time period:
					if BlockCase or RealCase:
						LiBlobEm = Li2p1.copy()						# copy Li-Matrix for BlobCount use
					#if DenCase:
						#LiBlobEm = ne.copy()						# copy density-Matrix for BlobCount use

					# Change matrix to 0s and 1s (0 below, 1 above threshold) for finding the position of blob occurance
					if not RealCase:
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
					if RealCase:
						# Change matrix to 0s and 1s (0 below, 1 above threshold) for finding the position of blob occurance
						for m in range (start_ind,end_ind):					# go through every time point
							if LiBlobEm[m,Refdec_ind]<LiCon:		# go in, if emission is below threshold
								LiBlobEm[m,Refdec_ind]=0		# set all values to zero, if it is below the threshold
							else:							# go in, if emission is above threshold
								LiBlobEm[m,Refdec_ind]=1		# set all values to 1, if it is above threshold
							m=m+1

						# Count Blobs at first time position, when threshold is passed, with condtion, that in LiBlobEm there should occur the event 0 1 1 or 0 1 0:
						BlobCount=0					# Counter counting the number of blobs
						remind=[np.NaN]*1000				# index used for time window in order to not double count an index
						for m in range (start_ind,end_ind):					# go through every time point
							if LiBlobEm[m,Refdec_ind]>0 and LiBlobEm[m+1,Refdec_ind]>0 and LiBlobEm[m-1,Refdec_ind]==0: # condition for counting an event as a blob (0 1 1)
								BlobCount=BlobCount+1							# raise blob count
								if (remind[BlobCount-1]>=m-WinBar):					# safty time intervall to not count one blob as two
									BlobCount=BlobCount-1
								remind[BlobCount]=m
							if LiBlobEm[m,Refdec_ind]>0 and LiBlobEm[m+1,Refdec_ind]==0 and LiBlobEm[m-1,Refdec_ind]==0:
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
						ax4.set_xlabel('time $t$ (ms)')
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

					if not RealCase:
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
					LiConAvrel = 0

					if (Standardblob==True and BlobCount>0):				# only occurs if switch is on and if there is certain number of blobs

					
					# 	Find emission values for standard blobs ###################################################################################

						Count2=0							# set blobcount to zero (now count again and find window around)
						LiConBlo = np.empty((BlobCount,int(TWind/avt),len(x)))		# Creating an empty, corectly shaped matrix for all blobs and their intensity values (later) [left time shift, right time shift, blobindex] and fill the seccond dim of the all-blob-matrix with zeros
						LiConBlo.fill(0)
						LiConAv = np.empty((int(TWind/avt), len(x)))			# array to hold values for standard blob later
						LiConAv.fill(0)
						print('LiConAv', np.shape(LiConAv))

						# find the time-windows of every blob event and put them into a matrix ######################################################
						if not RealCase:
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
						if RealCase:
							remind2 = [np.NaN]*1000			# reminding vector holding the indizes of blobs (maximum 1000 blobs)
							intermediatestore = np.NaN
							for m in range (start_ind,end_ind-int(TWind/avt)):					# go through every time point until one TWind before the edge --> a possible blob at the right edge-blob is not used
								if LiBlobEm[m,Refdec_ind]>0 and LiBlobEm[m+1,Refdec_ind]>0 and LiBlobEm[m-1,Refdec_ind]==0 or \
								LiBlobEm[m,Refdec_ind]>0 and LiBlobEm[m+1,Refdec_ind]==0 and LiBlobEm[m-1,Refdec_ind]==0: # condition for counting an event as a blob (0 1 1) or (0 1 0)
									Count2=Count2+1
									if (remind2[Count2-1]>=m-WinBar):								# safty time intervall to not count one blob as two
										intermediatestore = remind2[Count2-1]							# just for the entrance-condition to calculate emission time intervall
										Count2=Count2-1
									if not intermediatestore>=m-WinBar:								# enter only if that blob is not part of the first blob
										for l in range (0, int(TWind/avt)):							# go through time window
											for k in range (0, len(x)):							# go through all x-positions within time window
												LiConBlo[Count2-1,l,k]=Li2p1Z[(m-int(TWind/avt/2)+l),k]	# store emission values within time window in LiConBlo
										#print('Blob # {0} can be found at t={1:.2f}ms'.format(Count2,time[m]))
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
							if (Delt[m-1]<-1*10**(-6) and Delt[m]>-1*10**(-6)):
								NullPos_ind = m-1				# NullPos_ind is index, used for calculating blobwidth at HWHM#
							m=m+1
						if not Radial:
							BloX, BloT = np.meshgrid(x,Delt)			# new meshgrid for standard blob
						if (NullPos_ind==0):
							pdb.set_trace()
							sys.exit('The parameter for time window was not chosen to be even, NullPos_ind is not exact')
	# 	only for calculation for long time window for tau_B ###################################################################################

						if not RealCase:
							Count2=0							# set blobcount to zero (now count again and find window around)
							LiConBlo_long = np.empty((BlobCount,int(TWind_long/avt),len(x)))		# Creating an empty, corectly shaped matrix for all blobs and their intensity values (later) [left time shift, right time shift, blobindex] and fill the seccond dim of the all-blob-matrix with zeros
							LiConBlo_long.fill(0)
							LiConAv_long = np.empty((int(TWind_long/avt), len(x)))			# array to hold values for standard blob later
							LiConAv_long.fill(0)

							# find the time-windows of every blob event and put them into a matrix ######################################################
							remind2 = [np.NaN]*200			# reminding vector holding the indizes of blobs (maximum 200 blobs)
							intermediatestore = np.NaN
							for m in range (int(t_start/avt),int(t_end/avt)-1-int(TWind_long/avt)):					# go through every time point until one TWind before the edge --> a possible blob at the right edge-blob is not used
								if LiBlobEm[m,shift*2,Refdec_ind]>0 and LiBlobEm[m+1,shift*2,Refdec_ind]>0 and LiBlobEm[m-1,shift*2,Refdec_ind]==0 or \
								LiBlobEm[m,shift*2,Refdec_ind]>0 and LiBlobEm[m+1,shift*2,Refdec_ind]==0 and LiBlobEm[m-1,shift*2,Refdec_ind]==0: # condition for counting an event as a blob (0 1 1) or (0 1 0)
									Count2=Count2+1
									if (remind2[Count2-1]>=m-WinBar):								# safty time intervall to not count one blob as two
										intermediatestore = remind2[Count2-1]							# just for the entrance-condition to calculate emission time intervall
										Count2=Count2-1
									if not intermediatestore>=m-WinBar:								# enter only if that blob is not part of the first blob
										for l in range (0, int(TWind_long/avt)):							# go through time window
											for k in range (0, len(x)):							# go through all x-positions within time window
												LiConBlo_long[Count2-1,l,k]=Li2p1Z[(m-int(TWind_long/avt/2)+l),int(shift*2),k]	# store emission values within time window in LiConBlo
						#				print('Blob # {0} can be found at t={1:.2f}ms'.format(Count2,time[m]))
									remind2[Count2]=m

							# calculate average over all blob events --> Standard Blob in LiConAv #######################################################
							for m in range (0,len(x)):
								for l in range (0, int(TWind_long/avt)):
									LiConAv_long[l,m] = np.mean(LiConBlo_long[:,l,m])	 # calculate the values for the standard blob and put in array
									l=l+1
								m=m+1

							Delt_long=[None]*int(TWind_long/avt)				# create emty list for standard blob time window

							NullPos_ind_long = 0						# initialize to be zero (just for checking afterwards, that it is not any more
							for m in range (0,int(TWind_long/avt)):			# create a vector containing the times of the time window for standard blob
								Delt_long[m]=-TWind_long/2+avt*m
								if (Delt_long[m]==0.0):
									NullPos_ind_long = m				# NullPos_ind is index, used for calculating blobwidth at HWHM#
								m=m+1
						if RealCase:
							TWind_long = TWind
							Delt_long = Delt
							NullPos_ind_long = NullPos_ind
							LiConBlo_long = LiConBlo
							LiConAv_long = LiConAv


					# 	Define Time Window of Standard Blob Analysis ##############################################################################
						
						LookWinMin = Delt[0]
						LookWinMin_ind = 0
						LookWinMax = Delt[-1]
						LookWinMax_ind = len(Delt)-1

					# 	only for long time window for tau_B calculations
						LookWinMin_long = Delt_long[0]
						LookWinMin_ind_long = 0
						LookWinMax_long = Delt_long[-1]
						LookWinMax_ind_long = len(Delt_long)-1



					# 	interpolate data for block-case evaluation for calculation of FWHM (correlation time)
						if RealCase:
							a_an=a_an_real
							b_an=b_an_real
							c_an=c_an_real
						if BlockCase or RealCase:
							if x[-1]<x[0]:
								x_inter = np.linspace(x[-1],x[0],300)				# Attenation: axis has to go from low to high values!
							if x[-1]>x[0]:
								x_inter = np.linspace(x[0],x[-1],300)				# Attenation: axis has to go from low to high values!
							x_here=[np.NaN]*len(x)							# introduce new x-axis going from negative to positive values for spline interpolation
							LiConAv_here_a = [np.NaN]*len(x)
							LiConAv_here_b = [np.NaN]*len(x)
							LiConAv_here_c = [np.NaN]*len(x)
							LiConAv_here_Null = [np.NaN]*len(x)
							for i in range(len(x)):
								if x[-1]<x[0]:
									x_here[len(x)-i-1]=x[i]					# Attenation: axis has to go from low to high values!
									LiConAv_here_a[len(x)-i-1] = LiConAv[a_an,i]
									LiConAv_here_b[len(x)-i-1] = LiConAv[b_an,i]
									LiConAv_here_c[len(x)-i-1] = LiConAv[c_an,i]
									LiConAv_here_Null[len(x)-i-1] = LiConAv[NullPos_ind,i]
								if x[-1]>x[0]:
									x_here[i]=x[i]					# Attenation: axis has to go from low to high values!
									LiConAv_here_a[i] = LiConAv[a_an,i]
									LiConAv_here_b[i] = LiConAv[b_an,i]
									LiConAv_here_c[i] = LiConAv[c_an,i]
									LiConAv_here_Null[i] = LiConAv[NullPos_ind,i]
							LiConAv_smooth_a = spline(x_here,LiConAv_here_a,x_inter)
							LiConAv_smooth_b = spline(x_here,LiConAv_here_b,x_inter)
							LiConAv_smooth_c = spline(x_here,LiConAv_here_c,x_inter)
							LiConAv_smooth_Null = spline(x_here,LiConAv_here_Null,x_inter)

					# 	Calculate HWHM of standard blob (more work, since curve not necessarily gausian!) ##########################################

						if BlockCase or RealCase:
							LiConAv_max = max(LiConAv_smooth_Null[:])		# Maximum emission value
							for c in range (0,len(x_inter)):
								if (LiConAv_max==LiConAv_smooth_Null[c]):
									max_ind = c				# index on x-axis of maximum emission value
									
							# find closest index to half-max-position (lower index)	
							PeakC = [None]*100					# peak counter to determine for sure the position of half-max
							RandCount = 1
							for d in range (max_ind):
								if (LiConAv_smooth_Null[d]>LiConAv_max/2 and LiConAv_smooth_Null[d-1]<LiConAv_max/2):
									PeakC[RandCount]=int(d)			# store index in peak counter array (in best case only one entry peakC[1]!)
									RandCount = RandCount+1
							if max(PeakC)>=0:
								low_index = int(max(PeakC))					# lower index is maximum for last time it surpasses HM
							if max(PeakC)<0:
								low_index = np.nan 

							# find closest index to half-max-position (upper index)
							PeakC = [None]*100					# peak counter array to determine for sure the position of half-max
							PeakC = np.array(PeakC)
							PeakC.fill(10000)					# fill with high values (higher than len(x), since we need the lowest one
							RandCount = 1
							for e in range (max_ind+1,len(x_inter)):
								if (LiConAv_smooth_Null[e-1]>LiConAv_max/2 and LiConAv_smooth_Null[e]<LiConAv_max/2):
									PeakC[RandCount]=int(e)			# store index in peak counter array (in best case only one entry peakC[1]!)
									RandCount = RandCount+1	
							up_index = int(min(PeakC))					# up index is minimum index for first time it surpasses HM
							
							if (up_index>5000 or type(low_index)!=int or type(up_index)!=int):			# exit loop for this analysis if a proper calculation is not possible here

								blub=999
								print('blub in up/low_index in Delx BlockCase')
								if not BlockCase:
									shift_Block = 0
								Counterbad, blub, Delx, tau_B, Blobf,vr, rval, vdmax, vimax = write_nan(fu, shift_Block, Measurement,maxlen, NewData, Radial, WriteHeader, DenCase, EmCase, BlockCase, Counterbad, blub, t_start,t_end, x,y,Refdec_ind, shift, BlobCount, Delx, tau_B, Blobf,vr, rval, vdmax, vimax, B0, Lp, q95, Te0, Ti0, ne0, omegaCi, rhoS, endtime, LiConAvrel)
								continue
							
							Delx = abs(x_inter[low_index]-x_inter[up_index])/2.	# HWHM of intensity at Delt = 0.
							print('HWHM of standard blob (interpolated data): {0:.3f}cm'.format(Delx))
							
	#					if blub==999:							# exit analysis because of bad result
	#						print(blub, 'one level deeper')
	#						break

						if EmCase or DenCase:
							LiConAv_max = max(LiConAv[NullPos_ind,:])		# Maximum emission value
							for c in range (0,len(x)):
								if (LiConAv_max==LiConAv[NullPos_ind,c]):
									max_ind = c				# index on x-axis of maximum emission value

							# find closest index to half-max-position (lower index)	
							PeakC = [None]*100					# peak counter to determine for sure the position of half-max
							RandCount = 1
							for d in range (max_ind):
								if (LiConAv[NullPos_ind,d]>LiConAv_max/2 and LiConAv[NullPos_ind,d-1]<LiConAv_max/2):
									PeakC[RandCount]=int(d)			# store index in peak counter array (in best case only one entry peakC[1]!)
									RandCount = RandCount+1
							if (max(PeakC)>=0):			# go only in if there is a real value to calculate (otherwise PeakC is empty)
								low_index = int(max(PeakC))						# lower index is maximum for last time it surpasses HM
							else:
								low_index=np.nan

							# find closest index to half-max-position (upper index)
							PeakC = [None]*100					# peak counter array to determine for sure the position of half-max
							PeakC = np.array(PeakC)
							PeakC.fill(10000)					# fill with high values (higher than len(x), since we need the lowest one
							RandCount = 1
							for e in range (max_ind+1,len(x)):
								if (LiConAv[NullPos_ind,e-1]>LiConAv_max/2 and LiConAv[NullPos_ind,e]<LiConAv_max/2):
									PeakC[RandCount]=int(e)			# store index in peak counter array (in best case only one entry peakC[1]!)
									RandCount = RandCount+1	
							up_index = int(min(PeakC))					# up index is minimum index for first time it surpasses HM
							print('low index', low_index, 'up_index', up_index)
							
							if (up_index>5000 or low_index<1 or type(low_index)!=int):			# exit loop for this analysis if a proper calculation is not possible here
								blub=999
								print('blub in up/low_index in Delx EmCase')
							#	pdb.set_trace()
								if not BlockCase:
									shift_Block = 0
								Counterbad, blub, Delx, tau_B, Blobf,vr, rval, vdmax, vimax = write_nan(fu, shift_Block, Measurement,maxlen, NewData, Radial, WriteHeader, DenCase, EmCase, BlockCase, Counterbad, blub, t_start,t_end, x,y,Refdec_ind, shift, BlobCount, Delx, tau_B, Blobf,vr, rval, vdmax, vimax, B0, Lp, q95, Te0, Ti0, ne0, omegaCi, rhoS, endtime, LiConAvrel)
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

					# 	only for long time window for tau_B calculation
						TWind_inter_long = np.linspace(LookWinMin_long, LookWinMax_long,900)
						LiConAv_self_smooth_long = spline(Delt_long, LiConAv_long[:,Refdec_ind],TWind_inter_long)

						LiConAv_self_long = max(LiConAv_self_smooth_long[:])		# maximum of emission for longer time window
						LiConAv_self = max(LiConAv_self_smooth[:])
						RandCount=1 
						for c in range (0,len(TWind_inter_long)):
							if (LiConAv_self_long==LiConAv_self_smooth_long[c]):
								self_ind = c				# index on time-window-axis of maximum emission value

						# find closest index to half-max-position (lower index)	
						PeakC = [None]*100					# peak counter to determine for sure the position of half-max
						RandCount = 1
						for d in range (self_ind):
							if (LiConAv_self_smooth_long[d]>LiConAv_self_long/2 and LiConAv_self_smooth_long[d-1]<LiConAv_self_long/2):
								PeakC[RandCount]=int(d)			# store index in peak counter array (in best case only one entry peakC[1]!)
								RandCount = RandCount+1
						if max(PeakC)>=0:
							left_index = int(max(PeakC))					# lower index is maximum for last time it surpasses HM
						if max(PeakC)<0:
							left_index = 0					# lower index is maximum for last time it surpasses HM


						# find closest index to half-max-position (upper index)
						PeakC = [None]*100					# peak counter array to determine for sure the position of half-max
						PeakC = np.array(PeakC)
						PeakC.fill(10000)					# fill with high values (higher than len(x), since we need the lowest one
						RandCount = 1
						for e in range (self_ind+1,len(TWind_inter_long)):
							if (LiConAv_self_smooth_long[e-1]>LiConAv_self_long/2 and LiConAv_self_smooth_long[e]<LiConAv_self_long/2):
								PeakC[RandCount]=int(e)			# store index in peak counter array (in best case only one entry peakC[1]!)
								RandCount = RandCount+1	
						right_index = int(min(PeakC))				# up index is minimum index for first time it surpasses HM
						if (right_index>5000 or left_index==0):
							blub==999	
							print('blub in left/right index of tauB')			# exit analysis because of bad result
						#	pdb.set_trace()
							if not BlockCase:
								shift_Block = 0
							Counterbad, blub, Delx, tau_B, Blobf,vr, rval, vdmax, vimax = write_nan(fu, shift_Block, Measurement,maxlen, NewData, Radial, WriteHeader, DenCase, EmCase, BlockCase, Counterbad, blub, t_start,t_end, x,y,Refdec_ind, shift, BlobCount, Delx, tau_B, Blobf,vr, rval, vdmax, vimax, B0, Lp, q95, Te0, Ti0, ne0, omegaCi, rhoS, endtime, LiConAvrel)
							tau_B = 0
							continue
						if (right_index==10000 or left_index==0):
							sys.exit('Please change time-window, because HWHM-calculation could not take place')
						tau_B = (TWind_inter_long[right_index]-TWind_inter_long[left_index])/2.	# HWHM of intensity at Delt = 0.
			#			print('Self correlation time of standard blob: {0:.3f}ms'.format(tau_B))
						self_TWind_inter = TWind_inter
						self_TWind_inter_long = TWind_inter_long


					#	Calculate Blobfrequency #########################################################################################################
						Blobf = BlobCount/(time[-1]-time[0])*1000
						print('The Blobfrequency is: {0:d}'.format(int(Blobf)))


					# 	Calculate center of mass (COM) for blobs at different t in TWind ######################################################################

						# set all negative values of LiConAv to zero -- it should not count for calculation of COM
						for l in range (0,len(LiConAv[:,1])):
							for k in range (0,len(LiConAv[1,:])):
								if (LiConAv[l,k]<0):
									LiConAv[l,k]=0
						COM = [None]*len(LiConAv[:,1])		# create COM-vector with appropriate shape
						for m in range (0,int(TWind/avt)):
							#Blubern = ndimage.measurements.center_of_mass(LiConAv[m,:])		# Provides tuple (x and y) --> Blubern is intermediate to hold tuple and not used, while COM[m] holds the index of the COM on the x-axis
							#if np.isnan(Blubern[0]) or Blubern<0 or Blubern[0]+1==len(x):

							#	blub=999
							#	print('blub in left/right_index in COM')
						#	#	pdb.set_trace()
							#	if not BlockCase:
							#		shift_Block = 0	
							#	Counterbad, blub, Delx, tau_B, Blobf,vr, rval, vdmax, vimax = write_nan(fu, shift_Block, Measurement,maxlen, NewData, Radial, WriteHeader, DenCase, EmCase, BlockCase, Counterbad, blub, t_start,t_end, x,y,Refdec_ind, shift, BlobCount, Delx, tau_B, Blobf,vr, rval, vdmax, vimax, B0, Lp, q95, Te0, Ti0, ne0, omegaCi, rhoS, endtime, LiConAvrel)
							#	blub=999
							#	break
							#if m == 0:
							#	startblub = int(Blubern[0])
							# wrong: COM[m] = x[int(Blubern[0])]-abs(x[int(Blubern[0])]-x[int(Blubern[0])+1])*(Blubern[0]-int(Blubern[0]))						# use only first part of tuple
							#if x[0]>x[-1]:
							#	COM[m] = x[startblub]-abs(x[int(Blubern[0])]-x[int(Blubern[0])+1])*(Blubern[0]-startblub)		# use only first part of tuple and calculate the effective COM-Position
							#if x[0]<x[-1]:
							#	COM[m] = x[startblub]+abs(x[int(Blubern[0])+1]-x[int(Blubern[0])])*(Blubern[0]-startblub)		# use only first part of tuple and calculate the effective COM-Position
							COM[m] = np.sum(x*LiConAv[m,:])/np.sum(LiConAv[m,:])
						#if np.isnan(Blubern[0]) or blub == 999:
						#	blub = 0
						#	continue
						COM = np.array(COM)								# convert list to array

					# 	Calculation of relative fluction of intensity in spectrum or density ######################################################

						LiConAv_maxrel = max(LiConAv[:,Refdec_ind])		# Maximum emission value
						LiConAv_mean = np.mean(LiConAv[:Refdec_ind])		# Mean of emission value at this reference detector
						LiConAvrel = LiConAv_maxrel/LiAv			# relative flucutation of intensity/density
						print('Relative fluctuation of intensity/density for standard blob: ', LiConAvrel)
						print('Mean fluctuation of intensity/density for standard blob: ', LiConAv_mean)
						print('Max fluctuation of intensity/density for standard blob: ', LiConAv_maxrel)



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
							ax21.plot([self_TWind_inter_long[left_index],self_TWind_inter_long[right_index]],[x[Refdec_ind],x[Refdec_ind]],linewidth=3, color = 'w')
							ax21.set_xlabel('time $\Delta$t')			
							ax21.set_ylabel('beam axis $x$ (cm)')
							ax21.set_title('Spectrum for standard blob')
							#time points for representatives blobs are drawn into plot
							ax21.axvline(Delt[a_an],color='r', linestyle='--')
							ax21.axvline(Delt[b_an],color='g', linestyle='--')
							ax21.axvline(Delt[c_an],color='b', linestyle='--')
							ax21.axvline(Delt[NullPos_ind], color='w', linestyle='-.')						# Null pos v-oline
							# ax1.set_xlim(xmin,xmax)
							ax21.set_ylim(LookWinMin_long,LookWinMax_ind)
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
							ax22.set_xlabel('beam axis $x$ (cm)')			
							ax22.set_ylabel('I(x,$\Delta$t) (a.u.)',labelpad = -10)
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

						#	calculate speed only for area, where blob really is!
						if RealCase:
							#for better velocity calculation:
							LookWinMin	= -0.02
							LookWinMax	=  0.08
							for i in range(len(Delt)):
								if Delt[i]>LookWinMin:
									LookWinMin_ind = i-1
									break
							for i in range(len(Delt)):
								if Delt[i]>LookWinMax:
									LookWinMax_ind = i-1
									break
								

						if not RealCase:
							#for better velocity calculation:
							LookWinMin	= -0.01
							LookWinMax	=  0.03
							for i in range(len(Delt)):
								if Delt[i]>LookWinMin:
									LookWinMin_ind = i-1
									break
							for i in range(len(Delt)):
								if Delt[i]>LookWinMax:
									LookWinMax_ind = i-1
									break


						# calculate lienar regression statistics:
						v_r, v_intercept, r_value, p_value, std_err = stats.linregress(Delt[LookWinMin_ind:LookWinMax_ind],COM[LookWinMin_ind:LookWinMax_ind])

						#create plot information of linear regression
						(a_lgr,b_lgr) =polyfit(Delt[LookWinMin_ind:LookWinMax_ind],COM[LookWinMin_ind:LookWinMax_ind],1) 
						COM_lgr = np.polyval([a_lgr,b_lgr],Delt[LookWinMin_ind:LookWinMax_ind])

						# interpolate data for better calculation of slope
						if not RealCase:
							TWind_inter = np.linspace(LookWinMin,Delt[LookWinMax_ind+1],300)		
							XC_COM_smooth = spline(Delt[LookWinMin_ind:LookWinMax_ind+1],COM[LookWinMin_ind:LookWinMax_ind+1],TWind_inter)
						if RealCase:
							TWind_inter = np.linspace(LookWinMin,Delt[LookWinMax_ind+1],300)		
							XC_COM_smooth = spline(Delt[LookWinMin_ind:LookWinMax_ind+4],COM[LookWinMin_ind:LookWinMax_ind+4],TWind_inter)

						#calculate maximum slopes of real data and interpolated data:
						slope_data=[None]*(LookWinMax_ind-LookWinMin_ind)
						slope_inter=[None]*(len(TWind_inter))
						slope_data=np.array(slope_data)						# convert list to array
						slope_inter=np.array(slope_inter)					# convert list to array
						slope_data[0]=0
						slope_inter[0]=0
						for q in range (LookWinMin_ind+1,LookWinMax_ind):			# Calculate slopes between every point
							if x[0]>x[-1]:
								slope_data[q-LookWinMin_ind] = (COM[q]-COM[q-1])/abs(Delt[q]-Delt[q-1])
								if (slope_data[q-LookWinMin_ind]<0):
									slope_data[q-LookWinMin_ind]=0
							if x[0]<x[-1]:
								slope_data[q-LookWinMin_ind] = -(COM[q]-COM[q-1])/abs(Delt[q]-Delt[q-1])
								if (slope_data[q-LookWinMin_ind]<0):
									slope_data[q-LookWinMin_ind]=0

						v_dmax=np.amax(slope_data)						# find maximum slope of data 
						for p in range (0,len(TWind_inter)):					# Calculate slopes between every point
							if x[0]>x[-1]:
								slope_inter[p] = (XC_COM_smooth[p]-XC_COM_smooth[p-1])/abs(TWind_inter[p]-TWind_inter[p-1])
								if (slope_inter[p]<0):
									slope_inter[p]=0
							if x[0]<x[-1]:
								slope_inter[p] = -(XC_COM_smooth[p]-XC_COM_smooth[p-1])/abs(TWind_inter[p]-TWind_inter[p-1])
								if (slope_inter[p]<0):
									slope_inter[p]=0
						v_imax=np.amax(slope_inter)						# find maximum slope interpolation

						#convert into [m/s]
						if x[0]>x[-1]:
							vr = v_r*0.01*pow(10,3)
							rval = r_value*0.01*pow(10,3)
						if x[0]<x[-1]:
							vr = -v_r*0.01*pow(10,3)
							rval = -r_value*0.01*pow(10,3)
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
							ax23.set_xlabel('time $\Delta$t')			
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

							Em_Stan = ax24.plot(Delt_long,LiConAv_long[:,Refdec_ind], color = 'k')				# plot original data
							Em_Stan_inter = ax24.plot(self_TWind_inter_long, LiConAv_self_smooth_long, color = 'm')	# plot interpolated data
							ax24.plot([self_TWind_inter_long[left_index-1],self_TWind_inter_long[right_index]],[LiConAv_self_smooth_long[left_index-1],LiConAv_self_smooth[right_index]],linewidth=3, color = 'm')
							ax24.set_xlabel('time $\Delta$t')			
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
						if EmCase or DenCase and BlobCount>=0:
							data = [("00"+str(Measurement), ""), (t_start,".2f"), (t_end-t_start, ".2f"), (x[Refdec_ind], ".2f"), (y[shift*2],".2f"), (BlobCount, "d"), (Delx, ".3f"), (tau_B,".4f"), (Blobf,".1f"), (vr, ".1f"), (rval,".1f"), (vdmax,".1f"), (vimax,".1f"), (LiConAvrel,".2f"), (B0,".2f"), (Lp,".2f"), (q95,".2f"), (Te0,".2f"), (Ti0,".2f"), (ne0,".2e"), (omegaCi,".2e"), (rhoS,".6f"), (endtime,".2f")]
			#			if not BlockCase and np.isnan(BlobCount):
			#				data = [("00"+str(Measurement), ""), (t_start,".2f"), (t_end-t_start, ".2f"), (x[Refdec_ind], ".2f"), (y[shift*2],".2f"), (BlobCount, ".1f"), (Delx, ".1f"), (tau_B,".1f"), (Blobf,".1f"), (vr, ".1f"), (rval,".1f"), (vdmax,".1f"), (vimax,".1f"), (LiConAvrel,".2f"), (B0,".2f"), (Lp,".2f"), (q95,".2f"), (Te0,".2f"), (Ti0,".2f"), (ne0,".2e"), (omegaCi,".2e"), (rhoS,".6f"), (endtime,".2f")]
						if BlockCase and BlobCount>=0:
							data = [("00"+str(Measurement), ""), (t_start,".2f"), (t_end-t_start, ".2f"), (x[Refdec_ind], ".2f"), (y[shift_Block*2],".2f"), (BlobCount, "d"), (Delx, ".3f"), (tau_B,".4f"), (Blobf,".1f"), (vr, ".1f"), (rval,".1f"), (vdmax,".1f"), (vimax,".1f"), (LiConAvrel,".2f"), (B0,".2f"), (Lp,".2f"), (q95,".2f"), (Te0,".2f"), (Ti0,".2f"), (ne0,".2e"), (omegaCi,".2e"), (rhoS,".6f"), (endtime,".2f")]
						if RealCase and BlobCount>=0:
							data = [(str(shot), ""), (t_start,".2f"), (t_end-t_start, ".2f"), (x[Refdec_ind], ".2f"), (t_end,".2f"), (BlobCount, "d"), (Delx, ".3f"), (tau_B,".4f"), (Blobf,".1f"), (vr, ".1f"), (rval,".1f"), (vdmax,".1f"), (vimax,".1f"), (LiConAvrel,".2f"), (B0,".2f"), (Lp,".2f"), (q95,".2f"), (Te0,".2f"), (Ti0,".2f"), (ne0,".2e"), (omegaCi,".2e"), (rhoS,".6f"), (endtime,".2f")]
			#			
			#			if np.isnan(BlobCount) and BlockCase:
			#				data = [("00"+str(Measurement), ""), (t_start,".2f"), (t_end-t_start, ".2f"), (x[Refdec_ind], ".2f"), (y[shift_Block*2],".2f"), (BlobCount, ".1f"), (Delx, ".1f"), (tau_B,".1f"), (Blobf,".1f"), (vr, ".1f"), (rval,".1f"), (vdmax,".1f"), (vimax,".1f"), (LiConAvrel,".2f"), (B0,".2f"), (Lp,".2f"), (q95,".2f"), (Te0,".2f"), (Ti0,".2f"), (ne0,".2e"), (omegaCi,".2e"), (rhoS,".6f"), (endtime,".2f")]

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

