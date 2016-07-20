#!/usr/bin/env python
# -*- coding: utf-8 -*-


########### Importing modules ###########################################################################

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#from sympy import *
from pylab import plot
import matplotlib.gridspec as gridspec		# define layout of grid
from mpl_toolkits.axes_grid1 import make_axes_locatable # for colorbars
from scipy import ndimage			# calculate com
from scipy import signal			# for smoothing time signal
from scipy.interpolate import spline		# for interpolate COM position for good speed calculation of blob
from scipy import stats, polyfit		# for linear interpolation data
from scipy_cookbook_smooth import smooth	# for interpolating data
from hdf5_read import Li_read_hdf5		# read hdf5 file in
from block_func import block_beam_func		# for calculating the block data
from dect_func import dect_func
from Blob_params import *			# set a few parameters 
import timeit					# for timekeeping
from datetime import datetime			# for date and time printing in output file
import sys					# for exiting the program if errors occur
import pdb					# for error tracking (set_trace at right position)
from matplotlib.offsetbox import AnchoredText
import os.path

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
#
#
#
#
#
#
#

###############################################################################################################################
########### READING DATA AND PRELIMINARY settings #############################################################################
###############################################################################################################################

#use specific shots
Measurement=777
#shot = 29315
shot = 29306


# Select cases and switches ###################################################################################################

# Evaluates the 'infinite' resolved beam emission, the density evaluation or the Block evaluation case


EmCase = True 		# infinite emission case --> actually always true, since every different setting during analysis is made by if-statements with the other two cases --> can stay false
DenCase = True		# density analysis case
BlockCase = True			# block-averaging emission case
RealCase = False	# for reading in real emission data from ASDEX

realSNR = False			# reads in SNR profile for the Block Case

NewData = False				# Switch to write output (only data for this one run! -->  usually not needed, since the data is generated in radial run)

saving = True				# save images
if not saving:
	print('figures are not saved!!')

Noise = False			# switch on additive white Gaussian noise --> Very time consuming!
NoiseFast = False			# Enhances speed, if needed for quicker calculation
Smooth = False			# smooth noisy data
NoiseTime = False
SmoothTime = False

statist = True			# puts out a statistical analysis plot

if BlockCase:
	realSNR = True		# reads in SNR profile for the Block Case
	NoiseTime = True
	SmoothTime = True
	Noise = False			# switch on additive white Gaussian noise --> Very time consuming!
	NoiseFast = False			# Enhances speed, if needed for quicker calculation	
	Smooth = False			# smooth noisy data



# data reading  ################################################################################################################


if not RealCase:
						# will be used for file search and in output
	print('Data is going to be read from file: ASDEX.{0:03d}.h5'.format(Measurement))
	# read in data from function Li_read_hdf5: Specify correct hdf5-file
	#time, x, y, ne, LiTemp, Li2p1, B0, Lp, q95, Te0, Ti0, ne0, omegaCi, rhoS, endtime = Li_read_hdf5('/pfs/home/bschiess/public/hesel/data/ASDEX.{0:03d}.h5'.format(Measurement))
	time, x, y, ne, LiTemp, Li2p1, B0, Lp, q95, Te0, Ti0, ne0, omegaCi, rhoS, endtime = Li_read_hdf5('/pfs/home/bschiess/itmwork/hesel_data/ASDEX.{0:03d}.h5'.format(Measurement))
	# other settings and preparations  #############################################################################################

if not RealCase:
	# selected times, which are supposed to be investigated [us]
	t_start=0.2
	t_end = time[-1]
	print('time to be investigated betweem %.3gms and %.3gms:' % (t_start, t_end))


if not RealCase:
	# calculate time in [us] for snapshots:
	point=position*(time[11]-time[10])

# some lengthes used during the program:
mn=len(x)
x_axis = x					# just useful sometimes to avoid confusion
timestep=len(time)				# misleading name, sorry
avt=(time[timestep-1]-time[0])/timestep

if not RealCase:
	stepsize=len(y)					# misleading name, sorry
	avy=y[stepsize-1]/stepsize	

	#define center of beam and shift for correct beam position:
	shift=int(bcenter/y[-1]*stepsize/2)+1		# shift of beam center in order to place the center at correct position
	shift2 = int(b/2/y[-1]*stepsize)		# shift*2-shift2 is now left index of beam edge and shift*2+shift2 the right one
	print('Position of Beam: {0:.2f}cm'.format(y[shift*2]))
	print('Beamwidth: {0:.2f}cm'.format(y[shift2+shift2]))	
	
Li2p1raw = Li2p1.copy()	# raw Li2p1-signal for all other non-noise evaulations

if BlockCase:
	x_ori = x
#	Li2p1 = block_func(Li2p1,timestep,mn,b,avy,shift,shift2)			# if original block-function should be used
	Li2p1_beam, Li2p1, x, blur = dect_func(Li2p1,timestep,mn,b,avy,shift,shift2,x)	# if detector reduced block function should be used
	x_axis = x
	mn = len(x)
	shift_Block = shift
	shift = 0									# we are in the 1D-Case, so every time the entry at shift is called in matrix, this is set to the zero entry

	if NoiseTime:
		
		# use realistic profile for noise for BlockCase
		if realSNR:
			xSNR, rSNR = np.loadtxt('SNR_29302.txt', usecols = (0,1), unpack = True, skiprows=1)
			maxSNR = max(rSNR)
			for i in range(len(rSNR)):
				if maxSNR == rSNR[i]:
					xmaxSNR = xSNR[maxSNR]
					maxSNR_ind = i 

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

		Li2p1noise=Li2p1.copy()
		for x_ind in range(len(x)):
			if realSNR:
				noiLi = np.random.normal(0,np.mean(Li2p1[:,0,x_ind])/SNR[x_ind],timestep)
			if not realSNR:
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
		Li2p1blockraw=Li2p1noise.copy()
		if SmoothTime:
			Li2p1blocksmooth = Li2p1hn2.copy()

#Meshgrids
if not RealCase:
	X,Y=np.meshgrid(x_ori,y)						# transform x,y-axis in matrix for contourf-plot of 2D-plots
if BlockCase:
	X_Block,Y_Block = np.meshgrid(x,y)				# only needed, if BlockCase changes x-axis



###################################################################################################################################
########### FIGURE 1 ##############################################################################################################
###################################################################################################################################


# create fiugre with certain ratio and facecolor:
if not RealCase:
	f1=plt.figure(figsize=(16,16), facecolor = 'white')
	matplotlib.rcParams.update({'font.size': 27})

	#define limits for plots
	xmin=2
	xmax=4
	
	lwid = 2

if not RealCase:

	plot(time,ne[:,112,112]/np.mean(ne[:,112,112]), label = 'Density Case', linewidth = lwid)

	plot(time,Li2p1raw[:,112,112]/np.mean(Li2p1raw[:,112,112])+1, label = 'Emission Case', linewidth = lwid)

	plot(time,Li2p1[:,0,3]/np.mean(Li2p1[:,0,3])+2, label = 'Block Case', linewidth = lwid)
	
	plot(time,Li2p1blocksmooth[:,0,3]/np.mean(Li2p1blocksmooth[:,0,3])+3, label = 'Real SNR Block Case', linewidth = lwid)
	

	plt.ylabel('signal (a.u.)', fontsize = 45)
	plt.xlabel('time (ms)', fontsize = 45)
	#plt.title('normalized time signals of density and emission data', y = 1.03)
	plt.xlim(xmin,xmax)
	plt.ylim(0.5,5.8)
	plt.yticks(fontsize = 35)
	plt.xticks(fontsize = 35)


	plt.legend(loc='upper left',fancybox = True, numpoints = 1)



plt.tight_layout()

f1.savefig('All_in_one_1D')

plt.show()
