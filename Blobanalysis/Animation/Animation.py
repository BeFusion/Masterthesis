#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from pylab import plot
import matplotlib.gridspec as gridspec
from PIL import Image
import moviepy.editor as mpy				#creation of video
from moviepy.video.io.bindings import mplfig_to_npimage	#for conversion to RGB image of figure
from hdf5_read import Li_read_hdf5			# read hdf5 file in
from block_func import block_beam_func			# for calculating the block data
from dect_func import dect_func
from Blob_params import *				# set a few parameters 


#################################################################################################################################



########### READING DATA AND PRELIMINARY CALCULATIONS ###########################################################################
Measurement=4
print('Data is going to be read from file: /hesel/data/ASDEX.{0:03d}.h5'.format(Measurement))
# read in data from function Li_read_hdf5: Specify correct hdf5-file
time, x, y, ne, LiTemp, Li2p1, B0, Lp, q95, Te0, Ti0, ne0, omegaCi, rhoS, endtime = Li_read_hdf5('/pfs/home/bschiess/public/hesel/data/ASDEX.{0:03d}.h5'.format(Measurement))

# selected times, which are supposed to be investigated [us]
t_start=0.0
t_end = 0.6
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

#define center of beam and shift for correct beam position:
shift=int(bcenter/y[-1]*stepsize/2)+1		# shift of beam center in order to place the center at correct position
shift2 = int(b/2/y[-1]*stepsize)		# shift*2-shift2 is now left index of beam edge and shift*2+shift2 the right one
print('Position of Beam: {0:.2f}cm'.format(y[shift*2]))
print('Beamwidth: {0:.2f}cm'.format(y[shift2+shift2]))	

# use the block-emission data for the block case:
x_ori = x
#	Li2p1 = block_func(Li2p1,timestep,mn,b,avy,shift,shift2)			# if original block-function should be used
Li2p1_beam, Li2p1_block, x_block, blur = dect_func(Li2p1,timestep,mn,b,avy,shift,shift2,x)	# if detector reduced block function should be used
x_axis = x
mn = len(x)
shift_Block = shift
shift = 0									# we are in the 1D-Case, so every time the entry at shift is called in matrix, this is set to the zero entry


# Calculate index and position of Reference detector on xaxis (indices are counted backwards x[len(x)-1]<0!)
print('Reference detector has to be shifted to a position available for beam evaluation!')
countx=0
for g in range (0,len(x_block)):
	if (x_block[g]<0):				# take negative values into account for correct shift
		countx=countx+1			# number of elements, which are negative
Refdec_ind_block = len(x)-int(Refdec/blur)-1 - countx	# index for Reference detector position
print('Position along beam-axis for Reference Detector in Block-Case: {0:.2f}cm'.format(x[Refdec_ind])) # check, if Reference detector is at correct position
print('Position of Reference Detector set in Parameter-file: {0:.2f}cm'.format(Refdec)) # check, if Reference detector is at correct position
if (Refdec>x[0]):				# exit, if error with reference detector position occurs
	sys.exit('Reference Detector outside of observable region')

countx=0
for g in range (0,len(x)):
	if (x[g]<0):				# take negative values into account for correct shift
		countx=countx+1			# number of elements, which are negative
Refdec_ind = len(x)-Refdec/0.02-1 - countx	# index for Reference detector position
print('Position along beam-axis for Reference Detector: {0:.2f}cm'.format(x[Refdec_ind])) # check, if Reference detector is at correct position

#Meshgrids
X,Y=np.meshgrid(x,y)						# transform x,y-axis in matrix for contourf-plot of 2D-plots
X_Ori,Y_block = np.meshgrid(x_ori,y)				# only needed, if BlockCase changes x-axis
X2,T2=np.meshgrid(x,time[int(t_start/avt):int(t_end/avt)]) 	# Meshgrid for timeplots:



##################################################################################################################################




########### INITIALIZATION OF PLOTS AND CREATING SUBPLOTS 1 & 2 ##################################################################


for t_ind in range(timestep)
	if time[t_ind]>t_start
		t_start_ind = t_ind
		break
for t_ind in range(t_start_ind,timestep)
	if time[t_ind]>t_end
		t_end_ind = t_ind
		break

for position in range(t_start_ind,t_end_ind):

	# create fiugre with certain ratio and facecolor:
	f1=plt.figure(figsize=(20,20), facecolor = 'white')
	#f1.suptitle("t=%timestep [1/$\Omega_i$]".format(times[t]))	# set time step title
	gs = gridspec.GridSpec(1, 3,				# ratio of grid space (2 plots per collumn, 3 per row)
 	                      width_ratios=[1],		# width ratios of the 3 plots per row
 	                      height_ratios=[1,1,1]		# height ratios of the 2 polots per collumn
	                       )

	#define limits for plots
	xmin=0.6
	xmax=1.8
	ymin=0
	ymax=8


	# create axes with the ratios specified in gs above. ax5, ax6 are the axes for the colorbars:
	ax1 = plt.subplot(gs[0])
	ax2 = plt.subplot(gs[1])
	ax3 = plt.subplot(gs[2])


	# create a list in for simple modifications on all plots at the same time:
	allax=[ax1,ax2,ax3]
	font = {'family' :'normal', 'weight' : 'regular', 'size' : 16}
	matplotlib.rc('font', **font)
	# first subplot: density
	density = ax1.contourf(Y,X,ne[position,:,:],300)
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

	


	f1.savefig('image%04d.png' % (position))
	##################################################################################################################################

