# -*- coding: utf-8 -*-

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import *
import os
from pylab import plot
import matplotlib.gridspec as gridspec
from matplotlib import animation		#creation of video
from scipy_cookbook_smooth import smooth	# for interpolating data
from hdf5_read import Li_read_hdf5			# read hdf5 file in
from block_func import block_beam_func			# for calculating the block data
from dect_func import dect_func
from Blob_params import *				# set a few parameters 
from dect_func import dect_func
import pdb

#################################################################################################################################

ion()		# switch on interactive mode
Noise = False
Smooth = False

########### READING DATA AND PRELIMINARY CALCULATIONS ###########################################################################
Measurement=4
filename = '{0:03d}Animation_test'.format(Measurement)
print('Data is going to be read from file: /hesel/data/ASDEX.{0:03d}.h5'.format(Measurement))
# read in data from function Li_read_hdf5: Specify correct hdf5-file
time, x, y, ne, LiTemp, Li2p1, B0, Lp, q95, Te0, Ti0, ne0, omegaCi, rhoS, endtime = Li_read_hdf5('/pfs/home/bschiess/public/hesel/data/ASDEX.{0:03d}.h5'.format(Measurement))

# selected times, which are supposed to be investigated [us]
t_start=0.2
t_end = 1.5
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

if Noise:
	for x_ind in range(mn):
		for y_ind in range(stepsize):				# + 2 safty distance
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
		smoothlen = 71
		for t_ind in range(timestep):
			for y_ind in range(stepsize):
				Li2p1hn[t_ind,y_ind,:] = smooth(Li2p1[t_ind,y_ind,:],window_len=smoothlen,window='hanning')
			#	Li2p1hm[t_ind,shift*2,:] = smooth(Li2p1[t_ind,shift*2,:],window_len=smoothlen,window='hamming')
			#	Li2p1ba[t_ind,shift*2,:] = smooth(Li2p1[t_ind,shift*2,:],window_len=smoothlen,window='bartlett')
			#	Li2p1bl[t_ind,shift*2,:] = smooth(Li2p1[t_ind,shift*2,:],window_len=smoothlen,window='blackman')
		Li2p1 = Li2p1hn.copy()

# use the block-emission data for the block case:


#	Li2p1 = block_func(Li2p1,timestep,mn,b,avy,shift,shift2)			# if original block-function should be used
Li2p1_beam, Li2p1_block, x_block, blur = dect_func(Li2p1,timestep,mn,b,avy,shift,shift2,x)	# if detector reduced block function should be used

# Calculate index and position of Reference detector on xaxis (indices are counted backwards x[len(x)-1]<0!)
print('Reference detector has to be shifted to a position available for beam evaluation!')
countx=0
for g in range (0,len(x_block)):
	if (x_block[g]<0):				# take negative values into account for correct shift
		countx=countx+1			# number of elements, which are negative
Refdec_ind_block = len(x)-int(Refdec/blur)-1 - countx	# index for Reference detector position
print('Position along beam-axis for Reference Detector in Block-Case: {0:.2f}cm'.format(x[Refdec_ind_block])) # check, if Reference detector is at correct position
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
#X_Ori,Y_block = np.meshgrid(x_ori,y)				# only needed, if BlockCase changes x-axis
X2,T2=np.meshgrid(x,time[int(t_start/avt):int(t_end/avt)]) 	# Meshgrid for timeplots:



##################################################################################################################################

# Calculations for Block Case:
x_ori = x
#	Li2p1 = block_func(Li2p1,timestep,mn,b,avy,shift,shift2)			# if original block-function should be used
Li2p1_beam, Li2p1_block, x, blur = dect_func(Li2p1,timestep,mn,b,avy,shift,shift2,x)	# if detector reduced block function should be used
x_axis = x
mn = len(x)
shift_Block = shift
shift = 0									# we are in the 1D-Case, so every time the entry at shift is called in matrix, this is set to the zero entry



########### INITIALIZATION OF PLOTS AND CREATING SUBPLOTS ##################################################################

for t_ind in range(timestep):
	if time[t_ind]>t_start:
		t_start_ind = t_ind
		break
for t_ind in range(t_start_ind,timestep):
	if time[t_ind]>t_end:
		t_end_ind = t_ind
		break

files = []
#ax4pl = ax4.contourf(Y,X,[],300)
# start animation-
run =np.linspace(t_start_ind,t_end_ind, t_end_ind-t_start_ind)
print(t_start_ind, time[t_start_ind], t_end_ind, time[t_end_ind], t_end_ind-t_start_ind, len(run))
for run in range(len(run)):
	f1=plt.figure(figsize=(8,8), facecolor = 'white', dpi= 300)
	font = {'family' :'Helvetica', 'weight' : 'regular', 'size' : 12}
	matplotlib.rc('font', **font)
	gs = gridspec.GridSpec(3, 1, width_ratios=[1], height_ratios=[1,1,1])
	gs.update(hspace=0.55)
	# create axes with the ratios specified in gs above. ax5, ax6 are the axes for the colorbars:
	ax1 = plt.subplot(gs[0])
	ax2 = plt.subplot(gs[1])
	ax3 = plt.subplot(gs[2])
#	ax4 = plt.subplot(gs[3])

	xmax = 5
	xmin = min(x_ori)

	pady = 0
#	f1.suptitle("time {0:.2f}ms".format(time[t_start_ind+run]),x=0.8,y=1.00)	# set time step title
#	pdb.set_trace()
	ax1pl = ax1.contourf(Y,X,ne[run,:,:],300)
	ax1.set_xlabel('axis $y$ (cm)', labelpad = pady)			# switched of, since axis in row below
	ax1.set_ylim(xmin,xmax)
	ax1.set_ylabel('beam axis $x$ (cm)', labelpad = pady)
	ax1.set_title(r'Input Density $n_e$ (1/cm$^3$)                        time {0:.2f}ms'.format(time[t_start_ind+run]), y = 1)


	ax2pl = ax2.contourf(Y,X,Li2p1[run,:,:],300)
	ax2.set_ylabel('beam axis $x$ (cm)', labelpad = pady)
	ax2.set_xlabel('axis $y$ (cm)', labelpad = pady)
	ax2.set_ylim(xmin,xmax)
	ax2.set_title('High Resolution Beam Emission Spectrum', y = 1)


	ax3pl = contourf(Y,X,Li2p1_beam[run,:,:],300)
	ax3.set_ylabel('beam axis $x$ (cm)', labelpad = pady)
	ax3.set_xlabel('axis $y$ (cm)', labelpad = pady)
	ax3.set_ylim(xmin,xmax)
	ax3.set_title('Block Beam Emission Spectrum', y = 1)

	fname = '_tmp%04d.png'%run
	files.append(fname)
#	plt.tight_layout()
#	plt.show
#	pdb.set_trace()
	savefig(fname)
#	return ax1pl, ax2pl, ax3pl				# RGB image of the figure


#high quality:
os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=2 -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:mv0:trell:v4mv:cbp:last_pred=3:predia=2:dia=2:vmax_b_frames=2:vb_strategy=1:precmp=2:cmp=2:subcmp=2:preme=2:qns=2 -oac copy -o " + filename + ".mpeg")

# cleanup
for fname in files: os.remove(fname)
	
#anim = animation.FuncAnimation(f1,animate, frames = t_end_ind-t_start_ind, blit = False)
#anim.save('Animation2D_test.mp4',codec='mencoder',fps=30)
##################################################################################################################################




