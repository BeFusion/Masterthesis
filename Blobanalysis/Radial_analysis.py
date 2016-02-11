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
from scipy.interpolate import spline		# for interpolate COM position for good speed calculation of blob
from scipy import stats, polyfit		# for linear interpolation data
from Blob_params import *			# set a few parameters 
import timeit					# for timekeeping
from datetime import datetime			# for date and time printing in output file
import sys					# for exiting the program if errors occur
import pdb					# for error tracking (set_trace at right position)


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

filename = '_veltwind_SNR100_Smooth'				# specify fileextension to be read in and to be stored (e.g '_dens_test', if file is called '004Emresults_dens_test')

Measurement = 4				# specify number of measurement without leading zeros, which is the prefix of the file in form (will be extended later to e.g. 004) 

time, x_pos, y_pos, f_B, rho_x, tau_B, v_r, v_rmist, vdmax, vimax, Relem, B0, Lp,q95 = np.loadtxt('{0:03d}Emresults{1:}.txt'.format(Measurement,filename), usecols = (3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17), unpack = True, skiprows=2)
time2, x_pos2,y_pos2, f_B2, rho_x2, tau_B2, v_r2, v_rmist2, vdmax2, vimax2, Relem2, B02, Lp2,q952 = np.loadtxt('{0:03d}Denresults{1:}.txt'.format(Measurement,filename), usecols = (3, 4, 5,  6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17), unpack = True, skiprows=2)
time3, x_pos3, y_pos3, f_B3, rho_x3, tau_B3, v_r3, v_rmist3, vdmax3, vimax3, Relem3, B03, Lp3,q953 = np.loadtxt('{0:03d}Blockresults{1:}.txt'.format(Measurement,filename), usecols = (3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17), unpack = True, skiprows=2)


#	time, x_pos, y_pos, f_B, rho_x, tau_B, v_r, v_rmist, vdmax, vimax, B0, Lp,q95 = np.loadtxt('{0:03d}Emresults.txt'.format(Measurement), usecols = (3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16), unpack = True, skiprows=2)
#	time2, x_pos2,y_pos2, f_B2, rho_x2, tau_B2, v_r2, v_rmist2, vdmax2, vimax2, B02, Lp2,q952 = np.loadtxt('{0:03d}Denresults.txt'.format(Measurement), usecols = (3, 4, 5,  6, 7, 8, 10, 11, 12, 13, 14, 15, 16), unpack = True, skiprows=2)
#	time3, x_pos3, y_pos3, f_B3, rho_x3, tau_B3, v_r3, v_rmist3, vdmax3, vimax3, B03, Lp3,q953 = np.loadtxt('{0:03d}Blockresults.txt'.format(Measurement), usecols = (3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16), unpack = True, skiprows=2)


###################################################################################################################################
########### Reshaping axes ########################################################################################################
###################################################################################################################################

# for Em-Case #####################################################################################################################
y_min = min(y_pos)
lencount = 0
for elem in range(len(y_pos)):
	if y_pos[elem] == y_min:
		lencount = lencount +1
leny = len(y_pos)/lencount 			# number of y-steps
lenx = lencount 				# number of x-steps
x = [None]*lenx
for i in range (lenx):
	x[i]=x_pos[i]

y = [None]*leny
pos_ind = 0
for i in range(len(y)):
	y[i]= y_pos[pos_ind]
	pos_ind = (i+1)*lenx

# for Den-Case ####################################################################################################################
y_min2 = min(y_pos2)
lencount2 = 0
for elem in range(len(y_pos2)):
	if y_pos2[elem] == y_min2:
		lencount2 = lencount2 +1
leny2 = len(y_pos2)/lencount2			# number of y-steps
lenx2 = lencount2				# number of x-steps
x2 = [None]*lenx2
for i in range (lenx2):
	x2[i]=x_pos2[i]

y2 = [None]*leny2
pos_ind2 = 0
for i in range(len(y2)):
	y2[i]= y_pos2[pos_ind2]
	pos_ind2 = (i+1)*lenx2

# for Block-Case #####################################################################################################################

y_min3 = min(y_pos3)
lencount3 = 0
for elem in range(len(y_pos3)):
	if y_pos3[elem] == y_min3:
		lencount3 = lencount3 +1
leny3 = len(y_pos3)/lencount3			# number of y-steps
lenx3 = lencount3 				# number of x-steps
x3 = [None]*lenx3
for i in range (lenx3):
	x3[i]=x_pos3[i]

y3 = [None]*leny3
pos_ind3 = 0
for i in range(len(y3)):
	y3[i]= y_pos3[pos_ind3]
	pos_ind3 = (i+1)*lenx3

###################################################################################################################################
########### Reshaping axes ########################################################################################################
###################################################################################################################################


# for Em-Case #####################################################################################################################

f_Bm = np.reshape(f_B, (leny, lenx))
rho_xm =np.reshape(rho_x, (leny,lenx))
v_rm =np.reshape(v_r, (leny,lenx))
tau_Bm =np.reshape(tau_B, (leny,lenx))
v_rmistm =np.reshape(v_rmist, (leny,lenx))
vdmaxm =np.reshape(vdmax, (leny,lenx))
vimaxm =np.reshape(vimax, (leny,lenx))

Relemm = np.reshape(Relem, (leny,lenx))

print('shape of f_B_Mat', np.shape(f_Bm))
print('Number of beam positions over which has been averaged: ', lenx, leny, lenx*leny,  len(y_pos), len(f_B), f_B[0],f_B[-1])

# for Den-Case ####################################################################################################################

f_B2m = np.reshape(f_B2, (leny2, lenx2))
rho_x2m =np.reshape(rho_x2, (leny2,lenx2))
v_r2m =np.reshape(v_r2, (leny2,lenx2))
tau_B2m =np.reshape(tau_B2, (leny2,lenx2))
v_rmist2m =np.reshape(v_rmist2, (leny2,lenx2))
vdmax2m =np.reshape(vdmax2, (leny2,lenx2))
vimax2m =np.reshape(vimax2, (leny2,lenx2))

Relem2m = np.reshape(Relem2, (leny2,lenx2))

print('shape of f_B_Mat', np.shape(f_B2m))
print('Number of beam positions over which has been averaged: ', lenx2, leny2, lenx2*leny2,  len(y_pos2), len(f_B2), f_B2[0],f_B[2-1])

# for Block-Case ##################################################################################################################

print('Number of beam positions over which has been averaged: ', lenx3, leny3, lenx3*leny3,  len(y_pos3), len(f_B3), f_B3[0],f_B3[-1])


f_B3m = np.reshape(f_B3, (leny3, lenx3))
rho_x3m=np.reshape(rho_x3, (leny3,lenx3))
v_r3m =np.reshape(v_r3, (leny3,lenx3))
tau_B3m =np.reshape(tau_B3, (leny3,lenx3))
v_rmist3m =np.reshape(v_rmist3, (leny3,lenx3))
vdmax3m =np.reshape(vdmax3, (leny3,lenx3))
vimax3m =np.reshape(vimax3, (leny3,lenx3))

Relem3m = np.reshape(Relem3, (leny3,lenx3))

print('shape of f_B_Mat', np.shape(f_B3m))


###################################################################################################################################
########### Mean and error calculations ###########################################################################################
###################################################################################################################################

# for Em-Case #####################################################################################################################

f_B = [None]*lenx
f_Berr = [None]*lenx
for f in range (lenx):
	f_B[f] = stats.nanmean(f_Bm[:,f])
	f_Berr[f] = stats.nanstd(f_Bm[:,f])

rho_x = [None]*lenx
rho_xerr = [None]*lenx
for f in range (lenx):
	rho_x[f] = stats.nanmean(rho_xm[:,f])
	rho_xerr[f] = stats.nanstd(rho_xm[:,f])
		
tau_B = [None]*lenx
tau_Berr = [None]*lenx
for f in range (lenx):
	tau_B[f] = stats.nanmean(tau_Bm[:,f])
	tau_Berr[f] = stats.nanstd(tau_Bm[:,f])

v_r = [None]*lenx
v_rerr = [None]*lenx
for f in range (lenx):
	v_r[f] = stats.nanmean(v_rm[:,f])
	v_rerr[f] = stats.nanstd(v_rm[:,f])

vdmax = [None]*lenx
vdmaxerr = [None]*lenx
for f in range (lenx):
	vdmax[f] = stats.nanmean(vdmaxm[:,f])
	vdmaxerr[f] = stats.nanstd(vdmaxm[:,f])

vimax = [None]*lenx
vimaxerr = [None]*lenx
for f in range (lenx):
	vimax[f] = stats.nanmean(vimaxm[:,f])
	vimaxerr[f] = stats.nanstd(vimaxm[:,f])

Relem = [None]*lenx
Relemerr = [None]*lenx
for f in range (lenx):
	Relem[f] = stats.nanmean(Relemm[:,f])
	Relemerr[f] = stats.nanstd(Relemm[:,f])


# for Den-Case #####################################################################################################################
f_B2 = [None]*lenx2
f_B2err = [None]*lenx2
for f in range (lenx2):
	f_B2[f] = stats.nanmean(f_B2m[:,f])
	f_B2err[f] = stats.nanstd(f_B2m[:,f])
	
rho_x2 = [None]*lenx2
rho_x2err = [None]*lenx2
for f in range (lenx2):
	rho_x2[f] = stats.nanmean(rho_x2m[:,f])
	rho_x2err[f] = stats.nanstd(rho_x2m[:,f])

tau_B2 = [None]*lenx2
tau_B2err = [None]*lenx2
for f in range (lenx2):
	tau_B2[f] = stats.nanmean(tau_B2m[:,f])
	tau_B2err[f] = stats.nanstd(tau_B2m[:,f])

v_r2 = [None]*lenx2
v_r2err = [None]*lenx2
for f in range (lenx2):
	v_r2[f] = stats.nanmean(v_r2m[:,f])
	v_r2err[f] = stats.nanstd(v_r2m[:,f])

vdmax2 = [None]*lenx2
vdmax2err = [None]*lenx2
for f in range (lenx2):
	vdmax2[f] = stats.nanmean(vdmax2m[:,f])
	vdmax2err[f] = stats.nanstd(vdmax2m[:,f])

vimax2 = [None]*lenx2
vimax2err = [None]*lenx2
for f in range (lenx2):
	vimax2[f] = stats.nanmean(vimax2m[:,f])
	vimax2err[f] = stats.nanstd(vimax2m[:,f])

Relem2 = [None]*lenx2
Relem2err = [None]*lenx2
for f in range (lenx2):
	Relem2[f] = stats.nanmean(Relem2m[:,f])
	Relem2err[f] = stats.nanstd(Relem2m[:,f])



# for Block-Case #####################################################################################################################
f_B3 = [None]*lenx3
f_B3err = [None]*lenx3
for f in range (lenx3):
	f_B3[f] = stats.nanmean(f_B3m[:,f])
	f_B3err[f] = stats.nanstd(f_B3m[:,f])
	
rho_x3 = [None]*lenx3
rho_x3err = [None]*lenx3
for f in range (lenx3):
	rho_x3[f] = stats.nanmean(rho_x3m[:,f])
	rho_x3err[f] = stats.nanstd(rho_x3m[:,f])

tau_B3 = [None]*lenx3
tau_B3err = [None]*lenx3
for f in range (lenx3):
	tau_B3[f] = stats.nanmean(tau_B3m[:,f])
	tau_B3err[f] = stats.nanstd(tau_B3m[:,f])

v_r3 = [None]*lenx3
v_r3err = [None]*lenx3
for f in range (lenx3):
	v_r3[f] = stats.nanmean(v_r3m[:,f])
	v_r3err[f] = stats.nanstd(v_r3m[:,f])

vdmax3 = [None]*lenx3
vdmax3err = [None]*lenx3
for f in range (lenx3):
	vdmax3[f] = stats.nanmean(vdmax3m[:,f])
	vdmax3err[f] = stats.nanstd(vdmax3m[:,f])

vimax3 = [None]*lenx3
vimax3err = [None]*lenx3
for f in range (lenx3):
	vimax3[f] = stats.nanmean(vimax3m[:,f])
	vimax3err[f] = stats.nanstd(vimax3m[:,f])


Relem3 = [None]*lenx3
Relem3err = [None]*lenx3
for f in range (lenx3):
	Relem3[f] = stats.nanmean(Relem3m[:,f])
	Relem3err[f] = stats.nanstd(Relem3m[:,f])


# Masking nan-values for plots, if need ############################################################################################
masks = False
if masks:
# 	convert to array
	f_B3 = np.array(f_B3)
	f_Berr3 = np.array(f_Berr3)
	tau_B3 = np.array(tau_B3)
	tau_B3err = np.array(tau_B3err)
	rho_x3 = np.array(rho_x3)
	rho_x3err = np.array(rho_x3err)
	v_r3 = np.array(v_r3)
	v_r3err = np.array(v_r3err)
	vdmax3 = np.array(vdmax3)
	vdmax3err = np.array(vdmax3err)
	vimax3 = np.array(vimax3)
	vimax3err = np.array(vimax3err)

# 	built masks
	f_B3mask = np.isfinite(f_B3)
	f_B3errmask = np.isfinite(f_Berr3)
	tau_B3mask = np.isfinite(tau_B3)
	tau_B3errmask = np.isfinite(tau_B3err)
	rho_x3mask = np.isfinite(rho_x3)
	rho_x3errmask = np.isfinite(rho_x3err)
	v_r3mask = np.isfinite(v_r3)
	v_r3errmask = np.isfinite(v_r3err)
	vdmax3mask = np.isfinite(vdmax3)
	vdmax3errmask = np.isfinite(vdmax3err)
	vimax3mask = np.isfinite(vimax3)
	vimax3errmask = np.isfinite(vimax3err)

# plot(x[f_Bmask],y[f_Bmask])


###################################################################################################################################
########### FIGURE ################################################################################################################
###################################################################################################################################


# create fiugre with certain ratio and facecolor:
f1=plt.figure(figsize=(16,16), facecolor = 'white')
gs = gridspec.GridSpec(4, 2,				# ratio of grid space (2 plots per collumn, 3 per row)
                       width_ratios=[1,1],		# width ratios of the 3 plots per row
                       height_ratios=[1,1,1,1]		# height ratios of the 2 polots per collumn
                       )


# create axes with the ratios specified in gs above. ax5, ax6 are the axes for the colorbars:
ax1 = plt.subplot(gs[0])	# blob width
ax4 = plt.subplot(gs[1])	# vr
ax5 = plt.subplot(gs[3])	# vdmax
ax2 = plt.subplot(gs[2])	# tauB
ax7 = plt.subplot(gs[5])	# vimax
ax3 = plt.subplot(gs[4])	# fB
ax6 = plt.subplot(gs[6])	# fB freq
ax8 = plt.subplot(gs[7])	# rleative fluctuation


# create a list in for simple modifications on all plots at the same time:
allax=[ax1,ax2,ax3,ax4,ax5,ax6]
font = {'family' :'Arial', 'weight' : 'normal', 'size' : 14}
matplotlib.rc('font', **font)
# Plots ############################################################################################################

# first subplot: density

ax1.errorbar(x,rho_x,yerr = rho_xerr, marker='s', color = 'k')
ax1.errorbar(x2,rho_x2,yerr = rho_x2err,marker='o', color = 'b')
ax1.errorbar(x3,rho_x3,yerr = rho_x3err,marker='D', color = 'g')
ax1.axvline(0,color='k', linestyle='-.')
#ax1.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
ax1.set_ylabel(r'blob width $\rho_x$ (cm)')
ax1.get_yaxis().set_label_coords(-0.09,0.5)
ax1.set_title(r'blob width $\rho_x$ (cm)')
ax1.plot((0,4),(2.3,2.3), 'k',) # arrow line
ax1.plot((0,0),(2.3,2.3), 'k', marker='<',) # left arrowhead
ax1.plot((4,4),(2.3,2.3), 'k', marker='>',) # right arrowhead
ax1.text(2,2.1, 'SOL', font)
ax1.text(-0.5,2.1, 'LCFS', font)

tau_Bfig = ax2.errorbar(x,tau_B, yerr = tau_Berr, marker='s', color = 'k')
tau_Bfig2 = ax2.errorbar(x2,tau_B2, yerr = tau_B2err, marker='o', color = 'b')
tau_Bfig3 = ax2.errorbar(x3,tau_B3, yerr = tau_B3err,marker='D', color = 'g')
ax2.axvline(0,color='k', linestyle='-.')
#ax2.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
ax2.set_ylabel(r'Self-Correlation time $\tau_B$ (ms)')
ax2.get_yaxis().set_label_coords(-0.09,0.5)
ax2.set_title(r'Self-Correlation time $\tau_B$')

f_Bfig = ax3.errorbar(x,f_B, yerr = f_Berr, marker='s', color = 'k')
f_Bfig2 = ax3.errorbar(x2,f_B2, yerr = f_B2err, marker='o', color = 'b')
f_Bfig3 = ax3.errorbar(x3,f_B3, yerr = f_B3err,marker='D', color = 'g')
ax3.axvline(0,color='k', linestyle='-.')
#ax3.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
ax3.set_ylabel('Number of Blobs')
ax3.get_yaxis().set_label_coords(-0.09,0.5)
ax3.set_title('Number of Blobs in observed time intervall {0:.2f}ms'.format(time[10]))

v_rfig = ax4.errorbar(x,v_r, yerr = v_rerr,marker='s', color = 'k')
v_rfig2 = ax4.errorbar(x2,v_r2, yerr = v_r2err,marker='o', color = 'b')
v_rfig3 = ax4.errorbar(x3,v_r3, yerr = v_r3err,marker='s', color = 'g')
ax4.axvline(0,color='k', linestyle='-.')
#ax4.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
ax4.set_ylabel(r'Average velocity $v_r$ (m/s)')
ax4.get_yaxis().set_label_coords(-0.11,0.5)
ax4.set_title(r'Average velocity $v_r$')

v_dmaxfig = ax5.errorbar(x,vdmax,yerr = vdmaxerr, color = 'k',marker='s')
v_dmax2fig = ax5.errorbar(x2,vdmax2,yerr = vdmax2err,color = 'b',marker='o')
v_dmax3fig = ax5.errorbar(x3,vdmax3,yerr = vdmax3err,color =  'g',marker='D')
ax5.axvline(0,color='k', linestyle='-.')
#ax5.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
ax5.set_ylabel(r'Maximum velocity $v_{dmax}$ (m/s)')
ax5.get_yaxis().set_label_coords(-0.11,0.5)
ax5.set_title(r'Maximum velocity for real data')


f_B2fig = ax6.errorbar(x,f_B/time[10]*1000,yerr = f_Berr/time[10]*1000, marker='s', color = 'k')
f_B22fig = ax6.errorbar(x2,f_B2/time[10]*1000,yerr = f_B2err/time[10]*1000, marker='o', color = 'b')
f_B23fig = ax6.errorbar(x3,f_B3/time[10]*1000,yerr = f_B3err/time[10]*1000, marker='D', color = 'g')
ax6.axvline(0,color='k', linestyle='-.')
ax6.set_xlabel('beam axis $x$ (cm)', labelpad= 0)			# switched of, since axis in row below
ax6.set_ylabel('Blob frequency $f_B$ (1/s)')
ax6.get_yaxis().set_label_coords(-0.09,0.5)
ax6.set_title('Blob frequency $f_B$')

ax7.errorbar(x,vimax, yerr = vimaxerr,color =  'k',marker='s',label='Emission data')
ax7.errorbar(x2,vimax2, yerr = vimax2err, color = 'b',marker='o',label='Density data')
ax7.errorbar(x3,vimax3, yerr = vimax3err, color = 'g',marker='D',label='Block data')
ax7.axvline(0,color='k', linestyle='-.')
#ax7.set_xlabel('beam axis $x$ (cm)', labelpad= 0)			# switched of, since axis in row below
ax7.set_ylabel(r'Maximum velocities $v_{imax}$ (m/s)')
ax7.get_yaxis().set_label_coords(-0.11,0.5)
ax7.set_title(r'Maximum velocities for interpolated data')


ax8.errorbar(x,Relem, yerr = Relemerr,color =  'k',marker='s',label='Emission data')
ax8.errorbar(x2,Relem2, yerr = Relem2err, color = 'b',marker='o',label='Density data')
ax8.errorbar(x3,Relem3, yerr = Relem3err, color = 'g',marker='D',label='Block data')
ax8.axvline(0,color='k', linestyle='-.')
ax8.set_xlabel('beam axis $x$ (cm)', labelpad= 0)			# switched of, since axis in row below
ax8.set_ylabel(r'relative emission $\delta I/I$ or $\delta n/n$')
ax8.get_yaxis().set_label_coords(-0.11,0.5)
ax8.set_title(r'Blob amplitude')

leg = ax8.legend(loc='upper center', bbox_to_anchor = (-0.1,-0.05),fancybox = True, numpoints = 1)

# leg = ax7.legend(loc='upper center', bbox_to_anchor = (0.225,-0.18),fancybox = True, numpoints = 1)


f1.savefig('{0:03d}FigRadial{1:}'.format(Measurement,filename))

plt.show()






