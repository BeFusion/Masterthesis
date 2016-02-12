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

filename = '_veltwind_SNR'				# specify fileextension to be read in and to be stored (e.g '_dens_test', if file is called '004Emresults_dens_test')

Measurement = 4				# specify number of measurement without leading zeros, which is the prefix of the file in form (will be extended later to e.g. 004) 


EmCase		= False
BlockCase	= True

if EmCase:
	files=['{0:03d}Emresults{1:}010_Smooth.txt'.format(Measurement,filename),'{0:03d}Emresults{1:}020_Smooth.txt'.format(Measurement,filename),'{0:03d}Emresults{1:}030_Smooth.txt'.format(Measurement,filename),'{0:03d}Emresults{1:}035_Smooth.txt'.format(Measurement,filename),'{0:03d}Emresults{1:}050_Smooth.txt'.format(Measurement,filename), '{0:03d}Emresults{1:}100_Smooth.txt'.format(Measurement,filename), '{0:03d}Emresults_veltwind.txt'.format(Measurement)]
if BlockCase:
	files=['{0:03d}Blockresults{1:}010_Smooth.txt'.format(Measurement,filename),'{0:03d}Blockresults{1:}020_Smooth.txt'.format(Measurement,filename),'{0:03d}Blockresults{1:}030_Smooth.txt'.format(Measurement,filename),'{0:03d}Blockresults{1:}035_Smooth.txt'.format(Measurement,filename),'{0:03d}Blockresults{1:}050_Smooth.txt'.format(Measurement,filename), '{0:03d}Blockresults{1:}100_Smooth.txt'.format(Measurement,filename), '{0:03d}Blockresults_veltwind.txt'.format(Measurement)]

SNR = ['0.10','0.20','0.30','0.35','0.50','1.00','$\infty$']

numfil = len(files) 	# number of files to analyze

time=[None]*numfil
x_pos=[None]*numfil
y_pos=[None]*numfil
f_B=[None]*numfil
rho_x=[None]*numfil
tau_B=[None]*numfil
v_r=[None]*numfil
v_rmist=[None]*numfil
vdmax=[None]*numfil
vimax=[None]*numfil
Relem=[None]*numfil
B0=[None]*numfil
Lp=[None]*numfil
q95=[None]*numfil

#read in all emission files for all SNR-values
for i in range(numfil):
	time[i], x_pos[i], y_pos[i], f_B[i], rho_x[i], tau_B[i], v_r[i], v_rmist[i], vdmax[i], vimax[i], Relem[i], B0[i], Lp1,q95[i] = np.loadtxt(files[i], usecols = (3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17), unpack = True, skiprows=2)

#read in one density profile (should be the same for all SNR-results)
time2, x_pos2,y_pos2, f_B2, rho_x2, tau_B2, v_r2, v_rmist2, vdmax2, vimax2, Relem2, B02, Lp2,q952 = np.loadtxt('{0:03d}Denresults{1:}010_Smooth.txt'.format(Measurement,filename), usecols = (3, 4, 5,  6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17), unpack = True, skiprows=2)


###################################################################################################################################
########### Reshaping axes ########################################################################################################
###################################################################################################################################

# for Em-Cases #####################################################################################################################

y_min = min(y_pos[0])
lencount = 0
y_pos=np.array(y_pos[0])
x_pos=np.array(x_pos[0])
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

###################################################################################################################################
########### Reshaping axes ########################################################################################################
###################################################################################################################################



# for Den-Case ####################################################################################################################

f_B2m = np.reshape(f_B2, (leny2, lenx2))
rho_x2m =np.reshape(rho_x2, (leny2,lenx2))
v_r2m =np.reshape(v_r2, (leny2,lenx2))
tau_B2m =np.reshape(tau_B2, (leny2,lenx2))
v_rmist2m =np.reshape(v_rmist2, (leny2,lenx2))
vdmax2m =np.reshape(vdmax2, (leny2,lenx2))
vimax2m =np.reshape(vimax2, (leny2,lenx2))
Relem2m = np.reshape(Relem2, (leny2,lenx2))

###################################################################################################################################
########### Mean and error calculations ###########################################################################################
###################################################################################################################################

# for Em-Case #####################################################################################################################


f_Bcop = np.array(f_B).copy()
rho_xcop = np.array(rho_x).copy()
v_rcop = np.array(v_r).copy()
tau_Bcop = np.array(tau_B).copy()
v_rmistcop = np.array(v_rmist).copy()
vdmaxcop = np.array(vdmax).copy()
vimaxcop = np.array(vimax).copy()
Relemcop = np.array(Relem).copy()

for i in range(numfil):
	f_Bm = np.reshape(f_Bcop[i], (leny, lenx))
	rho_xm =np.reshape(rho_xcop[i], (leny,lenx))
	v_rm =np.reshape(v_rcop[i], (leny,lenx))
	tau_Bm =np.reshape(tau_Bcop[i], (leny,lenx))
	v_rmistm =np.reshape(v_rmistcop[i], (leny,lenx))
	vdmaxm =np.reshape(vdmaxcop[i], (leny,lenx))
	vimaxm =np.reshape(vimaxcop[i], (leny,lenx))
	Relemm = np.reshape(Relemcop[i], (leny,lenx))

	if i==0:
		f_B = [None]*lenx*numfil
		f_B = np.array(f_B)
		f_B = np.reshape(f_B,(numfil,lenx))
		f_Berr = [None]*lenx*numfil
		f_Berr = np.array(f_Berr)
		f_Berr = np.reshape(f_Berr,(numfil,lenx))

		rho_x = [None]*lenx*numfil
		rho_x = np.array(rho_x)
		rho_x = np.reshape(rho_x,(numfil,lenx))
		rho_xerr = [None]*lenx*numfil
		rho_xerr = np.array(rho_xerr)
		rho_xerr = np.reshape(rho_xerr,(numfil,lenx))

		tau_B = [None]*lenx*numfil
		tau_B = np.array(tau_B)
		tau_B = np.reshape(tau_B,(numfil,lenx))
		tau_Berr = [None]*lenx*numfil
		tau_Berr = np.array(tau_Berr)
		tau_Berr = np.reshape(tau_Berr,(numfil,lenx))

		v_r = [None]*lenx*numfil
		v_r = np.array(v_r)
		v_r = np.reshape(v_r,(numfil,lenx))
		v_rerr = [None]*lenx*numfil
		v_rerr = np.array(v_rerr)
		v_rerr = np.reshape(v_rerr,(numfil,lenx))

		vdmax = [None]*lenx*numfil
		vdmax = np.array(vdmax)
		vdmax = np.reshape(vdmax,(numfil,lenx))
		vdmaxerr = [None]*lenx*numfil
		vdmaxerr = np.array(vdmaxerr)
		vdmaxerr = np.reshape(vdmaxerr,(numfil,lenx))

		vimax = [None]*lenx*numfil
		vimax = np.array(vimax)
		vimax = np.reshape(vimax,(numfil,lenx))
		vimaxerr = [None]*lenx*numfil
		vimaxerr = np.array(vimaxerr)
		vimaxerr = np.reshape(vimaxerr,(numfil,lenx))

		Relem = [None]*lenx*numfil
		Relem = np.array(Relem)
		Relem = np.reshape(Relem,(numfil,lenx))
		Relemerr = [None]*lenx*numfil
		Relemerr = np.array(Relemerr)
		Relemerr = np.reshape(Relemerr,(numfil,lenx))

	for f in range (lenx):
		f_B[i,f] = stats.nanmean(f_Bm[:,f])
		f_Berr[i,f] = stats.nanstd(f_Bm[:,f])

		rho_x[i,f] = stats.nanmean(rho_xm[:,f])
		rho_xerr[i,f] = stats.nanstd(rho_xm[:,f])

		tau_B[i,f] = stats.nanmean(tau_Bm[:,f])
		tau_Berr[i,f] = stats.nanstd(tau_Bm[:,f])

		v_r[i,f] = stats.nanmean(v_rm[:,f])
		v_rerr[i,f] = stats.nanstd(v_rm[:,f])

		vdmax[i,f] = stats.nanmean(vdmaxm[:,f])
		vdmaxerr[i,f] = stats.nanstd(vdmaxm[:,f])

		vimax[i,f] = stats.nanmean(vimaxm[:,f])
		vimaxerr[i,f] = stats.nanstd(vimaxm[:,f])

		Relem[i,f] = stats.nanmean(Relemm[:,f])
		Relemerr[i,f] = stats.nanstd(Relemm[:,f])
	
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
markers=['v', 'H', 'D', 'd', '^', '<', '>', '8', 's', 'p', '*', 'h']
colors = ['g','r','c','m','orange','coral','pink','grey','darkred','peachpuff','goldenrod','plum','lavender']
font = {'family' :'Arial', 'weight' : 'normal', 'size' : 14}
matplotlib.rc('font', **font)
# Plots ############################################################################################################

# first subplot: density

ax1.errorbar(x2,rho_x2,yerr = rho_x2err, marker='o', color = 'b')
for i in range(numfil):
	ax1.errorbar(x,rho_x[i],yerr = rho_xerr[i],marker=markers[i], color = colors[i])
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

tau_Bfig2 = ax2.errorbar(x2,tau_B2, yerr = tau_B2err, marker='o', color = 'b')
for i in range(numfil):
	ax2.errorbar(x,tau_B[i], yerr = tau_Berr[i],marker=markers[i], color = colors[i])
ax2.axvline(0,color='k', linestyle='-.')
#ax2.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
ax2.set_ylabel(r'Self-Correlation time $\tau_B$ (ms)')
ax2.get_yaxis().set_label_coords(-0.09,0.5)
ax2.set_title(r'Self-Correlation time $\tau_B$')

f_Bfig2 = ax3.errorbar(x2,f_B2, yerr = f_B2err, marker='o', color = 'b')
for i in range(numfil):
	ax3.errorbar(x,f_B[i], yerr = f_Berr[i],marker=markers[i], color = colors[i])
ax3.axvline(0,color='k', linestyle='-.')
#ax3.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
ax3.set_ylabel('Number of Blobs')
ax3.get_yaxis().set_label_coords(-0.09,0.5)
ax3.set_title('Number of Blobs in observed time intervall {0:.2f}ms'.format(min(time[0])))


v_rfig2 = ax4.errorbar(x2,v_r2, yerr = v_r2err,marker='o', color = 'b')
for i in range(numfil):
	ax4.errorbar(x,v_r[i], yerr = v_rerr[i],marker=markers[i], color = colors[i])
ax4.axvline(0,color='k', linestyle='-.')
#ax4.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
ax4.set_ylabel(r'Average velocity $v_r$ (m/s)')
ax4.get_yaxis().set_label_coords(-0.11,0.5)
ax4.set_title(r'Average velocity $v_r$')


v_dmax2fig = ax5.errorbar(x2,vdmax2,yerr = vdmax2err,color = 'b',marker='o')
for i in range(numfil):
	ax5.errorbar(x,vdmax[i], yerr = vdmaxerr[i],marker=markers[i], color = colors[i])
ax5.axvline(0,color='k', linestyle='-.')
#ax5.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
ax5.set_ylabel(r'Maximum velocity $v_{dmax}$ (m/s)')
ax5.get_yaxis().set_label_coords(-0.11,0.5)
ax5.set_ylim(0,6000)
ax5.set_title(r'Maximum velocity for real data')


f_B22fig = ax6.errorbar(x2,np.array(f_B2)/min(time[0])*1000,yerr = np.array(f_B2err)/min(time[0])*1000, marker='o', color = 'b')
for i in range(numfil):
	ax6.errorbar(x,1/min(time[0])*1000*f_B[i], yerr = 1/min(time[0])*1000*f_Berr[i],marker=markers[i], color = colors[i])
ax6.axvline(0,color='k', linestyle='-.')
ax6.set_xlabel('beam axis $x$ (cm)', labelpad= 0)			# switched of, since axis in row below
ax6.set_ylabel('Blob frequency $f_B$ (1/s)')
ax6.get_yaxis().set_label_coords(-0.09,0.5)
ax6.set_title('Blob frequency $f_B$')

ax7.errorbar(x2,vimax2, yerr = vimax2err, color = 'b',marker='o',label='Density data')
for i in range(numfil):
	ax7.errorbar(x,vimax[i], yerr = vimaxerr[i],marker=markers[i], color = colors[i])
ax7.axvline(0,color='k', linestyle='-.')
#ax7.set_xlabel('beam axis $x$ (cm)', labelpad= 0)			# switched of, since axis in row below
ax7.set_ylabel(r'Maximum velocities $v_{imax}$ (m/s)')
ax7.get_yaxis().set_label_coords(-0.11,0.5)
ax7.set_ylim(0,8000)
ax7.set_title(r'Maximum velocities for interpolated data')



ax8.errorbar(x2,Relem2, yerr = Relem2err, color = 'b',marker='o',label='Density data')
for i in range(numfil):
	ax8.errorbar(x,Relem[i], yerr = Relemerr[i],marker=markers[i], color = colors[i], label='Emission for SNR = {0:}'.format(SNR[i]))
ax8.axvline(0,color='k', linestyle='-.')
ax8.set_xlabel('beam axis $x$ (cm)', labelpad= 0)			# switched of, since axis in row below
ax8.set_ylabel(r'relative emission $\delta I/I$ or $\delta n/n$')
#ax8.set_ylim(0,7)
ax8.get_yaxis().set_label_coords(-0.11,0.5)
ax8.set_title(r'Blob amplitude')

leg = ax8.legend(loc='upper center', bbox_to_anchor = (-0.16,-0.18), numpoints = 1,ncol=4)

# leg = ax7.legend(loc='upper center', bbox_to_anchor = (0.225,-0.18),fancybox = True, numpoints = 1)


if EmCase:
	f1.savefig('{0:03d}FigRadial_Em_data{1:}_Analysis'.format(Measurement,filename))
if BlockCase:
	f1.savefig('{0:03d}FigRadial_Block_data{1:}_Analysis'.format(Measurement,filename))
plt.show()






