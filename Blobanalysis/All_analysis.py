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
from Radial_analysis import radial_analysis
from Blob_params import *			# set a few parameters 
import timeit					# for timekeeping
from datetime import datetime			# for date and time printing in output file
import sys					# for exiting the program if errors occur
import pdb					# for error tracking (set_trace at right position)


# number of files to be analyzed:
#filename = '_veltwind_TEST_TIME_SNR200_Smooth'

NoiseSmooth = False
fourCase = True

if NoiseSmooth:
	#SNR = 1.00	# set for saving file correctly
	filename = '_veltwind_TIME_SNR100_Smooth'
if not NoiseSmooth:
	filename = '_veltwind'


# Give measurement numbers here:
Refdec = 1.7
# 111 is missing so far for Smooth-Case...

files = [ 102, 105, 108,109,110,111,112, 113,114,116]
le = len(files)

savenow = True 		# saves all images from the called function Radial_analysis.py
savenow2 = True		# saves the images created here

#define vectors, which are needed in the following
B0 = [None]*le
q95 = [None]*le
Lp = [None]*le

fB2 = [None]*le
fB2err = [None]*le
rho_x = [None]*le
rho_xerr = [None]*le
tau_B = [None]*le
tau_Berr = [None]*le
v_r = [None]*le
v_rerr = [None]*le
Relem = [None]*le
Relemerr = [None]*le

fB22 = [None]*le
fB22err = [None]*le
rho_x2 = [None]*le
rho_x2err = [None]*le
tau_B2 = [None]*le
tau_B2err = [None]*le
v_r2 = [None]*le
v_r2err = [None]*le
Relem2 = [None]*le
Relem2err = [None]*le

fB23 = [None]*le
fB23err = [None]*le
rho_x3 = [None]*le
rho_x3err = [None]*le
tau_B3 = [None]*le
tau_B3err = [None]*le
v_r3 = [None]*le
v_r3err = [None]*le
Relem3 = [None]*le
Relem3err = [None]*le

fB24 = [None]*le
fB24err = [None]*le
rho_x4 = [None]*le
rho_x4err = [None]*le
tau_B4 = [None]*le
tau_B4err = [None]*le
v_r4 = [None]*le
v_r4err = [None]*le
Relem4 = [None]*le
Relem4err = [None]*le



for i in range (0,len(files)):
	
	# BlockCase we need just the x-axis to find the right index
	x_pos3, y_pos3 = np.loadtxt('{0:03d}Blockresults{1:}.txt'.format(files[i],filename), usecols = (4,5), unpack = True, skiprows=2)

	y_min3 = min(y_pos3)
	lencount3 = 0
	for elem in range(len(y_pos3)):
		if y_pos3[elem] == y_min3:
			lencount3 = lencount3 +1

	leny3 = len(y_pos3)/lencount3 			# number of y-steps
	lenx3 = lencount3 				# number of x-steps
	x3 = [None]*lenx3
	for elem in range (lenx3):
		x3[elem]=x_pos3[elem]


	# find Reference detector index along beam 
	countx3=0
	for g in range (0,len(x3)):
		if (x3[g]<0):				# take negative values into account for correct shift
			countx3=countx3+1			# number of elements, which are negative
	Refdec_ind3 = len(x3)-int(Refdec/0.5)-1 - countx3	# index for Reference detector position

		
	# EmCase we need just the x-axis to find the right index
	x_pos, y_pos = np.loadtxt('{0:03d}Emresults{1:}.txt'.format(files[i],filename), usecols = (4,5), unpack = True, skiprows=2)

	y_min = min(y_pos)
	lencount = 0
	for elem in range(len(y_pos)):
		if y_pos[elem] == y_min:
			lencount = lencount +1

	leny = len(y_pos)/lencount 			# number of y-steps
	lenx = lencount 				# number of x-steps
	x = [None]*lenx
	for elem in range (lenx):
		x[elem]=x_pos[elem]

	countx=0
	for g in range (0,len(x)):		#500 random lenght bigger than len(x)
		if (x[g]<0):				# take negative values into account for correct shift
			countx=countx+1			# number of elements, which are negative
	Refdec_ind = int(len(x)-x3[Refdec_ind3]/0.02/SetResolution-1 - countx)	# index for Reference detector position

	# call function radial_analysis to get values at certain position
	B0[i], q95[i], Lp[i], fB2[i], fB2err[i], rho_x[i], rho_xerr[i], tau_B[i], tau_Berr[i], v_r[i], v_rerr[i], Relem[i], Relemerr[i], fB22[i], fB22err[i], rho_x2[i], rho_x2err[i], tau_B2[i], tau_B2err[i], v_r2[i], v_r2err[i], Relem2[i], Relem2err[i], fB23[i], fB23err[i], rho_x3[i], rho_x3err[i], tau_B3[i], tau_B3err[i], v_r3[i], v_r3err[i], Relem3[i], Relem3err[i],fB24[i], fB24err[i], rho_x4[i], rho_x4err[i], tau_B4[i], tau_B4err[i], v_r4[i], v_r4err[i], Relem4[i], Relem4err[i] = radial_analysis(filename,files[i],Refdec_ind, Refdec_ind3, savenow, SNR, NoiseSmooth, fourCase)
	plt.close()				#we do not want to see all the figures created in Radial_analysis

print('Refdec-Position in EmCase:', x[Refdec_ind], 'and was set to',  Refdec)
print('Refdec-Position in BlockCase:', x3[Refdec_ind3], 'and was set to',  Refdec)		
	

# create fiugre for frequency: figure 1
f1=plt.figure(figsize=(16,16), facecolor = 'white')
gs = gridspec.GridSpec(2, 2,				# ratio of grid space (2 plots per collumn, 3 per row)
                       width_ratios=[1,1],		# width ratios of the 3 plots per row
                       height_ratios=[1,1]		# height ratios of the 2 polots per collumn
                       )


font = {'family' :'Arial', 'weight' : 'normal', 'size' : 24}
matplotlib.rc('font', **font)

# create axes with the ratios specified in gs above. ax5, ax6 are the axes for the colorbars:
ax1 = plt.subplot(gs[0])	# vs B
ax2 = plt.subplot(gs[1])	# vs max fluctuation
ax3 = plt.subplot(gs[2])	# vs q95
ax4 = plt.subplot(gs[3])	# vs Lp


msize = 12
ewid = 4
csiz = 6
mwid = 0
# Plots ############################################################################################################

xshift = 0.05
yshift = 0.92
# first subplot: density


ax1.errorbar(B0,fB2/np.float64(1000), yerr = fB2err/np.float64(1000), marker='s', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'k', ls = 'None', label = 'Emission Case')
ax1.errorbar(B0,fB22/np.float64(1000), yerr = fB2err/np.float64(1000), marker='o', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'b', ls = 'None', label = 'Density Case')
ax1.errorbar(B0,fB23/np.float64(1000), yerr = fB23err/np.float64(1000), marker='D', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'g', ls = 'None', label = 'Block Case')
ax1.errorbar(B0,fB24/np.float64(1000), yerr = fB24err/np.float64(1000), marker='>', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'r', ls = 'None', label = 'Real SNR Block Case')
ax1.set_xlabel(r'$B_0$ (T)')			
ax1.set_ylabel(r'$f_B$ (1000/s)')
ax1.get_yaxis().set_label_coords(-0.12,0.5)


ax2.errorbar(Relem,fB2/np.float64(1000), xerr = Relemerr, yerr = fB23err/np.float64(1000), marker='s', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'k', ls = 'None')
ax2.errorbar(Relem2,fB22/np.float64(1000), xerr = Relem2err, yerr = fB23err/np.float64(1000), marker='o', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'b', ls = 'None')
ax2.errorbar(Relem3,fB23/np.float64(1000), xerr = Relem3err, yerr = fB23err/np.float64(1000), marker='D', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'g', ls = 'None')
ax2.errorbar(Relem4,fB24/np.float64(1000), xerr = Relem4err, yerr = fB24err/np.float64(1000), marker='>', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'r', ls = 'None')
ax2.set_xlabel(r'$\delta I/I$ or $\delta n/n$ ')			
ax2.set_ylabel(r'$f_B$ (1000/s)')
ax2.set_xlim(0.21,1.6)
ax2.get_yaxis().set_label_coords(-0.12,0.5)

ax3.errorbar(q95,fB2/np.float64(1000), yerr = fB2err/np.float64(1000), marker='s', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'k', ls = 'None')
ax3.errorbar(q95,fB22/np.float64(1000), yerr = fB22err/np.float64(1000), marker='o', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'b', ls = 'None')
ax3.errorbar(q95,fB23/np.float64(1000), yerr = fB23err/np.float64(1000), marker='D', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'g', ls = 'None')
ax3.errorbar(q95,fB24/np.float64(1000), yerr = fB24err/np.float64(1000), marker='>', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'r', ls = 'None')
ax3.set_xlabel(r'$q_{95}$')
ax3.set_ylabel(r'$f_B$ (1000/s)')
ax3.get_yaxis().set_label_coords(-0.12,0.5)
ax3.set_xlim(3.6,6.5)

ax4.errorbar(Lp,fB22/np.float64(1000), yerr = fB22err/np.float64(1000),marker='o', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'b', ls = '', label = 'Density Case')
ax4.errorbar(Lp,fB2/np.float64(1000), yerr = fB2err/np.float64(1000),marker='s', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'k', ls = '', label = 'Emission Case')
ax4.errorbar(Lp,fB23/np.float64(1000), yerr = fB23err/np.float64(1000),marker='D', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'g', ls = '', label = 'Block Case')
ax4.errorbar(Lp,fB24/np.float64(1000), yerr = fB24err/np.float64(1000),marker='>', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'r', ls = '', label = 'Real SNR Block Case')
ax4.set_xlabel(r'$L_{||}$')
ax4.set_ylabel(r'$f_B$ (1000/s)')
ax4.get_yaxis().set_label_coords(-0.12,0.5)
ax4.set_xlim(8,13.5)

ax1.annotate("a)", xycoords="axes fraction", xy = (xshift,yshift), va = 'center', ha = 'center', fontsize = 25)
ax2.annotate("b)", xycoords="axes fraction", xy = (xshift,yshift), va = 'center', ha = 'center', fontsize = 25)
ax3.annotate("c)", xycoords="axes fraction", xy = (xshift,yshift), va = 'center', ha = 'center', fontsize = 25)
ax4.annotate("d)", xycoords="axes fraction", xy = (xshift,yshift), va = 'center', ha = 'center', fontsize = 25)


plt.tight_layout()
#leg = ax4.legend(loc='upper center', bbox_to_anchor = (-0.1,-0.05),fancybox = True, numpoints = 1)
handles, labels = ax4.get_legend_handles_labels()
handles = [h[0] for h in handles]
leg1 = ax4.legend(handles, labels,loc='center left', bbox_to_anchor = (-0.9,-0.3), fancybox = True, numpoints = 1, ncol = 2, prop={'size':25})

#leg = ax4.legend(loc='upper center', bbox_to_anchor = (-0.1,-0.05),fancybox = True, numpoints = 1)


# create fiugre for dependence on B: figure 2
f2=plt.figure(figsize=(16,16), facecolor = 'white')
gs = gridspec.GridSpec(2, 2,				# ratio of grid space (2 plots per collumn, 3 per row)
                       width_ratios=[1,1],		# width ratios of the 3 plots per row
                       height_ratios=[1,1]		# height ratios of the 2 polots per collumn
                       )


# create axes with the ratios specified in gs above. ax5, ax6 are the axes for the colorbars:
ax5 = plt.subplot(gs[0])	# vs B
ax6 = plt.subplot(gs[1])	# vs max fluctuation
ax7 = plt.subplot(gs[2])	# vs q95
ax8 = plt.subplot(gs[3])	# vs Lp


# Plots ############################################################################################################


ax5.errorbar(B0,rho_x, yerr = rho_xerr, marker='s', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'k', ls = 'None', label = 'Emission Case')
ax5.errorbar(B0,rho_x2, yerr = rho_x2err, marker='o', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'b', ls = 'None', label = 'Density Case')
ax5.errorbar(B0,rho_x3, yerr = rho_x3err, marker='D', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'g', ls = 'None', label = 'Block Case')
ax5.errorbar(B0,rho_x4, yerr = rho_x4err, marker='>', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'r', ls = 'None', label = 'Real SNR Block Case')

#ax5.set_xlabel(r'$B_0$ (T)')			
ax5.set_ylabel(r'$\rho_x$ (cm)')
ax5.get_yaxis().set_label_coords(-0.12,0.5)
ax5.set_title(r'blob width $\rho_x$',y=1.02)

ax6.errorbar(B0,v_r, yerr = v_rerr, marker='s', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'k', ls = 'None')
ax6.errorbar(B0,v_r2, yerr = v_r2err, marker='o', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'b', ls = 'None')
ax6.errorbar(B0,v_r3, yerr = v_r3err, marker='D', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'g', ls = 'None')
ax6.errorbar(B0,v_r4, yerr = v_r4err, marker='>', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'r', ls = 'None')

#ax6.set_xlabel(r'$B_0$ (T)')			
ax6.set_ylabel(r'$v_r$ (m/s)')
ax6.set_ylim(0,300)
ax6.get_yaxis().set_label_coords(-0.12,0.5)
ax6.set_title(r'Average velocity $v_r$',y=1.02)

ax7.errorbar(B0,tau_B, yerr = tau_Berr, marker='s', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'k', ls = 'None')
ax7.errorbar(B0,tau_B2, yerr = tau_B2err, marker='o', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'b', ls = 'None')
ax7.errorbar(B0,tau_B3, yerr = tau_B3err, marker='D', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'g', ls = 'None')
ax7.errorbar(B0,tau_B4, yerr = tau_B4err, marker='>', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'r', ls = 'None')

ax7.set_xlabel(r'$B_0$ (T)')
ax7.set_ylabel(r'$\tau_B$ (ms)')
ax7.get_yaxis().set_label_coords(-0.12,0.5)
ax7.set_title(r'Self-Correlation time $\tau_B$',y=1.02)

ax8.errorbar(B0 ,Relem2, yerr = Relem2err,marker='o', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'b', ls = 'None', label = 'Density Case')
ax8.errorbar(B0 ,Relem, yerr = Relemerr,marker='s', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'k', ls = 'None', label = 'Emission Case')
ax8.errorbar(B0 ,Relem3, yerr = Relem3err,marker='D', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'g', ls = 'None', label = 'Block Case')
ax8.errorbar(B0 ,Relem4, yerr = Relem4err,marker='>', markersize = msize, capsize = csiz, markeredgewidth =mwid, elinewidth = ewid, color = 'r', ls = 'None', label = 'Real SNR BlockCase')

ax8.set_xlabel(r'$B_0$ (T)')
ax8.set_ylabel(r'$\delta I/I$ or $\delta n/n$')
ax8.get_yaxis().set_label_coords(-0.12,0.5)
ax8.set_title(r'Blob amplitude $\delta I/I$ or $\delta n/n$', y = 1.02)

#leg = ax4.legend(loc='upper center', bbox_to_anchor = (-0.1,-0.05),fancybox = True, numpoints = 1)
handles, labels = ax8.get_legend_handles_labels()
handles = [h[0] for h in handles]
leg2 = ax8.legend(handles, labels,loc='center left', bbox_to_anchor = (-0.9,-0.3), fancybox = True, numpoints = 1, ncol = 2, prop={'size':25})
Bax=[ax1, ax5,ax6,ax7,ax8]

for ax in Bax:
	ax.set_xlim(1.7,2.6)
	
ax5.annotate("a)", xycoords="axes fraction", xy = (xshift,yshift), va = 'center', ha = 'center', fontsize = 25)
ax6.annotate("b)", xycoords="axes fraction", xy = (xshift,yshift), va = 'center', ha = 'center', fontsize = 25)
ax7.annotate("c)", xycoords="axes fraction", xy = (xshift,yshift), va = 'center', ha = 'center', fontsize = 25)
ax8.annotate("d)", xycoords="axes fraction", xy = (xshift,yshift), va = 'center', ha = 'center', fontsize = 25)

plt.tight_layout()
if savenow2:
	if NoiseSmooth:
		f1.savefig('Fig_Blobfreq_analysis_TIME_SNR100_Smooth')
		f2.savefig('Fig_Magn_analysis_TIME_SNR100_Smooth')
	if not NoiseSmooth:
		f1.savefig('Fig_Blobfreq_analysis_r={0:}mm'.format(int(x[Refdec_ind]*100)),bbox_extra_artists=(leg1,), bbox_inches='tight')
		f2.savefig('Fig_Magn_analysis_r={0:}mm'.format(int(x[Refdec_ind]*100)),bbox_extra_artists=(leg2,), bbox_inches='tight')


plt.show()


	
