#!/usr/bin/env python
# -*- coding: utf-8 -*-

########### Importing modules ###########################################################################

import numpy as np
import matplotlib.pyplot as plt
#from sympy import *
from pylab import plot
import matplotlib
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

Measurement = 49

time, x, f_B, rho_x, tau_B, v_r, v_rmist, vdmax, vimax, B0, Lp,q95, v_rspe, v_rmistspe, vdmaxspe, vimaxspe = np.loadtxt('{0:03d}_speedtest_Emresults.txt'.format(Measurement), usecols = (3, 4, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 23,24,25,26 ), unpack = True, skiprows=2)
time2, x2, f_B2, rho_x2, tau_B2, v_r2, v_rmist2, vdmax2, vimax2, B02, Lp2,q952, v_rspe2, v_rmistspe2, vdmaxspe2, vimaxspe2  = np.loadtxt('{0:03d}_speedtest_Denresults.txt'.format(Measurement), usecols = (3, 4, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 23,24,25,26), unpack = True, skiprows=2)
time3, x3, f_B3, rho_x3, tau_B3, v_r3, v_rmist3, vdmax3, vimax3, B03, Lp3,q953, v_rspe3, v_rmistspe3, vdmaxspe3, vimaxspe3  = np.loadtxt('{0:03d}_speedtest_Blockresults.txt'.format(Measurement), usecols = (3, 4, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 23,24,25,26), unpack = True, skiprows=2)

###################################################################################################################################
########### FIGURE 1 ##############################################################################################################
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


# create a list in for simple modifications on all plots at the same time:
allax=[ax1,ax2,ax3,ax4,ax5,ax6,ax7]

font = {'family' :'normal', 'weight' : 'regular', 'size' : 16}
matplotlib.rc('font', **font)

# Plots ############################################################################################################

# first subplot: density

ax1.plot(x,rho_x,marker='s', color = 'k',label='Emission data')
ax1.plot(x2,rho_x2,marker='o', color = 'b',label='Density data')
ax1.plot(x3,rho_x3,marker='D', color = 'g',label='Block data')
ax1.axvline(0,color='k', linestyle='-')
ax1.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
ax1.set_ylabel(r'blob width $\rho_x$ (cm)')
ax1.set_title(r'blob width $\rho_x$ (cm)')
ax1.plot((0,4),(2.3,2.3), 'k',) # arrow line
ax1.plot((0,0),(2.3,2.3), 'k', marker='<',) # left arrowhead
ax1.plot((4,4),(2.3,2.3), 'k', marker='>',) # right arrowhead
ax1.text(2,2.1, 'SOL')


leg = ax1.legend(loc='best', fancybox = True, numpoints = 1)


tau_Bfig = ax2.plot(x,tau_B,marker='s', color = 'k')
tau_Bfig2 = ax2.plot(x2,tau_B2,marker='o', color = 'b')
tau_Bfig3 = ax2.plot(x3,tau_B3,marker='D', color = 'g')
ax2.axvline(0,color='k', linestyle='-')
ax2.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
ax2.set_ylabel(r'Self-Correlation time $\tau_B$ (ms)')
ax2.set_title(r'Self-Correlation time $\tau_B$')

f_Bfig = ax3.plot(x,f_B,marker='s', color = 'k')
f_Bfig2 = ax3.plot(x2,f_B2,marker='o', color = 'b')
f_Bfig3 = ax3.plot(x3,f_B3,marker='D', color = 'g')
ax3.axvline(0,color='k', linestyle='-')
ax3.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
ax3.set_ylabel('Number of Blobs')
ax3.set_title('Number of Blobs in observed time intervall {0:.2f}ms'.format(time[10]))

v_rfig = ax4.errorbar(x,v_r, yerr = v_rmist,marker='s', color = 'k', label = 'COM-data')
v_rfig2 = ax4.errorbar(x2,v_r2, yerr = v_rmist2,marker='o', color = 'b')
v_rfig3 = ax4.errorbar(x3,v_r3, yerr = v_rmist3,marker='s', color = 'g')

v_rfigspe = ax4.errorbar(x,v_rspe, yerr = v_rmist,marker='s', linestyle = '--', color = 'k', label = 'Max-data')
v_rfigspe2 = ax4.errorbar(x2,v_rspe2, yerr = v_rmist2,marker='o', linestyle = '--', color = 'b')
v_rfigspe3 = ax4.errorbar(x3,v_rspe3, yerr = v_rmist3,marker='s', linestyle = '--', color = 'g')
ax4.axvline(0,color='k', linestyle='-')
ax4.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
ax4.set_ylabel(r'Average velocity $v_r$ (m/s)')
ax4.set_title(r'Average velocity $v_r$')

leg = ax4.legend(loc='best', fancybox = True, numpoints = 1)

v_dmaxfig = ax5.plot(x,vdmax, 'k',marker='s', label = 'COM-data')
v_dmax2fig = ax5.plot(x2,vdmax2, 'b',marker='o')
v_dmax3fig = ax5.plot(x3,vdmax3, 'g',marker='D')
v_dmaxfigspe = ax5.plot(x,vdmaxspe, 'k',marker='s', linestyle = '--', label = 'Max-data')
v_dmax2figspe = ax5.plot(x2,vdmaxspe2, 'b',marker='o',linestyle = '--')
v_dmax3figspe = ax5.plot(x3,vdmaxspe3, 'g',marker='D', linestyle = '--')
ax5.axvline(0,color='k', linestyle='-')
ax5.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
ax5.set_ylabel(r'Maximum velocity $v_{dmax}$ (m/s)')
ax5.set_title(r'Maximum velocity for real data')

leg = ax5.legend(loc='best', fancybox = True, numpoints = 1)

f_B2fig = ax6.plot(x,f_B/time[10]*1000,marker='s', color = 'k')
f_B22fig = ax6.plot(x2,f_B2/time[10]*1000,marker='o', color = 'b')
f_B23fig = ax6.plot(x3,f_B3/time[10]*1000,marker='D', color = 'g')
ax6.axvline(0,color='k', linestyle='-')
ax6.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
ax6.set_ylabel('Blob frequency $f_B$ (1/s)')
ax6.set_title('Blob frequency $f_B$')

v_imaxfig = ax7.plot(x,vimax, 'k',marker='s', label = 'COM-data')
v_imax2fig = ax7.plot(x2,vimax2, 'b',marker='o')
v_imax3fig = ax7.plot(x3,vimax3, 'g',marker='D')

v_imaxfigspe = ax7.plot(x,vimaxspe, 'k',marker='s', linestyle = '--', label = 'Max-data')
v_imax2figspe = ax7.plot(x2,vimaxspe2, 'b',marker='o')
v_imax3figspe = ax7.plot(x3,vimaxspe3, 'g',marker='D')

ax7.axvline(0,color='k', linestyle='-')
ax7.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
ax7.set_ylabel(r'Maximum velocities $v_{imax}$ (m/s)')
ax7.set_title(r'Maximum velocities for interpolated data')

leg = ax7.legend(loc='best', fancybox = True, numpoints = 1)


f1.savefig('{0:03d}_speedtest_FigRadial_Analysis'.format(Measurement))
plt.show()






