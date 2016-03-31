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

RealCase = True
cmview = True
filename = ''

comp = True		# comparison of different thresholds etc


thressweep = False
shots = [29302,29303, 29306, 29307, 29308, 29309, 29310,29311,29312,29315]
if comp:
	shots = [29303,29303,29303,29303]
	suff = ['_thressweep','_threscon05','_threscon1','']
	label = [r'$\sigma$ varied in [1,2]', r'$\sigma$ = 0.5', r'$\sigma$ = 1',r'$\sigma$ = 2']
#29311,29312,29315]
		# specify number of measurement without leading zeros, which is the prefix of the file in form (will be extended later to e.g. 004) 
savenow = True

markers=['v', 'H', 'D', 'd', '^', '<', '>', '8', 's', 'p', '*', 'h']
colors = ['g','r','c','m','orange','coral','pink','grey','darkred','peachpuff','goldenrod','plum','lavender']

#def radial_analysis(filename,Measurement,Refdec_ind, Refdec_ind3, savenow, SNR, NoiseSmooth, fourCase):

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
allax=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
font = {'family' :'Arial', 'weight' : 'normal', 'size' : 14}
matplotlib.rc('font', **font)
# Plots ############################################################################################################


for run in range(len(shots)):
	shot = shots[run]
	if comp:
		suf = suff[run]
	if not thressweep and not comp:
		time, x_pos, t_end, f_B, rho_x, tau_B, v_r, v_rmist, vdmax, vimax, Relem, B0, Lp,q95 = np.loadtxt('{0:}RealresultsReal_data.txt'.format(shot), usecols = (3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17), unpack = True, skiprows=2)
	if thressweep:
		time, x_pos, t_end, f_B, rho_x, tau_B, v_r, v_rmist, vdmax, vimax, Relem, B0, Lp,q95 = np.loadtxt('{0:}Realresults_veltwind_thressweep.txt'.format(shot), usecols = (3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17), unpack = True, skiprows=2)
	if comp:
		time, x_pos, t_end, f_B, rho_x, tau_B, v_r, v_rmist, vdmax, vimax, Relem, B0, Lp,q95 = np.loadtxt('{0:}RealresultsReal_data{1:}.txt'.format(shot,suf), usecols = (3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17), unpack = True, skiprows=2)

	x1 = np.loadtxt('/afs/ipp-garching.mpg.de/u/bschiess/Analyzetools/Li-BES GUI/29302_LiBes_rhop_at2900ms.txt')
	if cmview:
		x_pos = -(x_pos-7.373)		#7.37 corresponds to the LCFS approximately

	#pdb.set_trace()
	#	time2, x_pos2,y_pos2, f_B2, rho_x2, tau_B2, v_r2, v_rmist2, vdmax2, vimax2, B02, Lp2,q952 = np.loadtxt('{0:03d}Denresults.txt'.format(Measurement), usecols = (3, 4, 5,  6, 7, 8, 10, 11, 12, 13, 14, 15, 16), unpack = True, skiprows=2)
	#	time3, x_pos3, y_pos3, f_B3, rho_x3, tau_B3, v_r3, v_rmist3, vdmax3, vimax3, B03, Lp3,q953 = np.loadtxt('{0:03d}Blockresults.txt'.format(Measurement), usecols = (3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16), unpack = True, skiprows=2)


	###################################################################################################################################
	########### Reshaping axes ########################################################################################################
	###################################################################################################################################

	# for Em-Case #####################################################################################################################
	t_min = min(t_end)
	lencount = 0
	for elem in range(len(t_end)):
		if t_end[elem] == t_min:
			lencount = lencount +1
	len_t_end = len(t_end)/lencount 			# number of y-steps

	lenx = lencount 				# number of x-steps
	x = [None]*lenx
	for i in range (lenx):
		x[i]=x_pos[i]

	if not cmview:
		#transform to rhop
		x= [None]*len(x)
		for x_ind in range(len(x)):
			x[x_ind] = x1[x_ind]

	###################################################################################################################################
	########### Reshaping axes ########################################################################################################
	###################################################################################################################################


	# for Em-Case #####################################################################################################################

	f_Bm = np.reshape(f_B, (len_t_end, lenx))
	rho_xm =np.reshape(rho_x, (len_t_end,lenx))
	v_rm =np.reshape(v_r, (len_t_end,lenx))
	tau_Bm =np.reshape(tau_B, (len_t_end,lenx))
	v_rmistm =np.reshape(v_rmist, (len_t_end,lenx))
	vdmaxm =np.reshape(vdmax, (len_t_end,lenx))
	vimaxm =np.reshape(vimax, (len_t_end,lenx))

	Relemm = np.reshape(Relem, (len_t_end,lenx))

#	print('shape of f_B_Mat', np.shape(f_Bm))
#	print('Number of beam positions over which has been averaged: ', lenx, len_t_end, lenx*len_t_end,  len(y_pos), len(f_B), f_B[0],f_B[-1])

	
	###################################################################################################################################
	########### Mean and error calculations ###########################################################################################
	###################################################################################################################################

	# for Em-Case #####################################################################################################################

	#pdb.set_trace()
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



	###################################################################################################################################
	########### FIGURE ################################################################################################################
	###################################################################################################################################



	# first subplot: density

	ax1.errorbar(x,rho_x,yerr = rho_xerr, marker=markers[run], color = colors[run])
	ax1.axvline(0,color='k', linestyle='-.')
	#ax1.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
	ax1.set_ylabel(r'blob width $\rho_x$ (cm)')
	ax1.get_yaxis().set_label_coords(-0.09,0.5)
	ax1.set_title(r'blob width $\rho_x$ (cm)')
	#ax1.plot((0,4),(2.3,2.3), 'k',) # arrow line
	#ax1.plot((0,0),(2.3,2.3), 'k', marker='<',) # left arrowhead
	#ax1.plot((4,4),(2.3,2.3), 'k', marker='>',) # right arrowhead
	#ax1.text(2,2.1, 'SOL', font)
	#ax1.text(-0.5,2.1, 'LCFS', font)

	tau_Bfig = ax2.errorbar(x,tau_B, yerr = tau_Berr, marker=markers[run], color = colors[run])
	ax2.axvline(0,color='k', linestyle='-.')
	#ax2.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
	ax2.set_ylabel(r'Self-Correlation time $\tau_B$ (ms)')
	ax2.get_yaxis().set_label_coords(-0.09,0.5)
	ax2.set_title(r'Self-Correlation time $\tau_B$')

	f_Bfig = ax3.errorbar(x,f_B, yerr = f_Berr, marker=markers[run], color = colors[run])
	ax3.axvline(0,color='k', linestyle='-.')
	#ax3.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
	ax3.set_ylabel('Number of Blobs')
	ax3.get_yaxis().set_label_coords(-0.09,0.5)
	ax3.set_title('Number of Blobs in observed time intervall {0:.2f}ms'.format(time[10]))

	v_rfig = ax4.errorbar(x,v_r, yerr = v_rerr,marker=markers[run], color = colors[run])
	ax4.axvline(0,color='k', linestyle='-.')
	#ax4.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
	ax4.set_ylabel(r'Average velocity $v_r$ (m/s)')
	ax4.set_ylim(-50,300)
	ax4.get_yaxis().set_label_coords(-0.11,0.5)
	ax4.set_title(r'Average velocity $v_r$')

	v_dmaxfig = ax5.errorbar(x,vdmax,yerr = vdmaxerr, color = colors[run],marker=markers[run])
	ax5.axvline(0,color='k', linestyle='-.')
	#ax5.set_xlabel('beam axis $x$ (cm)', labelpad= -10)			# switched of, since axis in row below
	ax5.set_ylabel(r'Maximum velocity $v_{dmax}$ (m/s)')
	ax5.set_ylim(0,4000)
	ax5.get_yaxis().set_label_coords(-0.11,0.5)
	ax5.set_title(r'Maximum velocity for real data')


	f_B2fig = ax6.errorbar(x,f_B/time[10]*1000,yerr = f_Berr/time[10]*1000, marker=markers[run], color = colors[run])
	ax6.axvline(0,color='k', linestyle='-.')
	if cmview:
		ax6.set_xlabel('beam axis $x$ (cm)', labelpad= 0)			# switched of, since axis in row below
	if not cmview:
		ax6.set_xlabel(r'$\rho_{pol}$', labelpad= 0)			# switched of, since axis in row below
	ax6.set_ylabel('Blob frequency $f_B$ (1/s)')
	ax6.get_yaxis().set_label_coords(-0.09,0.5)
	ax6.set_title('Blob frequency $f_B$')

	ax7.errorbar(x,vimax, yerr = vimaxerr,color =  colors[run],marker=markers[run],label='Real data')
	ax7.axvline(0,color='k', linestyle='-.')
	#ax7.set_xlabel('beam axis $x$ (cm)', labelpad= 0)			# switched of, since axis in row below
	ax7.set_ylabel(r'Maximum velocities $v_{imax}$ (m/s)')
	ax7.get_yaxis().set_label_coords(-0.11,0.5)
	ax7.set_ylim(0,4000)
	ax7.set_title(r'Maximum velocities for interpolated data')


	ax8.errorbar(x,Relem, yerr = Relemerr,color =  colors[run],marker=markers[run],label='#{0:}, {1:}'.format(shot,label[run]))
	ax8.axvline(0,color='k', linestyle='-.')
	if cmview:
		ax8.set_xlabel('beam axis $x$ (cm)', labelpad= 0)			# switched of, since axis in row below
	if not cmview:
		ax8.set_xlabel(r'$\rho_{pol}$', labelpad= 0)			# switched of, since axis in row below
	ax8.set_ylabel(r'relative emission $\delta I/I$ or $\delta n/n$')
	ax8.get_yaxis().set_label_coords(-0.11,0.5)
	ax8.set_title(r'Blob amplitude')


for ax in allax:
	if cmview:
		ax.set_xlim(x[-1]-0.3,x[0]+0.3)
	if not cmview:
		ax.set_xlim(x[-1]-0.01,x[0]+0.01)
	

leg = ax8.legend(loc='upper center', bbox_to_anchor = (-0.16,-0.18), numpoints = 1,ncol=4)
	# leg = ax7.legend(loc='upper center', bbox_to_anchor = (0.225,-0.18),fancybox = True, numpoints = 1)

if savenow:
	if RealCase and len(shots)>1 and not comp:
		f1.savefig('FigRadial_RealCase_Comparison')
	if RealCase and len(shots)>1 and thressweep and not comp:
		f1.savefig('FigRadial_RealCase_Comparison_thressweep')
	if RealCase and len(shots)<=1 and not comp:
		f1.savefig('{0:}FigRadial'.format(shot))
	if RealCase and comp:
		f1.savefig('{0:}FigRadial_threscomparison'.format(shot))

plt.show()

	#return B0[0], q95[0], Lp[0], f_B[Refdec_ind]/time[10]*1000, f_Berr[Refdec_ind]/time[10]*1000, rho_x[Refdec_ind], rho_xerr[Refdec_ind], tau_B[Refdec_ind], tau_Berr[Refdec_ind], v_r[Refdec_ind], v_rerr[Refdec_ind], Relem[Refdec_ind], Relemerr[Refdec_ind], f_B2[Refdec_ind]/time[10]*1000, f_B2err[Refdec_ind]/time[10]*1000, rho_x2[Refdec_ind], rho_x2err[Refdec_ind], tau_B2[Refdec_ind], tau_B2err[Refdec_ind], v_r2[Refdec_ind], v_r2err[Refdec_ind], Relem2[Refdec_ind], Relem2err[Refdec_ind],f_B3[Refdec_ind3]/time[10]*1000, f_B3err[Refdec_ind3]/time[10]*1000, rho_x3[Refdec_ind3], rho_x3err[Refdec_ind3], tau_B3[Refdec_ind3], tau_B3err[Refdec_ind3], v_r3[Refdec_ind3], v_r3err[Refdec_ind3], Relem3[Refdec_ind3], Relem3err[Refdec_ind3]
	#return 

'''if __name__ == "__main__":
	# call function from below:
	shot = 0
	Refdec_ind3 = 1
	Refdec_ind  = 30
	SNR = 1
	NoiseSmooth = False
	fourCase = False
	suff = 1

	radial_analysis(filename,shot,Refdec_ind, Refdec_ind3, savenow, SNR, NoiseSmooth, fourCase)

	plt.show()'''




