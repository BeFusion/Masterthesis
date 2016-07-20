# -*- coding: utf-8 -*-
#!/usr/bin/env python


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


#scalex = 2.4

cmview = True
filename = ''

shots = [29302,29303, 29306, 29307, 29308, 29309, 29310,29311,29312,29315]
Measurements = [108,109,004,110,112,113]
file = ''

savenow = True

markers=['v', 'H', 'D', 'd', '^', '<', '>', '8', 's', 'p', '*', 'h']
colors = ['g','r','c','m','orange','coral','pink','grey','darkred','peachpuff','goldenrod','plum','lavender']

#def radial_analysis(filename,Measurement,Refdec_ind, Refdec_ind3, savenow, SNR, NoiseSmooth, fourCase):

# Plots ############################################################################################################

if len(shots)>len(Measurements):
	length = len(shots)
else:
	length = len(Measurements)

for run in range(length):
	try:
		Measurement = Measurements[run]
		time1, x_pos1, y_pos1, f_B1, rho_x1, tau_B1, v_r1, v_rmist1, vdmax1, vimax1, Relem1, B01, Lp1,q951, wci1, rho_s1 = np.loadtxt('{0:03d}Emresults_scalingtest{1:}.txt'.format(Measurement, file), usecols = (3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 21, 22), unpack = True, skiprows=2)
		time2, x_pos2, y_pos2, f_B2, rho_x2, tau_B2, v_r2, v_rmist2, vdmax2, vimax2, Relem2, B02, Lp2, q952, wci2, rho_s2 = np.loadtxt('{0:03d}Denresults_scalingtest{1:}.txt'.format(Measurement, file), usecols = (3, 4, 5,  6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 21, 22), unpack = True, skiprows=2)
		time3, x_pos3, y_pos3, f_B3, rho_x3, tau_B3, v_r3, v_rmist3, vdmax3, vimax3, Relem3, B03, Lp3, q953,wci3, rho_s3 = np.loadtxt('{0:03d}Blockresults_scalingtest{1:}.txt'.format(Measurement, file), usecols = (3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 21, 22), unpack = True, skiprows=2)
		time4, x_pos4, y_pos4, f_B4, rho_x4, tau_B4, v_r4, v_rmist4, vdmax4, vimax4, Relem4, B04, Lp4, q954, wci4, rho_s4 = np.loadtxt('{0:03d}Blockresults_scalingtest_TIME_SNR_Comparison.txt'.format(Measurement), usecols = (3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 21,22), unpack = True, skiprows=2)
	except:
		time1, x_pos1, y_pos1, f_B1, rho_x1, tau_B1, v_r1, v_rmist1, vdmax1, vimax1, Relem1, B01, Lp1,q951, wci1, rho_s1 = np.loadtxt('{0:03d}Emresults_scalingtest{1:}.txt'.format(Measurements[0], file), usecols = (3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 21, 22), unpack = True, skiprows=2)
		time2, x_pos2, y_pos2, f_B2, rho_x2, tau_B2, v_r2, v_rmist2, vdmax2, vimax2, Relem2, B02, Lp2, q952, wci2, rho_s2 = np.loadtxt('{0:03d}Denresults_scalingtest{1:}.txt'.format(Measurements[0], file), usecols = (3, 4, 5,  6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 21, 22), unpack = True, skiprows=2)
		time3, x_pos3, y_pos3, f_B3, rho_x3, tau_B3, v_r3, v_rmist3, vdmax3, vimax3, Relem3, B03, Lp3, q953,wci3, rho_s3 = np.loadtxt('{0:03d}Blockresults_scalingtest{1:}.txt'.format(Measurements[0], file), usecols = (3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 21, 22), unpack = True, skiprows=2)
		time4, x_pos4, y_pos4, f_B4, rho_x4, tau_B4, v_r4, v_rmist4, vdmax4, vimax4, Relem4, B04, Lp4, q954, wci4, rho_s4 = np.loadtxt('{0:03d}Blockresults_scalingtest_TIME_SNR_Comparison.txt'.format(Measurements[0]), usecols = (3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 21,22), unpack = True, skiprows=2)

	try:
		shot = shots[run]
		time, x_pos, t_end, f_B, rho_x, tau_B, v_r, v_rmist, vdmax, vimax, Relem, B0, Lp,q95 = np.loadtxt('{0:03d}Realresults_scalingtest.txt'.format(shot), usecols = (3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17), unpack = True, skiprows=2)
	except:
		time, x_pos, t_end, f_B, rho_x, tau_B, v_r, v_rmist, vdmax, vimax, Relem, B0, Lp,q95 = np.loadtxt('{0:03d}Realresults_scalingtest.txt'.format(shots[0]), usecols = (3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17), unpack = True, skiprows=2)


	xread = np.loadtxt('/afs/ipp-garching.mpg.de/u/bschiess/Analyzetools/Li-BES GUI/29302_LiBes_rhop_at2900ms.txt')
	if cmview:
		x_pos = -(x_pos-7.373)		#7.37 corresponds to the LCFS approximately

	#pdb.set_trace()
	##	time2, x_pos2,y_pos2, f_B2, rho_x2, tau_B2, v_r2, v_rmist2, vdmax2, vimax2, B02, Lp2,q952 = np.loadtxt('{0:03d}Denresults.txt'.format(Measurement), usecols = (3, 4, 5,  6, 7, 8, 10, 11, 12, 13, 14, 15, 16), unpack = True, skiprows=2)
	#	time3, x_pos3, y_pos3, f_B3, rho_x3, tau_B3, v_r3, v_rmist3, vdmax3, vimax3, B03, Lp3,q953 = np.loadtxt('{0:03d}Blockresults.txt'.format(Measurement), usecols = (3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16), unpack = True, skiprows=2)


	###################################################################################################################################
	########### Reshaping axes ########################################################################################################
	###################################################################################################################################

	# for Em-Case #####################################################################################################################
	y_min1 = min(y_pos1)
	lencount1 = 0
	for elem in range(len(y_pos1)):
		if y_pos1[elem] == y_min1:
			lencount1 = lencount1 +1
	leny1 = len(y_pos1)/lencount1 			# number of y-steps
	lenx1 = lencount1 				# number of x-steps
	x1 = [None]*lenx1
	for i in range (lenx1):
		x1[i]=x_pos1[i]

	y1 = [None]*leny1
	pos_ind1 = 0
	for i in range(len(y1)):
		y1[i]= y_pos1[pos_ind1]
		pos_ind1 = (i+1)*lenx1

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

# for Block-Case for SNR-Comparison #####################################################################################################################

	y_min4 = min(y_pos4)
	lencount4 = 0
	for elem in range(len(y_pos4)):
		if y_pos4[elem] == y_min4:
			lencount4 = lencount4 +1
	leny4 = len(y_pos4)/lencount4			# number of y-steps
	lenx4 = lencount4 				# number of x-steps
	x4 = [None]*lenx4
	for i in range (lenx4):
		x4[i]=x_pos4[i]

	y4 = [None]*leny4
	pos_ind4 = 0
	for i in range(len(y4)):
		y4[i]= y_pos4[pos_ind4]
		pos_ind4 = (i+1)*lenx4

# for Real-Case #####################################################################################################################
	
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
			x[x_ind] = xread[x_ind]

	###################################################################################################################################
	########### Reshaping axes ########################################################################################################
	###################################################################################################################################


	# for Real-Case #####################################################################################################################

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

	# for Real-Case #####################################################################################################################

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


	# for Em-Case #####################################################################################################################

	f_B1m = np.reshape(f_B1, (leny1, lenx1))
	rho_x1m =np.reshape(rho_x1, (leny1,lenx1))
	v_r1m =np.reshape(v_r1, (leny1,lenx1))
	tau_B1m =np.reshape(tau_B1, (leny1,lenx1))
	v_rmist1m =np.reshape(v_rmist1, (leny1,lenx1))
	vdmax1m =np.reshape(vdmax1, (leny1,lenx1))
	vimax1m =np.reshape(vimax1, (leny1,lenx1))

	Relem1m = np.reshape(Relem1, (leny1,lenx1))

#	print('shape of f_B_Mat', np.shape(f_Bm))
#	print('Number of beam positions over which has been averaged: ', lenx, leny, lenx*leny,  len(y_pos), len(f_B), f_B[0],f_B[-1])

	# for Den-Case ####################################################################################################################

	f_B2m = np.reshape(f_B2, (leny2, lenx2))
	rho_x2m =np.reshape(rho_x2, (leny2,lenx2))
	v_r2m =np.reshape(v_r2, (leny2,lenx2))
	tau_B2m =np.reshape(tau_B2, (leny2,lenx2))
	v_rmist2m =np.reshape(v_rmist2, (leny2,lenx2))
	vdmax2m =np.reshape(vdmax2, (leny2,lenx2))
	vimax2m =np.reshape(vimax2, (leny2,lenx2))

	Relem2m = np.reshape(Relem2, (leny2,lenx2))

#	print('shape of f_B_Mat', np.shape(f_B2m))
#	print('Number of beam positions over which has been averaged: ', lenx2, leny2, lenx2*leny2,  len(y_pos2), len(f_B2), f_B2[0],f_B[2-1])

	# for Block-Case ##################################################################################################################

#	print('Number of beam positions over which has been averaged: ', lenx3, leny3, lenx3*leny3,  len(y_pos3), len(f_B3), f_B3[0],f_B3[-1])


	f_B3m = np.reshape(f_B3, (leny3, lenx3))
	rho_x3m=np.reshape(rho_x3, (leny3,lenx3))
	v_r3m =np.reshape(v_r3, (leny3,lenx3))
	tau_B3m =np.reshape(tau_B3, (leny3,lenx3))
	v_rmist3m =np.reshape(v_rmist3, (leny3,lenx3))
	vdmax3m =np.reshape(vdmax3, (leny3,lenx3))
	vimax3m =np.reshape(vimax3, (leny3,lenx3))

	Relem3m = np.reshape(Relem3, (leny3,lenx3))

#	print('shape of f_B_Mat', np.shape(f_B3m))


	f_B4m = np.reshape(f_B4, (leny4, lenx4))
	rho_x4m=np.reshape(rho_x4, (leny4,lenx4))
	v_r4m =np.reshape(v_r4, (leny4,lenx4))
	tau_B4m =np.reshape(tau_B4, (leny4,lenx4))
	v_rmist4m =np.reshape(v_rmist4, (leny4,lenx4))
	vdmax4m =np.reshape(vdmax4, (leny4,lenx4))
	vimax4m =np.reshape(vimax4, (leny4,lenx4))

	Relem4m = np.reshape(Relem4, (leny4,lenx4))


	###################################################################################################################################
	########### Mean and error calculations ###########################################################################################
	###################################################################################################################################

	# for Em-Case #####################################################################################################################

	#pdb.set_trace()
	f_B1 = [None]*lenx1
	f_B1err = [None]*lenx1
	for f in range (lenx1):
		f_B1[f] = stats.nanmean(f_B1m[:,f])
		f_B1err[f] = stats.nanstd(f_B1m[:,f])
		
	rho_x1 = [None]*lenx1
	rho_x1err = [None]*lenx1
	for f in range (lenx1):
		rho_x1[f] = stats.nanmean(rho_x1m[:,f])
		rho_x1err[f] = stats.nanstd(rho_x1m[:,f])

	tau_B1 = [None]*lenx1
	tau_B1err = [None]*lenx1
	for f in range (lenx1):
		tau_B1[f] = stats.nanmean(tau_B1m[:,f])
		tau_B1err[f] = stats.nanstd(tau_B1m[:,f])

	v_r1 = [None]*lenx1
	v_r1err = [None]*lenx1
	for f in range (lenx1):
		v_r1[f] = stats.nanmean(v_r1m[:,f])
		v_r1err[f] = stats.nanstd(v_r1m[:,f])

	vdmax1 = [None]*lenx1
	vdmax1err = [None]*lenx1
	for f in range (lenx1):
		vdmax1[f] = stats.nanmean(vdmax1m[:,f])
		vdmax1err[f] = stats.nanstd(vdmax1m[:,f])

	vimax1 = [None]*lenx1
	vimax1err = [None]*lenx1
	for f in range (lenx1):
		vimax1[f] = stats.nanmean(vimax1m[:,f])
		vimax1err[f] = stats.nanstd(vimax1m[:,f])

	Relem1 = [None]*lenx1
	Relem1err = [None]*lenx1
	for f in range (lenx1):
		Relem1[f] = stats.nanmean(Relem1m[:,f])
		Relem1err[f] = stats.nanstd(Relem1m[:,f])


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


	# for Block-Case for real-SNR-comparison #####################################################################################################################

	f_B4 = [None]*lenx4
	f_B4err = [None]*lenx4
	for f in range (lenx4):
		f_B4[f] = stats.nanmean(f_B4m[:,f])
		f_B4err[f] = stats.nanstd(f_B4m[:,f])
		
	rho_x4 = [None]*lenx4
	rho_x4err = [None]*lenx4
	for f in range (lenx4):
		rho_x4[f] = stats.nanmean(rho_x4m[:,f])
		rho_x4err[f] = stats.nanstd(rho_x4m[:,f])

	tau_B4 = [None]*lenx4
	tau_B4err = [None]*lenx4
	for f in range (lenx4):
		tau_B4[f] = stats.nanmean(tau_B4m[:,f])
		tau_B4err[f] = stats.nanstd(tau_B4m[:,f])

	v_r4 = [None]*lenx4
	v_r4err = [None]*lenx4
	for f in range (lenx4):
		v_r4[f] = stats.nanmean(v_r4m[:,f])
		v_r4err[f] = stats.nanstd(v_r4m[:,f])

	vdmax4 = [None]*lenx4
	vdmax4err = [None]*lenx4
	for f in range (lenx4):
		vdmax4[f] = stats.nanmean(vdmax4m[:,f])
		vdmax4err[f] = stats.nanstd(vdmax4m[:,f])

	vimax4 = [None]*lenx4
	vimax4err = [None]*lenx4
	for f in range (lenx4):
		vimax4[f] = stats.nanmean(vimax4m[:,f])
		vimax4err[f] = stats.nanstd(vimax4m[:,f])


	Relem4 = [None]*lenx4
	Relem4err = [None]*lenx4
	for f in range (lenx4):
		Relem4[f] = stats.nanmean(Relem4m[:,f])
		Relem4err[f] = stats.nanstd(Relem4m[:,f])


	###################################################################################################################################
	########### FIGURE ################################################################################################################
	###################################################################################################################################

	#find position, at which scaling should take place
	'''for i in range (1,len(x)-1):
		if x[i]<=scalex and x[i-1]>scalex:
			scalex_ind = i
			xis = x[i]
	for i in range (1,len(x1)-1):
		if x1[i]<=scalex and x1[i-1]>scalex:
			scalex_ind1 = i
			xis1 = x1[i]
	for i in range (1,len(x2)-1):
		if x2[i]<=scalex and x2[i-1]>scalex:
			scalex_ind2 = i
			xis2 = x2[i]
	for i in range (1,len(x3)-1):
		if x3[i]<=scalex and x3[i-1]>scalex:
			scalex_ind3 = i
			xis3 = x3[i]
	for i in range (1,len(x4)-1):
		if x4[i]<=scalex and x4[i-1]>scalex:
			scalex_ind4 = i
			xis4 = x4[i]'''
	R = 2.10
	Te = 12
	cs = (Te/(2.014*931.494*10**6))**(0.5)*299792458
	rho_s = (2.014*931.494*10**6*Te)**(0.5)/(B0[0]*299792458)

	Te1 = 20
	R1 = 2.10
	rho_s1 = rho_s1[0]
	cs1 = (Te1/(2.014*931.494*10**6))**(0.5)*299792458

	R2 = 2.10
	rho_s2 = rho_s2[0]
	cs2 = (Te1/(2.014*931.494*10**6))**(0.5)*299792458

	R3 = 2.10
	rho_s3 = rho_s3[0]
	cs3 = (Te1/(2.014*931.494*10**6))**(0.5)*299792458

	R4 = 2.10
	rho_s4 = rho_s4[0]
	cs4 = (Te1/(2.014*931.494*10**6))**(0.5)*299792458

	tau_i = 1
	
	rho_x = rho_x[0]*0.01
	rho_x1 = rho_x1[0]*0.01
	rho_x2 = rho_x2[0]*0.01
	rho_x3 = rho_x3[0]*0.01
	rho_x4 = rho_x4[0]*0.01

	v_r = v_r[0]
	v_r1 = v_r1[0]
	v_r2 = v_r2[0]
	v_r3 = v_r3[0]
	v_r4 = v_r4[0]

	v_rerr = v_rerr[0]
	v_r1err = v_r1err[0]
	v_r2err = v_r2err[0]
	v_r3err = v_r3err[0]
	v_r4err = v_r4err[0]

	rho_xerr = rho_xerr[0]*0.01
	rho_x1err = rho_x1err[0]*0.01
	rho_x2err = rho_x2err[0]*0.01
	rho_x3err = rho_x3err[0]*0.01
	rho_x4err = rho_x4err[0]*0.01

	vimax = vimax[0]
	vimax1 = vimax1[0]
	vimax2 = vimax2[0]
	vimax3 = vimax3[0]
	vimax4 = vimax4[0]

	vimaxerr = vimaxerr[0]
	vimax1err = vimax1err[0]
	vimax2err = vimax2err[0]
	vimax3err = vimax3err[0]
	vimax4err = vimax4err[0]
	
	xshift = 0.05
	yshift = 0.92
	

	# create fiugre with certain ratio and facecolor:
	f1=plt.figure(figsize=(16,16), facecolor = 'white')
	gs = gridspec.GridSpec(2, 1,				# ratio of grid space (2 plots per collumn, 3 per row)
			width_ratios=[1],		# width ratios of the 3 plots per row
			height_ratios=[1,1]		# height ratios of the 2 polots per collumn
			)
	# create axes with the ratios specified in gs above. ax5, ax6 are the axes for the colorbars:
	ax1 = plt.subplot(gs[0])
	ax2 = plt.subplot(gs[1])	
	# create a list in for simple modifications on all plots at the same time:
	font = {'family' :'Arial', 'weight' : 'normal', 'size' : 24}
	matplotlib.rc('font', **font)

	rho_x = np.linspace(0,30*rho_s1,1000)
	



	ax1.plot(rho_x/rho_s1,(2*rho_x*Relem1[0]/R1)**(0.5),linewidth = 4, color = 'coral', label = 'Garcia: intertial scaling')
	ax1.plot(rho_x/rho_s1,(1+tau_i)*(rho_s1/rho_x)**2*Lp1[0]*Relem1[0]/R1, color = 'g',linewidth = 4, label = 'Krasheninnikov: sheath-connected')
	ax1.plot(rho_x/rho_s1,(2*rho_x/R)**(0.5)*cs1/(1+1/((rho_s1**2)*Lp1[0])*((R/2)**(0.5))*rho_x**(2.5))*Relem1[0]/cs1, color = 'b',linewidth = 4, label = 'Theiler')
	
	g=(1+tau_i)*(2*rho_x/R)*Relem1[0]
	gerr = (1+tau_i)*(2*rho_x/R)*Relem1err[0]
	f=(tau_i*rho_s1*Relem1[0]/2/rho_x)**2
	ferr=(tau_i*rho_s1*Relem2[0]/2/rho_x)**2
	ax1.plot(rho_x/rho_s1,(0.5*(f*(1+g**2/f**2)**(0.5)-f))**(0.5), color = 'pink', linewidth = 4,label = 'Manz: warm ions inertial')

	f2=(tau_i*rho_s1*Relem1[0]/2/rho_x+rho_x**3/Lp1[0]/(rho_s1**2))**2
	#ax1.plot(rho_x/rho_s1,(1+tau_i)*(rho_s1/rho_x)**2*Lp1[0]*Relem1[0]/R1,  linewidth = 4, color = 'orange', label = 'Manz: warm ion sheath-connected'.format(tau_i))
	ax1.plot(rho_x/rho_s1,(0.5*(f2*(1+(g**2)/(f2**2))**(0.5)-f2))**(0.5), color = 'midnightblue', linewidth = 4,label = 'Manz: warm ion sheath-connected 2')


	ax1.set_title(r'Density Case ($\tau_i$=1,$\delta I/I$={0:.2f}, $\rho_s$={1:.2e}, $L$={2:.2f}m)'.format(Relem1[0],rho_s1,Lp1[0]), y=1.02)
	ax1.set_ylabel(r'blob velocity $v_r/c_s$', x=-0.15)
	ax1.set_ylim(0,0.1)
	leg = ax1.legend(loc='best', fancybox = True, numpoints = 1, ncol = 1, prop={'size':23})



	ax2.plot(rho_x/rho_s4,(2*rho_x*Relem4[0]/R1)**(0.5),linewidth = 4, color = 'coral', label = 'Garcia: intertial scaling')
	ax2.plot(rho_x/rho_s4,(1+tau_i)*(rho_s4/rho_x)**2*Lp4[0]*Relem4[0]/R1, color = 'g',linewidth = 4, label = 'Krasheninnikov: sheath-connected')
	ax2.plot(rho_x/rho_s4,(2*rho_x/R)**(0.5)*cs4/(1+1/((rho_s4**2)*Lp4[0])*((R/2)**(0.5))*rho_x**(2.5))*Relem4[0]/cs4, color = 'b',linewidth = 4, label = 'Theiler')

	g=(1+tau_i)*(2*rho_x/R)*Relem4[0]
	gerr = (1+tau_i)*(2*rho_x/R)*Relem4err[0]
	f=(tau_i*rho_s4*Relem4[0]/2/rho_x)**2
	ferr=(tau_i*rho_s4*Relem4[0]/2/rho_x)**2
	ax2.plot(rho_x/rho_s4,(0.5*((g**2+f**2)**(0.5)-f))**(0.5), color = 'pink', linewidth = 4,label = 'Manz warm ions inertial')
	f2=(tau_i*rho_s4*Relem4[0]/2/rho_x+rho_x**3/Lp1[0]/(rho_s1**2))**2
	#ax2.plot(rho_x/rho_s4,(1+tau_i)*(rho_s4/rho_x)**2*Lp4[0]*Relem4[0]/R1,  linewidth = 4, color = 'orange', label = 'Manz: warm ion \n' r'sheath-connected')
	ax2.plot(rho_x/rho_s1,(0.5*(f2*(1+g**2/f2**2)**(0.5)-f2))**(0.5), color = 'midnightblue', linewidth = 4,label = 'Manz: warm ions sheath-connected')
	
	ax2.set_xlabel(r'blob width $\rho_x/\rho_s$')
	ax2.set_ylabel(r'blob velocity $v_r/c_s$', x=-0.15)
	ax2.set_ylim(0,0.1)
	#pdb.set_trace()
	ax2.set_title(r'Real SNR Block Case ($\tau_i$=1,$\delta I/I$={0:.2f}, $\rho_s$={1:.2e}, $L$={2:.2f}m)'.format(Relem4[0],rho_s4,Lp4[0]), y=1.02)


	f1.savefig('Scalinglaws')
	break


#if savenow:
	#plt.tight_layout()
	

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




