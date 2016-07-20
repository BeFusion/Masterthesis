#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pylab as plt
import matplotlib
import pdb
from matplotlib import rcParams

x, ne, Li2s, Li2p, Li3s, Li3p,Li4s = np.loadtxt('sim45.out', delimiter=',', usecols=(0, 1, 4, 5,6,7,8), unpack=True,skiprows=2)

# x45, ne45, Li2p45 = np.loadtxt('sim45.out', delimiter=',', usecols=(0, 1, 6), unpack=True,skiprows=2)

font = {'family' :'Arial', 'weight' : 'normal', 'size' : 20}
matplotlib.rc('font', **font)

#############################################################################
#############################################################################
#Purpose: This is for generation of example plots in thesis for explanation of Li-Beam.
# --> density, occupation number and shot example pictures

num = 4

#############################################################################
#############################################################################
#############################################################################

#############################################################################
# plot for density example

if num == 1:
	fig = plt.figure(num=None, figsize=(8, 6), facecolor='w')
	ax = fig.add_subplot(111)
	ax.plot(x,ne/(10**13),'b-',linewidth=4)

	ax.set_xlabel('beam axis x (cm)')
	ax.set_ylabel('density $n_e$ (cm$^{-3}$)')
	ax.set_xlim(0,17)
	ax.set_ylim(0,5)
	#plt.annotate('a)', xy=(0.05,0.9), fontsize = 22, xycoords='axes fraction')

	fig.savefig('density_example', bbox_inches='tight', dpi = 300)

#############################################################################
# plot for occupation number example

if num == 2:
	fig, ax1 = plt.subplots(figsize=(8, 6))
	ab, = ax1.plot(x,Li2s,'b-',linewidth=4, label = 'Li(2s)')
	ax1.set_xlabel('beam axis x (cm)')
	ax1.set_xlim(0,17)
	ax1.set_ylim(0,1)
	ax1.set_ylabel('Occupation number N$_{Li2s}$', color = 'b')
	for tl in ax1.get_yticklabels():
		tl.set_color('b')

	ax2 = ax1.twinx()
	ax2.set_xlim(0,19)
	bb, = ax2.plot(x,Li2p,'r-',linewidth=4, label = 'Li(2p)')
	cb, = ax2.plot(x,Li3s,'r--',linewidth=2, label = r'Li($\geq$3s)')
	ax2.plot(x,Li3p,'r--',linewidth=2)
	ax2.plot(x,Li4s,'r--',linewidth=2)
	ax2.set_ylabel('Occupation number N$_{Li2p}$', color = 'r')
	ax2.set_ylim(0,0.2)
	for tl in ax2.get_yticklabels():
		tl.set_color('r')
	leg = plt.legend((ab,bb,cb), ('Li(2s)','Li(2p)', r'Li($\geq$3s)'), loc='center left', fancybox = True, numpoints = 1)

	#plt.annotate('a)', xy=(0.05,0.9), fontsize = 22, xycoords='axes fraction')

	fig.savefig('occupation_example', bbox_inches='tight', dpi = 300)


#############################################################################
# plot for shot example
if num == 3:
	fig = plt.figure(num=None, figsize=(8, 6), facecolor='w')
	ax = fig.add_subplot(111)
	shot = 29306 #(at t = t[2000])
	Li2p1 = np.array([ 0.03243181,  0.01785843,  0.0664832 ,  0.15242967,  0.27252736,
		0.30272112,  0.5270197 ,  1.0589739 ,  1.35997985,  1.98507982,
		2.34491021,  2.95644483,  3.45147427,  3.77090283,  3.3382132 ,
		3.60251751,  4.00786373,  2.85789071,  2.63264945,  2.4540447 ,
		2.51683776,  1.24410304,  1.86663311,  1.90314567,  1.57757274,
		2.4088046 ])/10

	x = np.array([0,1.2,1.7,2.4,3,3.6,4.22,5.45,6.1,6.75,\
		7.37,8.02,8.65,9.3,9.95,10.58,11.21,11.86,12.56,13.26,\
		13.99,14.66,15.39,16.13,16.89,17.61]) 

	ax.plot(x,Li2p1,'b.', marker = 's', markersize = 9, markeredgecolor = 'b',  linewidth=4, label = '#{0:} at t=2.91s'.format(shot))

	ax.set_xlabel('beam axis x (cm)')
	ax.set_ylabel('Li$_{2p}$ emission (a.u.)')
	ax.set_xlim(0,17)
	#plt.ylim(0,0.2)
	#plt.annotate('c)', xy=(0.05,0.9), fontsize = 22, xycoords='axes fraction')
	leg = ax.legend(loc='lower right', fancybox = True, numpoints = 1)

	fig.savefig('shot_example', bbox_inches='tight', dpi = 300)

#############################################################################
# plot for occupation number example with integtrated density profile

if num == 4:
	fig, ax1 = plt.subplots(figsize=(8, 6))
	ab, = ax1.plot(x,Li2p,'r-',linewidth=4, label = 'Li(2p)')
	cb, = ax1.plot(x,Li3s,'r--',linewidth=2, label = r'Li($\geq$3s)')
	ax1.plot(x,Li3p,'r--',linewidth=2)
	ax1.plot(x,Li4s,'r--',linewidth=2)
	ax1.set_xlabel('beam axis x (cm)')
	ax1.set_xlim(0,17)
	ax1.set_ylim(0,0.2)
	ax1.set_ylabel('Occupation number N$_{Li}$', color = 'r')
	for tl in ax1.get_yticklabels():
		tl.set_color('r')

	ax2 = ax1.twinx()
	ax2.set_xlim(0,19)	
	bb, = ax2.plot(x,ne/(10**13),'b-',linewidth=4)
	ax2.set_ylabel('density $n_e$ (10$^{13}$cm$^{-3}$)', color = 'b')
	ax2.set_ylim(0,3)
	for tl in ax2.get_yticklabels():
		tl.set_color('b')
	leg = plt.legend((bb,ab,cb), ('$n_e$','Li(2p)', r'Li($\geq$3s)'), loc='center left', fancybox = True, numpoints = 1)

	#plt.annotate('a)', xy=(0.05,0.9), fontsize = 22, xycoords='axes fraction')

	fig.savefig('occupation_density_example', bbox_inches='tight', dpi = 300)



plt.show()
