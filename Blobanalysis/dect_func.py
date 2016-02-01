#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

#Function for block beam calculations (returns new valued matrix for Li2p1 ##########################################################################
def dect_func(Li2p1,timestep,mn,b,avy,shift,shift2,x):

	Li2pcop=Li2p1.copy()
	Li2pcop.fill(0) 				# fill with zeros

	# s=int(b/avy)					# step forward and backward for average, b=beamwidth, avy=average stepsize in y-direction
	c=0.5						# bluriness in beam-axis direction (0.5cm
	stepnum=int(c/0.02)				# number of steps per block
	end2=int(mn/(c/0.02))				# number of blocks in x-direction
	Blo = end2-1
	
	dec_vec = [np.NaN]*Blo
	dec_ma = [np.NaN]*timestep*Blo
	dec_ma = np.array(dec_ma)
	dec_ma = dec_ma.reshape(timestep,1,Blo)	# shaped in order to cause no problems in the main program (time, y (is only one point now!), detector-x-position)

	dec_ind = 0
	for dec_ind in range(Blo):
		dec_vec[dec_ind] = x[dec_ind*stepnum + int(stepnum*3/2)]
		
	for t_index in range (0,timestep):
		for b_index in range (1,end2):		# parameter for number of blocks to produce in x-direction. Upper end
							# ..is given via number of steps in x-direction (mn) divided by c/0.02,
							# ..number of steps in beam width
			#calculating the rectangular entries for matrices:
			Li2pcop[t_index,shift*2-shift2:shift*2+shift2,mn-b_index*stepnum:mn-(b_index-1)*stepnum]=np.mean(Li2p1[t_index,shift*2-shift2:shift*2+shift2,mn-b_index*stepnum:mn-(b_index-1)*stepnum])
			dec_ma[t_index,0,end2-b_index-1]=np.mean(Li2p1[t_index,shift*2-shift2:shift*2+shift2,mn-b_index*stepnum:mn-(b_index-1)*stepnum])
	#		Li2p7[position,(((k-1)*s)+shift):((k*s)+shift),mn-(l*p):(mn-(l-1)*p)]= ( 	     # starting calculating the mean from the core
	#		np.mean(Li2p1[position,(((k-1)*s)+shift):((k*s)+shift),mn-(l*p):(mn-(l-1)*p)])    		     # ..outwards in order to hit the edge for sure
	#		)	



	Li2p1=Li2pcop.copy()

	return Li2p1, dec_ma, dec_vec, c


