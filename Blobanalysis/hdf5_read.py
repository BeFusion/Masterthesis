# -*- coding: utf-8 -*-
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pdb

def Li_read_hdf5(fileName):

	# Read HDF5 file into Python, define geometry of machine
	#fileName = 'ASDEX.001_test5.h5'
	file = h5py.File(fileName,'r')
	R0 = 0 # Major radius in m
	a = 0 # Minor radius in m

	# Conversion factors from normalized to SI-units
	data = file['params/structure_data']
	param = file['params/structure_param']
	B0 = param.attrs['B0']
	q95 = param.attrs['q']
	Lp = param.attrs['Lp']
	n0= 1e19*param.attrs['ne0']
	T0 = param.attrs['Te0']
	Ti0 = param.attrs['Ti0']
	rhoS = param.attrs['rho_s']
	omegaCi = param.attrs['w_ci']
	outTime = param.attrs['out_time']
	otMult = param.attrs['otmult']
	LCF_ind = param.attrs['edge']
	LCF    = LCF_ind*param.attrs['dx']*rhoS*100		# read in LCF, convert it to cm


	# Read density field as array
#	dens = file['data/var2d/Density']
#	density = np.empty(dens.shape)
#	dens.read_direct(density)
#	density *= n0

	# Read electron temperature field as array
#	electTemp = file['data/var2d/Temperature']
#	electronTemperature = np.empty(electTemp.shape)
#	electTemp.read_direct(electronTemperature)
#	electronTemperature *= T0

	# Read ion temperature field as array
#	ionTemp = file['data/var2d/Ion_temp']
#	ionTemperature = np.empty(ionTemp.shape)
#	ionTemp.read_direct(ionTemperature)
#	ionTemperature *= T0


#	Read original x-axis for conversion
	BeamAxShape= file['data/var2d/grid/x']
	x2 = np.empty(BeamAxShape.shape)
	BeamAxShape.read_direct(x2)

	#Read Li axes
	BeamAxisShape= file['data/Li-Profiles/Li_beam_coordinate_X']
	x = np.empty(BeamAxisShape.shape)
	BeamAxisShape.read_direct(x)
	NewAxisShape= file['data/Li-Profiles/Li_beam_coordinate_X']			# create a new axis with old shape, to convert to LCF-dimension
	xaxis = np.empty(NewAxisShape.shape)
	NewAxisShape.read_direct(xaxis)
	yAxisShape= file['data/Li-Profiles/Li_coordinate_Y']
	y = np.empty(yAxisShape.shape)
	yAxisShape.read_direct(y)
	tAxisShape= file['data/Li-Profiles/Li_time']
	t = np.empty(tAxisShape.shape)
	tAxisShape.read_direct(t)
	  
	#Convert t-axis to axis in ms:
	time=t/omegaCi*1000

	# Read Li density field as array
	#Shape of LiDensity: (time, beam-x-axis, y-axis)
	LiDensityShape = file['data/Li-Profiles/Li_2D_density']
	LiDensity = np.empty(LiDensityShape.shape)
	LiDensityShape.read_direct(LiDensity)

	# Read Li temperature field as array
	# Shape of LiTemp: (time, beam-x-axis, y-axis)
	LiTempShape = file['data/Li-Profiles/Li_2D_temperature']
	LiTemp = np.empty(LiTempShape.shape)
	LiTempShape.read_direct(LiTemp)

	# Read Li density field as array
	# Shape of LiDensity: (time, beam-x-axis, y-axis)
	Li2pShape = file['data/Li-Profiles/Li_2D_Li2p']
	Li2p = np.empty(Li2pShape.shape)
	Li2pShape.read_direct(Li2p)

	#Convert Li-x-Axis to axis with LCF as zero-position
	slope_dens=[None]*(len(x)-1)
	slope_dens=np.array(slope_dens)						# convert list to array
	slope_dens[0]=0
	for q in range (1,len(x)-1):			# Calculate slopes between every point
		 slope_dens[q] = (LiDensity[1,1,q]-LiDensity[1,1,q-1])/abs(xaxis[q]-xaxis[q-1])
	slope_dens_max=np.amax(slope_dens)	

	for q in range (1,len(x)-1):
		 if (slope_dens[q]==slope_dens_max):
			 slope_ind = q

	LCF_dens_max = xaxis[slope_ind]

	#x2[0] corresponds to xaxis[-1] and xaxis[0] corresponds to x2[0]+xaxis[-1]
#	LCF_x2pos=x2[-1,-1]-x2[LCF_ind,LCF_ind]			# distance of LCF from the zero point of x[2
	for m in range (0, len(x)):
#		x[m]=xaxis[len(x)-m-1]+LCF_x2pos-xaxis[-1]
		x[m] = -xaxis[-1]+LCF_dens_max+xaxis[len(x)-m-1]


# 	only use, if hesel run is incomplete (time-axis is too long and this will cause errors:	
	time2 = time[0:len(LiDensity[:,1,1])]
	

	endtime = time2[-1]

	return time2, x, y, LiDensity, LiTemp, Li2p, B0, Lp, q95, T0, Ti0, n0, omegaCi, rhoS, endtime

if __name__ == "__main__":
	t, x, y, LiDensity, LiTemp, Li2p = Li_read_hdf5('ASDEX.001_test5.h5')
	#meshgrid for plots
	X,Y= np.meshgrid(x,y)

	density = plt.contourf(Y,X,LiDensity[34,:,:],300)
	#plt.set_xlabel('axis $y$ (cm)')			# switched of, since axis in row below
	#plt.set_ylabel('beam axis $x$ (cm)')
	#set_title('input density $n_e$ (cm$^{-3}$)')
	#set_xlim(xmin,xmax)
	#set_ylim(ymin,ymax)
	chbar_den = plt.colorbar(density)
	chbar_den.set_label('density $n_e$ (cm$^{-3}$)', size=10)

	plt.show()
