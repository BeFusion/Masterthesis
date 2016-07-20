# -*- coding: utf-8 -*-

# Parameters to be used in files like the timefunction analysis and block_timefunction

#specify timepoint for snapshot [in multiples of timesteps used for simula calling]:
position=200						#for certain point in time

#specify reference detector position [cm]
Refdec=2.1					# Take care in Block-Case
Refdec_ind_real = 6 				# if real data is analyzed, the index of the Reference detector has to be set here

# specify beam width
b=1.2

# define beam center position on y-axis [cm]
bcenter=1

#settings for conditional average
ThresCon = 2				# Muliplier of Standard deviation for counting something as a blob
TWindMult = 34				# time window for conditional averaging: Multiple of one timestep (!! Always an even number --> clear Nullposition!!) for time window for standard blob
TWindMult_long = 60			# only for calculation of FWHM of tau_B --> does not affect the conditional averaging

# specify number of timesteps around a blob, in which a second event is not counted as another blob:
WinBar = 12

# Specify noise level as SNR:

SNR = 1.00
smoothlenx = 35					# for x-smoothing: must be uneven, gives number of steps, taken into account for smoothing, smaller value for x-smoothing
smoothlen = 25					# for t-smoothing: must be uneven, gives number of steps, taken into account for smoothing
#smoothlen = 25	#originial .

#only used, if radial scanning case is selected:
SetResolution = 5				# Multiples of steps in x-direction (usually multiples of 0.02) for going through all positions as Refdec-positions
SetStartrange = int(1.5/0.02)			# Where to start running through all Refdec-positions

#only used, if y-variation is swichted on:
y_starting = 1
y_ending = 6
yResolution = 0.5			# [cm] defines, in which distance a run should be done for multiple beam-center positions.

# Define positions for representation of different blobs:

#High-Resolution Cases
a_an = 18
b_an = 25
c_an = 32

# RealCase
a_an_real=26
b_an_real=33
c_an_real=39
