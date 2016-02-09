# -*- coding: utf-8 -*-

# Parameters to be used in files like the timefunction analysis and block_timefunction

#specify timepoint for snapshot [in multiples of timesteps used for simula calling]:
position=30						#for certain point in time

#specify reference detector position [cm]
Refdec=0.5				# Take care in Block-Case

# specify number of timesteps around a blob, in which a second event is not counted as another blob:
WinBar = 15

# specify beam width
b=1.2

# define beam center position on y-axis [cm]
bcenter=2

#settings for conditional average
ThresCon = 2				# Muliplier of Standard deviation for counting something as a blob
TWindMult = 40				# Multiple of one timestep (!! Always an even number --> clear Nullposition!!) for time window for standard blob

#only used, if radial scanning case is selected:
Resolution = 5				# Multiples of steps in x-direction (usually multiples of 0.02) for going through all positions as Refdec-positions

#only used, if y-variation is swichted on:
y_starting = 1
y_ending = 6
yResolution = 0.5			# [cm] defines, in which distance a run should be done for multiple beam-center positions.

# Define positions for representation of different blobs:
a_an=12
b_an=16
c_an=22