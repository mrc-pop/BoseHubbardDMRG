#!/usr/bin/julia

LL = collect(L for L in 10:10:70)		# Sizes
NN = LL									# Unitary filling
nmax = 4								# Max: four bosons per site

# Mott Insulator DMRG parameters
nsweeps = 10
maxlinkdim = [10,50,75,200,500]
cutoff = [1E-8]
DMRGParametersMI = [nsweeps, maxlinkdim, cutoff]

# Superfluid phase DMRG paramters
nsweeps = 20
maxlinkdim = [10,50,75,200,500]
cutoff = [1E-8]
DMRGParametersSF = [nsweeps, maxlinkdim, cutoff]

Resolution = Dict([
			 	("Low", 10),			# 10*10 steps
			 	("Medium", 30),			# 30*30 steps
			 	("High", 50)			# 50*50 steps
			 ])
			 
# Horizontal sweeps
HorizontalLL = [10, 20, 30, 40, 50, 60, 70] 
HorizontalJJ = collect(range(start=0.0, stop=0.35, length=Resolution["High"]))
Horizontalμμ = [0.2, 0.4, 0.6, 0.8]

# Rectangular sweeps
RectangularLL = [20] # [10, 15, 20] 
RectangularJJ = collect(range(start=0.0, stop=0.35, length=Resolution["Medium"]))
Rectangularμμ = collect(range(start=0.0, stop=1.0, length=Resolution["Medium"]))

# Rectangular selection sweeps
RectangularSelectionLL = [50] # [40, 50, 60]
# JJ and μμ defined inside the script
