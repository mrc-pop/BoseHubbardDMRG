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
