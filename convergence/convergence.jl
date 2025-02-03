#/usr/bin/julia

using Dates
using Plots
using LaTeXStrings

PROJECT_ROOT = @__DIR__ # Absloute path up to .../BoseHubbardDMRG/convergence/

include(PROJECT_ROOT * "/../src/dmrg.jl")
include(PROJECT_ROOT * "/../src/graphic_setup.jl")

"""
Parameters to be investigated: 
- nmax (maximum number of bosons per site, expected to be the weakest one); 
- nsweeps (for 1D chains it should be a weaker parameter as opposed to maxm due
  to rapid convergence); 
- maxm (maximum bond link, its optimal value is expected to increase with the 
  lattice size); 
- cutoff (dynamical updates along bonds, it can be lowered until the sweep 
  maximum link dimension saturates to maxm).
  
Adopted strategy: for a single phase space point (thus specific a specific point
(J,μ) on the plane) we perform many simulation on the same lattice at fixed
setting (L,N) and analyze eventual saturation of the observables.
"""

# TODO Import setting L, N and parameters J, μ

L = 10
N = 10
J = 0.2
μ = 0.6

# First test: nmax is fixed at 3, and cutoff is fixed at 1E-12

nmax = 3
Cutoff = [1E-12]

ModelParameters = [L, N, nmax, J, μ]
UserMode = "maxm" # "maxm" / "nsweeps"

# TODO Print also maxlinkdim (screen printed)

if UserMode == "maxm"

	# First subtest: nsweeps is fixed at 5, investigate maxm

	nSweeps = 5
	# MaxDims = [100 250 500 750 1000]
	MaxDims = [x for x in 10:10:100]	# Max dimensions to investigate
	
	# File write
	FilePathOut = PROJECT_ROOT * "/maxdim_data.txt"
	DataFile = open(FilePathOut,"w")
	write(DataFile,"# Hubbard model DMRG. L=$L, N=$N, J=$J, μ=$μ, nmax=$nmax\n")
    write(DataFile,"# maxm, E [calculated $(now()) @ nsweeps=$nSweeps]\n")

	plot(size=(600,400), 
	     xlabel=L"$\max_m$", ylabel=L"$E_g$",
	     title=L"Energy after $%$nSweeps$ DMRG sweeps")
	for SweepsMaxDim in MaxDims
		
		# TODO Improve: use increasing dim for each sweep (mamx locally defined)
		MaxM = SweepsMaxDim				# Equal dimension for all sweeps
		DMRGParameters = [nSweeps, MaxM, Cutoff]
		E, _, _, _, _, _ = RunDMRGAlgorithm(ModelParameters, DMRGParameters; verbose=true)
		
		write(DataFile,"$MaxM, $E\n")
		scatter!([MaxM], [E])

	end
	
	savefig(PROJECT_ROOT * "/maxdim_plot.pdf")

elseif UserMode == "nsweeps"

	# Second subtest: maxm is fixed at 200, investigate nsweeps

	nmax = 3
	cutoff = 1E-12
	ModelParameters = [L, N, nmax, J, μ]

	MaxM = 200
	# DMRGSweeps = [2 3 4 5 8 10 15 20]
	DMRGSweeps = [x for x in 2:2:20]
	
	# File write
	FilePathOut = PROJECT_ROOT * "/convergence/nsweeps_data.txt"
	DataFile = open(FilePathOut,"w")
	write(DataFile,"# Hubbard model DMRG. L=$L, N=$N, J=$J, μ=$μ, nmax=$nmax\n")
    write(DataFile,"# nsweeps, E [calculated $(now()) @ maxm=$MaxM]\n")

	plot(size=(600,400),
		 xlabel=L"$n_\mathrm{sweeps}$", ylabel=L"$E_g$",
		 title=L"Energy using maximum bond link $%MaxM$")
	for nSweeps in DMRGSweeps
		
		# TODO Improve: use increasing dims
		DMRGParameters = [nSweeps, MaxM, cutoff]
		E, _, _, _, _, _ = RunDMRGAlgorithm(ModelParameters, DMRGParameters; verbose=true)
		
		write(DataFile,"$nSweeps, $E\n")
		scatter!([nSweeps], [E])

	end
	
	savefig(PROJECT_ROOT * "/nsweeps_plot.pdf")
	
end

