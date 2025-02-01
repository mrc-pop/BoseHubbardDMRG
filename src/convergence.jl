#/usr/bin/julia

using Dates
using Plots

# TODO Import simulations

PROJECT_ROOT = @__DIR__ # Absloute path up to .../BoseHubbardDMRG/convergence/

# TODO Set PROJECT_ROOT as the absolute path to BoseHubbardDMRG

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

# First test: nmax is fixed at 3, and cutoff is fixed at 1E-12

nmax = 3
cutoff = 1E-12
ModelParameters = [L, N, nmax, J, μ]
UserMode = "maxm" # "maxm" / "nsweeps"

if UserMode == "maxm"

	# First subtest: nsweeps is fixed at 5, investigate maxm

	nsweeps = 5
	MaxDims = [100, 250, 500, 750, 1000]	# Max dimensions to investigate

	plot(size=(600,400), 
	     xlabel=L"$\max_m$", ylabel=L"$E_g$",
	     title=L"Energy after %$nsweeps DMRG sweeps")
	for SweepsMaxDim in MaxDims
		
		# TODO Improve: use increasing dim for each sweep (mamx locally defined)
		maxm = SweepsMaxDim		# Equal dimension for all sweeps
		DMRGParameters = [nsweeps, maxm, cutoff]
		E, _, _, _, _, _ = RunDMRGAlgorithm(ModelParameters, DMRGParameters; verbose=true)
		plot!(maxm, E)

	end
	
	savefig(PROJECT_ROOT * "/convergence/maxm_plot.pdf")

elseif UserMode == "nsweeps"

	# Second subtest: maxm is fixed at 200, investigate nsweeps

	nmax = 3
	cutoff = 1E-12
	ModelParameters = [L, N, nmax, J, μ]

	maxm = 200
	DMRGSweeps = 2:2:20

	plot(size=(600,400),
		 xlabel=L"$n_\mathrm{sweeps}$", ylabel=L"$E_g$",
		 title=L"Energy using maximum bond link $%maxm$")
	for nsweeps in DMRGSweeps
		
		# TODO Improve: use increasing dims
		DMRGParameters = [nsweeps, maxm, cutoff]
		E, _, _, _, _, _ = RunDMRGAlgorithm(ModelParameters, DMRGParameters; verbose=true)
		plot!(nsweeps, E)

	end
	
	savefig(PROJECT_ROOT * "/convergence/nsweeps_plot.pdf")	
	
end

