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

# TODO Print also maxlinkdim (screen printed)

function RunConvergenceMaxM(ModelParameters::Vector{Float64},
							nSweeps::Int64,
							MaxDims::Vector{Int64},
							Cutoff::Float64)

	# First subtest: nsweeps is fixed, investigate maxm
	FilePathOut = PROJECT_ROOT * "/maxdim_data.txt"
	
	DataFile = open(FilePathOut,"w")
	write(DataFile,"# Hubbard model DMRG. L=$L, N=$N, J=$J, μ=$μ, nmax=$nmax\n")
	write(DataFile,"# maxm, E [calculated $(now()) @ nsweeps=$nSweeps]\n")

	for SweepsMaxDim in MaxDims
		
		# TODO Improve: use increasing dim for each sweep (mamx locally defined)
		local MaxM = SweepsMaxDim				# Equal dimension for all sweeps
		DMRGParameters = [nSweeps, MaxM, Cutoff]
		E, _, _, _, _, _ = RunDMRGAlgorithm(ModelParameters, DMRGParameters; verbose=true)

		write(DataFile,"$MaxM, $E\n")

	end
	
	# TODO Separate plot from simulations
	
	plot(size=(600,400), 
	     xlabel=L"$\max_m$", ylabel=L"$E_g$",
	     formatter = :plain,
	     title=L"Energy after $%$nSweeps$ DMRG sweeps")
	    
	Data = readdlm(FilePathOut, ',', Float64, '\n'; comments=true)
	xData = Data[:,1]
	yData = Data[:,2]
	scatter!(xData, yData, label=L"$\max_m$")
	savefig(PROJECT_ROOT * "/maxdim_plot.pdf")
	
	return None
end

function RunConvergenceNSweeps(ModelParameters::Vector{Float64},
							   DMRGSweeps::Vector{Int64},
							   MaxM::Int64,
							   Cutoff::Float64)

	# Second subtest: maxm is fixed, investigate nsweeps
	
	# File write
	FilePathOut = PROJECT_ROOT * "/convergence/nsweeps_data.txt"
	
	DataFile = open(FilePathOut,"w")
	write(DataFile,"# Hubbard model DMRG. L=$L, N=$N, J=$J, μ=$μ, nmax=$nmax\n")
    write(DataFile,"# nsweeps, E [calculated $(now()) @ maxm=$MaxM]\n")

	for nSweeps in DMRGSweeps
	
		DMRGParameters = [nSweeps, MaxM, cutoff]
		E, _, _, _, _, _ = RunDMRGAlgorithm(ModelParameters, DMRGParameters; verbose=true)
		
		write(DataFile,"$nSweeps, $E\n")

	end
	
	# TODO Separate plot from simulations
	
	plot(size=(600,400), 
	     xlabel=L"$n_\mathrm{sw}$", ylabel=L"$E_g$",
	     formatter = :plain,
	     title=L"Energy using $%$MaxM$ max bond link")
	    
	Data = readdlm(FilePathOut, ',', Float64, '\n'; comments=true)
	xData = Data[:,1]
	yData = Data[:,2]
	scatter!(xData, yData, label=L"$n_\mathrm{sw}$")
	savefig(PROJECT_ROOT * "/nsweeps_plot.pdf")
	
	return None
end

function main()

	# TODO Import model parameters from user input
	L = 10
	N = 10
	nmax = 3
	J = 0.2
	μ = 0.6
	ModelParameters = [L, N, nmax, J, μ]
	
	# TODO import DMRG parameters from user input 
	Cutoff = [1E-12]

	error_msg = "Input error: use option --maxm or --nsweeps"
	
	if length(ARGS) != 1
		error(error_msg)
		exit()
	else
		UserMode = ARGS[1]
		
		if UserMode == "--maxm"
			nSweeps = 5
			MaxDims = [x for x in 10:10:100]	# Max dimensions to investigate
			println("Running convergence simulations on parameter \"maxm\" using the following model:
· L = $L
· N = $N
· nmax = $nmax
· J = $J
· μ = $μ
and the DMRG setup:
· nsweeps = $nSweeps
· cutoff = $Cutoff
The investigated values of maxm are: $MaxDims")
			
			RunConvergenceMaxM(ModelParameters,nSweeps,MaxDims,Cutoff)
			
		elseif UserMode == "--nsweeps"
			MaxM = 200
			DMRGSweeps = [x for x in 2:2:20]	# Number of sweeps to investigate
			println("Running convergence simulations on parameter \"nsweeps\" using the following model:
· L = $L
· N = $N
· nmax = $nmax
· J = $J
· μ = $μ
and the DMRG setup:
· maxm = $MaxM
· cutoff = $Cutoff
The investigated values of nsweeps are: $DMRGSweeps")

			RunConvergenceNSweeps(ModelParameters,DMRGSweeps,MaxM,Cutoff)
			
		else
			error(error_msg)
			exit()
		end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end

