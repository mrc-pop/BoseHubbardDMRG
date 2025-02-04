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

function RunConvergenceMaxM(MaxMFilePath::String,
						    ModelParameters::Vector{Float64},
							nSweeps::Int64,
							MaxDims::Vector{Int64},
							Cutoff::Vector{Float64})

	# First subtest: nsweeps is fixed, investigate maxm
	
	DataFile = open(MaxMFilePath,"a")
	write(DataFile,"# maxm, E [calculated $(now()) @ nsweeps=$nSweeps]\n")

	for SweepsMaxDim in MaxDims
		
		# TODO Improve: use increasing dim for each sweep (mamx locally defined)
		local MaxM = SweepsMaxDim				# Equal dimension for all sweeps
		DMRGParameters = [nSweeps, MaxM, Cutoff]
		E, _, _, _, _, _, psi = RunDMRGAlgorithm(ModelParameters, DMRGParameters; verbose=true)

		write(DataFile,"$MaxM, $E\n")

	end
end

function RunConvergenceNSweeps(NSweepsFilePath::String,
							   ModelParameters::Vector{Float64},
							   DMRGSweeps::Vector{Int64},
							   MaxM::Int64,
							   Cutoff::Vector{Float64})

	# Second subtest: maxm is fixed, investigate nsweeps
	
	DataFile = open(NSweepsFilePath,"a")
    write(DataFile,"# nsweeps, E [calculated $(now()) @ maxm=$MaxM]\n")

	for nSweeps in DMRGSweeps
	
		DMRGParameters = [nSweeps, MaxM, Cutoff]
		E, _, _, _, _, _, psi = RunDMRGAlgorithm(ModelParameters, DMRGParameters; verbose=true)
		
		write(DataFile,"$nSweeps, $E\n")

	end
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

	ModeErrorMsg = "Input error: use option --maxm or --nsweeps"
	RunErrorMsg = "Input error: set y or n"
	Header = "# Hubbard model DMRG. L=$L, N=$N, J=$J, μ=$μ, nmax=$nmax\n"
	
	if length(ARGS) != 1
		
		# If user does not specify the user mode
		error(ModeErrorMsg)
		exit()
		
	else
	
		UserMode = ARGS[1]
	
		# Let the user decide to run simulations or not
		print("Should I run simulations (already simulated data will be overwritten)? (y/n) ")
		UserRun = readline()
		if UserRun == "y"
			run = true
			println("Starting computation...")
		elseif UserRun == "n"
			run = false
			println("Reading data...")
		else
			error(RunErrorMsg)
			exit()
		end
		
		if UserMode == "--maxm"

			# First subtest: nsweeps is fixed, investigate maxm
		
			nSweeps = 5
			MaxDims = [x for x in 10:10:100]	# Max dimensions to investigate
			MaxMFilePath = PROJECT_ROOT * "/maxdim_data.txt"
			
			if run
				DataFile = open(MaxMFilePath,"w")
				write(DataFile,Header)
				println("Running convergence simulations on parameter \"maxm\". The investigated values of maxm are: $MaxDims")
				RunConvergenceMaxM(MaxMFilePath,ModelParameters,nSweeps,MaxDims,Cutoff)
			end
			
			println("Plotting...")
			plot(size=(600,400), 
				 xlabel=L"$\max_m$", ylabel=L"$E_g$",
				 formatter = :plain,
				 title=L"Energy after $%$nSweeps$ DMRG sweeps")

			Data = readdlm(MaxMFilePath, ',', Float64, '\n'; comments=true)
			xData = Data[:,1]
			yData = Data[:,2]
			scatter!(xData, yData, label=L"$\max_m$")
			savefig(PROJECT_ROOT * "/maxdim_plot.pdf")
			
		elseif UserMode == "--nsweeps"
		
			# Second subtest: maxm is fixed, investigate nsweeps	
		
			MaxM = 200
			DMRGSweeps = [x for x in 2:2:20]	# Number of sweeps to investigate
			NSweepsFilePath = PROJECT_ROOT * "/nsweeps_data.txt"
			
			if run
				DataFile = open(NSweepsFilePath,"w")
				write(DataFile,Header)
				println("Running convergence simulations on parameter \"nsweeps\". The investigated values of nsweeps are: $DMRGSweeps")
				RunConvergenceNSweeps(NSweepsFilePath,ModelParameters,DMRGSweeps,MaxM,Cutoff)
			end
		
			println("Plotting...")
			plot(size=(600,400), 
				 xlabel=L"$n_\mathrm{sw}$", ylabel=L"$E_g$",
				 formatter = :plain,
				 title=L"Energy using $%$MaxM$ max bond link")
				
			Data = readdlm(NSweepsFilePath, ',', Float64, '\n'; comments=true)
			xData = Data[:,1]
			yData = Data[:,2]
			scatter!(xData, yData, label=L"$n_\mathrm{sw}$")
			savefig(PROJECT_ROOT * "/nsweeps_plot.pdf")
			
		else
			error(ModeErrorMsg)
			exit()
		end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end

