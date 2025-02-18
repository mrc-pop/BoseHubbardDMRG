#!/usr/bin/julia

using DelimitedFiles

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/../setup/graphic_setup.jl")
include(PROJECT_ROOT * "/dmrg.jl")
include(PROJECT_ROOT * "/plots.jl")
PROJECT_ROOT *= "/../.."

function GetStateProperties(FilePathOut::String,
							ModelParameters::Vector{Float64},
							ConserveNumber::Bool)

	# ------------------------------- Simulation ------------------------------- 

	nsweep = 10
	maxlinkdim = [10,50,75,200,500]
	cutoff = [1E-8]
	DMRGParameters = [nsweep, maxlinkdim, cutoff]

	Observables = RunDMRGAlgorithm(ModelParameters,
								   DMRGParameters,
								   "Debug";
								   FixedN=ConserveNumber,
								   RandomPsi0=false)
									
	E, nMean, nVariance, LocalE, P, S = Observables

	# --------------------------------- Write ----------------------------------

	DataFile = open(FilePathOut, "a")
	write(DataFile, "$ConserveNumber; $E; $LocalE; $nMean; $nVariance; $P; $S\n")
	close(DataFile)
	
	return Observables
end

# ----------------------------------- Main -------------------------------------
	
function main()

	global Counter = 0
	L = 30
	N = 30
	nmax = 4
	
	print("Should I run the simulations? (y/n) ")
	Compute = readline()
	
	for SF in [true false]
	
		DirPathOut = PROJECT_ROOT * "/simulations/states_properties/"
		mkpath(DirPathOut)
		
		if SF			# Superfluid test state
			J = 0.33
			μ = 0.8
			FilePathOut = DirPathOut * "SF_J=$(J)_μ=$(μ).txt"
			
		elseif !SF		# Mott insulating test state
			J = 0.06
			μ = 0.4
			FilePathOut = DirPathOut * "MI_J=$(J)_μ=$(μ).txt"
		
		end
		
		ModelParameters = [L, N, nmax, J, μ]
		
		DataFile = open(FilePathOut, "w")
		write(DataFile, "# FixedN; E; LocalE; nMean; nVariance; Populations; Entropy\n")
		close(DataFile)
		
		DirPathOut = PROJECT_ROOT * "/analysis/states_properties/"
		mkpath(DirPathOut)
		
		for ConserveNumber in [true false]
	
			global Counter += 1
			print("Simulation ($Counter/4): SF=$SF and FixedN=$ConserveNumber.")
	
			if Compute=="y"
				Observables = GetStateProperties(FilePathOut, ModelParameters, ConserveNumber)
				_, _, _, _, P, S = Observables
			elseif Compute=="n"
				ObservablesData = readdlm(FilePathOut, ';', Any, comments=true)
				for Index in 1:2
					if ObservablesData[Index,1]==ConserveNumber
						_, _, _, _, P, S = ObservablesData[2:7]
					end
				end			
			end
	
			PlotPopulations(DirPathOut, P, ModelParameters, SF, ConserveNumber)
			PlotBipartiteEntropy(DirPathOut, S, ModelParameters, SF, ConserveNumber)
			println(" Plots ready!")
		end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
