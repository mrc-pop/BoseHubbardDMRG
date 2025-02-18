#!/usr/bin/julia

using DelimitedFiles

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/dmrg.jl")
include(PROJECT_ROOT * "/../setup/graphic_setup.jl")
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
	write(DataFile, "$ConserveNumber; $E; $LocalE; $nMean; $nVariance; $P; $S;\n")
	close(DataFile)
	
	return Observables
end

# ----------------------------------- Plots ------------------------------------

function PlotPopulations(DirPathOut::String,
						 P::Matrix{Float64},
						 ModelParameters::Vector{Float64},
						 SF::Bool,				# Redundant?
						 ConserveNumber::Bool)
	
	L, N, nmax = Int64.(ModelParameters[1:3])
	J, μ = ModelParameters[4:5]
	
	plot(grid=false,
		 xlabel=L"$i$ (site)",
		 ylabel=L"$|n\rangle$ (state)",
		 size = (440, 220),
		 minorticks = false)

	heatmap!(1:size(P,1), 0:(size(P,2)-1), P', clim=(0,1),
			 label=L"$\mathrm{Tr}\lbrace | i;n \rangle \langle i;n | \rho \rbrace$",
			 legend=true,
			 color=:coolwarm)

	# Draw by hand: grid
	m, n = size(P')
	vline!(0.5:(n+0.5), c=:black, alpha=0.2, linewidth=0.2, label=false)
	hline!(-0.5:(m-0.5), c=:black, alpha=0.2, linewidth=0.2, label=false)

	if SF
		if ConserveNumber
			title!(L"SF state population (fixed $N=%$N$, $J=%$J$, $\mu=%$μ$)")
			savefig(DirPathOut * "/SF_J=$(J)_μ=$(μ)_fixedN_populations.pdf")
		elseif !ConserveNumber
			title!(L"SF state population (optimal $N$, $J=%$J$, $\mu=%$μ$)")
			savefig(DirPathOut * "/SF_J=$(J)_μ=$(μ)_free_populations.pdf")
		end
	elseif !SF
		if ConserveNumber
			title!(L"MI state population (fixed $N=%$N$, $J=%$J$, $\mu=%$μ$)")
			savefig(DirPathOut * "/MI_J=$(J)_μ=$(μ)_fixedN_populations.pdf")
		elseif !ConserveNumber
			title!(L"MI state population (optimal $N$, $J=%$J$, $\mu=%$μ$)")
			savefig(DirPathOut * "/MI_J=$(J)_μ=$(μ)_free_populations.pdf")
		end
	end

end

function PlotBipartiteEntropy(DirPathOut::String,
						 	  S::Vector{Float64},
							  ModelParameters::Vector{Float64},
							  SF::Bool,				# Redundant?
							  ConserveNumber::Bool)
	
	L, N, nmax = Int64.(ModelParameters[1:3])
	J, μ = ModelParameters[4:5]
				  
	scatter([l for l in 1:length(S)], S,
		 	 xlabel=L"$\ell$ (link index)",
			 ylabel=L"$S$ (bipartition entropy)",
			 legend=false)

	if SF
		if ConserveNumber
			title!(L"SF bipartite entropy (fixed $N=%$N$, $J=%$J$, $\mu=%$μ$)")
			savefig(DirPathOut * "/SF_J=$(J)_μ=$(μ)_fixedN_entropy.pdf")
		elseif !ConserveNumber
			title!(L"SF bipartite entropy (optimal $N$, $J=%$J$, $\mu=%$μ$)")
			savefig(DirPathOut * "/SF_J=$(J)_μ=$(μ)_free_entropy.pdf")
		end
	elseif !SF
		if ConserveNumber
			title!(L"MI bipartite entropy (fixed $N=%$N$, $J=%$J$, $\mu=%$μ$)")
			savefig(DirPathOut * "/MI_J=$(J)_μ=$(μ)_fixedN_entropy.pdf")
		elseif !ConserveNumber
			title!(L"MI bipartite entropy (optimal $N$, $J=%$J$, $\mu=%$μ$)")
			savefig(DirPathOut * "/MI_J=$(J)_μ=$(μ)_free_entropy.pdf")
		end
	end

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
		
		DirPathOut = PROJECT_ROOT * "/analysis/states_properties/"
		mkpath(DirPathOut)
		
		for ConserveNumber in [true false]
	
			global Counter += 1
			print("Simulation ($Counter/4): SF=$SF and FixedN=$ConserveNumber.")
	
			if Compute=="y"
			DataFile = open(FilePathOut, "w")
			write(DataFile, "# FixedN; E; LocalE; nMean; nVariance; Populations; Entropy\n")
			close(DataFile)
				Observables = GetStateProperties(FilePathOut, ModelParameters, ConserveNumber)
				_, _, _, _, P, S = Observables
			elseif Compute=="n"
				ObservablesData = readdlm(FilePathOut, ';', Any)
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
