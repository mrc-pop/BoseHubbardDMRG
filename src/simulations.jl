#!/usr/bin/julia

PROJECT_ROOT = @__DIR__ # Absloute path up to .../BoseHubbardDMRG/src
include(PROJECT_ROOT * "/modules/dmrg.jl")
include(PROJECT_ROOT * "/modules/sweeps.jl")
PROJECT_ROOT *= "/.."	# Absloute path up to .../BoseHubbardDMRG/

function main()

    # TODO Remove 11-23, import from simulations/setup
    nmax = 4

    # Mott Insulator DMRG parameters
    nsweeps = 10
    maxlinkdim = [10,50,75,200,500]		# TODO Change with optimal values
    cutoff = [1E-8]					    # TODO Change with optimal values
    DMRGParametersMI = [nsweeps, maxlinkdim, cutoff]

	# Superfluid phase DMRG paramters
	nsweeps = 20
	maxlinkdim = [10,50,75,200,500]		# TODO Change with optimal values
	cutoff = [1E-8]					    # TODO Change with optimal values
	DMRGParametersSF = [nsweeps, maxlinkdim, cutoff]
    
	ModeErrorMsg = "Input error: use option --horizontal, --rectangular or --path"
	
	if length(ARGS) != 1
		# If user does not specify the user mode
		error(ModeErrorMsg)
		exit()
	else
		
		UserMode = ARGS[1]
		if UserMode=="--horizontal"

			# Horizontal sweep
	    	# TODO Import model parameters from user input
	    	JJ = collect(range(start=0.0, stop=0.35, length=50))
	    	LL = [50, 60, 70] # [10, 20, 30, 40] 
			μ0 = 0.0
	    	
	    	DirPathOut = PROJECT_ROOT * "/simulations/horizontal_sweep"
    		mkpath(DirPathOut)
			FilePathOut = DirPathOut * "/L=$LL.txt"
	    	
	    	DataFile = open(FilePathOut,"w")
			write(DataFile,"# Hubbard model DMRG. This file contains many sizes. nmax=$nmax, μ0=$μ0, nsweeps=$nsweeps, cutoff=$cutoff\n")
			write(DataFile,"# L; J; E; deltaE_g^+; deltaE_g^-; nVariance; Γ; eΓ; C; eC [calculated $(now())]\n")
			close(DataFile)
	    	
	    	for L in LL
        		println("Starting calculation of observables for L=$L...")
				
				# TODO Evaluate: use superfluid DMRG paramters?
				HorizontalSweep(L, nmax, JJ, DMRGParametersSF, FilePathOut; μ0=μ0)
			end
			
			println("Done!")

		elseif UserMode=="--rectangular"

			# Rectangular sweep
		    NN = [12]					# TODO Change
		    LL = [12] 					# TODO Change
		    JJ = [J for J in range(start=0.0, stop=0.35, length=20)]		# TODO Change, exclude J=0
		    μμ = [μ for μ in range(start=0.1, stop=1.0, length=20)]			# TODO Change

			# Uncomment here to use L-wise phase boundaries
			# FilePathIn =  PROJECT_ROOT * "/simulations/horizontal_sweep/L=$LL.txt"
			FilePathIn =  PROJECT_ROOT * "/analysis/phase_boundaries/fitted_phase_boundaries.txt"

			for i in 1:length(LL)
			
				L = LL[i]
				N = NN[i]
				i = ceil(Int64, L/2) # Site to calculate variance on
				
				DirPathOut = PROJECT_ROOT * "/simulations/rectangular_sweep"
	    		mkpath(DirPathOut)
				FilePathOut = DirPathOut * "/L=$(L)_site=$i.txt"
				
				RectangularSweep(i, L, N, nmax, JJ, μμ, DMRGParametersMI, DMRGParametersSF, FilePathIn, FilePathOut)
			end
			
			println("Done!")
			
		elseif UserMode=="--path"

			# Arbitrary path sweep
		    N = 10 								# TODO Change
		    LL = [10, 20, 30] 					# TODO Change
		    Path = GeneratePath()				# TODO From scratch

			write(DataFile,"# Hubbard model DMRG along path $Path.Name. nmax=$nmax.\n")
			write(DataFile,"# L, J, μ, D, K [calculated $(now())]\n")

			for L in LL
				FilePathIn =  PROJECT_ROOT * "/simulations/horizontal_sweep.txt"
				
				# TODO Cycle over sites, extend to different fillings
				i = ceil(Int64, L/2) 				# Site to calculate variance on
				
				DirPathOut = PROJECT_ROOT * "/simulations/path=$Path.Name_sweep"
	    		mkpath(DirPathOut)
				FilePathOut = DirPathOut * "/L=$LL.txt"
				
				PathSweep(Path, L, N, nmax, DMRGParametersMI, DMRGParametersSF, FilePathIn, FilePathOut)
			end
			
			println("Done!")
			
		else
			error(ModeErrorMsg)
			exit()
		end
		
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
