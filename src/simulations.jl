#!/usr/bin/julia

PROJECT_ROOT = @__DIR__ # Absloute path up to .../BoseHubbardDMRG/src
include(PROJECT_ROOT * "/modules/dmrg.jl")
include(PROJECT_ROOT * "/modules/sweeps.jl")
include(PROJECT_ROOT * "/modules/subdomain_selection.jl")
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
    
	ModeErrorMsg = "Input error: use option --horizontal, --rectangular or --rectangular-selection"
	
	if length(ARGS) != 1
		# If user does not specify the user mode
		error(ModeErrorMsg)
		exit()
	else
		
		UserMode = ARGS[1]
		if UserMode=="--horizontal"
		
			# TODO Import settings

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
				HorizontalSweep(L, nmax, JJ, DMRGParametersSF, FilePathOut; μ0=μ0)
			end
			
			println("Done!")

		elseif UserMode=="--rectangular"
		
			# TODO Import settings

			# Rectangular sweep
		    LL = [12] # [10, 20, 30, 40, 50, 60, 70]
			NN = LL
		    JJ = [J for J in range(start=0.0, stop=0.35, length=20)]		# TODO Change, exclude J=0
		    μμ = [μ for μ in range(start=0.1, stop=1.0, length=20)]			# TODO Change

			# Uncomment here to use L-wise phase boundaries
			# FilePathIn =  PROJECT_ROOT * "/simulations/horizontal_sweep/L=$LL.txt"
			FilePathIn =  PROJECT_ROOT * "/analysis/phase_boundaries/fitted_phase_boundaries.txt"

			for j in 1:length(LL)
			
				L = LL[j]
				N = NN[j]
				i = ceil(Int64, L/2) # Site to calculate variance on
				
				DirPathOut = PROJECT_ROOT * "/simulations/rectangular_sweep"
	    		mkpath(DirPathOut)
				FilePathOut = DirPathOut * "/L=$(L)_site=$i.txt"
				
				println("Starting calculation of observables for L=$L...")
				RectangularSweep(i, L, N, nmax, JJ, μμ, DMRGParametersMI, DMRGParametersSF, FilePathIn, FilePathOut)
			end
			
			println("Done!")
			
		elseif UserMode=="--rectangular-selection"
			
			println("Starting tip selection...")
			FilePathIn = PROJECT_ROOT * "/analysis/phase_boundaries/fitted_phase_boundaries.txt"
			Selections = FindMottTip(FilePathIn; verbose=false)
			PlotSelection(FilePathIn, Selections)
			
			print("Should I run a rectangular simulation inside the selected region? (y/n) ")
			Waiting=true
			PerformTipSweep = readline()
			while Waiting
				if PerformTipSweep=="y"
					global Waiting=false
					println("Selection accepted. Starting simulations around Mott lobe tip...")
				elseif PerformTipSweep=="n"
					global Waiting=false
					println("Selection rejected. Change the setting of this program to improve selection.")
					exit()
				else
					global Waiting=true
					print("Invalid input. Please use valid input. (y/n) ")
					global PerformTipSweep = readline()
				end
			end
			
			LL = [10, 20, 30, 40, 50, 60, 70]
			NN = LL
			
			Left, Right, Up, Down = Selections[i,:]
			JJ = collect(range(start=Left, stop=Right, length=50))
			μμ = collect(range(start=Down, stop=Up, length=50))

			DirPathOut = PROJECT_ROOT * "/simulations/rectangular_sweep_tip"
			mkpath(DirPathOut)
			FilePathOut = DirPathOut * "/L=$LL.txt"

			DataFile = open(FilePathOut,"w")
			write(DataFile,"# Hubbard model DMRG. This file contains many sizes. nmax=$nmax, μ0=$μ0, nsweeps=$nsweeps, cutoff=$cutoff\n")
			write(DataFile,"# L; J; E; nVariance; C; eC [calculated $(now())]\n")
			close(DataFile)

			for j in 1:length(LL)
				
				L = LL[j]
				N = NN[j]
				i = ceil(Int64, L/2) # Site to calculate variance on
				
				DirPathOut = PROJECT_ROOT * "/simulations/tip_sweep"
	    		mkpath(DirPathOut)
				FilePathOut = DirPathOut * "/L=$(L)_site=$i.txt"	
				
				println("Starting calculation of observables for L=$L...")
				RectangularSweep(i, L, N, nmax, JJ, μμ, DMRGParametersMI, DMRGParametersSF, FilePathIn, FilePathOut; ComputeCorrelators=true)
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
