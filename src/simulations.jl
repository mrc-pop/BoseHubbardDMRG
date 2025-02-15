#!/usr/bin/julia

PROJECT_ROOT = @__DIR__ # Absloute path up to .../BoseHubbardDMRG/src

# Include modules
include(PROJECT_ROOT * "/modules/dmrg.jl")
include(PROJECT_ROOT * "/modules/sweeps.jl")
include(PROJECT_ROOT * "/modules/subdomain_selection.jl")

# Include setup
include(PROJECT_ROOT * "/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/setup/simulations_setup.jl")

# Move back
PROJECT_ROOT *= "/.."	# Absloute path up to .../BoseHubbardDMRG/

function main()

	"""
	Note: the main() functions at this level are the only functions in the whole
	project referring to externally defined global variables. All simulations
	variables have been defined in ./setup/simulations_setup.jl and are here
	referred to with their names. Check the setup file to modify the
	definitions.
	"""
    
	ModeErrorMsg = "Input error: use option --horizontal, --rectangular or " *
				   "--rectangular-selection"
	
	if length(ARGS) != 1
		# If user does not specify the user mode
		error(ModeErrorMsg)
		exit()
	else
		
		UserMode = ARGS[1]

        # -------------------------- Horizontal sweep --------------------------
        
        if UserMode=="--horizontal"

			# Horizontal sweeps
			LL = HorizontalLL				# Imported from setup
	    	JJ = HorizontalJJ				# Imported from setup
	    	μ0μ0 = Horizontalμμ				# Imported from setup
	    	
	    	DirPathOut = PROJECT_ROOT * "/simulations/horizontal_sweep"
    		mkpath(DirPathOut)
	    	
	    	for μ0 in μ0μ0
	    		
	    		FilePathOut = DirPathOut * "/μ0=$(μ0)_L=$(LL).txt"
	    	
				DataFile = open(FilePathOut,"w")
				write(DataFile,"# Hubbard model DMRG. This file contains many sizes. nmax=$nmax, μ0=$μ0, nsweeps=$nsweeps, cutoff=$cutoff\n")
				write(DataFile,"# L; J; E; μ_c^+; μ_c^-; nVariance; Γ; eΓ; C; eC [calculated $(now())]\n")
				close(DataFile)
	    	
		    	for L in LL
    	    		println("Starting calculation of observables for L=$L...")
    	    		
    	    		# Note: here we use superfluid DMRG parameters everywhere
					HorizontalSweep(L, nmax, JJ, DMRGParametersSF, FilePathOut; μ0=μ0,)
				end
				
				DataFile = open(FilePathOut,"a")
				write(DataFile,"# [finished at $(now())]\n")
				close(DataFile)
				
			end
					
			println("Done!")

        # -------------------------- Rectangular sweep -------------------------
        
		elseif UserMode=="--rectangular"
		
			# Rectangular sweeps
	    	JJ = RectangularJJ				# Imported from setup
	    	LL = RectangularLL				# Imported from setup
	    	μμ = Rectangularμμ				# Imported from setup

			μ0 = 0.0 # μ0 from which to take phase boundaries

			NN = LL							# Unitary filling

			# Uncomment here to use L-wise phase boundaries
			# FilePathIn =  PROJECT_ROOT * "/simulations/horizontal_sweep/L=$LL.txt"
			FilePathIn =  PROJECT_ROOT * "/analysis/phase_boundaries/μ0=$μ0/fitted_phase_boundaries_μ0=$μ0.txt"

			for (j,L) in enumerate(LL)
				N = NN[j]
				i = ceil(Int64, L/2) # Site to calculate variance on
				
				DirPathOut = PROJECT_ROOT * "/simulations/rectangular_sweep"
	    		mkpath(DirPathOut)
				FilePathOut = DirPathOut * "/L=$(L)_site=$i.txt"
				
				println("Starting calculation of observables for L=$L...")
				RectangularSweep(i, L, N, nmax, JJ, μμ, DMRGParametersMI, DMRGParametersSF, FilePathIn, FilePathOut)
				
				DataFile = open(FilePathOut,"a")
				write(DataFile,"# [finished at $(now())]\n")
				close(DataFile)
			end
			
			println("Done!")

        # ----------------------------------------------------------------------
        # --------------------- Rectangular selection sweep --------------------
        # ----------------------------------------------------------------------
		elseif UserMode=="--rectangular-selection"
		
			# Rectangular selection sweep
			LL = RectangularSelectionLL		# Imported from setup
			
			NN = LL							# Unitary filling
			
			println("Starting tip selection...")
			FilePathIn = PROJECT_ROOT * "/analysis/phase_boundaries/fitted_phase_boundaries.txt"
			Selections = FindMottTip(FilePathIn; verbose=false)
			PlotSelection(FilePathIn, Selections)
			
			print("Should I run a rectangular simulation inside the selected region? (y/n) ")
			PerformTipSweep = readline()
			println("You have chosen ", PerformTipSweep)
			while true
				if PerformTipSweep == "y"
					println("Selection accepted. Starting simulations around Mott lobe tip...")
					break
				elseif PerformTipSweep == "n"
					println("Selection rejected. Change the setting of this program to improve selection.")
					exit()
				else
					print("Invalid input. Please use valid input. (y/n) ")
					PerformTipSweep = readline()
				end
			end
			
			if size(Selections, 1) > 1
				print("There are $(size(Selections,1)) intersection points. Which one do you select? ")
				UserTipSelection = parse(Int64, readline())
				if UserTipSelection>0 & UserTipSelection < size(Selections,1)
					Left, Right, Up, Down = Selections[UserTipSelection,:]
				else
					print("Invalid input. Please use valid input. (1, ..., $(size(Selections,1)))")
				end
			else
				Left, Right, Up, Down = Selections
			end

			JJ = collect(range(start=Left, stop=Right, length=Resolution["Medium"]))
			μμ = collect(range(start=Down, stop=Up, length=Resolution["Medium"]))

			DirPathOut = PROJECT_ROOT * "/simulations/rectangular_sweep_tip"
			mkpath(DirPathOut)
			FilePathOut = DirPathOut * "/L=$LL.txt"

			DataFile = open(FilePathOut,"w")
			write(DataFile,"# Hubbard model DMRG. This file contains many sizes. nmax=$nmax, nsweeps=$nsweeps, cutoff=$cutoff\n")
			write(DataFile,"# L; J; E; nVariance; C; eC [calculated $(now())]\n")
			close(DataFile)

			for (j, L) in enumerate(LL)
				N = NN[j]
				i = ceil(Int64, L/2) # Site to calculate variance on
				
				DirPathOut = PROJECT_ROOT * "/simulations/tip_sweep"
	    		mkpath(DirPathOut)
				FilePathOut = DirPathOut * "/L=$(L)_site=$i.txt"	
				
				println("Starting calculation of observables for L=$L...")
				RectangularSweep(i, L, N, nmax, JJ, μμ, DMRGParametersMI, DMRGParametersSF, FilePathIn, FilePathOut; ComputeCorrelators=true)
				
				DataFile = open(FilePathOut,"a")
				write(DataFile,"# [finished at $(now())]\n")
				close(DataFile)
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
