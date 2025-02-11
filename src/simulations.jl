#!/usr/bin/julia

using Base.Threads
using DelimitedFiles  # For writedlm

PROJECT_ROOT = @__DIR__ # Absloute path up to .../BoseHubbardDMRG/src

# Include the dmrg.jl file
include(PROJECT_ROOT * "/dmrg.jl")

PROJECT_ROOT *= "/.."	# Absloute path up to .../BoseHubbardDMRG/

# TODO Add """[overall introduction]"""

# ------------------------------------------------------------------------------
# ------------------------------ Sweep functions -------------------------------
# ------------------------------------------------------------------------------

# ------------------------------ Horizontal sweep ------------------------------

function HorizontalSweep(L::Int64,
                         nmax::Int64,
                         JJ::Array{Float64},
                         DMRGParameters::Vector{Any},                         
                         FilePathOut::String; 
                         μ0=0.0)
    """
    Calculate two relevant observables for the MI-SF transition.
        - Phase boundaries,
        - Correlation function.
    Use μ0=0 by default (SF phase).
    """
    
    DataFile = open(FilePathOut, "a")
    
    for (j, J) in enumerate(JJ)
        println("Running DMRG for J=$(round(J, digits=3)), μ=$μ0 (simulation ",
                "$j/$(length(JJ)) for L=$L)")
        
        # Use @sync to wait for all tasks to complete
        @sync begin
            # Create tasks for each DMRG call
            task1 = @spawn RunDMRGAlgorithm([L, L-1, nmax, J, μ0], 
                                             DMRGParameters; 
                                             FixedN = true)
            task2 = @spawn RunDMRGAlgorithm([L, L, nmax, J, μ0], 
                                             DMRGParameters; 
                                             ComputeGamma=true,
                                             ComputeC = true,
                                             FixedN = true)
            task3 = @spawn RunDMRGAlgorithm([L, L+1, nmax, J, μ0], 
                                             DMRGParameters; 
                                             FixedN = true)

            # Wait for all tasks to complete and collect results
            E1, _ = fetch(task1)
            E2, nVariance, Γ, eΓ, C, eC = fetch(task2)
            E3, _ = fetch(task3)

            # Store results in the E array
            E = [E1, E2, E3]

            # Calculate chemical potentials
            ΔEplus = E[3] - E[2]
            ΔEminus = E[2] - E[1]
            μUp = ΔEplus + μ0
            μDown = -ΔEminus + μ0

            # Write results to the file
            write(DataFile, "$L; $J; $(E[2]); $μUp; $μDown; $nVariance; $Γ; $eΓ; $C; $eC\n")
        end
    end
    
    close(DataFile)
    println("Data for L=$L saved on file!")
end

# ----------------------------- Rectangular sweep ------------------------------

function RectangularSweep(i::Int64,
    					  L::Int64,
    					  N::Int64,
    					  nmax::Int64,
    					  JJ::Vector{Float64},				# Horizontal sweep
						  μμ::Vector{Float64},				# Vertical sweep
    					  DMRGParametersMI::Vector{Any},
    					  DMRGParametersSF::Vector{Any},
    					  FilePathIn::String,				# To evaluate if a given point is MI or SF
    					  FilePathOut::String)
    """
    Calculate the variance of the number of particles on site i for a range of
    hopping J and chemical potential μ values. Results are saved to a file.
    """
    
    DataFile = open(FilePathOut, "w")
    write(DataFile,"# Hubbard model DMRG. L=$L, N=$N, nmax=$nmax\n")
    write(DataFile,"# NOTE: Different DMRG settings have been used for MI and SF. Check the code.\n")
    write(DataFile,"# J, μ, E, n_variance, <a_i> [calculated $(now()) @site  i=$i]\n")

	# Take data from fitting data of horizontal sweeps to separate MI from SF 
	# ( use ΔE^+(∞) and ΔE^-(∞) )
	BoundariesData = readdlm(FilePathIn, ',', Float64, '\n'; comments=true)

    for (j,J) in enumerate(JJ)
		# We take the best approximating J ∈ JJ, to assess whether we are in the MI or SF phase

		# Index = findall(==(J), BoundariesData[:,1]) # this would work if there is the EXACT J in the fit results
		Index = argmin(abs.(JJ .- J)) # this always works, gives the best approximation
		μUp = BoundariesData[Index,2][1]
		μDown = -BoundariesData[Index,3][1]
    	
		println("\nJ=$J, phase boundaries: μ^+=$μUp, μ^-=$μDown")

        for (m,μ) in enumerate(μμ)
        
            ModelParameters = [L, N, nmax, J, μ]
			inMottLobe = false
			
			if (μ>=μDown && μ<=μUp)
				inMottLobe=true
			end

			println("Running DMRG for L=$L, J=$(round(J, digits=3)), μ=$(round(μ, digits=3)) (simulation ",
			"$m/$(length(μμ)) in μ, $j/$(length(JJ)) in J, MI=$inMottLobe)")

			if inMottLobe
	            E, nVariance, aAvg = RunDMRGAlgorithm(ModelParameters,
	                                            DMRGParametersMI;
	                                            FixedN = false,
												RandomMPS = true)
	            write(DataFile,"$J, $μ, $E, $(nVariance[i]), $(aAvg[i]) # MI\n")
	        else
	        	E, nVariance, aAvg = RunDMRGAlgorithm(ModelParameters,
	                                            DMRGParametersSF;
  		                                        FixedN = false,
												RandomMPS = true)
	            write(DataFile,"$J, $μ, $E, $(nVariance[i]), $(aAvg[i]) # SF\n")
	        end
        end
    end

    close(DataFile)
end

# -------------------------------- Path sweep ----------------------------------

# function PathSweep(Path::NamedTuple(Name::String, Values::Vector{Tuple{Int64, Int64}}),
# 				   L::Int64,
# 				   N::Int64,
# 				   nmax::Int64,
# 				   DMRGParametersMI::Vector{Any},
# 				   DMRGParametersSF::Vector{Any},
# 				   FilePathIn::String,				# To evaluate if a given point is MI or SF
# 				   FilePathOut::String)
				   
# 	Sizes, Couplings, Energies, UpBoundaries, DownBoundaries, _, _ = readdlm(FilePathIn, ';', Float64, '\n'; comments=true)
# 	IndicesList = findall(==(L), Sizes)
				   
# 	for Point in Path.Values
# 		J, μ = Point
				                      
# 		# Possible improvement: we know how many simulations have been performed for each L
#     	Index = IndicesList[findall(==(J), Couplings[IndicesList])]
#     	μUp = UpBoundaries[Index]
#     	μDown = DownBoundaries[Index]
        
#         ModelParameters = [L, N, nmax, J, μ]
# 		inMottLobe = false
		
# 		if (μ>=μDown && μ<=μUp)
# 			inMottLobe=true
# 		end
		
# 		if inMottLobe
#             _, _, _, C = RunDMRGAlgorithm(ModelParameters,
#                                           DMRGParametersMI;
#                                           ComputeC = true,
#                                           FixedN = false)
#         else
#         	_, _, _, C = RunDMRGAlgorithm(ModelParameters,
#                                           DMRGParametersSF;
# 	                                      ComputeC = true,
# 	                                      FixedN = false)
#         end                
        
#         # Perform DFT
#         D = 0
#         q = 2*pi/L
#         for r in 1:length(C)
#         	D += (-1)^(2*r/L) * C[r] # Numerically smarter than using the imaginary unit
#         end
#         K = 1/(L*D)
        
#         DataFile = open(FilePathOut, "a")
#         write(DataFile, "$L, $J, $μ, $D, $K")
		
# 	end
# 	close(DataFile)
#     println("Data for L=$L saved on file!")
# end

# ------------------------------------------------------------------------------
# ------------------------------------ Main ------------------------------------
# ------------------------------------------------------------------------------

function main()

    # TODO Import DMRG parameters from user input
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

			#FilePathIn =  PROJECT_ROOT * "/simulations/horizontal_sweep/L=$LL.txt"
			FilePathIn =  PROJECT_ROOT * "/analysis/phase_boundaries/fitted_phase_boundaries.txt"

			for i in 1:length(LL)
			
				L = LL[i]
				N = NN[i]
			
				# TODO As for now, we are using for each L the local boundary.
				# TODO We may as well import the fitted function and locally
				# TODO calculate the thermodynamic boundary value.
				# ah, forse allora ho fatto proprio questo, usando direttamente i boundaries "veri" L → ∞
				
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
