#!/usr/bin/julia

using Base.Threads
using DelimitedFiles  # For writedlm

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
    
    l = length(JJ)
    for (j, J) in enumerate(JJ)
        println("Running DMRG for J=$(round(J, digits=3)), μ=$μ0 (simulation ",
                "$j/$l for L=$L)")
        
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
            E2, nVariance, aAvg, Γ, eΓ, C, eC = fetch(task2)
            E3, _ = fetch(task3)

            # Store results in the E array
            E = [E1, E2, E3]

            # Calculate chemical potentials
            ΔEplus = E[3] - E[2]
            ΔEminus = E[2] - E[1]
            μUp = ΔEplus + μ0
            μDown = -ΔEminus + μ0

            # Write results to the file
            # TODO Add aAvg, make it uniform in scripts
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

