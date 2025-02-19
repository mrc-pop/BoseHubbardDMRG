#!/usr/bin/julia

using Base.Threads
using DelimitedFiles  # For writedlm

# ------------------------------ Boundaries sweep ------------------------------

function BoundariesSweep(L::Int64,
                         nmax::Int64,
                         JJ::Array{Float64},
                         DMRGParameters::Vector{Any},                         
                         FilePathOut::String; 
                         μ0=0.0)
    """
    Calculate two relevant observables for the MI-SF transition.
        - Phase boundaries,
    Use μ0=0 by default (SF phase).
    """
    
    DataFile = open(FilePathOut, "a")
    println("Extracting boundaries...")
    
    l = length(JJ)
    for (j, J) in enumerate(JJ)
        println("Running DMRG for J=$(round(J, digits=3)), μ=$μ0 (simulation ",
                "$j/$l for L=$L)")
        
        # Use @sync to wait for all tasks to complete
        @sync begin
            # Create tasks for each DMRG call
            task1 = @spawn RunDMRGAlgorithm([L, L-1, nmax, J, μ0], 
                                             DMRGParameters,
                                             "Fast"; 
                                             FixedN=true,
                                             RandomPsi0=false)
            task2 = @spawn RunDMRGAlgorithm([L, L, nmax, J, μ0], 
                                             DMRGParameters,
                                             "Fast";
                                             FixedN=true,
                                             RandomPsi0=false)
            task3 = @spawn RunDMRGAlgorithm([L, L+1, nmax, J, μ0], 
                                             DMRGParameters,
                                             "Fast"; 
                                             FixedN=true,
                                             RandomPsi0=false)

            # Wait for all tasks to complete and collect results
            E1 = fetch(task1)
            E2 = fetch(task2)
            E3 = fetch(task3)

            # Calculate chemical potentials
            ΔEplus = E3 - E2
            ΔEminus = E2 - E1
            μUp = ΔEplus + μ0
            μDown = - ΔEminus - μ0 # Sign problem otherwise

            # Write results to the file
            write(DataFile, "$L; $J; $E2; $μUp; $μDown\n")
        end
    end
    
    close(DataFile)
    println("Data for L=$L saved on file!")
end

# ------------------------------ Horizontal sweep ------------------------------

function HorizontalSweep(L::Int64,
                         nmax::Int64,
                         JJ::Array{Float64},
                         μ0::Float64,
                         DMRGParameters::Vector{Any},                         
                         FilePathOut::String)
    """
    Run horizontal sweep (fixed μ0) to extract the Green's function behaviour
    at increasing J.
    """
    
    DataFile = open(FilePathOut, "a")
    println("Performing horizontal sweep at μ=$μ0...")
    
    l = length(JJ)
    for (j, J) in enumerate(JJ)
        println("Running DMRG for J=$(round(J, digits=3)), μ=$μ0 (simulation ",
                "$j/$l for L=$L)")
        
        # Use @sync to wait for all tasks to complete
        @sync begin
            # Create tasks for each DMRG call
            E, Γ, eΓ =  RunDMRGAlgorithm([L, L, nmax, J, μ0], 
                                 		 DMRGParameters,
                                 		 "Correlator";
							             FixedN=false,
							             RandomPsi0=false)

            # Write results to the file
            write(DataFile, "$L; $J; $E; $Γ; $eΓ;\n")
        end
    end
    
    close(DataFile)
    println("Data for L=$L saved on file!")
end


# ----------------------------- Rectangular sweep ------------------------------

function RectangularSweep(L::Int64,
    					  N::Int64,
    					  nmax::Int64,
    					  JJ::Vector{Float64},				# Horizontal domain
						  μμ::Vector{Float64},				# Vertical domain
    					  DMRGParametersMI::Vector{Any},
    					  DMRGParametersSF::Vector{Any},
    					  FilePathIn::String,				# To evaluate if a given point is MI or SF
    					  FilePathOut::String)				# Save computation time if correlators are not needed
    """
    Calculate the variance of the number of particles on site i and ⟨b_i⟩ for a
    range of hopping J and chemical potential μ values. Results are saved to a file.
    """
    
    DataFile = open(FilePathOut,"a")

	# Take data from fitting data of horizontal sweeps to separate MI from SF 
	# ( use ΔE^+(∞) and ΔE^-(∞) )
	BoundariesData = readdlm(FilePathIn, ',', Float64, '\n'; comments=true)
    JJFitted = BoundariesData[:,1]

    for (j,J) in enumerate(JJ)
		# We take the best approximating J ∈ JJ, to assess whether we are in the MI or SF phase

		# Index = findall(==(J), BoundariesData[:,1]) # this would work if there is the EXACT J in the fit results
		Index = argmin(abs.(JJFitted .- J)) # this always works, gives the best approximation
		μUp = BoundariesData[Index,2][1]
		μDown = -BoundariesData[Index,3][1]
    	
		println("\nJ=$J, phase boundaries: μ^+=$μUp, μ^-=$μDown")

		CachedRho = 1
		
        for (m,μ) in enumerate(μμ)
        
            ModelParameters = [L, N, nmax, J, μ]
			inMottLobe = false
			
			if (μ>=μDown && μ<=μUp)
				inMottLobe = true
			end

			println("Running DMRG for L=$L, J=$(round(J, digits=3)), ",
			"μ=$(round(μ, digits=3)) (simulation $m/$(length(μμ)) in μ, ",
			"$j/$(length(JJ)) in J, MI=$inMottLobe)")

			if inMottLobe
		        Results = RunDMRGAlgorithm(ModelParameters,
										   DMRGParametersMI,
										   "OrderParameters";
		                                   FixedN=false,
										   RandomPsi0=false)
				E, NTotAvg, nVariance, aAvg = Results
				
				if m>1
					Rho = NTotAvg/L
					K = (Rho - CachedRho) / (μ - μμ[m-1])	# Compressibility
					K /= (Rho^2)
					CachedRho = Rho							# Next step
		        elseif m==1
		        	K = NaN									# Avoid segmentation fault
		        end
		        
		        write(DataFile,"$J; $μ; $E; $nVariance; $aAvg; $K # MI\n")
		    else
		    	Results = RunDMRGAlgorithm(ModelParameters,
				                           DMRGParametersSF,
				                           "OrderParameters";
	  		                               FixedN=false,
										   RandomPsi0=true)
				E, NTotAvg, nVariance, aAvg = Results
				
				if m>1
					Rho = NTotAvg/L
					K = (Rho - CachedRho) / (μ - μμ[m-1])	# Compressibility
					K /= (Rho^2)
					CachedRho = Rho							# Next step
		        elseif m==1
		        	K = NaN									# Avoid segmentation fault
		        end
				
		        write(DataFile,"$J; $μ; $E; $nVariance; $aAvg; $K # SF\n")
		    end
        end
    end

    close(DataFile)
end

# -------------------------------- Path sweep ----------------------------------

# UNDER CONSTRUCTION!

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

