#!/usr/bin/julia

PROJECT_ROOT = @__DIR__ # Absloute path up to .../BoseHubbardDMRG/src

# Include the dmrg.jl file
include(PROJECT_ROOT * "/dmrg.jl")

PROJECT_ROOT *= "/.."	# Absloute path up to .../BoseHubbardDMRG/

# TODO Add """[overall introduction]"""

# -----------------
# --- Functions ---
# -----------------

function HorizontalSweep(LL::Array{Int64},
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
    
    # TODO Generalize this function for N != L

    DataFile = open(FilePathOut,"w")
    write(DataFile,"# Hubbard model DMRG. nmax=$nmax, μ0=$μ0.\n")
    write(DataFile,"# L; J; E; deltaE_g^+; deltaE_g^-; nVariance; Γ [calculated $(now())]\n")

    E = zeros(Float64, 3)

    for L in LL
        println("Starting calculation of observables for L=$L...")

        for (j, J) in enumerate(JJ)
            println("Running DMRG for J=$(round(J,digits=3)), μ=$μ0 (simulation ",
            		"$j/$(length(JJ)) for L=$L)")
            
           	# Different filling simulations; ComputeAllObservables=false by default
            E[1], _, _, = RunDMRGAlgorithm([L, L-1, nmax, J, μ0], 
                                    	   DMRGParameters;
                                    	   FixedN = true)
            E[2], nVariance, Γ, = RunDMRGAlgorithm([L, L, nmax, J, μ0], 
                                            	   DMRGParameters;
				                            	   ComputeGamma=true, 
						                           FixedN = true)
            E[3], _, _, = RunDMRGAlgorithm([L, L+1, nmax, J, μ0], 
                                    	   DMRGParameters; 
		                                   FixedN = true)
            ΔEplus = E[3] - E[2]
            ΔEminus = E[2] - E[1]

            μplus = ΔEplus + μ0
            μminus = -ΔEminus + μ0

            write(DataFile,"$L; $J; $(E[2]); $μplus; $μminus; $nVariance; $Γ\n")
        end
    end
    close(DataFile)
    println("Data saved on file!")
end

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
    
    DataFile = open(FilePathOut,"w")
    write(DataFile,"# Hubbard model DMRG. L=$L, N=$N, nmax=$nmax\n")
    write(DataFile,"# NOTE: Different DMRG settings have been used for MI and SF. Check the code.\n")
    write(DataFile,"# J, μ, E, n_variance [calculated $(now()) @site  i=$i]\n")

	# TODO Generate file from fitting data of horizontal sweeps to separate MI from SF
	Sizes, Couplings, Energies, UpBoundaries, DownBoundaries, _, _ = readdlm(FilePathIn, ';', Float64, '\n'; comments=true)
	IndicesList = findall(==(L), Sizes)

    for J in JJ
    	
    	# Possible improvement: we know how many simulations have been performed for each L
    	Index = IndicesList[findall(==(J), Couplings[IndicesList])]
    	μUp = UpBoundaries[Index]
    	μDown = DownBoundaries[Index]
    	
        for μ in μμ
        
            ModelParameters = [L, N, nmax, J, μ]
			inMottLobe = false
			
			if (μ>=μDown && μ<=μUp)
				inMottLobe=true
			end
			
			if inMottLobe
	            E, nVariance, _ = RunDMRGAlgorithm(ModelParameters,
	                                               DMRGParametersMI;
	                                               FixedN = false)
	            write(DataFile,"$J, $μ, $E, $(nVariance[i]) # MI\n")
	        else
	        	E, nVariance, _ = RunDMRGAlgorithm(ModelParameters,
	                                               DMRGParametersSF;
  		                                           FixedN = false)
	            write(DataFile,"$J, $μ, $E, $(nVariance[i]) # SF\n")
	        end
        end
    end

    close(DataFile)
    println("Data saved on file!")
end

function main()

    # TODO Import DMRG parameters from user input
    nmax = 3

    # Mott Insulator DMRG parameters
    nsweeps = 5
    maxlinkdim = [10,50,75,200,500]		# TODO Change with optimal values
    cutoff = [1E-12]					# TODO Change with optimal values
    DMRGParametersMI = [nsweeps, maxlinkdim, cutoff]

	# Superfluid phase DMRG paramters
	nsweeps = 5
	maxlinkdim = [10,50,75,200,500]		# TODO Change with optimal values
	cutoff = [1E-12]					# TODO Change with optimal values
	DMRGParametersSF = [nsweeps, maxlinkdim, cutoff]
    
	ModeErrorMsg = "Input error: use option --horizontal or --rectangular"
	
	if length(ARGS) != 1
		# If user does not specify the user mode
		error(ModeErrorMsg)
		exit()
	else
		
		UserMode = ARGS[1]
		if UserMode=="--horizontal"

			# Horizontal sweep
	    	# TODO Import model parameters from user input
	    	JJ = collect(range(start=0.0, stop=0.35, length=100))
	    	LL = [10, 20, 30]
	    
	    	FilePathOut = PROJECT_ROOT * "/simulations/horizontal_sweep.txt"
	    	# TODO Evaluate: use superfluid DMRG paramters?
	    	HorizontalSweep(LL, nmax, JJ, DMRGParametersSF, FilePathOut)
			

		elseif UserMode=="--rectangular"

			# Rectangular sweep
		    N = 10 								# TODO Change
		    L = 10 								# TODO Change
		    JJ = [J for J in 0.0:0.03:0.3]		# TODO Change, exclude J=0
		    μμ = [μ for μ in 0.0:0.1:1.0]		# TODO Change

			FilePathIn =  PROJECT_ROOT * "/simulations/horizontal_sweep.txt"
		    # TODO Cycle over sites, extend to different fillings
		    i = ceil(Int64, L/2) 				# Site to calculate variance on
		    mkpath(PROJECT_ROOT * "/simulations/L=$L")
			FilePathOut = PROJECT_ROOT * "/simulations/L=$L/rectangular_sweep_site=$i.txt"
		    RectangularSweep(i, L, N, nmax, JJ, μμ, DMRGParametersMI, DMRGParametersSF, FilePathIn, FilePathOut)
		
		else
			error(ModeErrorMsg)
			exit()
		end
		
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end