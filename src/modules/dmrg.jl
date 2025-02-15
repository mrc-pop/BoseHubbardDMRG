#!/usr/bin/julia

using ITensors, ITensorMPS
using Statistics
using Dates
using DelimitedFiles

# ------------------------------------------------------------------------------
# --------------------------------- Physics ------------------------------------
# ------------------------------------------------------------------------------

# Full hamiltonian

function GetHamiltonian(sites, J::Float64, U::Float64, μ::Float64; pbc=false)
    """
    Construct the 1D Bose-Hubbard Hamiltonian as an MPO for given sites,
    hopping `J`, interaction `U`, and chemical potential `μ`.
    If pbc=true, add periodic BC, otherwise open BC. (default: false)
    """
    os = OpSum()
    L = length(sites)
    for j=1:L
        if j < L
            os += -J,"Adag",j,"A",j+1
            os += -J,"Adag",j+1,"A",j
        end
        os += 0.5*U,"N",j,"N",j
        os += -0.5*U,"N",j
        os += -μ,"N",j
    end

    if pbc
        os += -J,"Adag",L,"A",1
        os += -J,"Adag",1,"A",L
    end
    return MPO(os,sites)
end

# Local hamiltonian
# TODO Evaluate if full hamiltonian can be obtained by calling GetLocalH

function GetLocalH(sites, i, J::Float64, U::Float64, μ::Float64)
    """
    Construct the local Hamiltonian term MPO for site `i` in the 1D Bose-Hubbard
    model with hopping `J`, interaction `U`, and chemical potential `μ`.
    """
    os = OpSum()
    if i > 1
        os += -J,"Adag",i-1,"A",i
        os += -J,"Adag",i,"A",i-1
    end
    if i < length(sites)
        os += -J,"Adag",i,"A",i+1
        os += -J,"Adag",i+1,"A",i
    end
    os += 0.5*U,"N",i,"N",i
    os += (-0.5*U-μ),"N",i
    return MPO(os,sites)
end

# Total number operator

function GetNumber(sites)
    """
    Calculate the total number operator `N` as a MPO for the given sites.
    """
    os = OpSum()
    for j=1:length(sites)
        os += "N",j
    end

    return MPO(os, sites)
end

# Two points correlator

function GetTwoPointCorrelator(sites, i::Int64, j::Int64)
    
    """
    Construct the two-point correlator operator `a_i^+ a_j` as an MPO
    for sites `i` and `j`.
    """
    
    os = OpSum()
    os += "Adag",i,"A",j
    return MPO(os,sites)
end

# Density-density correlator

function GetNumberCorrelator(psi, sites, i::Int64, j::Int64)
    
    """
    Extract the full density-density correlator beween sites i, j. The function
    does not return a MPO, but a scalar number.
    """
    
    if (i < 0) | (i > length(sites))
        error("First site index out of bounds! Use a proper index.")
        return
    elseif (j < 0) | (j > length(sites))
        error("Second site index out of bounds! Use a proper index.")
        return
    end
    
    os = 0
    os = OpSum()
    os += "N",i
    hatNi = MPO(os, sites)
    
    os = 0
    os = OpSum()
    os += "N",j
    hatNj = MPO(os, sites)
    
    NiNj = inner(hatNi, psi, hatNj, psi)
    Ni = inner(psi', hatNi, psi)
    Nj = inner(psi', hatNj, psi)
    
    return NiNj - (Ni*Nj)
end

# Number variance operator per site

function GetNumberVariance(psi, sites, i::Int64)
    """
    Calculate the variance of the number of particles on site `i` for the state
    `psi`.
    """
    if (i < 0) | (i > length(sites))
        error("Site index out of bounds! Use a proper index.")
        return
    end

    os = OpSum()
    os += "N",i
    hatNi = MPO(os, sites)

    Ni2 = inner(hatNi, psi, hatNi, psi)
    Ni1 = inner(psi', hatNi, psi)

    return Ni2 - (Ni1^2)
end

# Von Neumann entropy

function GetVonNeumannEntropy(psi::MPS, sites, i::Int64)
    """
    Calculate Von Neumann entropy of the bipartition `1, ..., i` and `i+1, ..., L`,
    tracing `psi` on all other sites,
    """
    # Perform SVD to all states except i, to prepare for calculating a local observable (end of Lect 4 page 4)
    orthogonalize!(psi, i)

    # First argument: the tensor on which to perform the SVD, i.e. psi[i]
    # Second argument: the indices on which to perform the SVD, i.e. linkind(psi, i-1) and sites[i] (left link and vertical link)
    _, S = svd(psi[i], (linkind(psi, i-1), sites[i]))

    # S is the diagonal matrix containing the singular values
    SvN = 0.0 # von Neumann entropy
    for n in 1:dim(S, 1)
      p = S[n,n]^2
      SvN -= p * log(p)
    end
    return SvN
end

# ------------------------------------------------------------------------------
# ------------------------------------ DMRG ------------------------------------ 
# ------------------------------------------------------------------------------

function SetStartingState(sites, N, d)
    """
    Create an initial MPS with `N` particles and local Hilbert space dimension
    `d` on the given sites.
    """
    L = length(sites)
    NumFullSites = floor(Int64, N/(d-1))
    k = NumFullSites + 1

    # For generic N, set |0,...,0,r,d-1,...,d-1,0,...,0>, where r=remainder
    states = ["0" for i in 1:L] # empty sites
    states[1] = string(floor(Int64, N%(d-1))) # remainder
    for i in 2:k # full sites
        states[i] = string(d-1)
    end
    circshift!(states, floor(Int64,N/2) - k + floor(Int64,k/2))

    # If N=L-1, L or L+1, set a smarmelled initial state
    if N == L
        states = ["1" for _ in 1:L]
    elseif N == (L-1)
        states = ["1" for _ in 1:L]
        states[ceil(Int64, L/2)] = "0"
    elseif N == (L+1)
        states = ["1" for _ in 1:L]
        states[ceil(Int64, L/2)] = "2"
    end

    # For debugging only
    # println("Starting state (N=$N) ", states)

    return MPS(sites, states)
end

function RunDMRGAlgorithm(ModelParameters::Vector{Float64},
						  DMRGParameters::Vector{Any},
						  UserMode::String;
						  verbose=false,				# Do not print at line
						  U=1.0,						# Simplify model
						  FixedN=false,					# Let N vary
						  pbc=false,                    # No periodic boundary conditions
                          RandomPsi0=true)				# Initialize a random state	
    
    """
    Run DMRG algorithm with chosen parameters, and return results.
    Truncate the local Hilbert space with a maximum occupation of nmax.
    The on-site interaction U is by default fixed to 1.
    Input:
        - ModelParameters: array of [L::Int64, N::Int64, nmax::Int64,
        							 J::Float64, μ::Float64]
        - DMRGParameters: array of [nsweeps::Int64, maxdim::Int64, 
        							cutoff::Vector{Float64}]
    Parametric input:
        - ComputeAllObservables: boolean variable, if false only E, 
          nVariance and Γ are computed (default: false)
        - ComputeGamma: boolean variable, if false all positional correalators
          are not extracted (default: false)
        - ComputeC: boolean variable, if false all positional correalators are
          not extracted (default: false)
    Output:
        - Results of DMRG optimization and relevant observables
    """
    
    L, N, nmax = Int64.(ModelParameters[1:3])
    J = ModelParameters[4]
    μ = ModelParameters[5]
    
    nsweeps, maxdim = DMRGParameters[1:2]
    cutoff = DMRGParameters[3]

	# Evaluate UserMode
	ModeErrorMsg = "Input error: use as argument \"OrderParameters\"
- \"Correlator\", \"Debug\" or \"Fast\":
- \"OrderParameters\": return E, nVariance, aAvg;
- \"Correlator\": return E, Γ, eΓ;
- \"Debug\": return E, LocalE, nMean, nVariance;
- \"Fast\": return E."
	
	OrderParameters=false
	Correlator=false
	Debug=false
	Fast=false
	
	if UserMode=="OrderParameters"
		OrderParameters=true
		nVariance = 0		 		# Variance on central site
   		aAvg = 0		 			# <a> on central site
	
	elseif UserMode=="Correlator"
		Correlator=true
		# Initialize later
		
	elseif UserMode=="Debug"
		Debug=true
		nMean = zeros(L)			# Mean number of particles per site
		nVariance = zeros(L) 		# Variance on n_i, for all sites i
    	LocalE = zeros(L)			# Local contribution to the energy
    
    elseif UserMode == "Fast"
    	Fast = true
    		
	else
		error(ModeErrorMsg)
	end

    if verbose
    	if !Fast
    		println("Selected mode: ", UserMode)
        end
        println("Model parameters: L=$L, N=$N, nmax=$nmax, J=$J, U=$U, μ=$μ.")
        println("Starting simulation...\n")
    end

    # Calculate hamiltonian and number operators
    d = nmax + 1
    sites = siteinds("Boson", L, dim=d; conserve_number=FixedN)
    H = GetHamiltonian(sites, J, U, μ; pbc)
    Ntot = GetNumber(sites)
    
    # Set starting state and print observables
    if RandomPsi0
        psi0 = randomMPS(sites)					# Initialize random state
    else
        psi0 = SetStartingState(sites, N, d)	# Initialize Mott state
    end

    if verbose
        N0 = inner(psi0', Ntot, psi0)
        E0 = inner(psi0', H, psi0)
        println("Expectation values on the initial state: N=$N0, E=$E0\n")
    end

    # Run DMRG algorithm and print results
    E, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff, outputlevel=verbose)

    # Sanity checks: calculate whether Ntot has been conserved, and 
    # the found ground state is actually an eigenstate of H.
    if verbose
        VarE = inner(H, psi, H, psi) - E^2
        NtotAvg = inner(psi', Ntot, psi)
        println("\nFinal N=$(round(NtotAvg, digits=3)), ", 
        		"E=$(round(E,digits=4)), ",
        		"VarE=$(round(VarE,digits=4))\n")
    end

    if OrderParameters
    
    	Index = ceil(Int64, L/2)
		nVariance = GetNumberVariance(psi, sites, Index)
	    aAvg = expect(psi, "a"; sites=Index)
		
		return E, nVariance, aAvg
    
    elseif Correlator

		"""
		Procedure: discard half of the lattice (the first quarter and the last);
		what remains is used as a domain to sweep across, collecting 
		symmetrically two-points and density-density correlators; these last
		are then avereaged; error is taken as the statistical standard 
		deviation.
		"""
		
		Trash = floor(Int64, L/4)	# How many to discard from each side?
		Start = Trash+1				# Start of useful segment
		Stop = L-Trash				# End of useful segment
		
		Segment = L-2*Trash			# Segment length
		Spacings = collect(Int64, 2:2:Segment)

		Γ = zeros(Float64, length(Spacings))
		eΓ = zeros(Float64, length(Spacings))
		
		# (Center) symmetric sweep for Γ
		for r in Spacings

		    ΓTmpArray = [] # stores the different Γ(r)

		    i = Start
		    j = Start+r-1

		    while j<=Stop
		        ΓOp = GetTwoPointCorrelator(sites, i, j)
		        ΓTmp = inner(psi', ΓOp, psi)
		        push!(ΓTmpArray, ΓTmp)	

		        i += 1
		        j += 1   
		    end

		    Γ[Int64(r/2)] = mean(ΓTmpArray)
		    eΓ[Int64(r/2)] = std(ΓTmpArray)
        end
        
        eΓ[end] = eΓ[end-1] # Correct NaN
		
        return E, Γ, eΓ
    
    elseif Debug
	    
	    for i in 1:L
			nMean[i] = expect(psi, "n"; sites=i)
			nVariance[i] = GetNumberVariance(psi, sites, i)
			LocalE[i] = inner(psi', GetLocalH(sites, i, J, U, μ), psi)
    	end
    	
    	return E, nMean, nVariance, LocalE
    	
    elseif Fast
    	return E
    end
end

# ------------------------------------------------------------------------------
# ------------------------------------ Main ------------------------------------ 
# ------------------------------------------------------------------------------

function main()

	"""
    If the script is called directly from terminal, run one full DMRG routineù
    with the following parameters.
    """
    
    L = 10
    N = 10
    nmax = 3
    J = 0.4
    μ = 1.0

    nsweep = 10
    maxlinkdim = [10,50,75,200,500]
    cutoff = [1E-8]
    DMRGParameters = [nsweep, maxlinkdim, cutoff]
    
    UserMode = "OrderParameters" # "OrderParameters" / "Correlator" / "Debug" / "Fast"

	if UserMode=="Fast"
		for i in 1:100
    		print("\r")
    		for j in 1:i
    			print(" ")
    		end
    		print(">>> Speedin' >>>")
    		sleep(0.05)
    	end
    	print("\n")
	end
	
    # ------------------------ Results for a single run ------------------------
    
    Observables = RunDMRGAlgorithm([L, N, nmax, J, μ],
    							    DMRGParameters,
    							    UserMode; 
    								verbose=true)
    								
	if UserMode=="OrderParameters"
		E, nVariance, aAvg = Observables
		println("Results of the simulation:
Energy of ground state: $(round.(E, digits=4))
Number variance on central site: $(round.(nVariance, digits=4))
Average of <a> on central site: $(round.(aAvg, digits=4))")
	
	elseif UserMode=="Correlator"
		E, Γ, eΓ = Observables
		println("Results of the simulation:
Energy of ground state: $(round.(E, digits=4))
Green's function: $(round.(Γ, digits=4))
Error on Green's function: $(round.(eΓ, digits=4))")
		
	elseif UserMode=="Debug"
		E, nMean, nVariance, LocalE = Observables
		println("Results of the simulation:
Energy of ground state: $(round.(E, digits=4))
\"Local\" part of energy: $(round.(LocalE, digits=4))
Mean number of particles: $(round.(nMean, digits=4))
Variance number of particles: $(round.(nVariance, digits=4))
Relative fluctuation: $(round.(sqrt.(nVariance)./nMean, digits=4))")
    	
	else
		E = Observables
		println("Results of the simulation:
Energy of ground state: $(round.(E, digits=4))")

		for i in 1:100
    		print("\r")
    		for j in 1:i
    			print(" ")
    		end
    		print(">>> I'm speed >>>")
    		sleep(0.05)
    	end
    	print("\n")
	end

end

if abspath(PROGRAM_FILE) == @__FILE__ # equivalent to if __name__ == "__main__"
	
	"""
	If this file is directly compiled, run a single DMRG simulation by the
    parameters defined in function main().
    """
    
    main()
end
