using ITensors, ITensorMPS
using Plots; pgfplotsx()
using LaTeXStrings
using Dates

function GetHamiltonian(sites, J::Float64, U::Float64, μ::Float64)
    """
    Calculate the 1d Bose-Hubbard hamiltonian with open boundary conditions
    as a MPO.
    """
    os = OpSum()
    N = length(sites)
    for j=1:(N-1)
        os += -J,"Adag",j,"A",j+1
        os += -J,"Adag",j+1,"A",j
        os += 0.5*U,"N",j,"N",j
        os += (-0.5*U-μ),"N",j
    end
    os += 0.5*U,"N",N,"N",N
    os += (-0.5*U-μ),"N",N

    # os += -J,"Adag",N,"A",1 # PERIODIC BC
    # os += -J,"Adag",1,"A",N # PERIODIC BC

    return MPO(os,sites)
end

function GetLocalH(sites, i, J::Float64, U::Float64, μ::Float64)
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

function GetNumber(sites)
    """
    Calculate the total number operator N as a MPO.
    """
    os = OpSum()
    for j=1:length(sites)
        os += "N",j
    end
    return MPO(os, sites)
end

function GetTwoPointCorrelator(sites,i::Int64,j::Int64)
    os = OpSum()
    os += "Adag",i,"A",j
    return MPO(os,sites)
end

function GetNumberVariance(psi, sites, i::Int64)
    """
    Calculate the variance of the number of particles on site i, on state psi.
    """
    if i < 0 | i > length(sites)
        error("Site index out of bounds! Use a proper index.")
        return
    end

    n = OpSum()
    n += "N",i
    n = MPO(n, sites)

    return inner(n, psi, n, psi) - inner(psi', n, psi)^2
end

function SetStartingState(sites, N, d)
    """
    Create a MPS on sites (system of size L), made of N particles, 
    and where the local Hilbert space has dimension d.
    """
    L = length(sites)
    NumFullSites = floor(Int64, N/(d-1))
    k = NumFullSites + 1

    states = ["0" for i in 1:L] # empty sites

    states[1] = string(floor(Int64, N%(d-1))) # remainder

    for i in 2:k # full sites
        states[i] = string(d-1)
    end

    circshift!(states, floor(Int64,N/2) - k + floor(Int64,k/2))

    println(states)

    return MPS(sites, "1") # TODO change
end

function RunDMRGAlgorithm(ModelParameters, DMRGParameters; 
    verbose=false, U=1.0, calculate_gamma=false)
    """
    Run DMRG algorithm with chosen parameters, and return results.
    Truncate the local Hilbert space with a maximum occupation of nmax.
    The on-site interaction U is by default fixed to 1.
    Input:
        - ModelParameters: array of L::Int64, N::Int64, nmax::Int64, J::Float64, 
            μ::Float64
        - DMRGParameters: array of nsweeps, maxdim, cutoff
    Output:
        - Results of DMRG optimization and relevant observables
    """
    L, N, nmax = Int64.(ModelParameters[1:3])
    J, μ = ModelParameters[4:5]
    nsweeps, maxdim = DMRGParameters[1:2]
    cutoff = DMRGParameters[3]

    if verbose
        println("Model parameters: L=$L, N=$N, nmax=$nmax, J=$J, U=$U, μ=$μ")
        println("Starting simulation...\n")
    end

    # Calculate hamiltonian and number operators
    d = nmax + 1
    sites = siteinds("Boson", L, dim=d; conserve_number=true)
    H = GetHamiltonian(sites, J, U, μ)
    Ntot = GetNumber(sites)
    
    # Set starting state and print observables
    psi0 = SetStartingState(sites, N, d)
    if verbose
        N0 = inner(psi0', Ntot, psi0)
        E0 = inner(psi0', H, psi0)
        println("Expectation values on the initial state: N=$N0, E=$E0\n")
    end

    # Run DMRG algorithm and print results
    E, psi = dmrg(H,psi0;nsweeps,maxdim,cutoff,outputlevel=verbose)

    # Sanity checks: calculate whether Ntot has been conserved, and 
    # the found ground state is actually an eigenstate of H.
    if verbose
        VarE = inner(H,psi,H,psi) - E^2
        NtotAvg = inner(psi', Ntot, psi)
        println("\nFinal N: $NtotAvg, variance on E: $VarE\n")
    end

    # Calculate relevant observables in the ground state.
    aAvg = zeros(L) # < a_j > (possible order parameter for SF)
    nMean = zeros(L)
    nVariance = zeros(L) # variance on n_i, for all sites i
    LocalE = zeros(L) # "local part" of contribution to the energy

    for i in 1:L
        aAvg[i] = expect(psi, "a"; sites=i)
        nMean[i] = expect(psi, "n"; sites=i)
        nVariance[i] = GetNumberVariance(psi, sites, i)
        LocalE[i] = inner(psi', GetLocalH(sites, i, J, U, μ), psi)
    end

    iCenter = floor(L/2)
    Γ = zeros(floor(Int64, iCenter/2)) # 2 point correlator Γ(r), r even, r < L/2
    if calculate_gamma
        for j in 1:floor(Int64,iCenter/2)
            Γop = GetTwoPointCorrelator(sites, iCenter-j, iCenter+j)
            Γ[j] = inner(psi', Γop, psi)
        end
    end

    return E, aAvg, nMean, nVariance, LocalE, Γ
end

function CalculatePlotVariance(JJ::Array{Float64}, μμ::Array{Float64}, i::Int64,
    L::Int64, N::Int64, nmax::Int64, DMRGParameters)
    """
    Given the arrays JJ and μμ, and the index site i, calculate, save on file 
    and plot the variance.
    """
    FilePath = "./data/data_variance.txt"
    DataFile = open(FilePath,"w")
    write(DataFile,"# Hubbard model DMRG. L=$L, N=$N, nmax=$nmax, i=$i.\n")
    write(DataFile,"# J, mu, n_variance [calculated $(now())]")

    vars = zeros(length(μμ), length(JJ))

    for (i,J) in enumerate(JJ)
        for (j,μ) in enumerate(μμ)
            ModelParameters = [L, N, nmax, J, μ]
            _, _, _, nVariance, _, _ = RunDMRGAlgorithm(ModelParameters,
                                            DMRGParameters)
            vars[j,i] = nVariance[i]
            write(DataFile,"$J, $μ, $(vars[j,i])\n")
        end
    end

    close(DataFile)

    heatmap(JJ, μμ, vars, xlabel=L"J", ylabel=L"μ", 
        title=L"Variance on $n_i$ ($L=%$L, n_\mathrm{max}=%$nmax, i=%$i$)", 
        size=(600, 400))
    
    savefig("./results/variance.pdf")
end

function main()
    # Run DMRG for different algorithm parameters, to check
    # how many sweeps are sufficient for convergence.
    N = 14
    L = 14
    nmax = 3
    J = 0.3
    μ = 0.4

    nsweep = 5
    maxlinkdim = [10,50,75,100,100]
    cutoff = [1E-10]
    DMRGParameters = [nsweep, maxlinkdim, cutoff]

    E, aAvg, nMean, nVariance, LocalE, Γ = RunDMRGAlgorithm([N, L, nmax, J, μ], DMRGParameters; verbose=true)

    println("\"Local\" part of energy: $(round.(LocalE, digits=4))\n")
    println("Mean number of particles: $(round.(nMean, digits=4))\n")
    println("Variance number of particles: $(round.(nVariance, digits=4))\n")
    println("Relative fluctuation: $(round.(sqrt.(nVariance)./nMean, digits=4))\n")

    # ---------------------
    # --- Plot Variance ---
    # ---------------------
    # # Calculate and plot variance on n_i, with i=3
    #JJ = collect(range(0.0,0.3,6))
    #mumu = collect(range(0.0,1,6))
    #CalculatePlotVariance(JJ, mumu, 3, L, N, nmax, DMRGParameters)
end

main()
