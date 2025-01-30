using ITensors, ITensorMPS
using Plots; pgfplotsx()
using Dates

function GetHamiltonian(sites, J::Float64, U::Float64, μ::Float64)
    # 
    os = OpSum()
    for j=1:length(sites)-1
        os += -J,"Adag",j,"A",j+1
        os += -J,"Adag",j+1,"A",j
        os += 0.5*U,"N",j,"N",j
        os += (-0.5*U-μ),"N",j
    end
    return MPO(os,sites)
end

function GetNumber(sites)
    os = OpSum()
    for j=1:length(sites)
        os += "N",j
    end
    return MPO(os, sites)
end

function GetNumberVariance(psi, sites, j)
    n = OpSum()
    n += "N",j
    n = MPO(n, sites)

    n2 = OpSum()
    n2 += "N",j,"N",j
    n2 = MPO(n2, sites)

    return inner(psi', n2, psi) - inner(psi', n, psi)^2
end

function SetStartingState(sites, N, d)
    L = length(sites)
    coefficients = zeros(d^L)
    # index = Int64((1-d^L)/(1-d)) + 1 # state |11...1> (Julia counts from 1)
    # vedi note tablet:
    n = d-1 # maximum number of bosons per site
    NumFullSites = floor(N/n)
    index = Int64( (N % n + 1) * d^NumFullSites ) # Julia counts from 1
    coefficients[index] = 1.0

    cutoff = 1E-8
    maxdim = 10
    return MPS(coefficients, sites; cutoff=cutoff, maxdim=maxdim)
end

function SweepVariance(J::Float64, μ::Float64) # TEMPORANEO!!
       # Set model parameters
       L = 6 # number of sites
       N = 6 # number of particles
       NTrunc = 3
       U = 1.0 # possibly fixed
   
       # Set DMRG parameters
       nsweeps = 5 # number of sweeps is 5
       maxdim = [10,20,100,100,200] # gradually increase states kept
       cutoff = [1E-10] # desired truncation error
   
       # ---------------------------
       # --- Start of simulation ---
       # ---------------------------
   
       # Calculate hamiltonian and number operators
       d = NTrunc + 1
       sites = siteinds("Boson", L, dim=d; conserve_number=true)
       H = GetHamiltonian(sites, J, U, μ)
       Ntot = GetNumber(sites)
   
       # Set starting state and print observables
       psi0 = SetStartingState(sites, N, d)
       N0 = inner(psi0', Ntot, psi0)
       E0 = inner(psi0', H, psi0)
   
       # Run DMRG algorithm and print results
       E, psi = dmrg(H,psi0;nsweeps,maxdim,cutoff,outputlevel=0)
       NtotAvg = inner(psi', Ntot, psi)
       aAvg = expect(psi, "a";sites=3) # one possible order parameter for the superfluid transition
   
       j = 3
       nVariance = GetNumberVariance(psi, sites, j)
 
       return nVariance
end

function CalculatePlotVariance()
    FilePath = "./data/data_variance.txt"
    DataFile = open(FilePath,"w")
    write(DataFile,"# J, mu, n_variance (Hubbard model DMRG)")

    JJ = range(0.0, 0.3, 20)
    mumu = range(0.1, 1, 20)
    vars = zeros(length(mumu), length(JJ))

    for (i,J) in enumerate(JJ)
        for (j,μ) in enumerate(mumu)
            vars[j,i] = SweepVariance(J, μ)
            var = SweepVariance(J, μ)
            write(DataFile,"$J, $μ, $var \n")
        end
    end

    close(DataFile)

    #display(JJ)
    #display(mumu)
    #display(vars)

    heatmap(JJ, mumu, vars, xlabel="J", ylabel="μ", 
    title="Variance as a function of J and mu", size=(600, 400))
    savefig("./results/variance.pdf")
end

function main()
    # Set model parameters
    L = 8 # number of sites
    N = 8 # number of particles
    NTrunc = 3
    J = 0.0
    U = 1.0 # possibly fixed
    μ = 0.5

    # Set DMRG parameters
    nsweeps = 5 # number of sweeps is 5
    maxdim = [10,20,100,100,200] # gradually increase states kept
    cutoff = [1E-10] # desired truncation error

    # ---------------------------
    # --- Start of simulation ---
    # ---------------------------

    println("Model parameters: L=$L, N=$N, NTrunc=$NTrunc, J=$J, U=$U, μ=$μ")
    println("Starting simulation...\n")

    # Calculate hamiltonian and number operators
    d = NTrunc + 1
    sites = siteinds("Boson", L, dim=d; conserve_number=true)
    H = GetHamiltonian(sites, J, U, μ)
    Ntot = GetNumber(sites)

    # Set starting state and print observables
    psi0 = SetStartingState(sites, N, d)
    N0 = inner(psi0', Ntot, psi0)
    E0 = inner(psi0', H, psi0)
    println("Expectation values on the initial state: N=$N0, E=$E0")

    # Run DMRG algorithm and print results
    E, psi = dmrg(H,psi0;nsweeps,maxdim,cutoff)
    NtotAvg = inner(psi', Ntot, psi)
    aAvg = expect(psi, "a") # one possible order parameter for the superfluid transition

    # Print variance of number in site j
    j = 3
    nVariance = GetNumberVariance(psi, sites, j)
    println("Final N=$NtotAvg, final E=$E, average a=$aAvg")
    println("nVariance = $nVariance")

    # Calculate and plot Variance
    CalculatePlotVariance()
end

main()
