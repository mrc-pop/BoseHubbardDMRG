#/usr/bin/julia

# Include the dmrg.jl file
include("./dmrg.jl")

function CalculateVariance(JJ::Array{Float64}, μμ::Array{Float64}, i::Int64,
    L::Int64, N::Int64, nmax::Int64, DMRGParameters, FilePathOut)
    """
    Calculate the variance of the number of particles on site `i` for a range of
    hopping `J` and chemical potential `μ` values. Results are saved to a file.
    """
    DataFile = open(FilePathOut,"w")
    write(DataFile,"# Hubbard model DMRG. L=$L, N=$N, nmax=$nmax, i=$i.\n")
    write(DataFile,"# J, mu, n_variance [calculated $(now())]\n")

    # vars = zeros(length(μμ), length(JJ))

    for J in JJ
        for μ in μμ
            ModelParameters = [L, N, nmax, J, μ]
            _, _, _, nVariance, _, _ = RunDMRGAlgorithm(ModelParameters,
                                            DMRGParameters; FixedN = false)
            write(DataFile,"$J, $μ, $(nVariance[i])\n")
        end
    end

    close(DataFile)
    println("Data saved on file!")
end

function CalculatePhaseBoundaries(LL::Array{Int64}, JJ::Array{Float64},
    nmax::Int64, μ0::Float64, DMRGParameters, FilePathOut)
    """
    Calculate the phase boundaries between the Mott Insulator (MI) and 
    Superfluid (SF) phases, for the first lobe (N=L, meaning ρ=1).
    The results are saved to a file.
    """

    DataFile = open(FilePathOut,"w")
    write(DataFile,"# Hubbard model DMRG. nmax=$nmax, μ0=$μ0.\n")
    write(DataFile,"# L, J, DeltaE_g^+, DeltaE_g^- [calculated $(now())]\n")

    E = zeros(Float64, 3)

    for L in LL
        println("Processing L = $L")
        for J in JJ
            E[1], _, _, _, _, _ = RunDMRGAlgorithm([L, L-1, nmax, J, μ0], 
                                                    DMRGParameters; FixedN = true)
            E[2], _, _, _, _, _ = RunDMRGAlgorithm([L, L, nmax, J, μ0], 
                                                    DMRGParameters; FixedN = true)
            E[3], _, _, _, _, _ = RunDMRGAlgorithm([L, L+1, nmax, J, μ0], 
                                                    DMRGParameters; FixedN = true)
            ΔEplus = E[3] - E[2]
            ΔEminus = E[2] - E[1]
            write(DataFile,"$L, $J, $ΔEplus, $ΔEminus\n")
        end
    end
    close(DataFile)
    println("Data saved on file!")
end

function CalculateCorrelationFunction(LL::Array{Int64}, JJ::Array{Float64},
    nmax::Int64, DMRGParameters, FilePathOut::String; pbc=false)
    """
    Calculate Γ(r) for various L and J. Save data on file.
    """

    μ = 1 # SF phase

    DataFile = open(FilePathOut,"w")
    write(DataFile,"# Hubbard model DMRG. nmax=$nmax, μ=$μ (SF phase).\n")
    write(DataFile,"# J, L, Γ(r) [r even, r<L/2] [calculated $(now())]\n")

    for J in JJ
        for L in LL
            N = L # arbitrary
            println("Calculating correlation functions for L=$L, J=$J.")
            _, _, _, _, _, Γ = RunDMRGAlgorithm([L, N, nmax, J, μ], DMRGParameters; 
            calculate_gamma=true, FixedN=true, pbc=false)
            write(DataFile,"$J; $L; $Γ\n")
        end
    end

    close(DataFile)
    println("Correlation functions saved to: ", FilePathOut)
end

function main()
    # Set "global" parameters
    nmax = 3
    nsweep = 5
    maxlinkdim = [10,50,75,200,500]
    cutoff = [1E-12]
    DMRGParameters = [nsweep, maxlinkdim, cutoff]

    # ----------------------
    # --- Study Variance ---
    # ----------------------
    N = 10
    L = 10
    JJ = collect(0.0:0.03:0.3)
    μμ = collect(range(0.0,1,5))
    i = ceil(Int64, L/2) # site on which to calculate variance
    #CalculateVariance(JJ, μμ, i, L, N, nmax, DMRGParameters,
    #    "../simulations/data_variance.txt")

    # ----------------------------------
    # --- Boundary between SF and MI ---
    # ----------------------------------
    # LL = [10, 15, 20, 25, 30]
    # μ0 = 0.0
    # CalculatePhaseBoundaries(LL, JJ, nmax, μ0, DMRGParameters, 
    #     string(@__DIR__, "/../simulations/phaseboundaries.txt"))

    # ------------------------------
    # --- Correlation function Γ ---
    # ------------------------------
    CalculateCorrelationFunction([10,20,30,40], collect(0.24:0.02:0.36), 3, 
        DMRGParameters, string(@__DIR__, "/../simulations/correlators.txt");
        pbc=false)
end

main()