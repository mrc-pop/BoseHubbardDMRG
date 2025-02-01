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

function CalculatePhaseBoundaries(L::Int64, N::Int64, nmax::Int64, 
    JJ::Array{Float64}, μ0::Float64, DMRGParameters, FilePathOut)
    """
    Calculate the phase boundaries between the Mott Insulator (MI) and 
    Superfluid (SF) phases. The results are saved to a file.
    """

    DataFile = open(FilePathOut,"w")
    write(DataFile,"# Hubbard model DMRG. L=$L, N=$N, nmax=$nmax, μ0=$μ0.\n")
    write(DataFile,"# J, DeltaE_g^+, DeltaE_g^- [calculated $(now())]\n")

    E = zeros(Float64, 3)

    for J in JJ
        E[1], _, _, _, _, _ = RunDMRGAlgorithm([L, N-1, nmax, J, μ0], 
                                                DMRGParameters; FixedN = true)
        E[2], _, _, _, _, _ = RunDMRGAlgorithm([L, N, nmax, J, μ0], 
                                                DMRGParameters; FixedN = true)
        E[3], _, _, _, _, _ = RunDMRGAlgorithm([L, N+1, nmax, J, μ0], 
                                                DMRGParameters; FixedN = true)
        ΔEplus = E[3] - E[2]
        ΔEminus = E[2] - E[1]
        write(DataFile,"$J, $ΔEplus, $ΔEminus\n")
    end

    close(DataFile)
    println("Data saved on file!")
end

function main()
    L = 10
    N = 10
    nmax = 3
    JJ = collect(range(0.0,0.3,5))
    μμ = collect(range(0.0,1,5))
    i = ceil(Int64, L/2) # site on which to calculate variance
    μ0 = 0.0

    nsweep = 5
    maxlinkdim = [10,50,75,200,500]
    cutoff = [1E-13]
    DMRGParameters = [nsweep, maxlinkdim, cutoff]

    CalculateVariance(JJ, μμ, i, L, N, nmax, DMRGParameters,
        "../simulations/data_variance.txt")

    CalculatePhaseBoundaries(L, N, nmax, JJ, μ0, DMRGParameters, 
         "../simulations/phaseboundaries.txt")
end

main()