#/usr/bin/julia

# Include the dmrg.jl file
include("./dmrg.jl")

function PlotVariance(FilePath, L, nmax, i)
    """
    Plot the variance of the number of particles on site `i` from data saved
    in `FilePath`.
    """
    VarianceData = readdlm(FilePath, ',', Float64, '\n'; comments=true)

    JJ = VarianceData[:,1]
    μμ = VarianceData[:,2]
    varvar = VarianceData[:,3]

    NumJ = length(unique(JJ))
    Numμ = length(unique(μμ))

    vars = zeros(Numμ, NumJ)

    for jj in 1:NumJ
        vars[:,jj] = varvar[ Numμ*(jj-1)+1 : Numμ*jj ]
    end

    heatmap(unique(JJ), unique(μμ), vars, 
        xlabel=L"J", ylabel=L"μ", 
        title=L"Variance on $n_i$ ($L=%$L, n_\mathrm{max}=%$nmax, i=%$i$)", 
        size=(600, 400))
    println("Variance plot saved on file!")
end

function PlotPhaseBoundaries(FilePath, μ0; gap=false)
    """
    Plot the phase boundaries between the Mott Insulator (MI) and Superfluid 
    (SF) phases. If gap=true, plot the charge gap instead of the phase 
    boundaries.
    """
    BoundariesData = readdlm(FilePath, ',', Float64, '\n'; comments=true)

    JJ = BoundariesData[:,1]
    ΔEplus = BoundariesData[:,2]
    ΔEminus = BoundariesData[:,3]

    if gap
        plot!(JJ, ΔEplus - ΔEminus, 
            size=(600, 400),
            label=false,
            xlabel=L"$J$", ylabel=L"$\Delta E_{\mathrm{gap}}$", 
            title=L"Charge gap  as a function of $J$ ($\mu_0=%$μ0$)")
    else
        plot(JJ, μ0 .+ ΔEminus , # TODO add ! if superposed=true
        xlabel=L"$J$", ylabel=L"$\mu$",
        size=(600, 400), 
        label=L"$\mu_c^-$",
        title=L"$\mu_c^\pm(J)$ as a function of $J$ ($\mu_0=%$μ0$)")

        plot!(JJ, μ0 .+ ΔEplus, label=L"\mu_c^+")
    end
    println("Phase boundaries plot saved on file!")
end


function main()    
    L = 10
    N = 10
    nmax = 3
    JJ = collect(range(0.0,0.3,5))
    μμ = collect(range(0.0,1,5))
    i = ceil(Int64, L/2) # site on which to calculate variance
    μ0 = 0.0

    # ---------------------
    # --- Plot Variance ---
    # ---------------------
    PlotVariance("../simulations/data_variance.txt", L, nmax, i)
    savefig("../analysis/variance.pdf")

    # ----------------------------------
    # --- Boundary between SF and MI ---
    # ----------------------------------
    PlotPhaseBoundaries("../simulations/phaseboundaries.txt", μ0)
    savefig("../analysis/phaseboundaries.pdf")
end

main()