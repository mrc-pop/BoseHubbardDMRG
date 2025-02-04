#/usr/bin/julia

using LsqFit

# Include the dmrg.jl file
include("./dmrg.jl")
include(string(@__DIR__, "/../src/graphic_setup.jl"))

function PlotVariance(FilePath::String, L, nmax, i)
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

function PlotPhaseBoundaries(FilePathIn::String; 
    FilePathOut="", gap=false, overwrite=true, CustomLL=[], μ0=0.0)
    """
    Plot the phase boundaries between the Mott Insulator (MI) and Superfluid 
    (SF) phases, calculated from CalculateObservables.
    If gap=true, plot the charge gap instead of the phase boundaries.
    If overwrite=true, clears previous plots. 
    If CustomLL specified, plot only those sizes.
    """
    BoundariesData = readdlm(FilePathIn, ';', '\n'; comments=true)

    if overwrite
        plot()
    end

    # Extract unique L values
    if CustomLL==[]
        LL = unique(BoundariesData[:, 1])
    else
        LL = CustomLL
    end

    # Get color scheme
    MyColors = ColorSchemes.tab10

    for (l, L) in enumerate(LL)
        indices = (BoundariesData[:, 1] .== L)
        JJ = BoundariesData[indices,2]
        ΔEplus = BoundariesData[indices,4]
        ΔEminus = BoundariesData[indices,5]
        
        if gap
            plot!(JJ, ΔEplus - ΔEminus, 
                xlabel=L"$J$", ylabel=L"$\Delta E_{\mathrm{gap}}$", 
                title=L"Charge gap  as a function of $J$ ($\mu_0=%$μ0$)",
                seriestype=:scatter,
                markersize=1.5,
                label=L"$L=%$L$",
                color=MyColors[l % length(MyColors)])
        else
            plot!(JJ, μ0 .+ ΔEminus, 
                xlabel=L"$J$", ylabel=L"$\mu$",
                label=L"$L=%$L$",
                title=L"Extrapolation of $\mu_c^\pm(J)$ ($\mu_0=%$μ0$)",
                seriestype=:scatter,
                markersize=1.5,
                color=MyColors[l % length(MyColors)])
            plot!(JJ, μ0 .+ ΔEplus, seriestype=:scatter,
                label="",
                markersize=1.5,
                color=MyColors[l % length(MyColors)])
        end
    end
    if !=(FilePathOut,"")
        savefig(FilePathOut)
        if gap
            println("Gap for L=$(Int.(LL)) plotted to ", FilePathOut)
        else
            println("Phase boundaries for L=$(Int.(LL)) plotted to ", FilePathOut)
        end
    end
end

function FitPhaseBoundaries(FilePathIn::String, FilePathOut::String;
    FilePathPlotOut="", FilePathSinglePlotOut="")
    #makeplot=false, makesinglefitplot=false)
    """
    Fit ΔEplus(L) and ΔEminus(L) as functions of 1/L for each J, and extract 
    their values at L -> ∞. Save the results to a file. 
    If makeplot=true, plots the result.

    The output file will have the following columns:
    # J, ΔEplus(∞), ΔEminus(∞), error_ΔEplus(∞), error_ΔEminus(∞), chi2_plus, chi2_minus
    """

    # Read the input data
    data = readdlm(FilePathIn, ';', '\n'; comments=true)

    # Extract unique J values
    JJ = unique(data[:, 2])

    ΔEplus_fit = []
    ΔEminus_fit = []
    error_ΔEplus = []
    error_ΔEminus = []
    chi2_plus = []
    chi2_minus = []
    
    # Define the fit function
    fit_func(x, p) = p[1] * x .+ p[2]
    p0 = [1.0, 0.0]

    # Perform the fit for each J
    for J in JJ
        # Filter data for the current J
        filter = (data[:, 2] .== J)
        L = data[filter, 1]
        ΔEplus = data[filter, 4]
        ΔEminus = data[filter, 5]

        inv_L = 1.0 ./ L

        # Fit ΔEplus
        fit_plus = curve_fit(fit_func, inv_L, ΔEplus, p0) 
        ΔEplus_inf = fit_plus.param[2]  # Intercept (value at L -> ∞)
        error_plus = stderror(fit_plus)[2]
        chi2p = sum(fit_plus.resid.^2)

        # Fit ΔEminus
        fit_minus = curve_fit(fit_func, inv_L, ΔEminus, p0)
        ΔEminus_inf = fit_minus.param[2]  # Intercept (value at L -> ∞)
        error_minus = stderror(fit_minus)[2]  
        chi2m = sum(fit_minus.resid.^2)

        # Save the results for this J
        push!(ΔEplus_fit, ΔEplus_inf)
        push!(ΔEminus_fit, ΔEminus_inf)
        push!(error_ΔEplus, error_plus)
        push!(error_ΔEminus, error_minus)
        push!(chi2_plus, chi2p)
        push!(chi2_minus, chi2m)

        # Get color scheme
        MyColors = ColorSchemes.tab10

        # Plot the fit results for one specific J if makesinglefitplot=true
        if FilePathSinglePlotOut != "" && J == maximum(JJ)
            # Define the best-fit lines
            fit_x = range(minimum(inv_L), maximum(inv_L), length=100)
            fit_y_plus = fit_func(fit_x, fit_plus.param)
            fit_y_minus = fit_func(fit_x, fit_minus.param)

            # Plot the data and best-fit lines
            scatter(inv_L, ΔEplus, label=L"$\Delta E^+$ data", xlabel=L"1/L", 
                ylabel=L"$\Delta E$", title="Fit results for J = $(J)",
                legend=:topleft, markersize=2, color=MyColors[1])
            plot!(fit_x, fit_y_plus, label=L"$\Delta E^+$ fit", color=MyColors[1],
                alpha=0.7)
            scatter!(inv_L, ΔEminus, label=L"$\Delta E^-$ data", markersize=2,
                color=MyColors[2])
            plot!(fit_x, fit_y_minus, label=L"$\Delta E^-$ fit", color=MyColors[2],
                alpha=0.7)

            savefig(FilePathSinglePlotOut)
        end
    end

    # Join the results for saving correctly on file
    results = hcat(JJ, ΔEplus_fit, ΔEminus_fit, error_ΔEplus, error_ΔEminus, 
    chi2_plus, chi2_minus)

    # Save the results
    open(FilePathOut, "w") do file
        write(file, "# J, ΔEplus(∞), ΔEminus(∞), error_ΔEplus(∞), ", 
        "error_ΔEminus(∞), chi2_plus, chi2_minus\n")
        writedlm(file, results, ',')
    end

    println("Fitted results saved to: ", FilePathOut)
    
    if FilePathPlotOut != ""
        # PlotPhaseBoundaries(FilePathIn; gap=false, overwrite=true)
        plot()
        plot!(JJ,
              [ΔEplus_fit ΔEminus_fit],
              label=[L"\mu_c^+ \, (L \rightarrow \infty)" L"\mu_c^- \, (L \rightarrow \infty)"], 
              xlabel=L"J", ylabel=L"$\mu$", 
              title="Fitted phase boundaries",
              alpha=1.0)
        savefig(FilePathPlotOut)
        println("Phase boundaries fit plotted to ", FilePathOut)
    end
end

function FitCorrelationFunction(FilePathIn::String, FilePathOut::String)
    """
    Read data from file. Then fit a power-law extracting the exponent K(J,L). 
    Then fit it against L to extract the thermodynamic limit K_∞(J).
    Finally, plot K_∞ against J.
    """
    # Read the input data
    data = readdlm(FilePathIn, ';', '\n'; comments=true)

    # Extract unique J, L values
    JJ = unique(data[:, 2])
    LL = unique(data[:, 1])

    println("JJ = $JJ")
    println(" LL = $LL")

    # TODO CAMBIA
    # JJ = JJ[50:end]

    # Mastruzzo to extract array of Γ
    function parse_array(str)
        return parse.(Float64, split(strip(str, ['[', ']', ' ']), ','))
    end
    Γall = [parse_array(row[7]) for row in eachrow(data)]
    display(Γall) # (TODO check, but it seems to work)

    # Define the fit function [power-law fit of Γ(r)]
    fit_func(x, p) = p[1]*x.^(-p[2]/2) # TODO disastri con 3 parametri
    p0 = [1.0, 1.0] # Γ(0) = 1, K = boh

    K = [] # K_∞ (thermodynamic limit) as a function of j
    e_K = [] # error on K_∞
    chi2n_K = []

    for J in JJ
        # Initialize {K(J,L)}_L at fixed J, its error and χ^2/ndof
        K_FSS = [] 
        e_K_FSS = []
        chi2n_FSS = []

        for L in LL
            # Filter data for the current J and L
            filter = (data[:, 2] .== J) .& (data[:, 1] .== L)

            # Extract the correlators {Γ(r,J,L)}_r for the selected J,L
            Γeven = Γall[filter] # this is an array, Γ(r even)
            r = range(start=2, step=2, length=length(Γeven)) # r even

            println("Fitting J=$J, L=$L, r=$r.")

            # Perform the fit of Γ(r) vs r, at fixed (J,L).
            # The argument w means weights, they should equal 1/σ^2. 
            # weights = 1.0 ./ e_Γeven.^2 # TODO estimate errors for Γeven!
            # fit = curve_fit(fit_func, r, Γeven, p0) # w=weights 

            # println("Best-fit parameters:", fit.param)

            # K_JL = fit.param[2] # K(J,L)
            # e_K_JL = stderror(fit)[2]
            # chi2n = sum(fit.resid.^2) / (length(r) - length(p0))
            # push!(K_FSS, K_JL)
            # push!(e_K_FSS, e_K_JL)
            # push!(chi2n_FSS, chi2n)

            if J == maximum(JJ) && L == minimum(LL)
                # Make a single plot, hopefully representative of the first
                # series of fits. Save the results.
                scatter(r, Γeven,
                    xlabel=L"$r$",
                    ylabel=L"$\Gamma(r)$",
                    title=L"Correlation function ($L=%$L, J=%$J, \mu=1.0$)",
                    label="DMRG data",
                    legend=:topright)
                #fit_x = range(1.0, maximum(r), length=80)
                #fit_y = fit_func(fit_x, fit.param)
                #K_JL_plot = round(K_JL, digits=2)
                #plot!(fit_x, fit_y, label=L"Best-fit ($K=%$K_JL_plot$)")
                savefig(string(@__DIR__,"/../analysis/correlators_singleplot.pdf"))
            end
        end

        inv_LL = 1 ./ LL

        fit_func_fss(x, p) = p[1] * x .+ p[2]
        p0_FSS = [-1.0, 1.0]

        fit_fss = curve_fit(fit_func_fss, inv_LL, K_FSS, p0_FSS)
        K_J = fit_fss.param[2]
        e_K_J = stderror(fit_fss)[2]
        chi2n_J = sum(fit_fss.resid.^2)/(length(LL)-2)
        push!(K, K_J)
        push!(e_K, e_K_J)
        push!(chi2n_K, chi2n_J)

        if J == maximum(JJ)
            # Make a single FSS plot, hopefully representative of the second
            # series of fits. Save the results.
            fit_x = range(0.0, maximum(inv_LL), length=50)
            fit_y = fit_func_fss(fit_x, fit_fss.param)
            scatter(inv_LL, K_FSS,
                xlabel=L"$1/L$",
                ylabel=L"$K(J,L)$",
                title=L"Exponent $K$ when $J=%$J$",
                label="Fits from DMRG data",
                legend=:topright)
                K_infty_plot = round(K_J, digits=2)
            plot!(fit_x, fit_y, label=L"Best-fit ($K_\infty = %$K_infty_plot$)")
            savefig(string(@__DIR__,"/../analysis/correlators_singleplot_FSS.pdf"))
        end
    end

    # Join the results for saving correctly on file
    results = hcat(JJ, K, e_K, chi2n_K)

    # Save the results
    open(FilePathOut, "w") do file
        write(file, "# J, K_∞, e_K_∞, chi2n\n")
        writedlm(file, results, ',')
    end

    scatter(JJ, K, xlabel=L"$J$", ylabel=L"$K_\infty$", yerr=e_K, 
        title=L"Power-law decay of correlation functions $\Gamma(r)$",
        markersize=2,
        label="Fitted data")
    plot!(JJ, 0.5*ones(length(JJ)), label=L"$K_\mathrm{th} = 1/2$")
    savefig(string(@__DIR__,"/../analysis/correlators_K.pdf"))
end


function main()    
    # L = 10
    # N = 10
    # nmax = 3
    # JJ = collect(range(0.0,0.3,5))
    # μμ = collect(range(0.0,1,5))
    # i = ceil(Int64, L/2) # site on which to calculate variance
    # μ0 = 0.0

    # ---------------------
    # --- Plot Variance ---
    # ---------------------
    #PlotVariance("../simulations/data_variance.txt", L, nmax, i)
    #savefig("../analysis/variance.pdf")

    # ----------------------------------
    # --- Boundary between SF and MI ---
    # ----------------------------------
    FilePathIn = string(@__DIR__, "/../simulations/simulations_240204.txt")

    PhaseBoundariesDir = string(@__DIR__, "/../analysis/phase_boundaries/")
    FilePathPlot = PhaseBoundariesDir*"phaseboundaries_240204.pdf"
    FilePathFit = PhaseBoundariesDir*"fitted_phase_boundaries.txt"

    PlotPhaseBoundaries(FilePathIn; gap=false, FilePathOut = FilePathPlot)

    FitPhaseBoundaries(FilePathIn, FilePathFit; 
        FilePathPlotOut=PhaseBoundariesDir*"phaseboundaries_240204_fit.pdf",
        FilePathSinglePlotOut=PhaseBoundariesDir*"phaseboundaries_240204_fit_single.pdf")
    
    # ------------------------------
    # --- Correlation function Γ ---
    # ------------------------------
    # FilePathIn = string(@__DIR__, "/../simulations/simulations_server_240203.txt")
    # FilePathFit = string(@__DIR__, "/../analysis/correlators_fitted_server.txt")
    # FitCorrelationFunction(FilePathIn, FilePathFit)
end

main()
