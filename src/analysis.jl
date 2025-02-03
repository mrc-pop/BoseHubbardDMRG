#/usr/bin/julia

using LsqFit

# Include the dmrg.jl file
include("./dmrg.jl")

# Set default options for plotting
Plots.default(
    size = (550, 400),
    bgcolor = :white,
    palette = :seaborn_colorblind,
    
    #gridlinewidth = 0.5,
    #gridstyle = :dot,
    #gridcolor = :gray,
    
    linewidth = 1.0,
    
    foreground_color = :black,
    foreground_color_axis = :black,
    
    tick_direction = :in,
    minorticks = true,
    
    # fontfamily = "Computer Modern",
    # legendfontfamily = "Computer Modern",
    # titlefontfamily = "Computer Modern",
    # guidefontfamily = "Computer Modern",
    # tickfontfamily = "Computer Modern",
    
    legend_background_color = nothing,
    legend=:topright,
    legend_foreground_color= nothing,
    legend_font_halign=:center,

    plot_titlefontsize = 12,
    legendfontsize = 9,
    guidefontsize = 10,
    tickfontsize = 10,

    markerstrokewidth=0,
    margin = 5Plots.mm,
    framestyle = :box,
)

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

function PlotPhaseBoundaries(FilePath::String, LL::Array{Int64}, μ0::Float64; 
    gap=false, overwrite=true)
    """
    Plot the phase boundaries between the Mott Insulator (MI) and Superfluid 
    (SF) phases. If gap=true, plot the charge gap instead of the phase 
    boundaries. If overwrite=true, clears previous plots. 
    """
    BoundariesData = readdlm(FilePath, ',', Float64, '\n'; comments=true)

    if overwrite
        plot()
    end

    for L in LL
        indices = (BoundariesData[:, 1] .== L)
        JJ = BoundariesData[indices,2]
        ΔEplus = BoundariesData[indices,3]
        ΔEminus = BoundariesData[indices,4]
        
        if gap
            plot!(JJ, ΔEplus - ΔEminus, 
                size=(600, 400),
                label=false,
                xlabel=L"$J$", ylabel=L"$\Delta E_{\mathrm{gap}}$", 
                title=L"Charge gap  as a function of $J$ ($\mu_0=%$μ0$)",
                linestyle=:dot)
        else
            plot!(JJ, [μ0 .+ ΔEminus, μ0 .+ ΔEplus], 
            xlabel=L"$J$", ylabel=L"$\mu$",
            size=(600, 400), 
            label=[L"$\mu_c^- (L=%$L)$" L"$\mu_c^+ (L=%$L)$"],
            title=L"Extrapolation of $\mu_c^\pm(J)$ ($\mu_0=%$μ0$)",
            seriestype=:scatter,
            markersize=3)
        end
    end
end

function FitPhaseBoundaries(FilePathIn::String, FilePathOut::String; 
    makeplot=false, makesinglefitplot=false)
    """
    Fit ΔEplus(L) and ΔEminus(L) as functions of 1/L for each J, and extract 
    their values at L -> ∞. Save the results to a file. 
    If makeplot=true, plots the result.

    The output file will have the following columns:
    # J, ΔEplus(∞), ΔEminus(∞), error_ΔEplus(∞), error_ΔEminus(∞), chi2_plus, chi2_minus
    """

    # Read the input data
    data = readdlm(FilePathIn, ',', Float64, '\n'; comments=true)

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
        ΔEplus = data[filter, 3]
        ΔEminus = data[filter, 4]

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

        # Plot the fit results for one specific J if makesinglefitplot=true
        if makesinglefitplot && J == JJ[8] # TODO change which one to plot
            # Define the best-fit lines
            fit_x = range(minimum(inv_L), maximum(inv_L), length=100)
            fit_y_plus = fit_func(fit_x, fit_plus.param)
            fit_y_minus = fit_func(fit_x, fit_minus.param)

            # Plot the data and best-fit lines
            scatter(inv_L, ΔEplus, label=L"$\Delta E^+$ data", xlabel=L"1/L", 
                ylabel=L"$\Delta E$", title="Fit results for J = $(J)")
            plot!(fit_x, fit_y_plus, label=L"$\Delta E^+$ fit")
            scatter!(inv_L, ΔEminus, label=L"$\Delta E^-$ data")
            plot!(fit_x, fit_y_minus, label=L"$\Delta E^-$ fit")

            savefig(string(@__DIR__,"/../analysis/phaseboundaries_singleplot.pdf"))
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
    
    if makeplot
        plot!(JJ,
              [ΔEplus_fit ΔEminus_fit],
              label=[L"$\mu_c^+$ (fitted)" L"$\mu_c^-$ (fitted)"], 
              xlabel=L"J", ylabel=L"$\mu$", 
              title="Fitted phase boundaries")
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
    JJ = unique(data[:, 1])
    LL = unique(data[:, 2])

    # Mastruzzo to extract array of Γ
    function parse_array(str)
        return parse.(Float64, split(strip(str, ['[', ']', ' ']), ','))
    end
    Γall = [parse_array(row[3]) for row in eachrow(data)]
    # display(Γall) # (TODO check, but it seems to work)

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
            filter = (data[:, 1] .== J) .& (data[:, 2] .== L)

            # Extract the correlators {Γ(r,J,L)}_r for the selected J,L
            Γeven = Γall[filter][1] # this is an array, Γ(r even)
            r = range(start=2, step=2, length=length(Γeven)) # r even

            println("Fitting J=$J, L=$L, r=$r.")

            # Perform the fit of Γ(r) vs r, at fixed (J,L).
            # The argument w means weights, they should equal 1/σ^2. 
            # weights = 1.0 ./ e_Γeven.^2 # TODO estimate errors for Γeven!
            fit = curve_fit(fit_func, r, Γeven, p0) # w=weights 

            println("Best-fit parameters:", fit.param)

            K_JL = fit.param[2] # K(J,L)
            e_K_JL = stderror(fit)[2]
            chi2n = sum(fit.resid.^2) / (length(r) - length(p0))
            push!(K_FSS, K_JL)
            push!(e_K_FSS, e_K_JL)
            push!(chi2n_FSS, chi2n)

            if J == JJ[5] && L == LL[4]
                # Make a single plot, hopefully representative of the first
                # series of fits. Save the results.
                scatter(r, Γeven,
                    xlabel=L"$r$",
                    ylabel=L"$\Gamma(r)$",
                    title=L"Correlation function ($L=%$L, J=%$J, \mu=1.0$)",
                    label="DMRG data",
                    legend=:topright)
                fit_x = range(1.0, maximum(r), length=80)
                fit_y = fit_func(fit_x, fit.param)
                K_JL_plot = round(K_JL, digits=2)
                plot!(fit_x, fit_y, label=L"Best-fit ($K=%$K_JL_plot$)")
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

        if J == JJ[5]
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
    #PlotVariance("../simulations/data_variance.txt", L, nmax, i)
    #savefig("../analysis/variance.pdf")

    # ----------------------------------
    # --- Boundary between SF and MI ---
    # ----------------------------------
    # LL = [10,20,30]
    # FilePathIn = string(@__DIR__, "/../simulations/phaseboundaries.txt")
    # # FilePathPlot = string(@__DIR__, "/../analysis/phaseboundaries.pdf")
    # FilePathFit = string(@__DIR__, "/../analysis/phaseboundaries_fitted.txt")

    # PlotPhaseBoundaries(FilePathIn, LL, μ0)
    # FitPhaseBoundaries(FilePathIn, FilePathFit; makeplot=true, makesinglefitplot=false)

    # savefig(string(@__DIR__,"/../analysis/phaseboundaries.pdf"))
    
    # ------------------------------
    # --- Correlation function Γ ---
    # ------------------------------
    FilePathIn = string(@__DIR__, "/../simulations/correlators.txt")
    FilePathFit = string(@__DIR__, "/../analysis/correlators_fitted.txt")
    FitCorrelationFunction(FilePathIn, FilePathFit)
end

main()