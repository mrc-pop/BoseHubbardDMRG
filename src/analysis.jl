#/usr/bin/julia

PROJECT_ROOT = @__DIR__ # Absloute path up to .../BoseHubbardDMRG/src

using LsqFit

# Include functions and graphic file
include(PROJECT_ROOT * "/dmrg.jl")
include(PROJECT_ROOT * "/graphic_setup.jl")

function PlotHeatmap(L::Int64,
                      FilePath::String;
                      FilePathHorizontalSweepIn="",
                      FilePathOut="",
                      FilePathOutA="")
    """
    Plot the variance of the number of particles from data saved
    in `FilePath`.
    """
    
    # Extract data coming from rectangular_sweep
    VarianceData = readdlm(FilePath, ',', Float64, '\n'; comments=true)

    JJ = VarianceData[:,1]
    μμ = VarianceData[:,2]
    EE = VarianceData[:,3]
    varvar = VarianceData[:,4]
    aa = VarianceData[:,5]

    NumJ = length(unique(JJ))
    Numμ = length(unique(μμ))

    # Store variance and order parameter <a_i>
    vars = zeros(Numμ, NumJ)
    orderparameters = zeros(Numμ, NumJ)

    for jj in 1:NumJ
        vars[:,jj] = varvar[ Numμ*(jj-1)+1 : Numμ*jj ]
        orderparameters[:,jj] = aa[ Numμ*(jj-1)+1 : Numμ*jj ]
    end

    i = (ceil(Int64, L/2)) # site index

    # Plot variance
    heatmap(unique(JJ), unique(μμ), vars, 
            xlabel=L"J",
            ylabel=L"μ",
            title=L"Variance $\delta n_i^2$ ($L=%$L, i=%$i$)")
    ylabel!(L"μ")

    if FilePathHorizontalSweepIn != ""
        # If this argument is specified, plot the finite-size
        # phase boundaries, for the closest L available
        BoundariesData = readdlm(FilePathHorizontalSweepIn, ';', '\n'; comments=true)
        LL = unique(BoundariesData[:, 1])

        L_PB_index = argmin(abs.(LL .- L)) # best approximation of L
        L_PhaseBoundaries = LL[L_PB_index]

        if L_PhaseBoundaries !== L
            println("L=$L is not among horizontal data. Using the closest available size, L=$L_PhaseBoundaries. ")
        end

        # Filter data corresponding to L_PhaseBoundaries
        indices = (BoundariesData[:, 1] .== L_PhaseBoundaries)
        JJ_PB = BoundariesData[indices, 2]
        ΔEplus = BoundariesData[indices, 4]
        ΔEminus = BoundariesData[indices, 5]
        
        plot!(JJ_PB, -ΔEminus, 
            #xlabel=L"$J$", ylabel=L"$\mu$",
            label=L"$L=%$L_PhaseBoundaries$",
            seriestype=:scatter,
            markersize=1.5,
            color="green",
            xlimits=(minimum(JJ), maximum(JJ)),
            ylimits=(minimum(μμ), maximum(μμ)))
        plot!(JJ_PB, ΔEplus, seriestype=:scatter,
            label="",
            markersize=1.5,
            color="green",
            xlimits=(minimum(JJ), maximum(JJ)),
            ylimits=(minimum(μμ), maximum(μμ)))        
    end

    if FilePathOut != ""
        # If this argument is specified, save plot.
        savefig(FilePathOut)
        println("Variance plot for L=$L saved on file!")
    end

    # Plot order parameter <a_i>
    heatmap(unique(JJ), unique(μμ), orderparameters, 
            xlabel=L"J",
            ylabel=L"μ",
            title=L"$\langle a_i \rangle$ ($L=%$L, i=%$i$)")

    if FilePathOutA != ""
        savefig(FilePathOutA)
        println("Order parameter plot for L=$L saved on file!")
    end
end

function PlotPhaseBoundaries(FilePathIn::String; 
    						 FilePathOut="",
    						 gap=false,
    						 overwrite=true, 
    						 CustomLL=[],
    						 μ0=0.0)
    
    """
    Plot the phase boundaries between the Mott Insulator (MI) and Superfluid 
    (SF) phases, calculated from RectangularSweep.
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
            plot!(JJ, ΔEplus + ΔEminus, 
                xlabel=L"$J$", ylabel=L"$\Delta E_{\mathrm{gap}}$", 
                title=L"Charge gap  as a function of $J$ ($\mu_0=%$μ0$)",
                seriestype=:scatter,
                markersize=1.5,
                label=L"$L=%$L$",
                color=MyColors[l % length(MyColors)])
        else
            plot!(JJ, μ0 .- ΔEminus, 
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
            println("\nGap for L=$(Int.(LL)) plotted to ", FilePathOut)
        else
            println("\nPhase boundaries for L=$(Int.(LL)) plotted to ", FilePathOut)
        end
    end
end

function FitPhaseBoundaries(FilePathIn::String,
							FilePathOut::String;
							FilePathPlotOut="",
							FilePathSinglePlotOut="")
    						# makeplot=false, 
    						# makesinglefitplot=false)
    						
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
        # Filter data for the current J # and only include rows where L > 20
        filtered_data = (data[:, 2] .== J) #.& (data[:,1] .> 20.0)
        L = data[filtered_data, 1]
        ΔEplus = data[filtered_data, 4]
        ΔEminus = data[filtered_data, 5]

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
        if FilePathSinglePlotOut != "" && J == JJ[end-20]
            # Define the best-fit lines
            fit_x = range(0.0, maximum(inv_L)*1.1, length=100)
            fit_y_plus = fit_func(fit_x, fit_plus.param)
            fit_y_minus = fit_func(fit_x, fit_minus.param)

            # Plot the data and best-fit lines
            round_J = round(J, digits=3)
            scatter(inv_L, ΔEplus, label=L"$\Delta E^+$ data", xlabel=L"1/L", 
                ylabel=L"$\Delta E$", title="FSS fit results for J = $round_J",
                legend=:topleft, markersize=2, color=MyColors[1],
                xlimits=(0.0,maximum(inv_L)*1.1))
            plot!(fit_x, fit_y_plus, label=L"$\Delta E^+$ fit", color=MyColors[1],
                alpha=0.7)
            #scatter!(inv_L, -ΔEminus, label=L"$\Delta E^-$ data", markersize=2,
            #    color=MyColors[2])
            #plot!(fit_x, -fit_y_minus, label=L"$\Delta E^-$ fit", color=MyColors[2],
            #    alpha=0.7)

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
              [ΔEplus_fit, -ΔEminus_fit],
              label=[L"\mu_c^+ \, (L \rightarrow \infty)" L"\mu_c^- \, (L \rightarrow \infty)"], 
              xlabel=L"J", ylabel=L"$\mu$", 
              title="Fitted phase boundaries",
              alpha=1.0)
        savefig(FilePathPlotOut)
        println("Phase boundaries fit plotted to ", FilePathPlotOut)
    end
end

function PlotCorrelationFunction(FilePathIn::String,
                                 j::Int64;
                                 FilePathOut="",
                                 overwrite=true)

    """ Plot the correlation function Γ(r) for the chosen J (the j-th)."""

    # Read the input data
    data = readdlm(FilePathIn, ';', '\n'; comments=true)

    # Extract unique J, L values
    JJ = unique(data[:, 2])
    LL = unique(data[:, 1])

    println("\nPlotting correlation function.")
    println("From input file, there are $(length(JJ)) possible values of J.")
    
    # Mastruzzo to extract array of Γ
    function parse_array(str)
        return parse.(Float64, split(strip(str, ['[', ']', ' ']), ','))
    end
    Γall = [parse_array(row[7]) for row in eachrow(data)]

    if overwrite
        plot()
    end

    J = JJ[j] # choose the j-th J

    println("Chosen value of J: $J")

    for L in LL
        # Filter data for the current J and L
        filter = (data[:, 2] .== J) .& (data[:, 1] .== L)

        # Extract the correlators {Γ(r,J,L)}_r for the selected J,L
        Γeven = Γall[filter][1] # this is an array, Γ(r even)
        r = range(start=2, step=2, length=length(Γeven)) # r even

        # Check if Γeven has at least two elements
        if length(Γeven) < 2
            println("Warning: Not enough data points for L = $L. Skipping.")
            continue
        end

        J_round = round(J, digits=3)

        scatter!(r, Γeven,
            xlabel=L"$r$",
            ylabel=L"$\Gamma(r)$",
            title=L"Correlation function ($J=%$J_round$)",
            label=L"L=%$L",
            markersize=2,
            xscale=:log10,
            yscale=:log10,  
            xticks=[1,10,100],
            xlimits=(1,100),
            legend=:topright)
    end

    if !=(FilePathOut,"")
        savefig(FilePathOut)
        println("Correlator vs r plotted to ", FilePathOut)
    end
end

function FitCorrelationFunction(FilePathIn::String,
								FilePathOut::String;
                                SingleFitPlotPathOut="",
                                PlotPathOut="",
                                FSSPlotPathOut="",
                                SingleFitPlotj::Int64)
								
    """
    Read data from file. Then fit a power-law extracting the exponent K(J,L). 
    Then fit it against L to extract the thermodynamic limit K_∞(J).
    Finally, plot K_∞ against J.
    """

    # CHANGE: PLOT & FIT PARAMETERS
    J_min = 0.25
    L_min = 20
    r_min_fit = 6
    r_max_fit = 20
    fit_func(x, p) = p[1]*x.^(-p[2]/2)  # power-law
    p0 = [1.0, 1.0]                     # power-law
    # fit_func(x, p) = p[1]*x.^(-p[2]/2) .* exp.(-x./p[3]) # power-law exponential
    # p0 = [1.0, 1.0, 100]                                 # power-law exponential

    # Read the input data
    data = readdlm(FilePathIn, ';', '\n'; comments=true)

    # Extract unique J, L values
    JJ = unique(data[:, 2])
    LL = unique(data[:, 1])

    println("Putting a cutoff of J>$J_min")
    println("Fit range restricted to $r_min_fit<=r<=$r_max_fit")

    # Select which J to plot the fits
    J_plot = JJ[SingleFitPlotj]

    # Apply the cutoff
    JJ = JJ[ JJ .> J_min]
    LL = LL[ LL .> L_min]

    # Mastruzzo to extract array of Γ
    function ParseArray(str)
        return parse.(Float64, split(strip(str, ['[', ']', ' ']), ','))
    end
    
    Γall = [ParseArray(row[7]) for row in eachrow(data)]
    eΓall = [ParseArray(row[8]) for row in eachrow(data)]
    Call = [ParseArray(row[9]) for row in eachrow(data)]
    eCall = [ParseArray(row[10]) for row in eachrow(data)]
    
    # display(Γall)

    K = [] # K_∞ (thermodynamic limit) as a function of j
    e_K = [] # error on K_∞
    chi2n_K = []

    for J in JJ
        # Initialize {K(J,L)}_L at fixed J, its error and χ^2/ndof
        K_FSS = [] 
        e_K_FSS = []
        chi2n_FSS = []

        plot()
        MyColors = ColorSchemes.tab10

        for (l,L) in enumerate(LL)
            # Filter data for the current J and L
            filter = (data[:, 2] .== J) .& (data[:, 1] .== L)

            # Extract the correlators {Γ(r,J,L)}_r for the selected J,L
            Γeven = Γall[filter][1] # this is an array, Γ(r even)
            # display(Γeven)
            e_Γeven = eΓall[filter][1]
            r = collect(range(start=2, step=2, length=length(Γeven))) # r even

            # Restrict to fit range
            fit_range_indices = findall(x -> x >= r_min_fit && x <= r_max_fit, r)
            r = r[fit_range_indices]
            Γeven = Γeven[fit_range_indices]
            e_Γeven = e_Γeven[fit_range_indices]

            println("\nFitting J=$J, L=$L, r=$r (set by r_min_fit, r_max_fit).")

            # Perform the fit of Γ(r) vs r, at fixed (J,L).
            # The argument w means weights, they should equal 1/σ^2. 
            weights = 1.0 ./ e_Γeven.^2
            fit = curve_fit(fit_func, r, Γeven, p0) # (..., weights, p0)

            println("  Done. Best-fit parameters: ", round.(fit.param, digits=4))

            K_JL = fit.param[2] # K(J,L)
            e_K_JL = stderror(fit)[2]
            chi2n = sum(fit.resid.^2) / (length(r) - length(p0))
            push!(K_FSS, K_JL)
            push!(e_K_FSS, e_K_JL)
            push!(chi2n_FSS, chi2n)

            if J == J_plot
                # Make a single plot, hopefully representative of the first
                # series of fits. Save the results.
                color = MyColors[l % length(MyColors)]  # color for current L
                K_JL_plot = round(K_JL, digits=2)

                scatter!(r, Γeven,
                    xlabel=L"$r$",
                    ylabel=L"$\Gamma(r)$",
                    #title=L"Correlation function (J=%$J, \mu=1.0$)",
                    label=L"$L=%$L$ ($K_\mathrm{fit}=%$K_JL_plot$)",
                    legend=:topright,
                    markersize=2,
                    xscale=:log10,       # for log plot
                    yscale=:log10,       # for log plot
                    ##yticks=[0.1, 1],     # for log plot
                    ##ylimits=(0.1, 1),    # for log plot
                    xticks=[1,10,100],   # for log plot
                    xlimits=(1,100),     # for log plot
                    color=color
                    )
                fit_x = range(1.5, maximum(r)*1.5, length=80)
                fit_y = fit_func(fit_x, fit.param)
                plot!(fit_x, fit_y, label=false, color=color, alpha=0.7)
                if SingleFitPlotPathOut!=""
                    J_round = round(J, digits=3)
                    title!(L"Correlation function power-law fit ($J=%$J_round$)")
                    savefig(SingleFitPlotPathOut)
                end
            end
        end

        if length(LL)>1
            # FSS fit.
            inv_LL = 1 ./ LL

            fit_func_fss(x, p) = p[1] * x .+ p[2]
            p0_FSS = [-1.0, 1.0]

            # If there is more than 1 size, perform the FSS fit 
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
                if FSSPlotPathOut!=""
                    savefig(FSSPlotPathOut)
                end
            end
        elseif length(LL)==1
            push!(K, K_FSS[1])
        end
    end # end loop over J

    if length(LL)==1 
        # If only 1 size, bypass the FSS fit, plot the results (K(L) vs J).
        L = LL[1]

        println("Note: number of sizes is $(length(LL)). Skipping FSS fit.")

        scatter(JJ, K, xlabel=L"$J$", ylabel=L"$K$",
            title=L"Power-law decay of $\Gamma(r)$ ($%$r_min_fit \le r \le %$r_max_fit$)",
            markersize=2,
            label="Data for L=$L")
        plot!(JJ, 0.5*ones(length(JJ)), label=L"$K_\mathrm{th} = 1/2$")
        
        if PlotPathOut!=""
            savefig(PlotPathOut)
            println("Plotted K(L=$L) vs J on file.")
        end
        
    elseif length(LL)>1
        # Join the results for saving correctly on file
        results = hcat(JJ, K, e_K, chi2n_K)

        # Save the results
        open(FilePathOut, "w") do file
            write(file, "# J, K_∞, e_K_∞, chi2n\n")
            writedlm(file, results, ',')
        end

        # Plot the results (K_∞ vs J)
        scatter(JJ, K, xlabel=L"$J$", ylabel=L"$K_\infty$", yerr=e_K, 
            title=L"Power-law decay of $\Gamma(r)$ ($%$r_min_fit \le r \le %$r_max_fit$)",
            markersize=2,
            label="Fitted data")
        plot!(JJ, 0.5*ones(length(JJ)), label=L"$K_\mathrm{th} = 1/2$")
        if PlotPathOut!=""
            savefig(PlotPathOut)
        end     
    end
end


function main()    
    # L = 10
    # N = 10
    # nmax = 3
    # JJ = collect(range(0.0,0.3,5))
    # μμ = collect(range(0.0,1,5))
    # i = ceil(Int64, L/2) # site on which to calculate variance
    # μ0 = 0.0


    ModeErrorMsg = "Input error: use option --heatmap, --boundaries or --gamma"
	
	if length(ARGS) != 1
		# If user does not specify the user mode
		error(ModeErrorMsg)
		exit()
	else
        UserMode = ARGS[1]

        # TODO add other --options

        # --------------------
        # --- Plot Heatmap ---
        # --------------------
        if UserMode=="--heatmap"

            L = 12 # TODO change!!

            FilePathIn = PROJECT_ROOT * "/../simulations/rectangular_sweep/L=$(L)_site=$(ceil(Int64, L/2)).txt"
            FilePathHorizontalSweepIn = PROJECT_ROOT * "/../simulations/horizontal_sweep/L=[10, 20, 30, 40, 50, 60, 70, 80].txt"

            FilePathOut = PROJECT_ROOT * "/../analysis/heatmap/variance_250211.pdf" # Variance plot
            FilePathOutA = PROJECT_ROOT * "/../analysis/heatmap/a_250211.pdf"       # <a_i> plot

            PlotHeatmap(L, FilePathIn; FilePathOut=FilePathOut, FilePathOutA=FilePathOutA,
                FilePathHorizontalSweepIn=FilePathHorizontalSweepIn)

        # ----------------------------------
        # --- Boundary between SF and MI ---
        # ----------------------------------
        elseif UserMode=="--boundaries"    
            FilePathIn = PROJECT_ROOT * "/../simulations/horizontal_sweep/L=[10, 20, 30, 40, 50, 60, 70, 80].txt"
            PhaseBoundariesDir = PROJECT_ROOT * "/../analysis/phase_boundaries/"
            
            FilePathPlot = PhaseBoundariesDir * "phaseboundaries_250211.pdf"
            FilePathFit = PhaseBoundariesDir * "fitted_phase_boundaries.txt"
            
            FilePathPlotOut = PhaseBoundariesDir * "phaseboundaries_250211_fit.pdf"
            FilePathSinglePlotOut = PhaseBoundariesDir * "phaseboundaries_250211_fit_single.pdf"

            PlotPhaseBoundaries(FilePathIn; gap=false, FilePathOut=FilePathPlot)

            FitPhaseBoundaries(FilePathIn, FilePathFit; FilePathPlotOut, FilePathSinglePlotOut)

        # ------------------------------
        # --- Correlation function Γ ---
        # ------------------------------
        elseif UserMode=="--gamma"
            FilePathIn = PROJECT_ROOT * "/../simulations/horizontal_sweep/L=[10, 20, 30, 40, 50, 60, 70, 80].txt"
            
            # PLOT
            FileGammaPlot = PROJECT_ROOT * "/../analysis/gamma/gamma_data_plot.pdf"
            j = 50 # CHANGE: which J to choose for the plot
            PlotCorrelationFunction(FilePathIn, j; FilePathOut=FileGammaPlot, overwrite=false)

            # FIT (& PLOT)
            PlotPathOut = PROJECT_ROOT * "/../analysis/gamma/gamma_final_plot.pdf"
            SinglePlotPathOut = PROJECT_ROOT * "/../analysis/gamma/gamma_power_law_fit_plot.pdf"
            FilePathFit = PROJECT_ROOT * "/../analysis/gamma/fit_correlation.txt"
            FitCorrelationFunction(FilePathIn, FilePathFit; PlotPathOut=PlotPathOut, SingleFitPlotPathOut=SinglePlotPathOut, SingleFitPlotj=j)
		else
			error(ModeErrorMsg)
			exit()
        end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
