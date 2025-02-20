#!/usr/bin/julia

using LsqFit

# ----------------------------- Phase boundaries -------------------------------

function FitPhaseBoundaries(FilePathIn::String,
							FilePathOut::String;
                            μ0=0.0,
							FilePathPlotOut="",
							FilePathSinglePlotOut="")
    						# makeplot=false, 
    						# makesinglefitplot=false)
    						
    """
    Fit μUp(L) and μDown(L) as functions of 1/L for each J, and extract 
    their values at L -> ∞. Save the results to a file. 
    If makeplot=true, plots the result.

    The output file will have the following columns:
    # J, μUp(∞), μDown(∞), error_μUp(∞), error_μDown(∞), chi2_plus, chi2_minus
    """

    # Read the input data
    data = readdlm(FilePathIn, ';', '\n'; comments=true)

    # Extract unique J values
    JJ = unique(data[:, 2])

    μUp_fit = []
    μDown_fit = []
    error_μUp = []
    error_μDown = []
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
        μUp = data[filtered_data, 4]
        μDown = data[filtered_data, 5]

        inv_L = 1.0 ./ L

        # Fit μUp
        fit_plus = curve_fit(fit_func, inv_L, μUp, p0) 
        μUp_inf = fit_plus.param[2]  # Intercept (value at L -> ∞)
        error_plus = stderror(fit_plus)[2]
        chi2p = sum(fit_plus.resid.^2)

        # Fit μDown
        fit_minus = curve_fit(fit_func, inv_L, μDown, p0)
        μDown_inf = fit_minus.param[2]  # Intercept (value at L -> ∞)
        error_minus = stderror(fit_minus)[2]  
        chi2m = sum(fit_minus.resid.^2)

        # Save the results for this J
        push!(μUp_fit, μUp_inf)
        push!(μDown_fit, μDown_inf)
        push!(error_μUp, error_plus)
        push!(error_μDown, error_minus)
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
            scatter(inv_L, μUp, label=L"$\Delta^+$ data", xlabel=L"1/L", 
                ylabel=L"$\Delta^+$", title=L"FSS fit results for $J = %$round_J$ ($\mu_0 = %$μ0$)",
                legend=:topleft, markersize=2, color=MyColors[1],
                xlimits=(0.0,maximum(inv_L)*1.1))
            plot!(fit_x, fit_y_plus, label=L"$\Delta^+$ fit", color=MyColors[1],
                alpha=0.7)
            #scatter!(inv_L, -μDown, label=L"$\Delta E^-$ data", markersize=2,
            #    color=MyColors[2])
            #plot!(fit_x, -fit_y_minus, label=L"$\Delta E^-$ fit", color=MyColors[2],
            #    alpha=0.7)

            savefig(FilePathSinglePlotOut)
        end
    end

    # Join the results for saving correctly on file
    results = hcat(JJ, μUp_fit, μDown_fit, error_μUp, error_μDown, 
    chi2_plus, chi2_minus)

    # Save the results
    open(FilePathOut, "w") do file
        write(file, "# J, μUp(∞), μDown(∞), error_μUp(∞), ", 
        "error_μDown(∞), chi2_plus, chi2_minus\n")
        writedlm(file, results, ',')
    end

    println("Fitted results saved to: ", FilePathOut)
    
    if FilePathPlotOut != ""
        # PlotPhaseBoundaries(FilePathIn; gap=false, overwrite=true)
        plot()
        plot!(JJ,
              [μUp_fit, -μDown_fit .+ 2 * μ0],
              label=[L"\mu_c^+ \, (L \rightarrow \infty)" L"\mu_c^- \, (L \rightarrow \infty)"], 
              xlabel=L"J", ylabel=L"$\mu$", 
              title=L"Fitted phase boundaries ($\mu_0 = %$μ0$)",
              alpha=1.0)
        savefig(FilePathPlotOut)
        println("Phase boundaries fit plotted to ", FilePathPlotOut)
    end
end

# ------------------------------------------------------------------------------
# --------------------------- Correlation functions ----------------------------
# ------------------------------------------------------------------------------

"""
Fit Gamma(r) vs r as a power law. 
If ThreeParameterFit=true, fit a three-parameter exponential and power law
function instead.
If rMin and rMax are specified, the fit range is reduced.

Input:
    - Γeven: array of Γ(r), r even, r≤L/2
    - e_Γeven: array of errors
"""
function FitPowerLawGamma(Γeven::Array{Float64},
                          e_Γeven::Array{Float64};
                          rMin=1,
                          rMax=100,
                          ThreeParameterFit=false)

    ThreeParamsFunction(x, p) = p[1]*x.^(-p[2]/2) .* exp.(-x./p[3]) # power-law * exp
    TwoParamsFunction(x, p) = p[1]*x.^(-p[2]/2)  # power-law

    # Define fit funcion and p0
    if ThreeParameterFit
        p0 = [1.0, 1.0, 100]
        fit_func = ThreeParamsFunction
    else
        p0 = [1.0, 1.0]
        fit_func = TwoParamsFunction
    end

    # Define variable r
    r = collect(range(start=2, step=2, length=length(Γeven))) # r even

    # Restrict to fit range
    FitRangeIndices = findall(x -> x >= rMin && x <= rMax, r)
    r = r[FitRangeIndices]
    Γeven = Γeven[FitRangeIndices]
    e_Γeven = e_Γeven[FitRangeIndices]

    # Perform the fit.
    # The argument w means weights, they should equal 1/σ^2. 
    Weights = 1.0 ./ e_Γeven.^2
    Fit = curve_fit(fit_func, r, Γeven, Weights, p0)
    
    # Extract the best-fit parameters
    K = Fit.param[2]
    e_K = stderror(Fit)[2]
    chi2n = sum(Fit.resid.^2) / (length(r) - length(p0))

    if ThreeParameterFit
        ξ = Fit.param[3]
        e_ξ = stderror(Fit)[3]
        Results = (K, e_K, chi2n, ξ, e_ξ)
    else
        Results = (K, e_K, chi2n)
    end
    
    return Results
end


"""
Fit K vs L (at fixed J) in order to extract K_∞.
Outputs (K_∞, e_K_∞, chi2n).
"""
function FitFiniteSizeScalingK(LL::Array{Int64},
                               KK::Array{Float64},
                               e_KK::Array{Float64})
        
        # Define the fit function
        fit_func_fss(x, p) = p[1] * x .+ p[2]
        p0 = [-1.0, minimum(KK)]

        # The fit is linear in 1/L
        inv_LL = 1 ./ LL

        Weights = 1.0 ./ e_KK.^2
        Fit = curve_fit(fit_func_fss, inv_LL, KK, Weights, p0)

        K_∞ = Fit.param[2]
        e_K_∞ = stderror(Fit)[2]
        chi2n = sum(Fit.resid.^2)/(length(LL)-2)

        Results = (K_∞, e_K_∞, chi2n)

    return Results
end

"""
Perform fits on Γ(r) vs r, varying r_min and r_max. Extrapolate K_∞ as a function
of J, and save it on file.
"""
function FitRoutineGamma(FilePathIn, 
                         FilePathOut;
                         rrMin=[2, 4, 6, 8],
                         rrMax=[12, 14, 16, 18],
                         JMin = 0.0,
                         LMin = 10)

    # Read the input data
    HorizontalData = readdlm(FilePathIn, ';', '\n'; comments=true)

    # Extract unique J, L values
    LL = Int64.(unique(HorizontalData[:, 1]))

    JJ = unique(HorizontalData[:, 2])

    # Select the couplings and sizes on which to perform the fits
    JJ = JJ[ JJ .>= JMin]
    LL = LL[ LL .>= LMin]

    # Mastruzzo to extract array of Γ
    function ParseArray(str)
        return parse.(Float64, split(strip(str, ['[', ']', ' ']), ','))
    end
    Γall = [ParseArray(row[4]) for row in eachrow(HorizontalData)]
    eΓall = [ParseArray(row[5]) for row in eachrow(HorizontalData)]

    # Open output file
    DataFile = open(FilePathOut, "w")
    write(DataFile,"# Fit results on Γ. Minimum size considered: LMin = $LMin.\n")
    write(DataFile,"# J, rMin, rMax, K_∞, e_K_∞, chi2n\n")

    # Perform the fits and save the results
    for J in JJ
        println("Making all the fits for J=$J...")
        for rMin in rrMin
            for rMax in rrMax
                KK_J = zeros(Float64, length(LL))
                e_KK_J = zeros(Float64, length(LL))

                for (l,L) in enumerate(LL)
                    # Select data for current (J,L)
                    filter = (HorizontalData[:, 2] .== J) .& (HorizontalData[:, 1] .== L)
                    Γeven = Γall[filter][1]
                    e_Γeven = eΓall[filter][1]

                    # Error on Γ(r=L/2) is NaN. Copy the error of the preceding one.
                    e_Γeven[end] = e_Γeven[end-1]

                    """
                    println("Fitting Γ as a power-law for J=$(round(J,digits=3)), L=$L, rMin=$rMin, rMax=$rMax...")
                    """

                    # Perform the power law fit on Γ for the current (J,L)
                    K_JL, e_K_JL, chi2n = FitPowerLawGamma(Γeven, e_Γeven; rMin, rMax, ThreeParameterFit=false)
                
                    # Save on the arrays
                    KK_J[l] = K_JL
                    e_KK_J[l] = e_K_JL
                    
                    """
                    println("   Best-fit: K=$(round(K_JL, digits=4)) +/- $(round(e_K_JL, digits=4))\n")
                    """
                end

                # Perform the FSS fit on K for the current J
                K_∞, e_K_∞, chi2n = FitFiniteSizeScalingK(LL, KK_J, e_KK_J)

                # Save the results
                write(DataFile,"$J, $rMin, $rMax, $K_∞, $e_K_∞, $chi2n\n")
            end
        end
        println("   Done! Results saved on file.\n")
    end

    close(DataFile)
end


"""
Estimate the fraction of points in a rectangular sweep for which |⟨b_i⟩| > cutoff
"""
function EstimateFractionOrderParameter(L::Int64,
                                        FilePathIn::String,
                                        PhaseBoundariesFilePath::String,
                                        Cutoff::Float64;
                                        verbose=false)

        println("\nPerforming the fraction estimate with L=$L, cutoff=$Cutoff")
        # Extract data coming from rectangular_sweep
        VarianceData = readdlm(FilePathIn, ';', '\n'; comments=true)
        JJ = VarianceData[:,1]
        μμ = VarianceData[:,2]
        aa = VarianceData[:,5]

        # Extract phase boundaries data
        BoundariesData = readdlm(PhaseBoundariesFilePath, ';', '\n'; comments=true)
        LL = unique(BoundariesData[:, 1])

        if L ∉ LL
            error("L=$L is not among horizontal data.")
        end

        # Filter boundaries data corresponding to L
        indices = (BoundariesData[:, 1] .== L)
        JJ_PB = BoundariesData[indices, 2]

        NumberMI = 0 # how many points are in the MI phase
        NumberSFAboveCutoff = 0 # how many points of the SF phase are above cutoff

        for (i,J) in enumerate(JJ) # run over each point of the rect. sweep
            μ = μμ[i]

            # Find closest J of PB data, and save the corresponding μUp and μDown
            Index = argmin(abs.(JJ_PB .- J))
		    μUp = BoundariesData[Index,4][1]
		    μDown = -BoundariesData[Index,5][1]

            inMottLobe = false
			
            if (μ>=μDown && μ<=μUp)
                inMottLobe = true
                NumberMI += 1
            end

            if verbose && i % 50 == 13 # debug
                println("J=$J, μ=$μ, μDown=$μDown, μUp=$μUp, MI=$inMottLobe")
            end

            if !inMottLobe # for the SF points
                if abs(aa[i]) > Cutoff
                    NumberSFAboveCutoff += 1
                end
            end
        end

        TotalNumber = length(aa)
        NumberSF = TotalNumber - NumberMI
        Fraction = NumberSFAboveCutoff / NumberSF

        println("Out of $TotalNumber points, $NumberMI were in the MI phase.")
        println("Out of the other $NumberSF, exactly $NumberSFAboveCutoff were above",
            " the cutoff of $Cutoff.")
        println("Estimated fraction: $Fraction")
end


# """
# Old function which did it all
# """
# function FitCorrelationFunction(FilePathIn::String,
# 								FilePathOut::String;
#                                 rMinFit=2,
#                                 rMaxFit=100,
#                                 JMin=0.0,
#                                 LMin=0,
#                                 SingleFitPlotPathOut="",
#                                 PlotPathOut="",
#                                 FSSPlotPathOut="",
#                                 SingleFitPlotj::Int64,
#                                 ThreeParameterFit=false)
								
#     """
#     Read data from file. Then fit a power-law extracting the exponent K(J,L). 
#     Then fit it against L to extract the thermodynamic limit K_∞(J).
#     Finally, plot K_∞ against J.
#     """

#     # CHANGE: PLOT & FIT PARAMETERS
#     J_min = 0.25
#     L_min = 20

#     #if ThreeParameterFit
#     #    fit_func(x, p) = p[1]*x.^(-p[2]/2) .* exp.(-x./p[3]) # power-law exponential
#     #    p0 = [1.0, 1.0, 100]                                 # power-law exponential
#     #else
#         fit_func(x, p) = p[1]*x.^(-p[2]/2)  # power-law
#         p0 = [1.0, 1.0]                     # power-law
#     #end

#     # Read the input data
#     data = readdlm(FilePathIn, ';', '\n'; comments=true)

#     # Extract unique J, L values
#     JJ = unique(data[:, 2])
#     LL = unique(data[:, 1])

#     println("Putting a cutoff of J>$J_min")
#     println("Fit range restricted to $rMinFit<=r<=$rMaxFit")

#     # Select which J to plot the fits
#     J_plot = JJ[SingleFitPlotj]

#     # Apply the cutoff
#     JJ = JJ[ JJ .> J_min]
#     LL = LL[ LL .> L_min]

#     # Mastruzzo to extract array of Γ
#     function ParseArray(str)
#         return parse.(Float64, split(strip(str, ['[', ']', ' ']), ','))
#     end
    
#     Γall = [ParseArray(row[7]) for row in eachrow(data)]
#     eΓall = [ParseArray(row[8]) for row in eachrow(data)]
#     Call = [ParseArray(row[9]) for row in eachrow(data)]
#     eCall = [ParseArray(row[10]) for row in eachrow(data)]
    
#     # display(Γall)

#     K = [] # K_∞ (thermodynamic limit) as a function of j
#     e_K = [] # error on K_∞
#     chi2n_K = []

#     for J in JJ
#         # Initialize {K(J,L)}_L at fixed J, its error and χ^2/ndof
#         K_FSS = [] 
#         e_K_FSS = []
#         chi2n_FSS = []

#         plot()
#         MyColors = ColorSchemes.tab10

#         for (l,L) in enumerate(LL)
#             # Filter data for the current J and L
#             filter = (data[:, 2] .== J) .& (data[:, 1] .== L)

#             # Extract the correlators {Γ(r,J,L)}_r for the selected J,L
#             Γeven = Γall[filter][1] # this is an array, Γ(r even)
#             # display(Γeven)
#             e_Γeven = eΓall[filter][1]
#             r = collect(range(start=2, step=2, length=length(Γeven))) # r even

#             # Restrict to fit range
#             fit_range_indices = findall(x -> x >= rMinFit && x <= rMaxFit, r)
#             r = r[fit_range_indices]
#             Γeven = Γeven[fit_range_indices]
#             e_Γeven = e_Γeven[fit_range_indices]

#             println("\nFitting J=$J, L=$L, r=$r (set by rMinFit, rMaxFit).")

#             # Perform the fit of Γ(r) vs r, at fixed (J,L).
#             # The argument w means weights, they should equal 1/σ^2. 
#             weights = 1.0 ./ e_Γeven.^2
#             fit = curve_fit(fit_func, r, Γeven, p0) # (..., weights, p0)

#             println("  Done. Best-fit parameters: ", round.(fit.param, digits=4))

#             K_JL = fit.param[2] # K(J,L)
#             e_K_JL = stderror(fit)[2]
#             chi2n = sum(fit.resid.^2) / (length(r) - length(p0))
#             push!(K_FSS, K_JL)
#             push!(e_K_FSS, e_K_JL)
#             push!(chi2n_FSS, chi2n)

#             if J == J_plot
#                 # Make a single plot, hopefully representative of the first
#                 # series of fits. Save the results.
#                 color = MyColors[l % length(MyColors)]  # color for current L
#                 K_JL_plot = round(K_JL, digits=2)

#                 scatter!(r, Γeven,
#                     xlabel=L"$r$",
#                     ylabel=L"$\Gamma(r)$",
#                     #title=L"Correlation function (J=%$J, \mu=1.0$)",
#                     label=L"$L=%$L$ ($K_\mathrm{fit}=%$K_JL_plot$)",
#                     legend=:topright,
#                     markersize=2,
#                     xscale=:log10,       # for log plot
#                     yscale=:log10,       # for log plot
#                     ##yticks=[0.1, 1],     # for log plot
#                     ##ylimits=(0.1, 1),    # for log plot
#                     xticks=[1,10,100],   # for log plot
#                     xlimits=(1,100),     # for log plot
#                     color=color
#                     )
#                 fit_x = range(1.5, maximum(r)*1.5, length=80)
#                 fit_y = fit_func(fit_x, fit.param)
#                 plot!(fit_x, fit_y, label=false, color=color, alpha=0.7)
#                 if SingleFitPlotPathOut!=""
#                     J_round = round(J, digits=3)
#                     title!(L"Correlation function power-law fit ($J=%$J_round$)")
#                     savefig(SingleFitPlotPathOut)
#                 end
#             end
#         end

#         if length(LL)>1
#             # FSS fit.
#             inv_LL = 1 ./ LL

#             fit_func_fss(x, p) = p[1] * x .+ p[2]
#             p0_FSS = [-1.0, 1.0]

#             # If there is more than 1 size, perform the FSS fit 
#             fit_fss = curve_fit(fit_func_fss, inv_LL, K_FSS, p0_FSS)
#             K_J = fit_fss.param[2]
#             e_K_J = stderror(fit_fss)[2]
#             chi2n_J = sum(fit_fss.resid.^2)/(length(LL)-2)
#             push!(K, K_J)
#             push!(e_K, e_K_J)
#             push!(chi2n_K, chi2n_J)

#             if J == maximum(JJ)
#                 # Make a single FSS plot, hopefully representative of the second
#                 # series of fits. Save the results.
#                 fit_x = range(0.0, maximum(inv_LL), length=50)
#                 fit_y = fit_func_fss(fit_x, fit_fss.param)
#                 scatter(inv_LL, K_FSS,
#                     xlabel=L"$1/L$",
#                     ylabel=L"$K(J,L)$",
#                     title=L"Exponent $K$ when $J=%$J$",
#                     label="Fits from DMRG data",
#                     legend=:topright)
#                     K_infty_plot = round(K_J, digits=2)
#                 plot!(fit_x, fit_y, label=L"Best-fit ($K_\infty = %$K_infty_plot$)")
#                 if FSSPlotPathOut!=""
#                     savefig(FSSPlotPathOut)
#                 end
#             end
#         elseif length(LL)==1
#             push!(K, K_FSS[1])
#         end
#     end # end loop over J

#     if length(LL)==1 
#         # If only 1 size, bypass the FSS fit, plot the results (K(L) vs J).
#         L = LL[1]

#         println("Note: number of sizes is $(length(LL)). Skipping FSS fit.")

#         scatter(JJ, K, xlabel=L"$J$", ylabel=L"$K$",
#             title=L"Power-law decay of $\Gamma(r)$ ($%$rMinFit \le r \le %$rMaxFit$)",
#             markersize=2,
#             label="Data for L=$L")
#         plot!(JJ, 0.5*ones(length(JJ)), label=L"$K_\mathrm{th} = 1/2$")
        
#         if PlotPathOut!=""
#             savefig(PlotPathOut)
#             println("Plotted K(L=$L) vs J on file.")
#         end
        
#     elseif length(LL)>1
#         # Join the results for saving correctly on file
#         results = hcat(JJ, K, e_K, chi2n_K)

#         # Save the results
#         open(FilePathOut, "w") do file
#             write(file, "# J, K_∞, e_K_∞, chi2n\n")
#             writedlm(file, results, ',')
#         end

#         # Plot the results (K_∞ vs J)
#         scatter(JJ, K, xlabel=L"$J$", ylabel=L"$K_\infty$", yerr=e_K, 
#             title=L"Power-law decay of $\Gamma(r)$ ($%$rMinFit \le r \le %$rMaxFit$)",
#             markersize=2,
#             label="Fitted data")
#         plot!(JJ, 0.5*ones(length(JJ)), label=L"$K_\mathrm{th} = 1/2$")
#         if PlotPathOut!=""
#             savefig(PlotPathOut)
#         end     
#     end
# end
