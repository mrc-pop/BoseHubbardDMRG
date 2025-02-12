#!/usr/bin/julia

# ---------------------------------- Heatmap -----------------------------------

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

# ----------------------------- Phase boundaries -------------------------------

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

# --------------------------- Correlation functions ----------------------------

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
