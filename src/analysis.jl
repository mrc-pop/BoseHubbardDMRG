#!/usr/bin/julia

PROJECT_ROOT = @__DIR__ # Absloute path up to .../BoseHubbardDMRG/src

# Include functions and graphic file
# include(PROJECT_ROOT * "/dmrg.jl")
include(PROJECT_ROOT * "/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/modules/fits.jl")
include(PROJECT_ROOT * "/modules/plots.jl")

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
