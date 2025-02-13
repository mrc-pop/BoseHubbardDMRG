#!/usr/bin/julia

PROJECT_ROOT = @__DIR__ # Absloute path up to .../BoseHubbardDMRG/src

# Include functions and graphic file
include(PROJECT_ROOT * "/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/modules/fits.jl")
include(PROJECT_ROOT * "/modules/plots.jl")
using DelimitedFiles

function main()    

    ModeErrorMsg = "Input error: use option --heatmap, --boundaries or --gamma"
	
	if length(ARGS) != 1
		# If user does not specify the user mode
		error(ModeErrorMsg)
		exit()
	else
        UserMode = ARGS[1]

        # TODO add other --options

        # ----------------------------------------------------------------------
        # ---------------------------- Plot Heatmap ----------------------------
        # ----------------------------------------------------------------------
        if UserMode=="--heatmap"

            L = 20 # TODO change!!

            FilePathIn = PROJECT_ROOT * "/../simulations/tip_sweep/L=$(L)_site=$(ceil(Int64, L/2)).txt"
            FilePathHorizontalSweepIn = PROJECT_ROOT * "/../simulations/horizontal_sweep/L=[10, 20, 30, 40, 50, 60, 70].txt"

            FilePathOut = PROJECT_ROOT * "/../analysis/heatmap/variance_250212.pdf" # Variance plot
            FilePathOutA = PROJECT_ROOT * "/../analysis/heatmap/a_250212.pdf"       # <a_i> plot
            FilePathOutK = PROJECT_ROOT * "/../analysis/heatmap/K_250212.pdf"       # K plot

            PlotHeatmap(L, FilePathIn; FilePathOut, FilePathOutA,
                FilePathHorizontalSweepIn, FilePathOutK)
   

        # ----------------------------------------------------------------------
        # --------------------- Boundary between SF and MI ---------------------
        # ----------------------------------------------------------------------
        elseif UserMode=="--boundaries"

            FilePathIn = PROJECT_ROOT * "/../simulations/horizontal_sweep/L=[10, 20, 30, 40, 50, 60, 70].txt"
            PhaseBoundariesDir = PROJECT_ROOT * "/../analysis/phase_boundaries/"
            
            FilePathPlot = PhaseBoundariesDir * "phaseboundaries_250212.pdf"
            FilePathFit = PhaseBoundariesDir * "fitted_phase_boundaries.txt"
            
            FilePathPlotOut = PhaseBoundariesDir * "phaseboundaries_250212_fit.pdf"
            FilePathSinglePlotOut = PhaseBoundariesDir * "phaseboundaries_250212_fit_single.pdf"

            PlotPhaseBoundaries(FilePathIn; gap=false, FilePathOut=FilePathPlot)

            FitPhaseBoundaries(FilePathIn, FilePathFit; FilePathPlotOut, FilePathSinglePlotOut)

        # ----------------------------------------------------------------------
        # ----------------------- Correlation function Γ -----------------------
        # ----------------------------------------------------------------------
        elseif UserMode=="--gamma"
            FilePathIn = PROJECT_ROOT * "/../simulations/horizontal_sweep/L=[10, 20, 30, 40, 50].txt"
            
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
