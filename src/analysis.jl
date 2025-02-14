#!/usr/bin/julia

PROJECT_ROOT = @__DIR__ # Absloute path up to .../BoseHubbardDMRG/src

# Include modules
include(PROJECT_ROOT * "/modules/fits.jl")
include(PROJECT_ROOT * "/modules/plots.jl")

# Include setup
include(PROJECT_ROOT * "/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/setup/simulations_setup.jl")

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
        
        # ---------------------------- Plot Heatmap ----------------------------
        
        if UserMode=="--heatmap"

			PhaseBoundariesLL = HorizontalLL
            LL = RectangularLL

            FilePathIn = PROJECT_ROOT * "/../simulations/tip_sweep/L=$(L)_site=$(ceil(Int64, L/2)).txt"
            PhaseBoundariesFilePath = PROJECT_ROOT * "/../simulations/horizontal_sweep/L=$PhaseBoundariesLL.txt"

			for L in LL
				# TODO Do we need today()?
				VarianceFilePathOut = PROJECT_ROOT * "/../analysis/heatmap/variance_L=$L_$today().pdf" # Variance plot
	            AFilePathOut = PROJECT_ROOT * "/../analysis/heatmap/a_L=$L_$today().pdf"       # <a_i> plot
    	        KFilePathOut = PROJECT_ROOT * "/../analysis/heatmap/K_L=$L_$today().pdf"       # K plot
				
	            PlotHeatmap(L, FilePathIn; PhaseBoundariesFilePath, VarianceFilePathOut, AFilePathOut, KFilePathOut)
   			end

        # --------------------- Boundary between SF and MI ---------------------

        elseif UserMode=="--boundaries"
        	
        	# TODO Use μ0 = 0
			μ0 = Horizontalμμ[1]

            FilePathIn = PROJECT_ROOT * "/../simulations/horizontal_sweep/μ0=$(μ0)_L=$HorizontalLL.txt"
            PhaseBoundariesDir = PROJECT_ROOT * "/../analysis/phase_boundaries/"
            
            # TODO Do we need today()?
            FilePathPlot = PhaseBoundariesDir * "phaseboundaries_today().pdf"
            FilePathFit = PhaseBoundariesDir * "fitted_phase_boundaries_today().txt"
            
            # TODO Do we need today()?
            FilePathPlotOut = PhaseBoundariesDir * "phaseboundaries_fit_today().pdf"
            FilePathSinglePlotOut = PhaseBoundariesDir * "phaseboundaries_fit_single_today().pdf"

            PlotPhaseBoundaries(FilePathIn; gap=false, FilePathOut=FilePathPlot)
            FitPhaseBoundaries(FilePathIn, FilePathFit; FilePathPlotOut, FilePathSinglePlotOut)


        # ----------------------- Correlation function Γ -----------------------

        elseif UserMode=="--gamma"
        
        	# TODO Use μ0 = 0
			μ0 = Horizontalμμ[1]
			
            FilePathIn = PROJECT_ROOT * "/../simulations/horizontal_sweep/μ0=$(μ0)_L=$HorizontalLL.txt"
            
            # Plot
            # TODO Do we need today()?
            FileGammaPlot = PROJECT_ROOT * "/../analysis/gamma/gamma_data_plot_today().pdf"
            j = 50 # CHANGE: which J to choose for the plot
            PlotCorrelationFunction(FilePathIn, j; FilePathOut=FileGammaPlot, overwrite=false)

            # Fit, Plot
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
