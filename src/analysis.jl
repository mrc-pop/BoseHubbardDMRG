#!/usr/bin/julia

PROJECT_ROOT = @__DIR__ # Absloute path up to .../BoseHubbardDMRG/src

# Include setup
include(PROJECT_ROOT * "/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/setup/simulations_setup.jl")

# Include modules
include(PROJECT_ROOT * "/modules/fits.jl")
include(PROJECT_ROOT * "/modules/plots.jl")

using DelimitedFiles
using Dates

function main()    

    ModeErrorMsg = "Input error: use option --heatmap, --boundaries, --gamma,
    or --gamma-MISF"
	
	if length(ARGS) != 1
		# If user does not specify the user mode
		error(ModeErrorMsg)
		exit()
	else
        UserMode = ARGS[1]

        # ----------------------------------------------------------------------
        # ---------------------------- Plot Heatmap ----------------------------
        # ----------------------------------------------------------------------

        if UserMode=="--heatmap"

            global HorizontalLL, RectangularLL # Imported from setup
            
			PhaseBoundariesLL = HorizontalLL
            LL = RectangularLL               

            μ0 = 0.0 # TODO CHANGE!

            PhaseBoundariesFilePath = PROJECT_ROOT * "/../simulations/horizontal_sweep/μ0=$(μ0)_L=$PhaseBoundariesLL.txt"
            HeatmapDir = PROJECT_ROOT * "/../analysis/heatmap/"

			for L in LL
				# TODO Do we need today()?
                FilePathIn = PROJECT_ROOT * "/../simulations/rectangular_sweep/L=$(L)_site=$(ceil(Int64, L/2)).txt"

				VarianceFilePathOut = PROJECT_ROOT * "/../analysis/heatmap/variance_L=$(L)_$(today()).pdf" # Variance plot
	            AFilePathOut = HeatmapDir * "a_L=$(L)_$(today()).pdf"       # <a_i> plot
    	        KFilePathOut = HeatmapDir * "K_L=$(L)_$(today()).pdf"       # K plot
				
	            PlotHeatmap(L, FilePathIn; PhaseBoundariesFilePath, VarianceFilePathOut, AFilePathOut)
   			end

        # ----------------------------------------------------------------------
        # --------------------- Boundary between SF and MI ---------------------
        # ----------------------------------------------------------------------

        elseif UserMode=="--boundaries"
        	
			μ0 = 0.0#Horizontalμμ[3] # CHANGE!

            FilePathIn = PROJECT_ROOT * "/../simulations/boundaries_sweep/μ0=$(μ0)_L=$HorizontalLL.txt"
            PhaseBoundariesDir = PROJECT_ROOT * "/../analysis/phase_boundaries/μ0=$(μ0)/"
            mkpath(PhaseBoundariesDir)

            FilePathPlot = PhaseBoundariesDir * "phaseboundaries_μ0=$(μ0).pdf"
            FilePathFit = PhaseBoundariesDir * "fitted_phase_boundaries_μ0=$(μ0).txt"
            
            FilePathPlotOut = PhaseBoundariesDir * "phaseboundaries_fit_μ0=$(μ0).pdf"
            FilePathSinglePlotOut = PhaseBoundariesDir * "phaseboundaries_fit_single_μ0=$(μ0).pdf"

            #PlotPhaseBoundaries(FilePathIn; gap=false, FilePathOut=FilePathPlot, μ0)
            #FitPhaseBoundaries(FilePathIn, FilePathFit; FilePathPlotOut, FilePathSinglePlotOut, μ0)
            PlotPhaseBoundaries(FilePathIn;
                FilePathOut=PhaseBoundariesDir*"colored_PB.pdf",
                HideDataPoints=true,
                DrawMottLobe=true,
                MottLobeFilePath=PROJECT_ROOT * "/../analysis/phase_boundaries/μ0=$μ0/fitted_phase_boundaries_μ0=$μ0.txt",
            )

        # ----------------------------------------------------------------------
        # ----------------------- Correlation function Γ -----------------------
        # ----------------------------------------------------------------------

        elseif UserMode=="--gamma"
        
            global HorizontalLL, RectangularLL # Imported from setup

            μ0 = 0.0 #Horizontalμμ[2] # CHANGE!
            rrMin = [2, 4, 6, 8]
            rrMax = [12, 14, 16, 18]
            # FitRange = 10
            JMin = 0.20
            LMin = 20

            GammaDir = PROJECT_ROOT * "/../analysis/gamma/μ0=$(μ0)/"
            mkpath(GammaDir)

            FilePathIn = PROJECT_ROOT * "/../simulations/horizontal_sweep/μ0=$(μ0)_L=$HorizontalLL.txt"
            FilePathFit = GammaDir * "fitted_Luttinger_parameter_μ0=$μ0.txt"

            # (Step 1) Perform all the fits for all J>J_min
            FitRoutineGamma(FilePathIn, FilePathFit; rrMin, rrMax, JMin, LMin)

            # (Step 2) Plot Γ(r) vs r for one  given (J, μ0) ( j ∈ [1, 50] )
            # j = 43 corresponds to J=0.30
            # !!! NOTE: already plotted for all the μ0 we have. !!!
            # for j in 12:2:50
            #     mkpath(GammaDir*"data_plot/")
            #     FileGammaPlot = GammaDir * "data_plot/data_gamma_j=$(j)_μ0=$μ0.pdf"
            #     PlotPowerLawGamma(FilePathIn, μ0, j; FilePathOut=FileGammaPlot, overwrite=false)
            # end

            # (Step 3) Plot K_∞ vs J, reading data from Step 1
            FilePathInK = GammaDir * "fitted_Luttinger_parameter_μ0=$μ0.txt"
            FilePathOutKInfty = GammaDir * "fitted_Luttinger_plot_μ0=$μ0.pdf"
            PlotFitResultsK(FilePathInK, μ0; FilePathOut=FilePathOutKInfty)
        
        elseif UserMode=="--gamma-MISF"
            global HorizontalLL
            μ0 = 0.6

            FilePathIn = PROJECT_ROOT * "/../simulations/horizontal_sweep/μ0=$(μ0)_L=$HorizontalLL.txt"
            PhaseBoundariesFilePath = PROJECT_ROOT * "/../simulations/boundaries_sweep/μ0=0.0_L=$HorizontalLL.txt"

            GammaDir = PROJECT_ROOT * "/../analysis/gamma/illustrative_plots/"
            mkpath(GammaDir)

            FilePathOut1 = GammaDir * "gamma_MI_vs_SF.pdf"
            FilePathOut2 = GammaDir * "gamma_MI_vs_SF_phaseboundaries.pdf"
            PlotCorrelationFunctionsMIvsSF(FilePathIn; μ0, JMin=0.10, JMax=0.15, L=70, PhaseBoundariesFilePath,
                FilePathOut1, FilePathOut2)
        
		else
			error(ModeErrorMsg)
			exit()
        end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
