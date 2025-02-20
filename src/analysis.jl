#!/usr/bin/julia

PROJECT_ROOT = @__DIR__ # Absloute path up to .../BoseHubbardDMRG/src

# Include setup
include(PROJECT_ROOT * "/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/setup/simulations_setup.jl")

# Include modules
include(PROJECT_ROOT * "/modules/fits.jl")
include(PROJECT_ROOT * "/modules/plots.jl")
include(PROJECT_ROOT * "/modules/subdomain_selection.jl")

using DelimitedFiles
using Dates

function main()    

    ModeErrorMsg = "Input error: use option --heatmap, --boundaries, --gamma, gamma-K
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
            L = 10 # TODO CHANGE!             
            μ0 = 0.0 # TODO CHANGE!
            Cutoff = 0.05 # TODO CHANGE!

            PhaseBoundariesFilePath = PROJECT_ROOT * "/../simulations/horizontal_sweep/μ0=$(μ0)_L=$PhaseBoundariesLL.txt"
            HeatmapDir = PROJECT_ROOT * "/../analysis/heatmap/"


			for L in LL
                FilePathIn = PROJECT_ROOT * "/../simulations/rectangular_sweep/L=$(L)_High.txt"

				VarianceFilePathOut = HeatmapDir * "variance_L=$(L)_High.pdf" # Variance plot
	            AFilePathOut = HeatmapDir * "b_L=$(L)_High.pdf"       # <a_i> plot
    	        KFilePathOut = HeatmapDir * "K_L=$(L)_High.pdf"       # K plot
				
	            PlotHeatmap(L, FilePathIn; PhaseBoundariesFilePath, VarianceFilePathOut, AFilePathOut, KFilePathOut)
   			end
            
            FilePathIn = PROJECT_ROOT * "/../simulations/rectangular_sweep/L=$(L)_High.txt"

            EstimateFractionOrderParameter(L, FilePathIn, PhaseBoundariesFilePath, Cutoff)


        # ----------------------------------------------------------------------
        # --------------------- Boundary between SF and MI ---------------------
        # ----------------------------------------------------------------------

        elseif UserMode=="--boundaries"
        	
			μ0 = 0.0 #Horizontalμμ[3] # CHANGE!

            FilePathIn = PROJECT_ROOT * "/../simulations/boundaries_sweep/μ0=$(μ0)_L=$HorizontalLL.txt"
            PhaseBoundariesDir = PROJECT_ROOT * "/../analysis/phase_boundaries/μ0=$(μ0)/"
            mkpath(PhaseBoundariesDir)

            FilePathPlot = PhaseBoundariesDir * "phaseboundaries_μ0=$(μ0).pdf"
            FilePathFit = PhaseBoundariesDir * "fitted_phase_boundaries_μ0=$(μ0).txt"
            
            FilePathPlotOut = PhaseBoundariesDir * "phaseboundaries_fit_μ0=$(μ0).pdf"
            FilePathSinglePlotOut = PhaseBoundariesDir * "phaseboundaries_fit_single_μ0=$(μ0).pdf"

            # Uncomment the needed analysis.
            #PlotPhaseBoundaries(FilePathIn; gap=false, FilePathOut=FilePathPlot, μ0)
            FitPhaseBoundaries(FilePathIn, FilePathFit; FilePathPlotOut, FilePathSinglePlotOut, μ0)
            PlotPhaseBoundaries(FilePathIn;
                FilePathOut=PhaseBoundariesDir*"colored_PB.pdf",
                HideDataPoints=true,
                DrawMottLobe=true,
                MottLobeFilePath=PROJECT_ROOT * "/../analysis/phase_boundaries/μ0=$μ0/fitted_phase_boundaries_μ0=$μ0.txt",
            )

            FindMottTip(FilePathFit, verbose=true)

        elseif UserMode=="--boundaries-point"
            μ0 = 0.0
            MottLobeFilePath=PROJECT_ROOT * "/../analysis/phase_boundaries/μ0=$μ0/fitted_phase_boundaries_μ0=$μ0.txt"
            PlotPointAndPhaseBoundaries(MottLobeFilePath)

        # ----------------------------------------------------------------------
        # ----------------------- Correlation function Γ -----------------------
        # ----------------------------------------------------------------------

        elseif UserMode=="--gamma"
        
            global HorizontalLL, RectangularLL # Imported from setup

            μ0 = 0.08 # CHANGE!
            rrMin = [2, 4, 6, 8]
            rrMax = [12, 14, 16, 18]
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
            # ! NOTE: already plotted for all the μ0 we have. !
            # for j in 12:2:50
            #     mkpath(GammaDir*"data_plot/")
            #     FileGammaPlot = GammaDir * "data_plot/data_gamma_j=$(j)_μ0=$μ0.pdf"
            #     PlotPowerLawGamma(FilePathIn, μ0, j; FilePathOut=FileGammaPlot, overwrite=false)
            # end

            # (Step 3) Plot K_∞ vs J, reading data from Step 1
            FilePathInK = GammaDir * "fitted_Luttinger_parameter_μ0=$μ0.txt"
            for rMin in [2,4,6,8]
                FilePathOutKInfty = GammaDir * "fitted_Luttinger_plot_μ0=$(μ0)_rMin=$rMin.pdf"
                PlotFitResultsK(FilePathInK, μ0, rMin; FilePathOut=FilePathOutKInfty)
            end
        
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
        
        elseif UserMode=="--gamma-K"
            # Make an illustrative plot of K_∞ vs J, for one
            # rMin and one rMax, chosen with the definitions below.

            global HorizontalLL
            μ0 = 0.0
            rMin = 8 
            rMax = 18

            GammaDir = PROJECT_ROOT * "/../analysis/gamma/μ0=$(μ0)/"
            FilePathInK = GammaDir * "fitted_Luttinger_parameter_μ0=$μ0.txt"
            FilePathOutKInfty = GammaDir * "illustrative_K_μ0=$μ0.pdf"
            PlotIllustrativeResultK(FilePathInK, μ0, rMin, rMax; FilePathOut=FilePathOutKInfty)

		else
			error(ModeErrorMsg)
			exit()
        end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
