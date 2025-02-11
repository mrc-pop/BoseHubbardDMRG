#/!usr/bin/julia

PROJECT_ROOT = @__DIR__ # Absloute path up to .../BoseHubbardDMRG/src
include(PROJECT_ROOT * "/definitions/dmrg.jl")
PROJECT_ROOT *= "/.."

function SelectSubdomain(FilePathIn::String;
						 gap=true,
    				     verboseWaiting=false,
    				     μ0=0.0)
    				   
	println("\nPerforming linear extrapolation of gap closing point...")
    				   
	BoundariesData = readdlm(FilePathIn, ',', '\n'; comments=true)
	JJ = BoundariesData[:,1]
	ΔEp = BoundariesData[:,2]
	ΔEm = BoundariesData[:,3]
	eΔEp = BoundariesData[:,4]
	eΔEm = BoundariesData[:,5]
	Gap = ΔEp .+ ΔEm
	
	CrossingIndices = []
	for i in 1:(length(Gap)-1)
		if Gap[i]*Gap[i+1]<0
			push!(CrossingIndices, i)
		end
	end
	
	Selections = zeros(Float64, length(CrossingIndices), 4)
	
	for i in 1:length(CrossingIndices)
		Index = CrossingIndices[i]
		
		# Linear extrapolation
		x1 = JJ[Index]
		y1 = Gap[Index]
		x2 = JJ[Index+1]
		y2 = Gap[Index+1]
		
		a = (y1-y2)/(x1-x2)
		b = ( (y1+y2)-a*(x1+x2) )/2
		
		x3 = -b/a
		
		ReducedXStep = (x2-x1)/2								# Arbitrary
		Left = round(x3-ReducedXStep, digits=3)
		Right = round(x3+ReducedXStep, digits=3)
		
		ReducedYStep = Gap[Index-10]
		Up = round(ΔEp[Index]+ReducedYStep, digits=3)
		Down = round(ΔEp[Index]-ReducedYStep, digits=3)	# Arbitrary
		
		if verbose
			println("Number of sign changing points: ", length(CrossingIndices))
			println("Point above: ", Gap[Index])
			println("Point below; ", Gap[Index+1])
			println("Gap closes at J: ", x3)
			println("Restrict simulation to $Left≤J≤$Right, $Down≤μ≤$Up")
		end
		
		Selections[i,:] = [Left Right Up Down]
	end
	
	return Selections
end

function PlotSelection(FilePathIn::String,
					   Selections::Matrix{Float64})
	
	BoundariesData = readdlm(FilePathIn, ',', '\n'; comments=true)
	JJ = BoundariesData[:,1]
	ΔEp = BoundariesData[:,2]
	ΔEm = BoundariesData[:,3]
	
	plot()
    plot!(JJ,
          [ΔEp, -ΔEm],
          label=[L"\mu_c^+ \, (L \rightarrow \infty)" L"\mu_c^- \, (L \rightarrow \infty)"], 
          xlabel=L"J", ylabel=L"$\mu$", 
          title="Fitted phase boundaries",
          alpha=1.0)
     
    rectangle(l, r, u, d) = Shape(l .+ [0,r-l,r-l,0], d .+ [0,0,u-d,u-d])
    for i in 1:size(Selections,1)
    	Left, Right, Up, Down = Selections[i,:] 
    	plot!(rectangle(Left,Right,Up,Down),
    		  label="Selection",
    		  linewidth=0,
    		  opacity=0.2)
	end
	gui()		# .pdf file saved on /tmp, erased on boot
end

"""
function GetKIntercept(FilePathIn::String)
end
"""

println("\nBeta: no errors")
Selections = SelectSubdomain(PROJECT_ROOT * "/analysis/phase_boundaries/fitted_phase_boundaries.txt"; verbose=false)
PlotSelection(PROJECT_ROOT * "/analysis/phase_boundaries/fitted_phase_boundaries.txt", Selections)

print("Should I run a rectangular simulation inside the selected region? (y/n) ")
Waiting=true
PerformTipSweep = readline()
while Waiting
	if PerformTipSweep=="y"
		global Waiting=false
		println("Selection accepted. Starting simulations around Mott lobe tip...")
	elseif PerformTipSweep=="n"
		global Waiting=false
		println("Selection rejected. Change the setting of this program to improve selection.")
		exit()
	else
		global Waiting=true
		print("Invalid input. Please use valid input. (y/n) ")
		global PerformTipSweep = readline()
	end
end

display(PerformTipSweep)

Left, Right, Up, Down = Selections[i,:]
JJ = collect(range(start=Left, stop=Right, length=50))
μμ = collect(range(start=Down, stop=Right, length=50))
LL = [10, 20, 30, 40, 50, 60, 70]

DirPathOut = PROJECT_ROOT * "/simulations/rectangular_sweep_tip"
mkpath(DirPathOut)
FilePathOut = DirPathOut * "/L=$LL.txt"

DataFile = open(FilePathOut,"w")
write(DataFile,"# Hubbard model DMRG. This file contains many sizes. nmax=$nmax, μ0=$μ0, nsweeps=$nsweeps, cutoff=$cutoff\n")
write(DataFile,"# L; J; E; nVariance; C; eC [calculated $(now())]\n")
close(DataFile)

for L in LL
	println("Starting calculation of observables for L=$L...")
	RectangularSweep(i, L, N, nmax, JJ, μμ, DMRGParametersMI, DMRGParametersSF, FilePathIn, FilePathOut)
end

println("Done!")
