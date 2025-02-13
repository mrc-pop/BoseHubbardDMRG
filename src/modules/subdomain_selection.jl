#!/usr/bin/julia

include(PROJECT_ROOT * "/setup/graphic_setup.jl")

# ------------------------------- Tip finding ----------------------------------

function FindMottTip(FilePathIn::String;
					 gap=true,
    				 verbose=false,
    				 μ0=0.0)
    				   
	println("Performing linear extrapolation of gap closing point...")
    				   
	# Import phase boundaries 
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
	
	for (i, Index) in enumerate(CrossingIndices)

		# Linear extrapolation
		x1 = JJ[Index]
		y1 = Gap[Index]
		x2 = JJ[Index+1]
		y2 = Gap[Index+1]
		
		a = (y1-y2)/(x1-x2)
		b = ( (y1+y2)-a*(x1+x2) )/2
		
		x3 = -b/a
		
		# TODO Change steps
		ReducedXStep = 5*(x2-x1)/2								# Arbitrary
		Left = round(x3-ReducedXStep, digits=3)
		Right = round(x3+ReducedXStep, digits=3)
		
		ReducedYStep = 2*Gap[Index-10]
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

# ------------------------ Selection plot for check ----------------------------

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
