#/!usr/bin/julia

PROJECT_ROOT = @__DIR__ # Absloute path up to .../BoseHubbardDMRG/src
include(PROJECT_ROOT * "/dmrg.jl")
PROJECT_ROOT *= "/.."

function GetClosingGap(FilePathIn::String;
    				   gap=true,
    				   μ0=0.0)
    				   
   println("\nLinear extrapolation of gap closing")
    				   
	BoundariesData = readdlm(FilePathIn, ',', '\n'; comments=true)
	JJ = BoundariesData[:,1]
	ΔEp = BoundariesData[:,2]
	ΔEm = BoundariesData[:,3]
	Gap = ΔEp .+ ΔEm
	
	CrossingIndices = []
	for i in 1:(length(Gap)-1)
		if Gap[i]*Gap[i+1]<0
			push!(CrossingIndices, i)
		end
	end
	println("Number of sign changing points: ", length(CrossingIndices))
	
	for i in 1:length(CrossingIndices)
		Index = CrossingIndices[i]
		println("Point above: ", Gap[Index])
		println("Point below; ", Gap[Index+1])
		
		# Linear extrapolation
		x1 = JJ[Index]
		y1 = Gap[Index]
		x2 = JJ[Index+1]
		y2 = Gap[Index+1]
		
		a = (y1-y2)/(x1-x2)
		b = ( (y1+y2)-a*(x1+x2) )/2
		
		x3 = -b/a
		println("Gap closes at J: ", x3)
	end
end

function GetKIntercept(FilePathIn::String)

	println("\nLinear extrapolation of K crossing 1/2")
    				   
	BoundariesData = readdlm(FilePathIn, ',', '\n'; comments=true)
	JJ = BoundariesData[:,1]
	KK = BoundariesData[:,2]
	
	CrossingIndices = []
	for i in 1:(length(KK)-1)
		if ((KK[i]-0.5)*(KK[i+1]-0.5))<0
			push!(CrossingIndices, i)
		end
	end
	println("Number of points crossing 1/2: ", length(CrossingIndices))
	
	for i in 1:length(CrossingIndices)
		Index = CrossingIndices[i]
		println("Point above: ", KK[Index])
		println("Point above: ", KK[Index+1])
		
		# Linear extrapolation
		x1 = JJ[Index]
		y1 = KK[Index]
		x2 = JJ[Index+1]
		y2 = KK[Index+1]
		
		a = (y1-y2)/(x1-x2)
		b = ( (y1+y2)-a*(x1+x2) )/2
		
		x3 = (0.5-b)/a
		println("K crosses 1/2 at J: ", x3)
	end
end

println("\n Beta: no errors \n")
GetClosingGap(PROJECT_ROOT * "/analysis/phase_boundaries/fitted_phase_boundaries.txt")
GetKIntercept(PROJECT_ROOT * "/analysis/gamma/fit_correlation.txt")
