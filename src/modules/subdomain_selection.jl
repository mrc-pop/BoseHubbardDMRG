#!/usr/bin/julia

# ------------------------------- Tip finding ----------------------------------

function FindMottTip(FilePathIn::String;
					 gap=true,
    				 verbose=false,
    				 μ0=0.0)
    				   
	println("Performing linear interpolation of gap closing point...")
    				   
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

		# Linear interpolation
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

# ----------------------------- 1/2 K crossing ---------------------------------

function FindCriticalK(FilePathIn::String;
					   verbose=false)
	
	println("Performing linear interpolation of K parameter crossing 1/2...")
	
	KData = readdlm(FilePathIn, ',', '\n'; comments=true)
	JJ = KData[:,1]
	
	mm = KData[:,2]
	MM = KData[:,3]
	rrMin = unique(mm)
	rrMax = unique(MM)
	
	KKtmp = KData[:,4]
	eKKtmp = KData[:,5]
	CriticalPoints = zeros(length(rrMin),4)
	
	for (i,rMin) in enumerate(rrMin)
		rMax = rrMax[i]
		CriticalPoints[i,1:2] = [rMin rMax]
		KK = KKtmp[i:length(rrMin):end]
		eKK = eKKtmp[i:length(rrMin):end]		
		
		Index=0
		for j in 1:(length(KK)-1)
			if ((KK[j]-0.5)*(KK[j+1]-0.5))<0
				Index = j
			end
		end

		# Linear interpolation with redundant error
		x3 = zeros(2)
		for k in 1:2
			x1 = JJ[(i-1) + length(rrMin)*Index]
			y1 = KK[Index] + (-1)^k * eKK[Index]
			x2 = JJ[(i-1) + length(rrMin)*(Index+1)]
			y2 = KK[Index+1] + (-1)^k * eKK[Index+1]

			a = (y1-y2)/(x1-x2)
			b = ( (y1+y2)-a*(x1+x2) )/2

			x3[k] = (0.5-b)/a

			if verbose
				println("\nAnalysis of data points for rmin=$(Int64(rMin)), rmax=$(Int64(rMax))")
				if k==1
					println("First interpolation, enhanced tilting.")
				elseif k==2
					println("Second interpolation, lowered tilting.")
				end
				println("Point above: ", round(KK[Index], digits=3))
				println("Point below: ", round(KK[Index+1], digits=3))
			end
		end
		
		Jc = sum(x3)/2
		eJc = diff(x3)[1]/2
		println("Gap closes at J: ", round(Jc, digits=5), " +/- ", round(eJc, digits=5))
		
		CriticalPoints[i,3:4] = [Jc eJc]
	end
	
	return CriticalPoints
end
