"""
	flowpayoffs!(u0,u1,m::DDCbasic)

Updates flowpayoffs in matrices `u0` and `u1` based on parameters in model `m`. 
"""
function flowpayoffs!(u0,u1,m::DDCbasic)
	# Fill in u0 
	u0[:,2] .= -m.δ[1] 
	# Fill in u1 
	u1 .= m.β[1] .+ m.β[2] .* m.supportX 
	u1[:,1] .-= m.δ[2] 
	return nothing 
end

# NOTES:
# - functions with "!" (e.g. f!(a)) mutate the input argument. The "!" 
#	is not required, but it's convention to use it for functions that 
# 	mutate its inputs. 
# - Using mutating functions in inner functions that are run many times 
#	avoids allocations, which can lead to efficiency gains.  
# 	(see https://docs.julialang.org/en/v1/manual/performance-tips/)

""" 
	bellman!(U0,U1,u0,u1,m::DDCbasic) 

Bellman update, i.e. update expected payoffs in `U0` and `U1` given flowpayoffs `u0` and `u1`, discount factor `m.ρ` and transition matrix `m.Π`.
"""
function bellman!(U0,U1,u0,u1,m::DDCbasic) 
	r0 = log.(exp.(U0[:,1])+exp.(U1[:,1]))
	U0 .= u0 .+ m.ρ .* m.Π * repeat(r0,1,2) 

	r1 = log.(exp.(U0[:,2])+exp.(U1[:,2]))
	U1 .= u1 .+ m.ρ .* m.Π * repeat(r1,1,2)  
	return nothing 
end


"""
	solve_fixedpoint!(U0,U1,u0,u1,m::DDCbasic)

Solve fixed point and update expected payoffs in `U0` and `U1`, given flowpayoffs `u0` and `u1`, and parameters in model `m`. 
"""
function solve_fixedpoint!(U0,U1,u0,u1,m::DDCbasic)
	tol = m.options.tol_fixed_point 
	_U0 = U0 .+ 2 .* tol 
	_U1 = U1 .+ 2 .* tol 
	while maximum(abs.(_U0 .- U0)) > tol || maximum(abs.(_U1 .- U1)) > tol
		_U0 .= U0
		_U1 .= U1 
		bellman!(U0,U1,u0,u1,m::DDCbasic) 
	end
	return nothing 
end

"""
	init_flowpayoffs(m::DDCbasic)

Initialize flowpayoff matrices `u0` and `u1` based on model `m`. 
"""
function init_flowpayoffs(m::DDCbasic)
	n = length(m.supportX)
	u0 = zeros(eltype(m.β),n,2)
	u1 = zeros(eltype(m.β),n,2) 
	flowpayoffs!(u0,u1,m)
	return u0,u1
end

"""
	calc_ΔU(m::DDCbasic)

Calculate and return matrix `ΔU` for model `m`. 
"""
function calc_ΔU(m::DDCbasic)

	u0,u1 = init_flowpayoffs(m)

	U0 = zeros(eltype(u0),length(m.supportX),2)
	U1 = copy(U0) 
	solve_fixedpoint!(U0,U1,u0,u1,m)
	ΔU = U1 .- U0 
	return ΔU 
end

"""
	calc_pinf(Π) 

Calculate and return P∞ from transition matrix `Π`. 
"""
function calc_pinf(Π) 
	nk = size(Π,1)
	A = vcat(Matrix(I,nk,nk)[1:nk-1,:] .- Π[:,1:nk-1]',ones(1,nk))
	return A \ vcat(zeros(nk-1),1) 
end


"""
	loglikelihood(m::DDCbasic,X,choices)

Calculates and returns loglikelihood for model `m`, given observed payoff states `X` and choices `choices`. 
"""
function loglikelihood(m::DDCbasic,X,choices)
	nperiod,nfirm = size(choices)
	ΔU = calc_ΔU(m::DDCbasic)
	pexit = 1 ./ (1 .+ exp.(ΔU))
	ll = zeros(eltype(ΔU),Threads.nthreads())
	Threads.@threads for threadid in 1:Threads.nthreads() 
		for f in getrange(nfirm)
			# Period 1 
			ll[threadid] += log(choices[1,f] + (1-2*choices[1,f]) * 
								pexit[getindX(m.supportX,X[1,f]),1])
			# Loop through other periods 
			for t in 2:nperiod 
				ll[threadid] += log(choices[t,f] + (1-2*choices[t,f]) * 
									pexit[getindX(m.supportX,X[t,f]),
											choices[t-1,f]+1])
			end
		end
	end
	return sum(ll)
end

"""
	gen_data(m::DDCbasic,nperiod,nfirm)

Generates and returns data on payoff states and firms choices `X,choices` for model `m`, `nperiod` periods and `nfirm` firms. 
"""
function gen_data(m::DDCbasic,nperiod,nfirm)
	P∞ = calc_pinf(m.Π)
	dε = GeneralizedExtremeValue(0,1,0)
	ΔU = calc_ΔU(m::DDCbasic)
	
	X = zeros(Int,nperiod,nfirm)

	Δε = zeros(nperiod,nfirm)
	choices = zeros(Bool,nperiod,nfirm)

	dX0 = DiscreteNonParametric(m.supportX,P∞) 
	dX = [DiscreteNonParametric(m.supportX,m.Π[i,:]) 
			for i in eachindex(m.supportX)]

	Threads.@threads for threadid in 1:Threads.nthreads() 
		Random.seed!(m.options.seed+threadid)
		for f in getrange(size(X,2))
			# Period 1 
			X[1,f] = rand(dX0) 
			Δε[1,f] = rand(dε) - rand(dε) 
			choices[1,f] = ΔU[getindX(m.supportX,X[1,f]),1] >
								Δε[1,f]		
			# Loop through other periods 
			for t in 2:size(X,1)
				X[t,f] = rand(dX[getindX(m.supportX,X[t-1,f])]) 
				Δε[t,f] = rand(dε) - rand(dε) 
				choices[t,f] = ΔU[getindX(m.supportX,X[t,f]),
									choices[t-1,f] + 1 ] > Δε[t,f]		
			end

		end
	end

	return X,choices 
end


"""
	estimate_model(m::DDCbasic,X,choices; 
						startvals = zeros(3),
						estimateΠ=false) 

Estimates parameters in model `m`, given observed payoff states `X` and choices `choices`. Returns new model with estimate parameter values. 

## Optional keyword arguments 
`startvals`: Set 3x1 vector of starting values
`estimateΠ`: Currently not implemented, would first estimate `Π` and then use that in estimation. 
"""
function estimate_model(m::DDCbasic,X,choices; 
							startvals = zeros(3),
							estimateΠ=false) 

	# For parameter vector θ, return loglikelihood by creating new model struct 
	# with these values 
	function f(θ,m::DDCbasic,X,choices) 
		# Get parameters from parameter vector 
		β = θ[1:2]
		δ = [m.δ[1];θ[3]]
		# Construct model with parameters 
		_m = typeof(m)(;supportX=m.supportX,Π=m.Π,β,δ,ρ=m.ρ)
		# Return neg loglikelihood 
		return - loglikelihood(_m,X,choices)
	end
	
	# Wrap into objective function. Note, p is unused, 
	# but required by GalacticOptim.jl
	function f(θ,p) 
		f(θ,m,X,choices)
	end

	# Define objective function and problem for GalactiOptim.jl
	of = OptimizationFunction(f,m.options.AD)        
	prob = OptimizationProblem(of,startvals)

	# Find optimum 
	sol = solve(prob,m.options.algo) 

	# Put together in new estimated model 
	mout = deepcopy(m) 
	mout.β .= sol.u[1:2]
	mout.δ[2] = sol.u[3] 

	return mout 
end

"""
	calc_stderr(m::DDCbasic,X,choices)

For model `m`, calculate standard errors for parameters given data `X,choices`. 
"""
function calc_stderr(m::DDCbasic,X,choices)

	# Define function for parameter values in vector θ
	function f(θ,m::DDCbasic,X,choices) 
		β = θ[1:2]
		δ = [m.δ[1];θ[3]]
		_m = typeof(m)(;supportX=m.supportX,Π=m.Π,β,δ,ρ=m.ρ)
		return - loglikelihood(_m,X,choices)
	end

	# Wrap into function depending only on parameter vector θ	
	function f(θ) 
		f(θ,m,X,choices)
	end

	# Return standard errors (diagonal elements of 
	# inverse hessian are variance estimates)
   sqrt.(diag(inv(ForwardDiff.hessian(f,vcat(m.β,m.δ[2])))))
end