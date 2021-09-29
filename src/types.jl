#############################################################
# Define types
#############################################################
abstract type Model end  # General type for any model

# DDCbasic model type as described by Abbring and Klein 
# Type gathers model parameters, and offers additional options (see below )
# Uses Parameters.jl, to define default values  

@with_kw mutable struct DDCbasic <: Model
	supportX::Vector 		# Support points of profit state X_t
	Π::Matrix				# Transition matrix between profit states 
	β::Vector				# Profit parameters
	δ::Vector				# Entry/exit cost parameters 
	ρ::Float64 				# discount factor
	options = Options()		# Additional options 
end

# Options structure gathering additional model/estimation options
@with_kw mutable struct Options 
    # tolerance level for inner fixed point
	tol_fixed_point::Float64 = 1e-10    
    # seed for data generation, and simulation
	seed::Int = 42                  
    # Optimization algorithm (default LBFGS() from Optim.jl)
    # For alternatives see GalacticOptim.jl, which is used to 
    # implement optimization routine 
	algo = LBFGS()  
    # Automatic differentiation package (see GalacticOptim.jl for options)
	AD = GalacticOptim.AutoForwardDiff() 
end