module DDCModelsExample

# Types and functions package exports 
export DDCbasic 
export gen_data, estimateΠ
export estimate_model, calc_stderr

# Load additional packages  
using Parameters 
using LinearAlgebra, Distributions, Random
using GalacticOptim, Optim, ForwardDiff 

# Include  the source files 
include("types.jl")
include("support.jl")
include("models/DDCbasic.jl")
end
