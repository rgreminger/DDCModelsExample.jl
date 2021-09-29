using DDCModelsExample
using Test

@testset "DDCbasic.jl" begin
        
    # Define number of periods and number of firms 
    nperiod = 100
    nfirm = 10000
   
    # Define support of X 
    supportX = collect(1:5) 
    nsupportX = length(supportX)
    
    # Define Π
    Π = 1 ./ (1 .+ abs.(repeat(supportX',nsupportX,1) .- repeat(supportX,1,nsupportX)))
    Π ./= sum(Π;dims=2)

    # Define parameters 
    β = [-0.1*nsupportX;0.2]
    δ = [0.0;1.0] 
    ρ = 0.95 

    # Construct model 
    m = DDCbasic(;supportX,Π,β,δ,ρ) 

    # Simulate dataset 
    X,choices = gen_data(m,nperiod,nfirm) 

    # Estimate Π
    Πestim = estimateΠ(X,m.supportX)

    # Test whether estimated Π is approximately the same as 
    # true Π with which data was simulated 
    @test maximum(abs.(Πestim .- m.Π)) < 0.01

    # Estimate model 
    mestim = estimate_model(m,X,choices; startvals = [-1;-0.1;0.5])
    # Calculate standard errors
    stderr = calc_stderr(mestim,X,choices)
    
    # Test whether coefficient estimates within 2 std 
    # errors of true value 
    @test maximum(abs.(m.β .- mestim.β) .< 2 .* stderr[1:2])
    @test abs(m.δ[2] - m.δ[2])  < 2*stderr[3] 

end

