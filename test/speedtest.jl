using DDCModelsExample

# Define test_speed function 
function test_speed(β,δ,ρ,nperiod,nfirm,supportX,nrepeat)
    nsupportX = length(supportX)

    # Construct Π matrix 
    Π = 1 ./ (1 .+ abs.(repeat(supportX',nsupportX,1) .- repeat(supportX,1,nsupportX)))
    Π ./= sum(Π;dims=2)

    m = DDCbasic(;supportX,Π,β,δ,ρ) 

    t = 0 

    for r in 1:nrepeat 
        m = DDCbasic(;supportX,Π,β,δ,ρ) 
        m.options.seed = 098471 + r 
        X,choices = gen_data(m,nperiod,nfirm) 
        t += @elapsed begin 
            estimate_model(m,X,choices; startvals = [-1;-0.1;0.5])
            calc_stderr(m,X,choices)
        end
    end
    return t 
end

# Define arguments
nperiod = 100
nfirm = 10000
supportX = collect(1:5) 
β = [-0.1*length(supportX);0.2]
δ = [0.0;1.0] 
ρ = 0.95 
nrepeat = 10

# Run first time: includes precompilation time 
tpc = test_speed(β,δ,ρ,nperiod,nfirm,supportX,nrepeat)

# Run second time: no precompilation time 
tnopc = test_speed(β,δ,ρ,nperiod,nfirm,supportX,nrepeat)


        