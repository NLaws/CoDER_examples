include("CoDER_examples.jl")
using Test
using Gurobi

T = 24
LDFinputs, bigM, smlM = get_LDFinputs(T; v_lolim=0.0);


## check UL power flow feasibility
# model = Model(Gurobi.Optimizer)
# LDF.build_ldf!(model, LDFinputs)
# optimize!(model)

#= these voltages should be the same after changing Sbase, which they are.
julia> maximum(sqrt.(value.(model[:vsqrd])))
1.008374586764746

julia> minimum(sqrt.(value.(model[:vsqrd])))
0.9845817662226841
=#

model_linearized = linearized_problem_bess_bigM(Gurobi.Optimizer, T, LDFinputs, bigM, smlM);

model_bileveljump = bileveljump_bess_bigM(Gurobi.Optimizer, T, LDFinputs, bigM, smlM);

@test objective_value(model_linearized) ≈ objective_value(model_bileveljump)


T = 8760
LDFinputs, bigM, smlM = get_LDFinputs(T; v_lolim=0.0);

model_linearized = linearized_problem_bess_bigM(Gurobi.Optimizer, T, LDFinputs, bigM, smlM);

 # TODO bileveljump_bess_bigM is taking too much time to build problem due to loops over variables and constraints (including before linearization process, while making complementary constraints in BilevelJuMP.jl)

 model_bileveljump = bileveljump_bess_bigM(Gurobi.Optimizer, T, LDFinputs, bigM, smlM);

 @test objective_value(model_linearized) ≈ objective_value(model_bileveljump)