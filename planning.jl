#=
The planning example assumes that:
- there is at least one over loaded line in the network
- the power flow is approximated with the LinDistFlow model (Baran & Wu 1989)
- the system owner is considering battery energy storage, compensating customers for exported power, and upgrading the overloaded components
- the system owner pays the 2019 ERCOT average hourly real-time market price for imports at the feeder head
- customers have the option to purchase PV systems to reduce their costs
- there is a large refrigerated warehouse in the network with a price responsive HVAC system
- the base retail tariff is TODO
- the network model is based on the 38 node model from Andrianesis et al. 2019
- location is Austin, TX
=#
using Gurobi
using BilevelJuMP, JuMP
using REoptLite
import MathOptInterface
import LinDistFlow as LDF  # REoptLite exports LinDistFlow, maybe should remove that
using DelimitedFiles
using Logging
include("extend_lindistflow.jl")
global_logger(ConsoleLogger(stderr, Logging.Debug))
const MOI = MathOptInterface


#=
        PARAMETERS
=#

# scalars
T = 24 # time horizon
lat, lon = 30.2672, -97.7431  # Austin, TX
cpv = 1500
cbkW = 700
cbkWh = 350
η = 0.95
years = 10
discountrate = 0.05
M = 10_000  # default bound for decision variables
CHILLER_COP = 4.55

# arrays
LDFinputs = LDF.singlephase38linesInputs(Sbase=1);  # TODO this method is not released yet
# NOTE Sbase is 1, using absolute power units in model, kW to align with DoE profiles
loadnodes = collect(keys(LDFinputs.Pload))
# NOTE LDF using strings for node keys to allow for arbitrary names (should make them symbols)

# relax voltage limits while setting up scenario
# LDFinputs.v_uplim = 1.2
# LDFinputs.v_lolim = 0.8
# TODO current (squared) limits?


LLnodes = ["33", "34"]  # all nodes in LL model (that have decisions)

# TODO implement constraints/variables for varying nodes with PV and BESS
LLnodes_withPV = ["34"]
LLnodes_warehouse = ["33"]




# fill in loads, all nodes have some uncontrolled load, defined using DoE Commercial Reference Buildings
profile_names = ["FastFoodRest", "FullServiceRest", "Hospital", "LargeHotel", "LargeOffice", "MediumOffice", "MidriseApartment", "Outpatient", "PrimarySchool", "RetailStore", "SecondarySchool", "SmallHotel", "SmallOffice", "StripMall", "Supermarket", "Warehouse"]

#= 
    Define known loads, 
    - nodes w/o DER are set in the LDF model using build_ldf! 
    - nodes w/DER use the known load in the LL model load balance constraint
    - in the LDF model, nodes w/DER have decision variables yi and ye (which are complementary)

=#
# knownloads_nodestrings = setdiff(loadnodestrings, LLnodes)
rand_names = rand(profile_names, length(loadnodes))

for (i, node) in enumerate(loadnodes)
    doe_profile = REoptLite.BuiltInElectricLoad("", rand_names[i], lat, lon, 2017)
    LDFinputs.Pload[node] = doe_profile[1:T]
    LDFinputs.Qload[node] = doe_profile[1:T] * 0.1
end

# check power flow feasibility w/o DER
model = Model(Gurobi.Optimizer)
LDF.build_ldf!(model, LDFinputs)
optimize!(model)

minimum(sqrt.(value.(model[:vsqrd])))

# pepper some pv into the system
pvpf = REoptLite.get_pvwatts_prodfactor(lat, lon);  # TODO this function is only in flex branch
pvpf_complement = map(x -> x ≈ 0.0 ? M : 0.0, pvpf)
PVkW = 1e3   # TODO more baseline PV ?
LDFinputs.Pload["3"] .-= PVkW * pvpf[1:T]
LDFinputs.Qload["3"] .-= PVkW * pvpf[1:T] * 0.1
LDFinputs.Pload["30"] .-= PVkW * pvpf[1:T]
LDFinputs.Qload["30"] .-= PVkW * pvpf[1:T] * 0.1
LDFinputs.Pload["18"] .-= PVkW * pvpf[1:T]
LDFinputs.Qload["18"] .-= PVkW * pvpf[1:T] * 0.1

# check power flow feasibility w/baseline PV
model = Model(Gurobi.Optimizer)
LDF.build_ldf!(model, LDFinputs)
optimize!(model)

maximum(sqrt.(value.(model[:vsqrd])))


ci = repeat([0.1], T)
pwf = REoptLite.annuity(years, 0.0, discountrate)
clmp = vec(readdlm("./data/cleaned_ercot2019_RTprices.csv", ',', Float64, '\n'));
clmp = abs.(clmp) / 1e3;  # problem is unbounded with negative prices, convert from $/MWh to $/kWh
tamb = REoptLite.get_ambient_temperature(lat, lon);


# instantiate optimizer and empty model
optimizer = Gurobi.Optimizer()
model = BilevelModel(()->optimizer, linearize_bilinear_upper_terms=true)

#= 
        LOWER MODEL VARIABLES
=#

# TODO add more node indices to LL
@variable(Lower(model), 0 <= ye[LLnodes_withPV, 1:T] <= M) # TODO? start values from REopt
@variable(Lower(model), 0 <= yi[LLnodes_withPV, 1:T] <= M)
@variable(Lower(model), 0 <= ypv[LLnodes_withPV] <= M)
@variable(Lower(model), 0 <= ypvprod[LLnodes_withPV, 1:T] <= M)
@variable(Lower(model), 0 <= spvprod[LLnodes_withPV, 1:T] <= M) 
@variable(Lower(model), 0 <= dummyslack[LLnodes_withPV, 1:T] <= M) 
# @variable(Lower(model), -10 <= dvTemperature[LLnodes_warehouse, 1:T] <= 0.0)
# @variable(Lower(model), 0 <= dvThermalProduction[LLnodes_warehouse, 1:T] <= M)

# LL load balances
@constraint(Lower(model), loadbalance_withPV[n in LLnodes_withPV, t in 1:T],
    yi[n, t] + ypvprod[n, t] == LDFinputs.Pload[n][t] + ye[n, t]
)
# @constraint(Lower(model), loadbalance_warehouse[n in LLnodes_warehouse, t in 1:T],
#     yi[n, t] == LDFinputs.Pload[n][t] + dvThermalProduction[n, t] / CHILLER_COP
# )

# PV production
@constraint(Lower(model), pvprod[n in LLnodes_withPV, t in 1:T],
    ypvprod[n, t] + spvprod[n, t] == ypv[n] * pvpf[t]
)
@constraint(Lower(model), [n in LLnodes_withPV, t in 1:T],
    ypvprod[n, t] + dummyslack[n, t] == ypv[n] * pvpf_complement[t]
)

#=
        LL REFRIGERATED WAREHOUSE MODEL
=#
# R = 0.00025  # K/kW
# C = 1e5   # kJ/K
# A = reshape([-1/(R*C)], 1,1)
# B = [1/(R*C) 1/C]
# u = [tamb zeros(8760)]';  # could replace the zeros vector with endogenous heat input
# J = size(B,2)
# # state space temperature evolution
# @constraint(Lower(model), sstemperature[n in LLnodes_warehouse, t in 2:T],
#     dvTemperature[n, t] == dvTemperature[n, t-1] + A[1, 1] * dvTemperature[n, t-1] + 
#     sum(B[1, j] * u[j, t-1] for j=1:J) + B[1, 1] * (-dvThermalProduction[n, t-1])
# )
# @constraint(Lower(model), init_temperature[n in LLnodes_warehouse], dvTemperature[n, 1] == -1.0)  # initial temperature


#= 
        UPPER MODEL VARIABLES
=#

# duals of LL load balances
@variable(Upper(model), 0 <= lambdaPV <= M, DualOf(loadbalance_withPV))
# @variable(Upper(model), 0 <= lambda_warehouse <= M, DualOf(loadbalance_warehouse))

# UL variables in LL objective
@variable(Upper(model), 0 <= xe[LLnodes_withPV, 1:T] <= M)  # PV export compensation
# @variable(Upper(model), 0 <= xi[LLnodes_warehouse, 1:T] <= M)  # time varying tariff for Warehouses

@variable(Upper(model), 0 <= x0[1:T] <= M)  # positive power import at feeder head

# TODO which nodes can have BESS? Assume the unloaded nodes
LDF.build_ldf!(Upper(model), LDFinputs, LLnodes_withPV, ye, yi)
@constraint(Upper(model), [t in 1:T], x0[t] >= model[:Pⱼ]["0", t] )

# yi ⟂ ye
for (i, e) in zip(yi, ye)
    @constraint(Upper(model),
        [i, e] in MOI.SOS1([1.0, 2.0])
    )
end


#= 
        LL objective
=#
@objective(Lower(model), Min,
    pwf * sum(ci[t] * yi[n, t] for n in LLnodes_withPV, t in 1:T) - 
    # pwf * sum( -xi[n, t] * yi[n, t] for n in LLnodes_warehouse, t in 1:T) + 
    pwf * sum( xe[n, t] * ye[n, t] for n in LLnodes_withPV, t in 1:T) + 
    sum(cpv * ypv[n] for n in LLnodes_withPV)
)

#= 
        UL objective
=#

# TODO voltage cost? wires alternative? 
# TODO UL BESS cost (w/o it pwf has no effect)
@objective(Upper(model), Min,
    pwf * sum( clmp[t] * x0[t] for t in 1:T) +
    # pwf * sum( yi[n, t] * lambda_warehouse[n, t] for n in LLnodes_warehouse, t in 1:T) +
    pwf * sum( ye[n, t] * lambdaPV[n, t] for n in LLnodes_withPV, t in 1:T)
)


optimize!(model)
# add flexHVAC back in then create problems that lead to non-zero ye (xe) and xi


maximum(sqrt.(value.(model[:ye])))


maximum(sqrt.(value.(model[:xi])))