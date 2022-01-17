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


lat, lon = 30.2672, -97.7431  # Austin, TX
η = 0.95
years = 20
discountrate = 0.05
M = 10_000  # default bound for decision variables
CHILLER_COP = 4.55
#***********************#
Sbase=1e3  # best to use Sbase of 1kW because thermal model assumes kW and kJ
#***********************#

#   cost components
pwf = REoptLite.annuity(years, 0.0, discountrate)
T = 48

clmp = vec(readdlm("./data/cleaned_ercot2019_RTprices.csv", ',', Float64, '\n'));
clmp = abs.(clmp) / 1e3;  # problem is unbounded with negative prices, convert from $/MWh to $/kWh
cpv = 1400
cbkW = 700
cbkWh = 350
ci = repeat([0.25], T);

# outdoor temperature and PV production factor
tamb = REoptLite.get_ambient_temperature(lat, lon);
prod_factor = REoptLite.get_pvwatts_prodfactor(lat, lon);

# power flow inputs
LDFinputs = LDF.singlephase38linesInputs(Sbase=Sbase);
loadnodes = collect(keys(LDFinputs.Pload))
# remove some loadnodes to make problem smaller and keep voltage further from lower limit
# loadnodes = loadnodes[3:end-3]  # removes "32", "12", "25",  "18", "30", "3"
# for n in ["32", "12", "25",  "18", "30", "3"]
#     delete!(LDFinputs.Pload, n)
#     delete!(LDFinputs.Qload, n)
# end
LLnodes_withPV = ["9", "22", "31", "34", "17"]
LLnodes_warehouse = ["5", "10", "15"]  # (price responsive refrigeration), can add any of the loadnodes
LLnodes = union(LLnodes_withPV, LLnodes_warehouse)  # all nodes in LL model (that have decisions)

ULnodes_withBESS = ["2", "7", "24"]

profile_names = ["FastFoodRest", "FullServiceRest", "Hospital", "LargeHotel", "LargeOffice", 
"MediumOffice", "MidriseApartment", "Outpatient", "PrimarySchool", "RetailStore", "SecondarySchool", 
"SmallHotel", "SmallOffice", "StripMall", "Supermarket", "Warehouse"]
doe_profiles = Dict{String, Vector{Float64}}()
for name in profile_names
    doe_profiles[name] = REoptLite.BuiltInElectricLoad("", name, Real(lat), Real(lon), 2017, nothing, Real[])  
end
# TODO address type issues in REoptLite
rand_names = rand(profile_names, length(loadnodes))

# fill in net uncontrolled loads
LDFinputs.Ntimesteps = T
load_scaler = 0.7
for (i, node) in enumerate(loadnodes)
    LDFinputs.Pload[node] = load_scaler * doe_profiles[rand_names[i]][1:T] / Sbase;
    LDFinputs.Qload[node] = load_scaler * doe_profiles[rand_names[i]][1:T] / Sbase * 0.1;
end
# pepper some pv into the system
PVkW = 1e3   # TODO more baseline PV ?
LDFinputs.Pload["5"] .-= PVkW * prod_factor[1:T] / Sbase;
LDFinputs.Qload["5"] .-= PVkW * prod_factor[1:T] / Sbase * 0.1;
LDFinputs.Pload["31"] .-= PVkW * prod_factor[1:T] / Sbase;
LDFinputs.Qload["31"] .-= PVkW * prod_factor[1:T] / Sbase * 0.1;
LDFinputs.Pload["28"] .-= PVkW * prod_factor[1:T] / Sbase;
LDFinputs.Qload["28"] .-= PVkW * prod_factor[1:T] / Sbase * 0.1;
LDFinputs.Pload["27"] .-= PVkW * prod_factor[1:T] / Sbase;
LDFinputs.Qload["27"] .-= PVkW * prod_factor[1:T] / Sbase * 0.1;

LDFinputs.v_lolim = 0.0

peak_load = maximum(sum(values(LDFinputs.Pload)))
peak_single_load = maximum(maximum(values(LDFinputs.Pload)))
LDFinputs.P_up_bound =  peak_load * 100
LDFinputs.Q_up_bound =  peak_load * 10
LDFinputs.P_lo_bound = -peak_single_load * 100
LDFinputs.Q_lo_bound = -peak_single_load * 10

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

optimizer = Gurobi.Optimizer()
model = BilevelModel(()->optimizer, linearize_bilinear_upper_terms=true)
#= 
        LOWER MODEL VARIABLES
=#
@variable(Lower(model), 0 <= ye[LLnodes_withPV, 1:T] <= M);
@variable(Lower(model), 0 <= yi[LLnodes, 1:T] <= M);
@variable(Lower(model), 0 <= ypv[LLnodes_withPV] <= M);
@variable(Lower(model), 0 <= ypvprod[LLnodes_withPV, 1:T] <= M);
@variable(Lower(model), 0 <= spvprod[LLnodes_withPV, 1:T] <= M);

@variable(Lower(model), -10 <= dvTemperature[LLnodes_warehouse, 1:T] <= 0.0);
@variable(Lower(model), 0 <= dvThermalProduction[LLnodes_warehouse, 1:T] <= M);
# LL load balances
@constraint(Lower(model), loadbalance_withPV[n in LLnodes_withPV, t in 1:T],
    yi[n, t] + ypvprod[n, t] == LDFinputs.Pload[n][t] + ye[n, t]
);
@constraint(Lower(model), loadbalance_warehouse[n in LLnodes_warehouse, t in 1:T],
    yi[n, t] == LDFinputs.Pload[n][t] + dvThermalProduction[n, t] / CHILLER_COP
);
# PV production
@constraint(Lower(model), pvprod[n in LLnodes_withPV, t in 1:T],
    ypvprod[n, t] + spvprod[n, t] == ypv[n] * prod_factor[t]  # need curtailment to make UL pay for ye ?
);

#=
        LL REFRIGERATED WAREHOUSE MODEL
=#
R = 0.00025  # K/kW
C = 1e5   # kJ/K
A = reshape([-1/(R*C)], 1,1)
B = [1/(R*C) 1/C]
u = [tamb zeros(8760)]';  # could replace the zeros vector with endogenous heat input
J = size(B,2)
# state space temperature evolution
@constraint(Lower(model), sstemperature[n in LLnodes_warehouse, t in 2:T],
    dvTemperature[n, t] == dvTemperature[n, t-1] + A[1, 1] * dvTemperature[n, t-1] + 
    sum(B[1, j] * u[j, t-1] for j=1:J) + B[1, 1] * (-dvThermalProduction[n, t-1])
);
@constraint(Lower(model), init_temperature[n in LLnodes_warehouse], dvTemperature[n, 1] == -1.0);  # initial temperature

#= 
        UPPER MODEL VARIABLES
=#
# duals of LL load balances
@variable(Upper(model), 0 <= lambdaPV <= M, DualOf(loadbalance_withPV));
@variable(Upper(model), 0 <= lambda_warehouse <= M, DualOf(loadbalance_warehouse))

# UL variables in LL objective
@variable(Upper(model), 0 <= xe[LLnodes_withPV, 1:T] <= M);  # PV export compensation
@variable(Upper(model), 0 <= xi[LLnodes, 1:T] <= M)  # time varying tariff for Warehouses

@variable(Upper(model), 0 <= x0[1:T] <= M);  # positive power import at feeder head

# TODO which nodes can have BESS? Assume the unloaded nodes
LDF.build_ldf!(Upper(model), LDFinputs, LLnodes, ye, yi);
@constraint(Upper(model), [t in 1:T], x0[t] >= model[:Pⱼ]["0", t] );

# yi ⟂ ye
for nd in LLnodes_withPV, (i, e) in zip(yi[nd, :], ye[nd, :])
    @constraint(Upper(model),
        [i, e] in MOI.SOS1([1.0, 2.0])
    );
end
#= 
        LL objective
=#
@objective(Lower(model), Min,
    pwf * sum(ci[t] * yi[n, t] for n in LLnodes, t in 1:T) + 
    pwf * sum( -xe[n, t] * ye[n, t] for n in LLnodes_withPV, t in 1:T) +
    pwf * sum(  xi[n, t] * yi[n, t] for n in LLnodes_warehouse, t in 1:T) +
    sum(cpv * ypv[n] for n in LLnodes_withPV)
);
#= 
        UL objective
=#


@objective(Upper(model), Min,
    pwf * sum( clmp[t] * x0[t] for t in 1:T) +
    pwf * sum( ye[n, t] * lambdaPV[n, t] for n in LLnodes_withPV, t in 1:T) +
    pwf * sum( yi[n, t] * lambda_warehouse[n, t] for n in LLnodes_warehouse, t in 1:T)
    # UL should not benefit from charging LL w/o some cost
);

optimize!(model)

#= use model.var_lower and MOI.get(upper, MOI.ObjectiveFunction{MOI.get(upper, MOI.ObjectiveFunctionType())}()) to look up linearization
upperobj = MOI.get(model.solver, MOI.ObjectiveFunction{MOI.get(model.solver, MOI.ObjectiveFunctionType())}())

A = [(j,n), (1,1), (2,2)] since ye is the first variable in LL and its load balance is the first two constraints in LL

upper = JuMP.backend(model.upper)
lower = JuMP.backend(model.lower)

lower_var_indices_of_upper_vars = JuMP.index.(
    collect(values(model.upper_to_lower_link)))
upper_to_lower_var_indices = BilevelJuMP.convert_indices(model.link)
upper_var_lower_ctr = BilevelJuMP.index2(model.upper_var_to_lower_ctr_link)

# moi.jl
U, V, w = BilevelJuMP.standard_form(lower, upper_var_indices=lower_var_indices_of_upper_vars);
rows1, cols1, redundant_vals = BilevelJuMP.find_connected_rows_cols(V, 1, 1; skip_1st_col_check=true)
([1, 5], [4, 8, 10], false)
rows2, cols2, redundant_vals = BilevelJuMP.find_connected_rows_cols(V, 2, 2; skip_1st_col_check=true)
([2, 6], [6, 9, 11], false)

julia> for c in cols1 println(model.var_lower[c]) end
yi[34,1]
ypvprod[34,1]
spvprod[34,1]

julia> for c in cols2 println(model.var_lower[c]) end
yi[34,2]
ypvprod[34,2]
spvprod[34,2]

p = Ajn / Vjn = -1
MathOptInterface.ScalarAffineFunction{Float64}(MathOptInterface.ScalarAffineTerm{Float64}[
    MathOptInterface.ScalarAffineTerm{Float64}(115.82595, MathOptInterface.VariableIndex(1)), pwf*ci * yi[1]
    MathOptInterface.ScalarAffineTerm{Float64}(115.82595, MathOptInterface.VariableIndex(3)), pwf*ci * yi[2]
    MathOptInterface.ScalarAffineTerm{Float64}(-10000.0, MathOptInterface.VariableIndex(18)), M * mu_upper
    MathOptInterface.ScalarAffineTerm{Float64}(-10000.0, MathOptInterface.VariableIndex(24)), 
    MathOptInterface.ScalarAffineTerm{Float64}(-10000.0, MathOptInterface.VariableIndex(30)), 
    MathOptInterface.ScalarAffineTerm{Float64}(-10000.0, MathOptInterface.VariableIndex(42)), 
    MathOptInterface.ScalarAffineTerm{Float64}(-10000.0, MathOptInterface.VariableIndex(54)), 
    MathOptInterface.ScalarAffineTerm{Float64}(-10000.0, MathOptInterface.VariableIndex(66)), 
    MathOptInterface.ScalarAffineTerm{Float64}(0.01419, MathOptInterface.VariableIndex(75)),   x0
    MathOptInterface.ScalarAffineTerm{Float64}(0.0151825, MathOptInterface.VariableIndex(76)), x0
    MathOptInterface.ScalarAffineTerm{Float64}(-122.08466166373908, MathOptInterface.VariableIndex(466)),  LDFinputs.Pload["34"][1] (w1 * λ1)
    MathOptInterface.ScalarAffineTerm{Float64}(-105.96647268103634, MathOptInterface.VariableIndex(467))], 0.0) LDFinputs.Pload["34"][2] (w2 * λ2)

for j = 5 and 6 w_j = 0 so no terms needed (would be lambda_j * w_j)
=#
# create problems that lead to non-zero ye (xe) and xi
# xi = ci s.t. lambda_warehouse = 0 and both levels benefit
#=
with yi for flexHVAC and ye for PV the model violates condition 5 (p = A_jn/V_jn) b/c V_jn = +/-1 for yi/ye
    but are the two sets of bilinear terms separable and therefore the conditions as well? Yes
    dealt with this in BilevelJuMP

Try yi only s.t. UL can use thermal storage to absorb excess PV?

Need scenario with PV desirable to lower total cost of power sometimes and other time PV export is undesirable due to voltage rise.
xe will encourage export when desired. When it is not desired UL will use BESS to absorb it (w/o compensation?)
=#
value.(model[:ypv])
maximum(sqrt.(value.(model[:vsqrd])))


maximum(sqrt.(value.(model[:xi])))




# # epigraph for norm1 of voltage excursion
# @constraint(Upper(model), [t in 1:T],
#     epi[t] >= v1[t] - 1.0
# )
# @constraint(Upper(model), [t in 1:T],
#     -epi[t] <= v1[t] - 1.0
# )
# @constraint(Upper(model), [t in T+1:2T],
#     epi[t] >= v2[t-T] - 1.0
# )
# @constraint(Upper(model), [t in T+1:2T],
#     -epi[t] <= v2[t-T] - 1.0
# )
# @constraint(Upper(model), [t in 1:2T],
#     max_epi >= epi[t]
# )
# instantiate optimizer and empty model
optimizer = Gurobi.Optimizer()
model = BilevelModel(()->optimizer, linearize_bilinear_upper_terms=true)

#= 
        LOWER MODEL VARIABLES
=#
# TODO add more node indices to LL
@variable(Lower(model), 0 <= ye[LLnodes, 1:T] <= M) # TODO? start values from REopt
@variable(Lower(model), 0 <= yi[LLnodes, 1:T] <= M)
@variable(Lower(model), 0 <= ypv[LLnodes] <= M)
@variable(Lower(model), 0 <= ypvprod[LLnodes, 1:T] <= M)
@variable(Lower(model), 0 <= spvprod[LLnodes, 1:T] <= M) 
@variable(Lower(model), 0 <= dummyslack[LLnodes, 1:T] <= M) 
@variable(Lower(model), -10 <= dvTemperature[LLnodes, 1:T] <= 0.0)
@variable(Lower(model), 0 <= dvThermalProduction[LLnodes, 1:T] <= M)
# LL load balances
@constraint(Lower(model), loadbalance[n in LLnodes, t in 1:T],
    yi[n, t] + ypv[n] * prod_factor[t] == LDFinputs.Pload[n][t] + dvThermalProduction[n, t] / CHILLER_COP + ye[n, t]
)
#=
        LL REFRIGERATED WAREHOUSE MODEL
=#
R = 0.00025  # K/kW
C = 1e5   # kJ/K
A = reshape([-1/(R*C)], 1,1)
B = [1/(R*C) 1/C]
u = [tamb zeros(8760)]';  # could replace the zeros vector with endogenous heat input
J = size(B,2)
# state space temperature evolution
@constraint(Lower(model), sstemperature[n in LLnodes, t in 2:T],
    dvTemperature[n, t] == dvTemperature[n, t-1] + A[1, 1] * dvTemperature[n, t-1] + 
    sum(B[1, j] * u[j, t-1] for j=1:J) + B[1, 1] * (-dvThermalProduction[n, t-1])
)
@constraint(Lower(model), init_temperature[n in LLnodes], dvTemperature[n, 1] == -1.0)  # initial temperature

#= 
        UPPER MODEL VARIABLES
=#
# duals of LL load balances
@variable(Upper(model), 0 <= lambda <= M, DualOf(loadbalance))

# UL variables in LL objective
@variable(Upper(model), 0 <= xe[LLnodes, 1:T] <= M)  # PV export compensation

@variable(Upper(model), 0 <= x0[1:T] <= M)  # positive power import at feeder head

# TODO which nodes can have BESS? Assume the unloaded nodes
LDF.build_ldf!(Upper(model), LDFinputs, LLnodes, ye, yi)
@constraint(Upper(model), [t in 1:T], x0[t] >= model[:Pⱼ]["0", t] )

# yi ⟂ ye
for nd in LLnodes, (i, e) in zip(yi[nd, :], ye[nd, :])
    @constraint(Upper(model),
        [i, e] in MOI.SOS1([1.0, 2.0])
    )
end
#= 
        LL objective
=#
@objective(Lower(model), Min,
    pwf * sum(ci[t] * yi[n, t] for n in LLnodes, t in 1:T) + 
    pwf * sum( -xe[n, t] * ye[n, t] for n in LLnodes, t in 1:T) +
    sum(cpv * ypv[n] for n in LLnodes)
)
#= 
        UL objective
=#
# TODO voltage cost? wires alternative? 
# TODO UL BESS cost (w/o it pwf has no effect)
@objective(Upper(model), Min,
    pwf * sum( clmp[t] * x0[t] for t in 1:T) +
    pwf * sum( ye[n, t] * lambda[n, t] for n in LLnodes, t in 1:T)
)

optimize!(model)








# this model works w/o node index: (so try adding node index)

    T=2
    optimizer = Gurobi.Optimizer()
    model = BilevelModel(()->optimizer, linearize_bilinear_upper_terms=true)

    EXISTING_CHILLER_COP = 4.55
    R = 0.00025  # K/kW
    C = 1e5   # kJ/K

    A = reshape([-1/(R*C)], 1,1)
    B = [1/(R*C) 1/C]
    u = [tamb zeros(8760)]';
    J = 2

    M = 10_000

    nodes = 1:2
    ld = repeat([1], T)

    @variable(Lower(model), 19.8 <= dvTemperature[nodes, 1:T] <= 22.0)  # NEED HIGH-COST COMFORT VIOLATION?
    @variable(Lower(model), 0 <= dvThermalProduction[nodes, 1:T] <= M)

    @constraint(Lower(model), [n in nodes, ts in 2:T],
        dvTemperature[n,ts] == dvTemperature[n,ts-1] + A[1, 1] * dvTemperature[n,ts-1] + 
        sum(B[1, j] * u[j, ts-1] for j=1:J) + B[1, 1] * (-dvThermalProduction[n,ts-1])
    )

    # TODO initial_temperatures
    @variable(Lower(model), 0 <= ye[nodes, 1:T] <= M) # TODO? start values from REopt
    @variable(Lower(model), 0 <= yi[nodes, 1:T] <= M)
    @variable(Lower(model), 0 <= ypv[nodes] <= M) # TODO add more node indices to LL
    @variable(Lower(model), 0 <= ypvprod[nodes, 1:T] <= M) # TODO add more node indices to LL
    @variable(Lower(model), 0 <= spvprod[nodes, 1:T] <= M) # TODO add more node indices to LL
    @variable(Lower(model), 0 <= dummyslack[nodes, 1:T] <= M) # TODO add more node indices to LL


    for (i, e) in zip(yi, ye)
        @constraint(Upper(model),
            [i, e] in MOI.SOS1([1.0, 2.0])
        )
    end

    @constraint(Lower(model), loadbalance[n in nodes, t in 1:T],
        yi[n,t] + ypvprod[n,t] == ld[t] + ye[n,t] + dvThermalProduction[n,t] / EXISTING_CHILLER_COP
    )

    @constraint(Lower(model), [n in nodes, t in 1:T],
        ypvprod[n,t] + spvprod[n,t] == ypv[n] * prod_factor[t]
    )
    @constraint(Lower(model), [n in nodes, t in 1:T],
        ypvprod[n,t] + dummyslack[n,t] == ypv[n]
    )

    @variable(Upper(model), 0 <= λ <= M, DualOf(loadbalance))

    @variable(Upper(model), 0 <= xe[nodes, 1:T] <= M)
    #= 
            LL objective
    =#
    @objective(Lower(model), Min,
        pwf * sum(ci[t] * yi[n,t] - xe[n,t] * ye[n,t] for n in nodes, t in 1:T) + 
        sum(cpv * ypv[n] for n in nodes)
    )
    #= 
            UL objective
    =#
    @variable(Upper(model), 0 <= x0[1:T] <= M)
    @constraint(Upper(model), [t in 1:T],
        x0[t] == sum(yi[n,t] - ye[n,t] for n in nodes)
    )

    # TODO voltage cost? wires alternative? UL BESS ?
    @objective(Upper(model), Min,
        sum(
            clmp[t] * x0[t] +
            ye[n,t] * λ[n,t]
        for n in nodes, t in 1:T)
    )

    optimize!(model)

#=
Condition 3 is not met with node index, but what if treated the LL problems separately in linearization process? (along with UL bilinear terms)
Seems right that Condition 3 does not pass because it is for all of the y_n variables being connected, but the separate node y_n variables are not connected.
can compare quadratic solution to linearized solution with hack for not enforcing Condition 3

TODO add above model to tests in BilevelJuMP for testing separable problems

dev /Users/nlaws/Projects/BilevelJuMP.jl
NOTE does not work with "add" !
=#


    T=2
    optimizer = Gurobi.Optimizer()
    quadmodel = BilevelModel(()->optimizer, linearize_bilinear_upper_terms=false)

    @variable(Lower(quadmodel), 19.8 <= dvTemperature[nodes, 1:T] <= 22.0)  # NEED HIGH-COST COMFORT VIOLATION?
    @variable(Lower(quadmodel), 0 <= dvThermalProduction[nodes, 1:T] <= M)

    @constraint(Lower(quadmodel), [n in nodes, ts in 2:T],
        dvTemperature[n,ts] == dvTemperature[n,ts-1] + A[1, 1] * dvTemperature[n,ts-1] + 
        sum(B[1, j] * u[j, ts-1] for j=1:J) + B[1, 1] * (-dvThermalProduction[n,ts-1])
    )

    # TODO initial_temperatures
    @variable(Lower(quadmodel), 0 <= ye[nodes, 1:T] <= M) # TODO? start values from REopt
    @variable(Lower(quadmodel), 0 <= yi[nodes, 1:T] <= M)
    @variable(Lower(quadmodel), 0 <= ypv[nodes] <= M) # TODO add more node indices to LL
    @variable(Lower(quadmodel), 0 <= ypvprod[nodes, 1:T] <= M) # TODO add more node indices to LL
    @variable(Lower(quadmodel), 0 <= spvprod[nodes, 1:T] <= M) # TODO add more node indices to LL
    @variable(Lower(quadmodel), 0 <= dummyslack[nodes, 1:T] <= M) # TODO add more node indices to LL


    for (i, e) in zip(yi, ye)
        @constraint(Upper(quadmodel),
            [i, e] in MOI.SOS1([1.0, 2.0])
        )
    end

    @constraint(Lower(quadmodel), loadbalance[n in nodes, t in 1:T],
        yi[n,t] + ypvprod[n,t] == ld[t] + ye[n,t] + dvThermalProduction[n,t] / EXISTING_CHILLER_COP
    )

    @constraint(Lower(quadmodel), [n in nodes, t in 1:T],
        ypvprod[n,t] + spvprod[n,t] == ypv[n] * prod_factor[t]
    )
    @constraint(Lower(quadmodel), [n in nodes, t in 1:T],
        ypvprod[n,t] + dummyslack[n,t] == ypv[n]
    )

    @variable(Upper(quadmodel), 0 <= λ <= M, DualOf(loadbalance))

    @variable(Upper(quadmodel), 0 <= xe[nodes, 1:T] <= M)
    #= 
            LL objective
    =#
    @objective(Lower(quadmodel), Min,
        pwf * sum(ci[t] * yi[n,t] - xe[n,t] * ye[n,t] for n in nodes, t in 1:T) + 
        sum(cpv * ypv[n] for n in nodes)
    )
    #= 
            UL objective
    =#
    @variable(Upper(quadmodel), 0 <= x0[1:T] <= M)
    @constraint(Upper(quadmodel), [t in 1:T],
        x0[t] == sum(yi[n,t] - ye[n,t] for n in nodes)
    )

    # TODO voltage cost? wires alternative? UL BESS ?
    @objective(Upper(quadmodel), Min,
        sum(
            clmp[t] * x0[t] +
            ye[n,t] * λ[n,t]
        for n in nodes, t in 1:T)
    )

    optimize!(quadmodel)

    objective_value(quadmodel) ≈ objective_value(model) # true!

    #=
    Looks like both conditions 3 and 5 should be applied seperately to seperable LL problems, which requries
    seperate A sets (A_N and A respectively)
    =# 

    # jump.jl
    upper = JuMP.backend(model.upper)
    lower = JuMP.backend(model.lower)

    lower_var_indices_of_upper_vars = JuMP.index.(
        collect(values(model.upper_to_lower_link)))
    upper_to_lower_var_indices = BilevelJuMP.convert_indices(model.link)
    upper_var_lower_ctr = BilevelJuMP.index2(model.upper_var_to_lower_ctr_link)

# moi.jl
U, V, w = BilevelJuMP.standard_form(lower, upper_var_indices=lower_var_indices_of_upper_vars);
V

writedlm( "V.csv",  V, ',')
# print keys for labeling columns in Excel (inspecting V connections)
ks = sort!(collect(keys(model.var_lower)))
for k in ks
    println(model.var_lower[k])
end
