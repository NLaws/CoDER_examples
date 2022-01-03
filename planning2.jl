
using Gurobi
using JuMP
using REoptLite
using DelimitedFiles
import MathOptInterface
import LinDistFlow as LDF  # REoptLite exports LinDistFlow, maybe should remove that

include("extend_lindistflow.jl")
const MOI = MathOptInterface



lat, lon = 30.2672, -97.7431  # Austin, TX
cpv = 1500
cbkW = 700
cbkWh = 350
η = 0.95
years = 10
discountrate = 0.05
M = 10_000  # default bound for decision variables
CHILLER_COP = 4.55

pwf = REoptLite.annuity(years, 0.0, discountrate)
clmp = vec(readdlm("./data/cleaned_ercot2019_RTprices.csv", ',', Float64, '\n'));
clmp = abs.(clmp) / 1e3;  # problem is unbounded with negative prices, convert from $/MWh to $/kWh
tamb = REoptLite.get_ambient_temperature(lat, lon);
prod_factor = REoptLite.get_pvwatts_prodfactor(lat, lon);  # TODO this function is only in flex branch

LDFinputs = LDF.singlephase38linesInputs(Sbase=1);  # TODO this method is not released yet
loadnodes = collect(keys(LDFinputs.Pload))
LLnodes = ["33", "34"]  # all nodes in LL model (that have decisions)
LLnodes_withPV = ["34"]
LLnodes_warehouse = ["33"]

profile_names = ["FastFoodRest", "FullServiceRest", "Hospital", "LargeHotel", "LargeOffice", 
"MediumOffice", "MidriseApartment", "Outpatient", "PrimarySchool", "RetailStore", "SecondarySchool", 
"SmallHotel", "SmallOffice", "StripMall", "Supermarket", "Warehouse"]
rand_names = rand(profile_names, length(loadnodes))
doe_profiles = Any[]
for (i, node) in enumerate(loadnodes)
    push!(doe_profiles, REoptLite.BuiltInElectricLoad("", rand_names[i], lat, lon, 2017))
end


# fill in net uncontrolled loads
T = 10
LDFinputs.Ntimesteps = T
ci = repeat([0.25], T)
for (i, node) in enumerate(loadnodes)
    LDFinputs.Pload[node] = doe_profiles[i][1:T];
    LDFinputs.Qload[node] = doe_profiles[i][1:T] * 0.1;
end
# pepper some pv into the system
PVkW = 2e3   # TODO more baseline PV ?
LDFinputs.Pload["3"] .-= PVkW * prod_factor[1:T];
LDFinputs.Qload["3"] .-= PVkW * prod_factor[1:T] * 0.1;
LDFinputs.Pload["30"] .-= PVkW * prod_factor[1:T];
LDFinputs.Qload["30"] .-= PVkW * prod_factor[1:T] * 0.1;
LDFinputs.Pload["28"] .-= PVkW * prod_factor[1:T];
LDFinputs.Qload["28"] .-= PVkW * prod_factor[1:T] * 0.1;
LDFinputs.Pload["27"] .-= PVkW * prod_factor[1:T];
LDFinputs.Qload["27"] .-= PVkW * prod_factor[1:T] * 0.1;
LDFinputs.Pload["18"] .-= PVkW * prod_factor[1:T];
LDFinputs.Qload["18"] .-= PVkW * prod_factor[1:T] * 0.1;



model = Model(Gurobi.Optimizer)
LDF.build_ldf!(model, LDFinputs)
optimize!(model)


# compare to BilevelJuMP result with T = 2, 10?

function linearized_problem(cpv, ci, clmp, LLnodes, LLnodes_withPV, LLnodes_warehouse, LDFinputs; 
    T=8760)

    R = 0.00025  # K/kW
    C = 1e5   # kJ/K
    A = reshape([-1/(R*C)], 1,1)
    B = [1/(R*C) 1/C]
    u = [tamb zeros(8760)]';  # could replace the zeros vector with endogenous heat input
    J = size(B,2)
    M = 1e8
    T_hi = 0
    T_lo = -20

    model = JuMP.Model(Gurobi.Optimizer)

    @variables model begin
        M >= yi[LLnodes, 1:T] >= 0
        M >= ye[LLnodes_withPV, 1:T] >= 0
        M >= ypv[LLnodes_withPV] >=0
        M >= ypvprod[LLnodes_withPV, 1:T] >= 0
        T_hi >= ytemperature[LLnodes_warehouse, 1:T] >= T_lo
        ytherm[LLnodes_warehouse, 1:T] >= 0

        M >= xe[LLnodes_withPV, 1:T] >= 0
        M >= x0[1:T] >= 0

        M >= lambda[LLnodes_withPV, 1:T] >= -M
        M >= lambda_warehouse[LLnodes_warehouse, 1:T] >= -M
        M >= lambda_ss[LLnodes_warehouse, 2:T] >= -M
        M >= lambda_initTemperature >= -M
        M >= mu_i[LLnodes, 1:T] >= 0
        M >= mu_e[LLnodes_withPV, 1:T] >= 0
        M >= mu_pv[LLnodes_withPV] >= 0
        M >= mu_pvprod[LLnodes_withPV, 1:T] >= 0
        M >= mu_dd[LLnodes_withPV, 1:T] >= 0
        M >= mu_therm_lo[LLnodes_warehouse, 1:T] >= 0
        M >= mu_therm_hi[LLnodes_warehouse, 1:T] >= 0
        M >= mu_temperature_lo[LLnodes_warehouse, 1:T] >= 0
        M >= mu_temperature_hi[LLnodes_warehouse, 1:T] >= 0
    end

    # UL does not allow simultaneous export/import
    for n in LLnodes_withPV, t in 1:T
        @constraint(model,
            [ye[n,t], yi[n,t]] in MathOptInterface.SOS1([1.0, 2.0])
        )
    end
    
    # LinDistFlow (single phase, real power only)
    LDF.build_ldf!(model, LDFinputs, LLnodes, ye, yi);

    # Complementary slackness of KKT
    @constraint(model, [n in LLnodes_withPV],
        [mu_pv[n], ypv[n]] in MathOptInterface.SOS1([1.0, 2.0])
    )
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        [mu_e[n,t],  ye[n,t]] in MathOptInterface.SOS1([1.0, 2.0])
    )
    @constraint(model, [n in LLnodes, t in 1:T], 
        [mu_i[n,t],  yi[n,t]] in MathOptInterface.SOS1([1.0, 2.0])
    )
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        [mu_therm_lo[n,t],  ytherm[n,t]] in MathOptInterface.SOS1([1.0, 2.0])
    )
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        [mu_therm_hi[n,t],  ytherm[n,t] - M] in MathOptInterface.SOS1([1.0, 2.0])
    )
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        [mu_temperature_lo[n,t],  T_lo - ytemperature[n,t]] in MathOptInterface.SOS1([1.0, 2.0])
    )
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        [mu_temperature_hi[n,t], ytemperature[n,t] - T_hi] in MathOptInterface.SOS1([1.0, 2.0])
    )
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        [mu_pvprod[n,t],  ypvprod[n,t]] in MathOptInterface.SOS1([1.0, 2.0])
    )
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        [mu_dd[n,t], ypvprod[n,t] - ypv[n] * prod_factor[t]] in MathOptInterface.SOS1([1.0, 2.0])
    )
    # NEED -M in above constraints if upper bound variables by M

    # LL load balance with PV (lambda)
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        -ye[n,t] + yi[n,t] + ypvprod[n,t] - LDFinputs.Pload[n][t] == 0
    )

    # LL load balance for warehouse (lambda_warehouse)
    @constraint(model, [n in LLnodes_warehouse, t in 1:T],
        yi[n, t] - LDFinputs.Pload[n][t] - ytherm[n, t] / CHILLER_COP == 0
    );

    # state space temperature evolution (lambda_ss)
    @constraint(model, [n in LLnodes_warehouse, t in 2:T],
        ytemperature[n, t] - ytemperature[n, t-1] - A[1, 1] * ytemperature[n, t-1] -
        sum(B[1, j] * u[j, t-1] for j=1:J) + B[1, 1] * ytherm[n, t-1] ==  0
    );
    @constraint(model, [n in LLnodes_warehouse], ytemperature[n, 1] == -1.0);  # initial temperature (lambda_initTemperature)

    
    # LL operational
    @constraint(model, [n in LLnodes_withPV,t in 1:T], ypvprod[n,t] ≤ ypv[n] * prod_factor[t] )  # mu_dd[t]
    
    ## LL duals
    @constraint(model, [n in LLnodes_withPV], 
        cpv - mu_pv[n] - sum(mu_dd[n,t] * prod_factor[t] for t in 1:T) == 0)  # dual constraint of ypv
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        -xe[n,t] + lambda[n,t] - mu_e[n,t] == 0)  # dual constraint of ye[t]
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        ci[t] - lambda[n, t] - mu_i[n,t] == 0)  # dual constraint of yi[t] for LLnodes_withPV
    @constraint(model, [n in LLnodes_withPV, t in 1:T],       
        -lambda[n, t] - mu_pvprod[n,t] + mu_dd[n, t] == 0)  # dual constraint of ypvprod[t]
    @constraint(model, [n in LLnodes_warehouse, t in 1:T-1], 
        lambda_warehouse[n, t]/CHILLER_COP - B[1,1]*lambda_ss[n,t+1] - mu_therm_lo[n,t] + mu_therm_hi[n,t] == 0)  # dual constraint of ytherm
    @constraint(model, [n in LLnodes_warehouse], 
        lambda_warehouse[n, T]/CHILLER_COP - mu_therm_lo[n,T] + mu_therm_hi[n,T] == 0)  # dual constraint of ytherm[T]
    @constraint(model, [n in LLnodes_warehouse], (1+A[1,1]) * lambda_ss[n,2] - mu_temperature_lo[n,1] + mu_temperature_hi[n,1] - lambda_initTemperature == 0)  # dual constraint of ytemperature[1]
    @constraint(model, [n in LLnodes_warehouse, t in 2:T-1],       
        -lambda_ss[n, t] + (1+A[1,1])*lambda_ss[n,t+1] - mu_temperature_lo[n,t] + mu_temperature_hi[n,t] == 0)  # dual constraint of ytemperature
    @constraint(model, [n in LLnodes_warehouse], -lambda_ss[n,T] - mu_temperature_lo[n,T] + mu_temperature_hi[n,T] == 0)  # dual constraint of ytemperature[T]
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        ci[t] - lambda_warehouse[n, t] - mu_i[n,t] == 0)  # dual constraint of yi[t] for LLnodes_warehouse


    @constraint(model, [t in 1:T], x0[t] >= model[:Pⱼ]["0", t] );
    @objective(model, Min, 
        pwf * sum(x0[t] * clmp[t] for t in 1:T)
        + pwf * sum(
            ci[t] * yi[n,t] - lambda[n,t] * LDFinputs.Pload[n][t]
            for n in LLnodes_withPV, t in 1:T
        )
        + sum(ypv[n] * cpv for n in LLnodes_withPV)
    )

    optimize!(model)

    # objective_value should be 3.988362294337e+02 with T=2 Got it! (and match with T=5 1.018037704829e+03)
    return model
end
