
using Gurobi
using JuMP
using REoptLite
using DelimitedFiles
import MathOptInterface
import LinDistFlow as LDF  # REoptLite exports LinDistFlow, maybe should remove that
using Random
Random.seed!(42)
# next try adding one PV in LL
include("extend_lindistflow.jl")
const MOI = MathOptInterface

# TODO get PV and BESS in solution with lower T using high pwf? might need to raise LMP too
# maybe start values will help? previous commit solved with 10% gap, had v_lolim = 0 (and v=0 in solution)

lat, lon = 30.2672, -97.7431  # Austin, TX
cpv = 1400
cbkW = 700
cbkWh = 350
η = 0.95
years = 20
discountrate = 0.05
M = 10_000  # default bound for decision variables
CHILLER_COP = 4.55
#***********************#
Sbase=1e3
#***********************#
pwf = REoptLite.annuity(years, 0.0, discountrate)
clmp = vec(readdlm("./data/cleaned_ercot2019_RTprices.csv", ',', Float64, '\n'));
clmp = abs.(clmp) / Sbase;  # problem is unbounded with negative prices, convert from $/MWh to $/kWh
tamb = REoptLite.get_ambient_temperature(lat, lon);
prod_factor = REoptLite.get_pvwatts_prodfactor(lat, lon);  # TODO this function is only in flex branch
LDFinputs = LDF.singlephase38linesInputs(Sbase=Sbase);  # TODO this method is not released yet
loadnodes = collect(keys(LDFinputs.Pload))
LLnodes_withPV = ["25"]
LLnodes_warehouse = ["33"]
LLnodes = union(LLnodes_withPV, LLnodes_warehouse)  # all nodes in LL model (that have decisions)


ULnodes_withBESS = ["2", "7", "24"]

profile_names = ["FastFoodRest", "FullServiceRest", "Hospital", "LargeHotel", "LargeOffice", 
"MediumOffice", "MidriseApartment", "Outpatient", "PrimarySchool", "RetailStore", "SecondarySchool", 
"SmallHotel", "SmallOffice", "StripMall", "Supermarket", "Warehouse"]
doe_profiles = Dict{String, Vector{Float64}}()
for name in profile_names
    doe_profiles[name] = REoptLite.BuiltInElectricLoad("", name, lat, lon, 2017)
end


# fill in net uncontrolled loads
T = 8760
LDFinputs.Ntimesteps = T
ci = repeat([0.25], T)
rand_names = rand(profile_names, length(loadnodes))
for (i, node) in enumerate(loadnodes)
    LDFinputs.Pload[node] = doe_profiles[rand_names[i]][1:T] / Sbase;
    LDFinputs.Qload[node] = doe_profiles[rand_names[i]][1:T] / Sbase * 0.1;
end
# pepper some pv into the system
PVkW = 2e3   # TODO more baseline PV ?
LDFinputs.Pload["3"] .-= PVkW * prod_factor[1:T] / Sbase;
LDFinputs.Qload["3"] .-= PVkW * prod_factor[1:T] / Sbase * 0.1;
LDFinputs.Pload["30"] .-= PVkW * prod_factor[1:T] / Sbase;
LDFinputs.Qload["30"] .-= PVkW * prod_factor[1:T] / Sbase * 0.1;
LDFinputs.Pload["28"] .-= PVkW * prod_factor[1:T] / Sbase;
LDFinputs.Qload["28"] .-= PVkW * prod_factor[1:T] / Sbase * 0.1;
LDFinputs.Pload["27"] .-= PVkW * prod_factor[1:T] / Sbase;
LDFinputs.Qload["27"] .-= PVkW * prod_factor[1:T] / Sbase * 0.1;
LDFinputs.Pload["18"] .-= PVkW * prod_factor[1:T] / Sbase;
LDFinputs.Qload["18"] .-= PVkW * prod_factor[1:T] / Sbase * 0.1;

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


# compare to BilevelJuMP result with T = 2, 10?

function linearized_problem(cpv, ci, clmp, LLnodes, LLnodes_withPV, LLnodes_warehouse, LDFinputs; 
    T=8760)

    R = 0.00025  # K/kW
    C = 1e5   # kJ/K
    A = reshape([-1/(R*C)], 1,1)
    B = [1/(R*C) 1/C]
    u = [tamb zeros(8760)]';  # could replace the zeros vector with endogenous heat input
    J = size(B,2)
    M = 1e6
    T_hi = 0
    T_lo = -20

    model = JuMP.Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "MIPGap", 1e-1)

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

    # Complementary slackness of KKT, only modeling lower bound in most cases
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


function linearized_problem_bess(cpv, ci, clmp, LLnodes, LLnodes_withPV, LLnodes_warehouse, LDFinputs, ULnodes_withBESS; 
    T=8760)

    R = 0.00025  # K/kW
    C = 1e5   # kJ/K
    A = reshape([-1/(R*C)], 1,1)
    B = [1/(R*C) 1/C]
    u = [tamb zeros(8760)]';  # could replace the zeros vector with endogenous heat input
    J = size(B,2)
    Mbig = peak_load * 10
    Msml = peak_single_load * 10
    T_hi = 0
    T_lo = -20

    model = JuMP.Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "MIPGap", 5e-2)

    @variables model begin
        Msml >= yi[LLnodes, 1:T] >= 0
        Msml >= ye[LLnodes_withPV, 1:T] >= 0
        Mbig >= ypv[LLnodes_withPV] >=0
        Mbig >= ypvprod[LLnodes_withPV, 1:T] >= 0
        T_hi >= ytemperature[LLnodes_warehouse, 1:T] >= T_lo
        ytherm[LLnodes_warehouse, 1:T] >= 0

        Mbig >= xe[LLnodes_withPV, 1:T] >= 0
        Mbig >= x0[1:T] >= 0
        Mbig >= xbkW[ULnodes_withBESS] >= 0
        Mbig >= xbkWh[ULnodes_withBESS] >= 0
        Mbig >= xsoc[ULnodes_withBESS, 0:T] >= 0
        Mbig >= xbplus[ULnodes_withBESS, 1:T] >= 0
        Mbig >= xbminus[ULnodes_withBESS, 1:T] >= 0

        Mbig >= lambda[LLnodes_withPV, 1:T] >= -Mbig
        Mbig >= lambda_warehouse[LLnodes_warehouse, 1:T] >= -Mbig
        Msml >= lambda_ss[LLnodes_warehouse, 2:T] >= -Msml
        Msml >= lambda_initTemperature >= -Msml
        Msml >= mu_i[LLnodes, 1:T] >= 0
        Mbig >= mu_e[LLnodes_withPV, 1:T] >= 0
        Mbig >= mu_pv[LLnodes_withPV] >= 0
        Mbig >= mu_pvprod[LLnodes_withPV, 1:T] >= 0
        Mbig >= mu_dd[LLnodes_withPV, 1:T] >= 0
        Msml >= mu_therm_lo[LLnodes_warehouse, 1:T] >= 0
        Msml >= mu_therm_hi[LLnodes_warehouse, 1:T] >= 0
        Msml >= mu_temperature_lo[LLnodes_warehouse, 1:T] >= 0
        Msml >= mu_temperature_hi[LLnodes_warehouse, 1:T] >= 0
    end

    # UL does not allow simultaneous export/import
    for n in LLnodes_withPV, t in 1:T
        @constraint(model,
            [ye[n,t], yi[n,t]] in MathOptInterface.SOS1([1.0, 2.0])
        )
    end
    
    # LinDistFlow (single phase, real power only)
    LDF.build_ldf!(model, LDFinputs, LLnodes, ye, yi, xbplus, xbminus);

    # Complementary slackness of KKT, only modeling lower bound in most cases
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
        [mu_therm_hi[n,t],  ytherm[n,t] - Msml] in MathOptInterface.SOS1([1.0, 2.0])
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

    # UL constraints
    @constraint(model, [t in 1:T], x0[t] >= model[:Pⱼ]["0", t] );

    @constraint(model, [n in ULnodes_withBESS],
        xsoc[n,0] == 0.5 * xbkWh[n]
    )
    @constraint(model, [n in ULnodes_withBESS],
        xsoc[n,T] == 0.5 * xbkWh[n]
    )
    @constraint(model, [n in ULnodes_withBESS, t in 1:T],
        xsoc[n,t] == xsoc[n,t-1] + xbplus[n,t] * η - xbminus[n,t] / η
    )
    @constraint(model, [n in ULnodes_withBESS, t in 1:T],
        xbkW[n] >= xbminus[n,t]
    )
    @constraint(model, [n in ULnodes_withBESS, t in 1:T],
        xbkW[n] >= xbplus[n,t]
    )
    @constraint(model, [n in ULnodes_withBESS, t in 1:T],
        xbkW[n] >= xbplus[n,t] + xbminus[n,t]
    )
    @constraint(model, [n in ULnodes_withBESS, t in 1:T],
        xbkWh[n] >= xsoc[n,t]
    )

    @objective(model, Min, 
        pwf * sum(x0[t] * clmp[t] for t in 1:T)
        + sum(cbkW * xbkW[n] + cbkWh * xbkWh[n] for n in ULnodes_withBESS)
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


function linearized_problem_bess_bigM(cpv, ci, clmp, LLnodes, LLnodes_withPV, LLnodes_warehouse, LDFinputs, ULnodes_withBESS; 
    T=8760)

    R = 0.00025  # K/kW
    C = 1e5   # kJ/K
    A = reshape([-1/(R*C)], 1,1)
    B = [1/(R*C) 1/C]
    u = [tamb zeros(8760)]';  # could replace the zeros vector with endogenous heat input
    J = size(B,2)
    Mbig = peak_load * 100
    Msml = peak_single_load * 100
    T_hi = 0
    T_lo = -20

    model = JuMP.Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "MIPGap", 5e-2)

    @variables model begin
        Msml >= yi[LLnodes, 1:T] >= 0
        Mbig >= ye[LLnodes_withPV, 1:T] >= 0
        Mbig >= ypv[LLnodes_withPV] >=0
        Mbig >= ypvprod[LLnodes_withPV, 1:T] >= 0
        T_hi >= ytemperature[LLnodes_warehouse, 1:T] >= T_lo
        ytherm[LLnodes_warehouse, 1:T] >= 0

        Mbig >= xe[LLnodes_withPV, 1:T] >= 0
        Mbig >= x0[1:T] >= 0
        Mbig >= xbkW[ULnodes_withBESS] >= 0
        Mbig >= xbkWh[ULnodes_withBESS] >= 0
        Mbig >= xsoc[ULnodes_withBESS, 0:T] >= 0
        Mbig >= xbplus[ULnodes_withBESS, 1:T] >= 0
        Mbig >= xbminus[ULnodes_withBESS, 1:T] >= 0

        Mbig >= lambda[LLnodes_withPV, 1:T] >= -Mbig
        Mbig >= lambda_warehouse[LLnodes_warehouse, 1:T] >= -Mbig
        Mbig >= lambda_ss[LLnodes_warehouse, 2:T] >= -Msml
        Mbig >= lambda_initTemperature >= -Msml
        Mbig >= mu_i[LLnodes, 1:T] >= 0
        Mbig >= mu_e[LLnodes_withPV, 1:T] >= 0
        Mbig >= mu_pv[LLnodes_withPV] >= 0
        Mbig >= mu_pvprod[LLnodes_withPV, 1:T] >= 0
        Mbig >= mu_dd[LLnodes_withPV, 1:T] >= 0
        Mbig >= mu_therm_lo[LLnodes_warehouse, 1:T] >= 0
        Mbig >= mu_therm_hi[LLnodes_warehouse, 1:T] >= 0
        Mbig >= mu_temperature_lo[LLnodes_warehouse, 1:T] >= 0
        Mbig >= mu_temperature_hi[LLnodes_warehouse, 1:T] >= 0
    end

    # UL does not allow simultaneous export/import
    @variable(model, byeyi[LLnodes_withPV, t in 1:T], Bin);
    @constraint(model, [n in LLnodes_withPV, t in 1:T],
        ye[n,t] <= Mbig * byeyi[n,t]
    );
    @constraint(model, [n in LLnodes_withPV, t in 1:T],
        yi[n,t] <= Mbig * (1-byeyi[n,t])
    );
    
    # LinDistFlow (single phase, real power only)
    LDF.build_ldf!(model, LDFinputs, LLnodes, ye, yi, xbplus, xbminus);
    # TODO if !isempty(LLnodes_warehouse) then add cons, variables

    # Complementary slackness of KKT, only modeling lower bound in most cases
    @variable(model, bypv[n in LLnodes_withPV], Bin);
    @constraint(model, [n in LLnodes_withPV],
        ypv[n] <= Mbig * bypv[n]
    );
    @constraint(model, [n in LLnodes_withPV],
        mu_pv[n] <= Mbig * (1-bypv[n])
    );
    @variable(model, bye[n in LLnodes_withPV, t in 1:T], Bin);
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        ye[n,t] <= Msml * bye[n,t] 
    );
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        mu_e[n,t] <= Msml * (1- bye[n,t])
    );
    @variable(model, byi[n in LLnodes, t in 1:T], Bin);
    @constraint(model, [n in LLnodes, t in 1:T], 
        yi[n,t] <= Msml * byi[n,t]
    );
    @constraint(model, [n in LLnodes, t in 1:T], 
        mu_i[n,t] <= Msml * (1 - byi[n,t])
    );
    @variable(model, bytherm_lo[n in LLnodes_warehouse, t in 1:T], Bin);
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        ytherm[n,t] <= Msml * bytherm_lo[n,t]
    );
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        mu_therm_lo[n,t] <= Msml * (1-bytherm_lo[n,t])
    );
    @variable(model, bytherm_hi[n in LLnodes_warehouse, t in 1:T], Bin);
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        ytherm[n,t] - Msml <= Msml * bytherm_hi[n,t]
    );
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        mu_therm_hi[n,t] <=  Msml * (1 - bytherm_hi[n,t])
    );
    @variable(model, bytemperature_lo[n in LLnodes_warehouse, t in 1:T], Bin);
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        T_lo - ytemperature[n,t] <= Msml * bytemperature_lo[n,t]
    );
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        mu_temperature_lo[n,t] <= Msml * (1 - bytemperature_lo[n,t])
    );
    @variable(model, bytemperature_hi[n in LLnodes_warehouse, t in 1:T], Bin);
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        ytemperature[n,t] - T_hi <= Msml * bytemperature_hi[n,t]
    );
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        mu_temperature_hi[n,t] <= Msml * (1 - bytemperature_hi[n,t])
    );
    @variable(model, bypvprod[n in LLnodes_withPV, t in 1:T], Bin)
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        ypvprod[n,t] <= Mbig * bypvprod[n,t]
    );
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        mu_pvprod[n,t] <= Mbig * (1 - bypvprod[n,t])
    )
    @variable(model, bydd[n in LLnodes_withPV, t in 1:T], Bin)
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        ypvprod[n,t] - ypv[n] * prod_factor[t] <= Mbig * bydd[n,t]
    )
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        mu_dd[n,t] <= Mbig * (1 - bydd[n,t])
    )

    # LL load balance with PV (lambda)
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        -ye[n,t] + yi[n,t] + ypvprod[n,t] - LDFinputs.Pload[n][t] == 0
    );

    # LL load balance for warehouse (lambda_warehouse)
    @constraint(model, [n in LLnodes_warehouse, t in 1:T],
        yi[n, t] - LDFinputs.Pload[n][t] - ytherm[n, t] / CHILLER_COP == 0
    );

    # state space temperature evolution (lambda_ss)
    @constraint(model, [n in LLnodes_warehouse, t in 2:T],
        ytemperature[n, t] - ytemperature[n, t-1] - A[1, 1] * ytemperature[n, t-1] -
        sum(B[1, j] * u[j, t-1] for j=1:J) + B[1, 1] * ytherm[n, t-1] ==  0
    );
    # initial temperature (lambda_initTemperature)
    @constraint(model, [n in LLnodes_warehouse], ytemperature[n, 1] == -1.0);  

    
    # LL operational
    @constraint(model, [n in LLnodes_withPV,t in 1:T], ypvprod[n,t] ≤ ypv[n] * prod_factor[t] )  # mu_dd[t]
    
    ## LL duals
    @constraint(model, [n in LLnodes_withPV], 
        cpv - mu_pv[n] - sum(mu_dd[n,t] * prod_factor[t] for t in 1:T) == 0);  # dual constraint of ypv
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        -xe[n,t] + lambda[n,t] - mu_e[n,t] == 0);  # dual constraint of ye[t]
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        ci[t] - lambda[n, t] - mu_i[n,t] == 0);  # dual constraint of yi[t] for LLnodes_withPV
    @constraint(model, [n in LLnodes_withPV, t in 1:T],       
        -lambda[n, t] - mu_pvprod[n,t] + mu_dd[n, t] == 0);  # dual constraint of ypvprod[t]
    @constraint(model, [n in LLnodes_warehouse, t in 1:T-1], 
        lambda_warehouse[n, t]/CHILLER_COP - B[1,1]*lambda_ss[n,t+1] - mu_therm_lo[n,t] + mu_therm_hi[n,t] == 0);  # dual constraint of ytherm
    @constraint(model, [n in LLnodes_warehouse], 
        lambda_warehouse[n, T]/CHILLER_COP - mu_therm_lo[n,T] + mu_therm_hi[n,T] == 0);  # dual constraint of ytherm[T]
    @constraint(model, [n in LLnodes_warehouse], 
        (1+A[1,1]) * lambda_ss[n,2] - mu_temperature_lo[n,1] + mu_temperature_hi[n,1] - lambda_initTemperature == 0);  # dual constraint of ytemperature[1]
    @constraint(model, [n in LLnodes_warehouse, t in 2:T-1],       
        -lambda_ss[n, t] + (1+A[1,1])*lambda_ss[n,t+1] - mu_temperature_lo[n,t] + mu_temperature_hi[n,t] == 0);  # dual constraint of ytemperature
    @constraint(model, [n in LLnodes_warehouse], 
        -lambda_ss[n,T] - mu_temperature_lo[n,T] + mu_temperature_hi[n,T] == 0);  # dual constraint of ytemperature[T]
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        ci[t] - lambda_warehouse[n, t] - mu_i[n,t] == 0);  # dual constraint of yi[t] for LLnodes_warehouse

    # UL constraints
    @constraint(model, [t in 1:T], x0[t] >= model[:Pⱼ]["0", t] );

    @constraint(model, [n in ULnodes_withBESS],
        xsoc[n,0] == 0.5 * xbkWh[n]
    );
    @constraint(model, [n in ULnodes_withBESS],
        xsoc[n,T] == 0.5 * xbkWh[n]
    );
    @constraint(model, [n in ULnodes_withBESS, t in 1:T],
        xsoc[n,t] == xsoc[n,t-1] + xbplus[n,t] * η - xbminus[n,t] / η
    );
    @constraint(model, [n in ULnodes_withBESS, t in 1:T],
        xbkW[n] >= xbminus[n,t]
    );
    @constraint(model, [n in ULnodes_withBESS, t in 1:T],
        xbkW[n] >= xbplus[n,t]
    );
    @constraint(model, [n in ULnodes_withBESS, t in 1:T],
        xbkW[n] >= xbplus[n,t] + xbminus[n,t]
    );
    @constraint(model, [n in ULnodes_withBESS, t in 1:T],
        xbkWh[n] >= xsoc[n,t]
    );

    @objective(model, Min, 
        pwf * sum(x0[t] * clmp[t] for t in 1:T)
        + sum(cbkW * xbkW[n] + cbkWh * xbkWh[n] for n in ULnodes_withBESS)
        + pwf * sum(
            ci[t] * yi[n,t] - lambda[n,t] * LDFinputs.Pload[n][t]
            for n in LLnodes_withPV, t in 1:T
        )
        + sum(ypv[n] * cpv for n in LLnodes_withPV)
    );

    optimize!(model)

    return model
end


function upper_only_with_bess(clmp, LDFinputs, ULnodes_withBESS; 
    T=8760)

    model = JuMP.Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "MIPGap", 1e-2)
    M = 1e3

    @variables model begin
        M >= yi[["99"], 1:0] >= 0 # dummy for building LinDistFlow model
        M >= ye[["99"], 1:0] >= 0 # dummy for building LinDistFlow model
        
        M >= xbkW[ULnodes_withBESS] >= 0
        M >= xbkWh[ULnodes_withBESS] >= 0
        M >= xsoc[ULnodes_withBESS, 0:T] >= 0
        M >= xbplus[ULnodes_withBESS, 1:T] >= 0
        M >= xbminus[ULnodes_withBESS, 1:T] >= 0
        M >= x0[1:T] >= 0
    end
    
    # LinDistFlow (single phase, real power only)
    LDF.build_ldf!(model, LDFinputs, ULnodes_withBESS, ye, yi, xbplus, xbminus);

    # UL constraints
    @constraint(model, [t in 1:T], x0[t] >= model[:Pⱼ]["0", t] );

    @constraint(model, [n in ULnodes_withBESS],
        xsoc[n,0] == 0.5 * xbkWh[n]
    );
    # @constraint(model, [n in ULnodes_withBESS],
    #     xsoc[n,T] == 0.5 * xbkWh[n]
    # )
    @constraint(model, [n in ULnodes_withBESS, t in 1:T],
        xsoc[n,t] == xsoc[n,t-1] + xbplus[n,t] * η - xbminus[n,t] / η
    );
    @constraint(model, [n in ULnodes_withBESS, t in 1:T],
        xbkW[n] >= xbminus[n,t]
    );
    @constraint(model, [n in ULnodes_withBESS, t in 1:T],
        xbkW[n] >= xbplus[n,t]
    );
    @constraint(model, [n in ULnodes_withBESS, t in 1:T],
        xbkW[n] >= xbplus[n,t] + xbminus[n,t]
    );
    @constraint(model, [n in ULnodes_withBESS, t in 1:T],
        xbkWh[n] >= xsoc[n,t]
    );

    @objective(model, Min, 
        pwf * sum(x0[t] * clmp[t] for t in 1:T)
        + sum(cbkW * xbkW[n] + cbkWh * xbkWh[n] for n in ULnodes_withBESS)
    );

    optimize!(model)

    # pwf *=5 results in batteries
    return model
end

# bkW = value.(model[:xbkW])
# bkWh = value.(model[:xbkWh])
# xbminus = value.(model[:xbminus]);
# xbmplus = value.(model[:xbplus]);
# sum(xbminus.data - xbmplus.data) ≈ 0
