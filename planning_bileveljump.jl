#=
The planning example assumes that:
- the power flow is approximated with the LinDistFlow model (Baran & Wu 1989)
- the system owner is considering battery energy storage, compensating customers for exported power, and sending price signal to refrigerated warehouses
- the system owner pays the 2019 ERCOT average hourly real-time market price for imports at the feeder head
- customers have the option to purchase PV systems to reduce their costs
- there are large refrigerated warehouses in the network with price responsive HVAC systems
- the base retail tariff is defined in parameters.jl
- the network model is based on the 38 node model from Andrianesis et al. 2019
- location is Austin, TX
=#


function bileveljump_bess_bigM(optimizer, T, LDFinputs, bigM; 
    setgap=false, MIPGap=1e-2)

    model = BilevelModel(()->optimizer(), linearize_bilinear_upper_terms=true,
        mode=BilevelJuMP.FortunyAmatMcCarlMode(primal_big_M = bigM, dual_big_M = bigM)
    )
    if setgap  # depends on solver parameter names, this is for Gurobi
        set_optimizer_attribute(model, "MIPGap", MIPGap)
    end

    @variables Lower(model) begin
        bigM >= yi[LLnodes, 1:T] >= 0
        bigM >= ye[LLnodes_withPV, 1:T] >= 0
        bigM >= ypv[LLnodes_withPV] >=0
        bigM >= ypvprod[LLnodes_withPV, 1:T] >= 0
        bigM >= spvprod[LLnodes_withPV, 1:T] >= 0
        T_hi >= ytemperature[LLnodes_warehouse, 1:T] >= T_lo
        ytherm[LLnodes_warehouse, 1:T] >= 0
    end
    @variables Upper(model) begin
        bigM >= xe[LLnodes_withPV, 1:T] >= 0        
        bigM >= xi[LLnodes_warehouse, 1:T] >= 0
        bigM >= x0[1:T] >= 0
        bigM >= xbkW[ULnodes_withBESS] >= 0
        bigM >= xbkWh[ULnodes_withBESS] >= 0
        bigM >= xsoc[ULnodes_withBESS, 0:T] >= 0
        bigM >= xbplus[ULnodes_withBESS, 1:T] >= 0
        bigM >= xbminus[ULnodes_withBESS, 1:T] >= 0
    end

    # UL does not allow simultaneous export/import
    @variable(Upper(model), byeyi[LLnodes_withPV, t in 1:T], Bin);
    @constraint(Upper(model), [n in LLnodes_withPV, t in 1:T],
        ye[n,t] <= bigM * byeyi[n,t]
    );
    @constraint(Upper(model), [n in LLnodes_withPV, t in 1:T],
        yi[n,t] <= bigM * (1-byeyi[n,t])
    );
    
    # LinDistFlow (single phase, real power only)
    LDF.build_ldf!(Upper(model), LDFinputs, LLnodes, ye, yi, xbplus, xbminus);
    # TODO if !isempty(LLnodes_warehouse) then add cons, variables

    # LL load balance with PV (lambda)
    @constraint(Lower(model), loadbalance[n in LLnodes_withPV, t in 1:T], 
        -ye[n,t] + yi[n,t] + ypvprod[n,t] - LDFinputs.Pload[n][t] == 0
    );
    
    # duals of LL load balances
    @variable(Upper(model), 0 <= lambda <= M, DualOf(loadbalance));

    # LL load balance for warehouse (lambda_warehouse)
    @constraint(Lower(model), [n in LLnodes_warehouse, t in 1:T],
        yi[n, t] - LDFinputs.Pload[n][t] - ytherm[n, t] / CHILLER_COP == 0
    );

    # state space temperature evolution (lambda_ss)
    @constraint(Lower(model), [n in LLnodes_warehouse, t in 2:T],
        ytemperature[n, t] - ytemperature[n, t-1] - A[1, 1] * ytemperature[n, t-1] -
        sum(B[1, j] * u[j, t-1] for j=1:J) + B[1, 1] * ytherm[n, t-1] ==  0
    );
    # initial temperature (lambda_initTemperature)
    @constraint(Lower(model), [n in LLnodes_warehouse], ytemperature[n, 1] == -1.0);  

    
    # LL operational
    @constraint(Lower(model), [n in LLnodes_withPV,t in 1:T], 
        ypvprod[n,t] + spvprod[n,t] == ypv[n] * prod_factor[t] 
    )

    # UL constraints
    @constraint(Upper(model), [t in 1:T], x0[t] >= model[:Pⱼ]["0", t] );

    @constraint(Upper(model), [n in ULnodes_withBESS],
        xsoc[n,0] == 0.5 * xbkWh[n]
    );
    @constraint(Upper(model), [n in ULnodes_withBESS],
        xsoc[n,T] == 0.5 * xbkWh[n]
    );
    @constraint(Upper(model), [n in ULnodes_withBESS, t in 1:T],
        xsoc[n,t] == xsoc[n,t-1] + xbplus[n,t] * η - xbminus[n,t] / η
    );
    @constraint(Upper(model), [n in ULnodes_withBESS, t in 1:T],
        xbkW[n] >= xbminus[n,t]
    );
    @constraint(Upper(model), [n in ULnodes_withBESS, t in 1:T],
        xbkW[n] >= xbplus[n,t]
    );
    @constraint(Upper(model), [n in ULnodes_withBESS, t in 1:T],
        xbkW[n] >= xbplus[n,t] + xbminus[n,t]
    );
    @constraint(Upper(model), [n in ULnodes_withBESS, t in 1:T],
        xbkWh[n] >= xsoc[n,t]
    );

    @objective(Lower(model), Min, 
        pwf * sum(
            ci[t] * yi[n,t] - xe[n,t] * ye[n,t]
            for n in LLnodes_withPV, t in 1:T
        )
        + sum(ypv[n] * cpv for n in LLnodes_withPV)
        + pwf * sum(
            (ci[t] + xi[n,t]) * yi[n,t]
            for n in LLnodes_warehouse, t in 1:T
        )
    );

    @objective(Upper(model), Min, 
        5*pwf * sum(x0[t] * clmp[t] for t in 1:T)
        + sum(cbkW * xbkW[n] + cbkWh * xbkWh[n] for n in ULnodes_withBESS)
        + pwf * sum( lambda[n,t] * ye[n,t] for n in LLnodes_withPV, t in 1:T)
    );

    optimize!(model)

    return model
end
