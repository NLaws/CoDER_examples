"""
    bilinear_problem_bess_bigM(optimizer, T, LDFinputs, bigM)

Bilevel planning example without linearization of bilinear terms in upper objective.
Lower level problem is replaced with the KKT conditions.
"""
function bilinear_problem_bess_bigM(optimizer, T, LDFinputs, bigM)

    model = JuMP.Model(optimizer)
    set_optimizer_attribute(model, "MIPGap", 5e-2)
    set_optimizer_attribute(model, "NonConvex", 2)

    @variables model begin
        bigM >= yi[LLnodes, 1:T] >= 0
        bigM >= ye[LLnodes_withPV, 1:T] >= 0
        bigM >= ypv[LLnodes_withPV] >=0
        bigM >= ypvprod[LLnodes_withPV, 1:T] >= 0
        T_hi >= ytemperature[LLnodes_warehouse, 1:T] >= T_lo
        ytherm[LLnodes_warehouse, 1:T] >= 0

        bigM >= xe[LLnodes_withPV, 1:T] >= 0        
        bigM >= xi[LLnodes_warehouse, 1:T] >= 0
        bigM >= x0[1:T] >= 0
        bigM >= xbkW[ULnodes_withBESS] >= 0
        bigM >= xbkWh[ULnodes_withBESS] >= 0
        bigM >= xsoc[ULnodes_withBESS, 0:T] >= 0
        bigM >= xbplus[ULnodes_withBESS, 1:T] >= 0
        bigM >= xbminus[ULnodes_withBESS, 1:T] >= 0

        bigM >= lambda[LLnodes_withPV, 1:T] >= -bigM
        bigM >= lambda_warehouse[LLnodes_warehouse, 1:T] >= -bigM
        bigM >= lambda_ss[LLnodes_warehouse, 2:T] >= -bigM
        bigM >= lambda_initTemperature >= -bigM
        bigM >= mu_i[LLnodes, 1:T] >= 0
        bigM >= mu_e[LLnodes_withPV, 1:T] >= 0
        bigM >= mu_pv[LLnodes_withPV] >= 0
        bigM >= mu_pvprod[LLnodes_withPV, 1:T] >= 0
        bigM >= mu_dd[LLnodes_withPV, 1:T] >= 0
        bigM >= mu_therm_lo[LLnodes_warehouse, 1:T] >= 0
        bigM >= mu_therm_hi[LLnodes_warehouse, 1:T] >= 0
        bigM >= mu_temperature_lo[LLnodes_warehouse, 1:T] >= 0
        bigM >= mu_temperature_hi[LLnodes_warehouse, 1:T] >= 0
    end

    # UL does not allow simultaneous export/import
    @variable(model, byeyi[LLnodes_withPV, t in 1:T], Bin);
    @constraint(model, [n in LLnodes_withPV, t in 1:T],
        ye[n,t] <= bigM * byeyi[n,t]
    );
    @constraint(model, [n in LLnodes_withPV, t in 1:T],
        yi[n,t] <= bigM * (1-byeyi[n,t])
    );
    
    # LinDistFlow (single phase, real power only)
    LDF.build_ldf!(model, LDFinputs, LLnodes, ye, yi, xbplus, xbminus);
    # TODO if !isempty(LLnodes_warehouse) then add cons, variables

    # Complementary slackness of KKT, only modeling lower bound in most cases
    @variable(model, bypv[n in LLnodes_withPV], Bin);
    @constraint(model, [n in LLnodes_withPV],
        ypv[n] <= bigM * bypv[n]
    );
    @constraint(model, [n in LLnodes_withPV],
        mu_pv[n] <= bigM * (1-bypv[n])
    );
    @variable(model, bye[n in LLnodes_withPV, t in 1:T], Bin);
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        ye[n,t] <= bigM * bye[n,t] 
    );
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        mu_e[n,t] <= bigM * (1- bye[n,t])
    );
    @variable(model, byi[n in LLnodes, t in 1:T], Bin);
    @constraint(model, [n in LLnodes, t in 1:T], 
        yi[n,t] <= bigM * byi[n,t]
    );
    @constraint(model, [n in LLnodes, t in 1:T], 
        mu_i[n,t] <= bigM * (1 - byi[n,t])
    );
    @variable(model, bytherm_lo[n in LLnodes_warehouse, t in 1:T], Bin);
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        ytherm[n,t] <= bigM * bytherm_lo[n,t]
    );
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        mu_therm_lo[n,t] <= bigM * (1-bytherm_lo[n,t])
    );
    @variable(model, bytherm_hi[n in LLnodes_warehouse, t in 1:T], Bin);
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        ytherm[n,t] - bigM <= bigM * bytherm_hi[n,t]
    );
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        mu_therm_hi[n,t] <=  bigM * (1 - bytherm_hi[n,t])
    );
    @variable(model, bytemperature_lo[n in LLnodes_warehouse, t in 1:T], Bin);
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        T_lo - ytemperature[n,t] <= bigM * bytemperature_lo[n,t]
    );
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        mu_temperature_lo[n,t] <= bigM * (1 - bytemperature_lo[n,t])
    );
    @variable(model, bytemperature_hi[n in LLnodes_warehouse, t in 1:T], Bin);
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        ytemperature[n,t] - T_hi <= bigM * bytemperature_hi[n,t]
    );
    @constraint(model, [n in LLnodes_warehouse, t in 1:T], 
        mu_temperature_hi[n,t] <= bigM * (1 - bytemperature_hi[n,t])
    );
    @variable(model, bypvprod[n in LLnodes_withPV, t in 1:T], Bin)
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        ypvprod[n,t] <= bigM * bypvprod[n,t]
    );
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        mu_pvprod[n,t] <= bigM * (1 - bypvprod[n,t])
    )
    @variable(model, bydd[n in LLnodes_withPV, t in 1:T], Bin)
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        ypvprod[n,t] - ypv[n] * prod_factor[t] <= bigM * bydd[n,t]
    )
    @constraint(model, [n in LLnodes_withPV, t in 1:T], 
        mu_dd[n,t] <= bigM * (1 - bydd[n,t])
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
        ci[t] + xi[n,t] - lambda_warehouse[n, t] - mu_i[n,t] == 0);  # dual constraint of yi[t] for LLnodes_warehouse

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
        5*pwf * sum(x0[t] * clmp[t] for t in 1:T)
        + sum(cbkW * xbkW[n] + cbkWh * xbkWh[n] for n in ULnodes_withBESS)
        + pwf * sum(
            lambda[n,t] * ye[n,t]
            for n in LLnodes_withPV, t in 1:T
        )
    );

    optimize!(model)

    return model
end
