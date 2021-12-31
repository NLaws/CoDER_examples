#=
    Methods to extend LinDistFlow. Allows for defining import/export decision variables for 
    node loads/injections.
=#

function LDF.build_ldf!(m::JuMP.AbstractModel, p::LDF.Inputs, ps::Vector{Int},
    exportvarrefs::JuMP.Containers.DenseAxisArray, 
    importvarrefs::JuMP.Containers.DenseAxisArray
    )
    
    LDF.add_variables(m, p)
    LDF.constrain_power_balance(m, p)
    LDF.constrain_substation_voltage(m, p)
    LDF.constrain_KVL(m, p)
    LDF.constrain_loads(m, p, ps, exportvarrefs, importvarrefs)

end
# TODO confirm all equality constraints in LDF (or make equality with slack and method extensions)


function LDF.constrain_loads(m::JuMP.AbstractModel, p::LDF.Inputs, ps::Vector{Int}, 
    exportvarrefs::JuMP.Containers.DenseAxisArray, 
    importvarrefs::JuMP.Containers.DenseAxisArray,
    )

    Pⱼ = m[:Pⱼ]
    Qⱼ = m[:Qⱼ]
    # positive values are injections

    for j in p.busses
        if j in keys(p.Pload)
            if parse(Int, j) in ps
                @constraint(m, [t in 1:p.Ntimesteps],
                    Pⱼ[j,t] == 1e3/p.Sbase * (  # 1e3 b/c REopt values in kW
                        exportvarrefs[parse(Int, j), t]
                        - importvarrefs[parse(Int, j), t]
                    )
                )
            else
                @constraint(m, [t in 1:p.Ntimesteps],
                    Pⱼ[j,t] == -p.Pload[j][t]
                )
            end
        elseif j != p.substation_bus
            @constraint(m, [t in 1:p.Ntimesteps],
                Pⱼ[j,t] == 0
            )
        end
        
        if j in keys(p.Qload)
            if parse(Int, j) in ps
                @constraint(m, [t in 1:p.Ntimesteps],
                    Qⱼ[j,t] == 1e3/p.Sbase * p.pf * (  # 1e3 b/c REopt values in kW
                        exportvarrefs[parse(Int, j), t]
                        - importvarrefs[parse(Int, j), t]
                    )
                )
            else
                @constraint(m, [t in 1:p.Ntimesteps],
                    Qⱼ[j,t] == -p.Qload[j][t]
                )
            end
        elseif j != p.substation_bus
            @constraint(m, [t in 1:p.Ntimesteps],
                Qⱼ[j,t] == 0
            )
        end
    end
    p.Nequality_cons += 2 * (p.Nnodes - 1) * p.Ntimesteps
end
