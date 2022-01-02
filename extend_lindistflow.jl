#=
    Methods to extend LinDistFlow. Allows for defining import/export decision variables for 
    node loads/injections.
=#

function LDF.build_ldf!(m::JuMP.AbstractModel, p::LDF.Inputs, ps::AbstractVector,
    exportvarrefs::JuMP.Containers.DenseAxisArray, 
    importvarrefs::JuMP.Containers.DenseAxisArray
    )
    
    LDF.add_variables(m, p)
    LDF.constrain_power_balance(m, p)
    LDF.constrain_substation_voltage(m, p)
    LDF.constrain_KVL(m, p)
    LDF.constrain_loads(m, p, ps, exportvarrefs, importvarrefs)

end


"""
    constrain_loads(m, p::Inputs)

- set loads to negative of Inputs.Pload, which are normalized by Sbase when creating Inputs
- keys of Pload must match Inputs.busses. Any missing keys have load set to zero.
- Inputs.substation_bus is unconstrained, slack bus
"""
function LDF.constrain_loads(m::JuMP.AbstractModel, p::LDF.Inputs, ps::AbstractVector, 
    exportvarrefs::JuMP.Containers.DenseAxisArray, 
    importvarrefs::JuMP.Containers.DenseAxisArray,
    )

    Pⱼ = m[:Pⱼ]
    Qⱼ = m[:Qⱼ]
    # positive values are injections

    for j in p.busses
        if j in keys(p.Pload)
            if j in ps
                if j in exportvarrefs.axes[1]
                    @constraint(m, [t in 1:p.Ntimesteps],
                        Pⱼ[j,t] == exportvarrefs[j, t] - importvarrefs[j, t]
                    )
                else
                    @constraint(m, [t in 1:p.Ntimesteps],
                        Pⱼ[j,t] == -importvarrefs[j, t]
                    )
                end
            else
                @constraint(m, [t in 1:p.Ntimesteps],
                    Pⱼ[j,t] == -p.Pload[j][t]
                )
            end
        elseif j != p.substation_bus # TODO add BESS decisions here
            @constraint(m, [t in 1:p.Ntimesteps],
                Pⱼ[j,t] == 0
            )
        end
        
        if j in keys(p.Qload)
            if j in ps
                if j in exportvarrefs.axes[1]
                    @constraint(m, [t in 1:p.Ntimesteps],
                        Qⱼ[j,t] == exportvarrefs[j, t] - importvarrefs[j, t]
                    )
                else
                    @constraint(m, [t in 1:p.Ntimesteps],
                        Qⱼ[j,t] == -importvarrefs[j, t]
                    )
                end
            else
                @constraint(m, [t in 1:p.Ntimesteps],
                    Qⱼ[j,t] == -p.Qload[j][t]
                )
            end
        elseif j != p.substation_bus  # TODO add BESS decisions here
            @constraint(m, [t in 1:p.Ntimesteps],
                Qⱼ[j,t] == 0
            )
        end
    end
    p.Nequality_cons += 2 * (p.Nnodes - 1) * p.Ntimesteps
end
