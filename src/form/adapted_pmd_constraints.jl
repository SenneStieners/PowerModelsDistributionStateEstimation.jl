################################################################################
#  Copyright 2020, Marta Vanin, Tom Van Acker                                  #
################################################################################
# PowerModelsDistributionStateEstimation.jl                                    #
# An extention package of PowerModelsDistribution.jl for Static Power System   #
# State Estimation.                                                            #
################################################################################

@enum ConnConfig WYE DELTA

function constraint_mc_generator_power_se(pm::_PMD.AbstractUnbalancedIVRModel, id::Int; nw::Int=_IM.nw_id_default, report::Bool=true, bounded::Bool=true)
    generator = _PMD.ref(pm, nw, :gen, id)
    bus =  _PMD.ref(pm, nw,:bus, generator["gen_bus"])

    N = length(generator["connections"])
    pmin = get(generator, "pmin", fill(-Inf, N))
    pmax = get(generator, "pmax", fill( Inf, N))
    qmin = get(generator, "qmin", fill(-Inf, N))
    qmax = get(generator, "qmax", fill( Inf, N))

    if get(generator, "configuration", WYE) == _PMD.WYE
        constraint_mc_generator_power_wye_se(pm, nw, id, bus["index"], generator["connections"], pmin, pmax, qmin, qmax; report=report, bounded=bounded)
    else
        constraint_mc_generator_power_delta_se(pm, nw, id, bus["index"], generator["connections"], pmin, pmax, qmin, qmax; report=report, bounded=bounded)
    end
end

"wye connected generator setpoint constraint for IVR formulation - SE adaptation"
function constraint_mc_generator_power_wye_se(pm::_PMD.AbstractUnbalancedIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector, pmax::Vector, qmin::Vector, qmax::Vector; report::Bool=true, bounded::Bool=true)
    vr =  _PMD.var(pm, nw, :vr, bus_id)
    vi =  _PMD.var(pm, nw, :vi, bus_id)
    crg =  _PMD.var(pm, nw, :crg, id)
    cig =  _PMD.var(pm, nw, :cig, id)

    if bounded
        for (idx, c) in enumerate(connections)
            if pmin[c]>-Inf
                JuMP.@constraint(pm.model, pmin[idx] .<= vr[c]*crg[c]  + vi[c]*cig[c])
            end
            if pmax[c]< Inf
                JuMP.@constraint(pm.model, pmax[idx] .>= vr[c]*crg[c]  + vi[c]*cig[c])
            end
            if qmin[c]>-Inf
                JuMP.@constraint(pm.model, qmin[idx] .<= vi[c]*crg[c]  - vr[c]*cig[c])
            end
            if qmax[c]< Inf
                JuMP.@constraint(pm.model, qmax[idx] .>= vi[c]*crg[c]  - vr[c]*cig[c])
            end
        end
    end
end

"delta connected generator setpoint constraint for IVR formulation - adapted for SE"
function constraint_mc_generator_power_delta_se(pm::_PMD.AbstractUnbalancedIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector, pmax::Vector, qmin::Vector, qmax::Vector; report::Bool=true, bounded::Bool=true)
    vr =  _PMD.var(pm, nw, :vr, bus_id)
    vi =  _PMD.var(pm, nw, :vi, bus_id)
    crg =  _PMD.var(pm, nw, :crg, id)
    cig =  _PMD.var(pm, nw, :cig, id)

    nph = length(pmin)

    prev = Dict(c=>connections[(idx+nph-2)%nph+1] for (idx,c) in enumerate(connections))
    next = Dict(c=>connections[idx%nph+1] for (idx,c) in enumerate(connections))

    vrg = Dict()
    vig = Dict()
    for c in connections
        vrg[c] = JuMP.@expression(pm.model, vr[c]-vr[next[c]])
        vig[c] = JuMP.@expression(pm.model, vi[c]-vi[next[c]])
    end

    if bounded
        JuMP.@constraint(pm.model, [i in 1:nph], pmin[i] <= vrg[i]*crg[i]+vig[i]*cig[i])
        JuMP.@constraint(pm.model, [i in 1:nph], pmax[i] >= vrg[i]*crg[i]+vig[i]*cig[i])
        JuMP.@constraint(pm.model, [i in 1:nph], qmin[i] <= -vrg[i]*cig[i]+vig[i]*crg[i])
        JuMP.@constraint(pm.model, [i in 1:nph], qmax[i] >= -vrg[i]*cig[i]+vig[i]*crg[i])
    end
end

function constraint_mc_current_balance_se(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=_IM.nw_id_default)
    bus = _PMD.ref(pm, nw, :bus, i)
    bus_arcs = _PMD.ref(pm, nw, :bus_arcs_conns_branch, i)
    bus_arcs_sw = _PMD.ref(pm, nw, :bus_arcs_conns_switch, i)
    bus_arcs_trans = _PMD.ref(pm, nw, :bus_arcs_conns_transformer, i)
    bus_gens = _PMD.ref(pm, nw, :bus_conns_gen, i)
    bus_storage = _PMD.ref(pm, nw, :bus_conns_storage, i)
    bus_loads = _PMD.ref(pm, nw, :bus_conns_load, i)
    bus_shunts = _PMD.ref(pm, nw, :bus_conns_shunt, i)

    constraint_mc_current_balance_se(pm, nw, i, bus["terminals"], bus["grounded"], bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_shunts)
end

function constraint_mc_current_balance_se(pm::_PMD.IVRENPowerModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    #NB only difference with pmd is crd_bus replaced by crd, and same with cid
    vr = _PMD.var(pm, nw, :vr, i)
    vi = _PMD.var(pm, nw, :vi, i)

    cr    = get(_PMD.var(pm, nw),    :cr_bus, Dict()); _PMD._check_var_keys(cr, bus_arcs, "real current", "branch")
    ci    = get(_PMD.var(pm, nw),    :ci_bus, Dict()); _PMD._check_var_keys(ci, bus_arcs, "imaginary current", "branch")
    crd   = get(_PMD.var(pm, nw),   :crd, Dict()); _PMD._check_var_keys(crd, bus_loads, "real current", "load")
    cid   = get(_PMD.var(pm, nw),   :cid, Dict()); _PMD._check_var_keys(cid, bus_loads, "imaginary current", "load")
    crg   = get(_PMD.var(pm, nw),   :crg, Dict()); _PMD._check_var_keys(crg, bus_gens, "real current", "generator")
    cig   = get(_PMD.var(pm, nw),   :cig, Dict()); _PMD._check_var_keys(cig, bus_gens, "imaginary current", "generator")
    crs   = get(_PMD.var(pm, nw),   :crs, Dict()); _PMD._check_var_keys(crs, bus_storage, "real current", "storage")
    cis   = get(_PMD.var(pm, nw),   :cis, Dict()); _PMD._check_var_keys(cis, bus_storage, "imaginary current", "storage")
    crsw  = get(_PMD.var(pm, nw),  :crsw, Dict()); _PMD._check_var_keys(crsw, bus_arcs_sw, "real current", "switch")
    cisw  = get(_PMD.var(pm, nw),  :cisw, Dict()); _PMD._check_var_keys(cisw, bus_arcs_sw, "imaginary current", "switch")
    crt   = get(_PMD.var(pm, nw),   :crt_bus, Dict()); _PMD._check_var_keys(crt, bus_arcs_trans, "real current", "transformer")
    cit   = get(_PMD.var(pm, nw),   :cit_bus, Dict()); _PMD._check_var_keys(cit, bus_arcs_trans, "imaginary current", "transformer")

    Gs, Bs = _PMD._build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx, t) in ungrounded_terminals
        JuMP.@constraint(pm.model,  sum(cr[a][t] for (a, conns) in bus_arcs if t in conns)
                                    + sum(crsw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
                                    + sum(crt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
                                    ==
                                      sum(crg[g][t]         for (g, conns) in bus_gens if t in conns)
                                    - sum(crs[s][t]         for (s, conns) in bus_storage if t in conns)
                                    - sum(crd[d][t]         for (d, conns) in bus_loads if t in conns)
                                    - sum( Gs[idx,jdx]*vr[u] -Bs[idx,jdx]*vi[u] for (jdx,u) in ungrounded_terminals) # shunts
                                    )
        JuMP.@constraint(pm.model, sum(ci[a][t] for (a, conns) in bus_arcs if t in conns)
                                    + sum(cisw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
                                    + sum(cit[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
                                    ==
                                      sum(cig[g][t]         for (g, conns) in bus_gens if t in conns)
                                    - sum(cis[s][t]         for (s, conns) in bus_storage if t in conns)
                                    - sum(cid[d][t]         for (d, conns) in bus_loads if t in conns)
                                    - sum( Gs[idx,jdx]*vi[u] +Bs[idx,jdx]*vr[u] for (jdx,u) in ungrounded_terminals) # shunts
                                    )
    end
end
#####
# "KCL including transformer arcs and load variables."
function constraint_mc_power_balance_se(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=_IM.nw_id_default)
    bus = _PMD.ref(pm, nw, :bus, i)
    bus_arcs = _PMD.ref(pm, nw, :bus_arcs_conns_branch, i)
    bus_arcs_sw = _PMD.ref(pm, nw, :bus_arcs_conns_switch, i)
    bus_arcs_trans = _PMD.ref(pm, nw, :bus_arcs_conns_transformer, i)
    bus_gens = _PMD.ref(pm, nw, :bus_conns_gen, i)
    bus_storage = _PMD.ref(pm, nw, :bus_conns_storage, i)
    bus_loads = _PMD.ref(pm, nw, :bus_conns_load, i)
    bus_shunts = _PMD.ref(pm, nw, :bus_conns_shunt, i)

    constraint_mc_power_balance_se(pm, nw, i, bus["terminals"], bus["grounded"], bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_shunts)
end

function constraint_mc_power_balance_se(pm::_PMD.AbstractUnbalancedACRModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    #NB only diffeerence is in pd and qd we refer to :qd, :pd instead of :pd_bus, :qd_bus
    vr = _PMD.var(pm, nw, :vr, i)
    vi = _PMD.var(pm, nw, :vi, i)

    p    = get(_PMD.var(pm, nw), :p,      Dict()); _PMD._check_var_keys(p,   bus_arcs,       "active power",   "branch")
    q    = get(_PMD.var(pm, nw), :q,      Dict()); _PMD._check_var_keys(q,   bus_arcs,       "reactive power", "branch")
    pg   = get(_PMD.var(pm, nw), :pg, Dict()); _PMD._check_var_keys(pg,  bus_gens,       "active power",   "generator")
    qg   = get(_PMD.var(pm, nw), :qg, Dict()); _PMD._check_var_keys(qg,  bus_gens,       "reactive power", "generator")
    ps   = get(_PMD.var(pm, nw), :ps,     Dict()); _PMD._check_var_keys(ps,  bus_storage,    "active power",   "storage")
    qs   = get(_PMD.var(pm, nw), :qs,     Dict()); _PMD._check_var_keys(qs,  bus_storage,    "reactive power", "storage")
    psw  = get(_PMD.var(pm, nw), :psw,    Dict()); _PMD._check_var_keys(psw, bus_arcs_sw,    "active power",   "switch")
    qsw  = get(_PMD.var(pm, nw), :qsw,    Dict()); _PMD._check_var_keys(qsw, bus_arcs_sw,    "reactive power", "switch")
    pt   = get(_PMD.var(pm, nw), :pt,     Dict()); _PMD._check_var_keys(pt,  bus_arcs_trans, "active power",   "transformer")
    qt   = get(_PMD.var(pm, nw), :qt,     Dict()); _PMD._check_var_keys(qt,  bus_arcs_trans, "reactive power", "transformer")
    pd   = get(_PMD.var(pm, nw), :pd, Dict()); _PMD._check_var_keys(pd,  bus_loads,      "active power",   "load")
    qd   = get(_PMD.var(pm, nw), :qd, Dict()); _PMD._check_var_keys(pd,  bus_loads,      "reactive power", "load")

    Gs, Bs = _PMD._build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    # pd/qd can be NLexpressions, so cannot be vectorized
    for (idx, t) in ungrounded_terminals
        cp = JuMP.@constraint(pm.model, [p, q, pg, qg, ps, qs, psw, qsw, pt, qt, pd, qd, vr, vi],
              sum(  p[arc][t] for (arc, conns) in bus_arcs if t in conns)
            + sum(psw[arc][t] for (arc, conns) in bus_arcs_sw if t in conns)
            + sum( pt[arc][t] for (arc, conns) in bus_arcs_trans if t in conns)
            ==
              sum(pg[gen][t] for (gen, conns) in bus_gens if t in conns)
            - sum(ps[strg][t] for (strg, conns) in bus_storage if t in conns)
            - sum(pd[load][t] for (load, conns) in bus_loads if t in conns)
            + ( -vr[t] * sum(Gs[idx,jdx]*vr[u]-Bs[idx,jdx]*vi[u] for (jdx,u) in ungrounded_terminals)
                -vi[t] * sum(Gs[idx,jdx]*vi[u]+Bs[idx,jdx]*vr[u] for (jdx,u) in ungrounded_terminals)
            )
        )
        push!(cstr_p, cp)

        cq = JuMP.@constraint(pm.model, [p, q, pg, qg, ps, qs, psw, qsw, pt, qt, pd, qd, vr, vi],
              sum(  q[arc][t] for (arc, conns) in bus_arcs if t in conns)
            + sum(qsw[arc][t] for (arc, conns) in bus_arcs_sw if t in conns)
            + sum( qt[arc][t] for (arc, conns) in bus_arcs_trans if t in conns)
            ==
              sum(qg[gen][t] for (gen, conns) in bus_gens if t in conns)
            - sum(qd[load][t] for (load, conns) in bus_loads if t in conns)
            - sum(qs[strg][t] for (strg, conns) in bus_storage if t in conns)
            + ( vr[t] * sum(Gs[idx,jdx]*vi[u]+Bs[idx,jdx]*vr[u] for (jdx,u) in ungrounded_terminals)
               -vi[t] * sum(Gs[idx,jdx]*vr[u]-Bs[idx,jdx]*vi[u] for (jdx,u) in ungrounded_terminals)
            )
        )
        push!(cstr_q, cq)
    end
end

function constraint_mc_power_balance_se(pm::_PMD.AbstractUnbalancedACPModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vm   = _PMD.var(pm, nw, :vm, i)
    va   = _PMD.var(pm, nw, :va, i)

    p    = get(_PMD.var(pm, nw), :p,      Dict()); _PMD._check_var_keys(p,   bus_arcs,       "active power",   "branch")
    q    = get(_PMD.var(pm, nw), :q,      Dict()); _PMD._check_var_keys(q,   bus_arcs,       "reactive power", "branch")
    pg   = get(_PMD.var(pm, nw), :pg, Dict()); _PMD._check_var_keys(pg,  bus_gens,       "active power",   "generator")
    qg   = get(_PMD.var(pm, nw), :qg, Dict()); _PMD._check_var_keys(qg,  bus_gens,       "reactive power", "generator")
    ps   = get(_PMD.var(pm, nw), :ps,     Dict()); _PMD._check_var_keys(ps,  bus_storage,    "active power",   "storage")
    qs   = get(_PMD.var(pm, nw), :qs,     Dict()); _PMD._check_var_keys(qs,  bus_storage,    "reactive power", "storage")
    psw  = get(_PMD.var(pm, nw), :psw,    Dict()); _PMD._check_var_keys(psw, bus_arcs_sw,    "active power",   "switch")
    qsw  = get(_PMD.var(pm, nw), :qsw,    Dict()); _PMD._check_var_keys(qsw, bus_arcs_sw,    "reactive power", "switch")
    pt   = get(_PMD.var(pm, nw), :pt,     Dict()); _PMD._check_var_keys(pt,  bus_arcs_trans, "active power",   "transformer")
    qt   = get(_PMD.var(pm, nw), :qt,     Dict()); _PMD._check_var_keys(qt,  bus_arcs_trans, "reactive power", "transformer")
    pd   = get(_PMD.var(pm, nw), :pd, Dict()); _PMD._check_var_keys(pd,  bus_loads,      "active power",   "load")
    qd   = get(_PMD.var(pm, nw), :qd, Dict()); _PMD._check_var_keys(pd,  bus_loads,      "reactive power", "load")

    Gs, Bs = _PMD._build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx,t) in ungrounded_terminals
        if any(Bs[idx,jdx] != 0 for (jdx, u) in ungrounded_terminals if idx != jdx) || any(Gs[idx,jdx] != 0 for (jdx, u) in ungrounded_terminals if idx != jdx)
            cp = JuMP.@constraint(pm.model,
                  sum(  p[a][t] for (a, conns) in bus_arcs if t in conns)
                + sum(psw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                + sum( pt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                - sum( pg[g][t] for (g, conns) in bus_gens if t in conns)
                + sum( ps[s][t] for (s, conns) in bus_storage if t in conns)
                + sum( pd[l][t] for (l, conns) in bus_loads if t in conns)
                + ( # shunt
                    +Gs[idx,idx] * vm[t]^2
                    +sum( Gs[idx,jdx] * vm[t]*vm[u] * cos(va[t]-va[u])
                         +Bs[idx,jdx] * vm[t]*vm[u] * sin(va[t]-va[u])
                        for (jdx,u) in ungrounded_terminals if idx != jdx)
                )
                ==
                0.0
            )
            push!(cstr_p, cp)

            cq = JuMP.@constraint(pm.model,
                  sum(  q[a][t] for (a, conns) in bus_arcs if t in conns)
                + sum(qsw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                + sum( qt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                - sum( qg[g][t] for (g, conns) in bus_gens if t in conns)
                + sum( qs[s][t] for (s, conns) in bus_storage if t in conns)
                + sum( qd[l][t] for (l, conns) in bus_loads if t in conns)
                + ( # shunt
                    -Bs[idx,idx] * vm[t]^2
                    -sum( Bs[idx,jdx] * vm[t]*vm[u] * cos(va[t]-va[u])
                         -Gs[idx,jdx] * vm[t]*vm[u] * sin(va[t]-va[u])
                         for (jdx,u) in ungrounded_terminals if idx != jdx)
                )
                ==
                0.0
            )
            push!(cstr_q, cq)
        else
            cp = JuMP.@constraint(pm.model, [p, pg, ps, psw, pt, pd, vm],
                  sum(  p[a][t] for (a, conns) in bus_arcs if t in conns)
                + sum(psw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                + sum( pt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                - sum( pg[g][t] for (g, conns) in bus_gens if t in conns)
                + sum( ps[s][t] for (s, conns) in bus_storage if t in conns)
                + sum( pd[l][t] for (l, conns) in bus_loads if t in conns)
                + Gs[idx,idx] * vm[t]^2
                ==
                0.0
            )
            push!(cstr_p, cp)

            cq = JuMP.@constraint(pm.model, [q, qg, qs, qsw, qt, qd, vm],
                  sum(  q[a][t] for (a, conns) in bus_arcs if t in conns)
                + sum(qsw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                + sum( qt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                - sum( qg[g][t] for (g, conns) in bus_gens if t in conns)
                + sum( qs[s][t] for (s, conns) in bus_storage if t in conns)
                + sum( qd[l][t] for (l, conns) in bus_loads if t in conns)
                - Bs[idx,idx] * vm[t]^2
                ==
                0.0
            )
            push!(cstr_q, cq)
        end
    end
end

function constraint_mc_power_balance_se(pm::_PMD.SDPUBFPowerModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    Wr = _PMD.var(pm, nw, :Wr, i)
    Wi = _PMD.var(pm, nw, :Wi, i)
    P = get(_PMD.var(pm, nw), :P, Dict()); _PMD._check_var_keys(P, bus_arcs, "active power", "branch")
    Q = get(_PMD.var(pm, nw), :Q, Dict()); _PMD._check_var_keys(Q, bus_arcs, "reactive power", "branch")
    Psw  = get(_PMD.var(pm, nw),  :Psw, Dict()); _PMD._check_var_keys(Psw, bus_arcs_sw, "active power", "switch")
    Qsw  = get(_PMD.var(pm, nw),  :Qsw, Dict()); _PMD._check_var_keys(Qsw, bus_arcs_sw, "reactive power", "switch")
    Pt   = get(_PMD.var(pm, nw),   :Pt, Dict()); _PMD._check_var_keys(Pt, bus_arcs_trans, "active power", "transformer")
    Qt   = get(_PMD.var(pm, nw),   :Qt, Dict()); _PMD._check_var_keys(Qt, bus_arcs_trans, "reactive power", "transformer")

    pd = get(_PMD.var(pm, nw), :pd, Dict()); _PMD._check_var_keys(pd, bus_loads, "active power", "load")
    qd = get(_PMD.var(pm, nw), :qd, Dict()); _PMD._check_var_keys(qd, bus_loads, "reactive power", "load")
    pg = get(_PMD.var(pm, nw), :pg, Dict()); _PMD._check_var_keys(pg, bus_gens, "active power", "generator")
    qg = get(_PMD.var(pm, nw), :qg, Dict()); _PMD._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(_PMD.var(pm, nw),   :ps, Dict()); _PMD._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(_PMD.var(pm, nw),   :qs, Dict()); _PMD._check_var_keys(qs, bus_storage, "reactive power", "storage")

    Gs, Bs = _PMD._build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx,t) in ungrounded_terminals
        cp = JuMP.@constraint(pm.model,
            sum(diag(P[a])[findfirst(isequal(t), conns)] for (a, conns) in bus_arcs if t in conns)
            + sum(diag(Psw[a_sw])[findfirst(isequal(t), conns)] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(diag(Pt[a_trans])[findfirst(isequal(t), conns)] for (a_trans, conns) in bus_arcs_trans if t in conns)
            ==
            sum(pg[g][t] for (g, conns) in bus_gens if t in conns)
            - sum(ps[s][t] for (s, conns) in bus_storage if t in conns)
            - sum(pd[d][t] for (d, conns) in bus_loads if t in conns)
            - diag(Wr*Gs'+Wi*Bs')[idx]
        )
        push!(cstr_p, cp)

        cq = JuMP.@constraint(pm.model,
            sum(diag(Q[a])[findfirst(isequal(t), conns)] for (a, conns) in bus_arcs if t in conns)
            + sum(diag(Qsw[a_sw])[findfirst(isequal(t), conns)] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(diag(Qt[a_trans])[findfirst(isequal(t), conns)] for (a_trans, conns) in bus_arcs_trans if t in conns)
            ==
            sum(qg[g][t] for (g, conns) in bus_gens if t in conns)
            - sum(qs[s][t] for (s, conns) in bus_storage if t in conns)
            - sum(qd[d][t] for (d, conns) in bus_loads if t in conns)
            - diag(-Wr*Bs'+Wi*Gs')[idx]
        )
        push!(cstr_q, cq)
    end
end

function constraint_mc_power_balance_se(pm::_PMD.LPUBFDiagModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    w = _PMD.var(pm, nw, :w, i)
    p   = get(_PMD.var(pm, nw),      :p,   Dict()); _PMD._check_var_keys(p,   bus_arcs, "active power", "branch")
    q   = get(_PMD.var(pm, nw),      :q,   Dict()); _PMD._check_var_keys(q,   bus_arcs, "reactive power", "branch")
    psw = get(_PMD.var(pm, nw),    :psw, Dict()); _PMD._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw = get(_PMD.var(pm, nw),    :qsw, Dict()); _PMD._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt  = get(_PMD.var(pm, nw),     :pt,  Dict()); _PMD._check_var_keys(pt,  bus_arcs_trans, "active power", "transformer")
    qt  = get(_PMD.var(pm, nw),     :qt,  Dict()); _PMD._check_var_keys(qt,  bus_arcs_trans, "reactive power", "transformer")
    pg  = get(_PMD.var(pm, nw),     :pg,  Dict()); _PMD._check_var_keys(pg,  bus_gens, "active power", "generator")
    qg  = get(_PMD.var(pm, nw),     :qg,  Dict()); _PMD._check_var_keys(qg,  bus_gens, "reactive power", "generator")
    ps  = get(_PMD.var(pm, nw),     :ps,  Dict()); _PMD._check_var_keys(ps,  bus_storage, "active power", "storage")
    qs  = get(_PMD.var(pm, nw),     :qs,  Dict()); _PMD._check_var_keys(qs,  bus_storage, "reactive power", "storage")
    pd  = get(_PMD.var(pm, nw),     :pd,  Dict()); _PMD._check_var_keys(pd,  bus_loads, "active power", "load")
    qd  = get(_PMD.var(pm, nw),     :qd,  Dict()); _PMD._check_var_keys(qd,  bus_loads, "reactive power", "load")

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx,t) in ungrounded_terminals
        cp = JuMP.@constraint(pm.model,
              sum(  p[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(psw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
            + sum( pt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
            - sum( pg[g][t] for (g, conns) in bus_gens if t in conns)
            + sum( ps[s][t] for (s, conns) in bus_storage if t in conns)
            + sum( pd[d][t] for (d, conns) in bus_loads if t in conns)
            + sum(diag(ref(pm, nw, :shunt, sh, "gs"))[findfirst(isequal(t), conns)]*w[t] for (sh, conns) in bus_shunts if t in conns)
            ==
            0.0
        )
        push!(cstr_p, cp)

        cq = JuMP.@constraint(pm.model,
              sum(  q[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(qsw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
            + sum( qt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
            - sum( qg[g][t] for (g, conns) in bus_gens if t in conns)
            + sum( qs[s][t] for (s, conns) in bus_storage if t in conns)
            + sum( qd[d][t] for (d, conns) in bus_loads if t in conns)
            - sum(diag(_PMD.ref(pm, nw, :shunt, sh, "bs"))[findfirst(isequal(t), conns)]*w[t] for (sh, conns) in bus_shunts if t in conns)
            ==
            0.0
        )
        push!(cstr_q, cq)
   end
end

function variable_mc_transformer_tap(pm::_PMD.AbstractUnbalancedPowerModel;
    nw::Int=_IM.nw_id_default, bounded::Bool=true, report::Bool=true
)
    #p_oltc_ids = [id for (id, tr) in _PMD.ref(pm, nw, :transformer)
     #            if endswith(string(get(tr, "source_id", "")), ".2")]

    # enkel trafos waarvoor NIET alles fixed is (dus Bool[0,0,0] -> variabelen; Bool[1,1,1] -> skip)
    #p_oltc_var_ids = [i for i in p_oltc_ids
    #                  if !all(_PMD.ref(pm, nw, :transformer, i, "tm_fix"))]

    p_oltc_var_ids = [i for (i, tr) in _PMD.ref(pm, nw, :transformer)
                  if !all(get(tr, "tm_fix", Bool[1,1,1]))]
                    
    tap = _PMD.var(pm, nw)[:tap] = Dict(i => JuMP.@variable(pm.model,
        [p in 1:length(_PMD.ref(pm, nw, :transformer, i, "tm_set"))],
        base_name="$(nw)_tm_$(i)",
        start = _PMD.ref(pm, nw, :transformer, i, "tm_set")[p],
    ) for i in p_oltc_var_ids)

    if bounded
        for tr_id in p_oltc_var_ids, p in 1:length(_PMD.ref(pm, nw, :transformer, tr_id, "tm_set"))
            _PMD.set_lower_bound(_PMD.var(pm, nw)[:tap][tr_id][p], _PMD.ref(pm, nw, :transformer, tr_id, "tm_lb")[p])
            _PMD.set_upper_bound(_PMD.var(pm, nw)[:tap][tr_id][p], _PMD.ref(pm, nw, :transformer, tr_id, "tm_ub")[p])
        end
    end

    report && _IM.sol_component_value(pm, :pmd, nw, :transformer, :tap, p_oltc_var_ids, tap)
end

function constraint_mc_transformer_tap_time_invariant(pm::_PMD.AbstractUnbalancedPowerModel; nw_ref::Int=1)
    nws = [n for (n, _) in _PMD.nws(pm)]
    length(nws) <= 1 && return

    tap_ref_dict = get(_PMD.var(pm, nw_ref), :tap, nothing)
    tap_ref_dict === nothing && return

    for n in nws
        n == nw_ref && continue
        tap_n_dict = get(_PMD.var(pm, n), :tap, nothing)
        tap_n_dict === nothing && continue

        for (tr, tap_ref) in tap_ref_dict
            haskey(tap_n_dict, tr) || continue
            tap_n = tap_n_dict[tr]
            for p in 1:length(tap_ref)
                JuMP.@constraint(pm.model, tap_n[p] == tap_ref[p])
            end
        end
    end
end




"Calculates the tap scale factor for the non-dimensionalized equations."
function calculate_tm_scale(trans::Dict{String,Any}, bus_fr::Dict{String,Any}, bus_to::Dict{String,Any})
    tm_nom = trans["tm_nom"]

    f_vbase = haskey(bus_fr, "vbase") ? bus_fr["vbase"] : bus_fr["base_kv"]
    t_vbase = haskey(bus_to, "vbase") ? bus_to["vbase"] : bus_to["base_kv"]
    config = trans["configuration"]

    tm_scale = tm_nom
    if config == _PMD.DELTA
        #TODO is this still needed?
        #tm_scale *= sqrt(3)
    elseif config == "zig-zag"
        error("Zig-zag not yet supported.")
    end

    return tm_scale
end

"Enforces equal tap across phases for transformer i (only if a tap variable exists)"
function constraint_mc_transformer_tap_equal_phase(
    pm::_PMD.AbstractUnbalancedPowerModel,
    i::Int;
    nw::Int=_IM.nw_id_default
)
    tap_dict = get(_PMD.var(pm, nw), :tap, nothing)
    (tap_dict === nothing || !haskey(tap_dict, i)) && return

    tap = tap_dict[i]  # vector (per phase)
    for p in 2:length(tap)
        JuMP.@constraint(pm.model, tap[p] == tap[1])
    end
end

"Enforces tap of phase 1 to be fixed to a given value (only if a tap variable exists)"
function constraint_mc_transformer_tap_test(
    pm::_PMD.AbstractUnbalancedPowerModel,
    i::Int;
    nw::Int = _IM.nw_id_default
)
    tap_dict = get(_PMD.var(pm, nw), :tap, nothing)
    (tap_dict === nothing || !haskey(tap_dict, i)) && return

    tap = tap_dict[i]  # vector (per phase)

    JuMP.@constraint(pm.model, tap[1] == 1.0000007878)
end



function constraint_mc_transformer_voltage(pm::_PMD.ExplicitNeutralModels, i::Int; nw::Int=_IM.nw_id_default, fix_taps::Bool=true)
    transformer = _PMD.ref(pm, nw, :transformer, i)
    f_bus = transformer["f_bus"]
    t_bus = transformer["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    configuration = transformer["configuration"]
    f_connections = transformer["f_connections"]
    t_connections = transformer["t_connections"]
    tm_set = transformer["tm_set"]
    tm_fixed = fix_taps ? ones(Bool, length(tm_set)) : transformer["tm_fix"]
    tm_scale = calculate_tm_scale(transformer, _PMD.ref(pm, nw, :bus, f_bus), _PMD.ref(pm, nw, :bus, t_bus))

    #TODO change data model
    # there is redundancy in specifying polarity seperately on from and to side
    #TODO change this once migrated to new data model
    pol = transformer["polarity"]

    if configuration == _PMD.WYE
        constraint_mc_transformer_voltage_yy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_connections, t_connections, pol, tm_set, tm_fixed, tm_scale)
    elseif configuration == _PMD.DELTA
        constraint_mc_transformer_voltage_dy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_connections, t_connections, pol, tm_set, tm_fixed, tm_scale)
    elseif configuration == "zig-zag"
        error("Zig-zag not yet supported.")
    end
end

function constraint_mc_transformer_voltage_yy(pm::_PMD.RectangularVoltageExplicitNeutralModels, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    vr_fr_P = [_PMD.var(pm, nw, :vr, f_bus)[c] for c in f_connections[1:end-1]]
    vi_fr_P = [_PMD.var(pm, nw, :vi, f_bus)[c] for c in f_connections[1:end-1]]
    vr_fr_n = _PMD.var(pm, nw, :vr, f_bus)[f_connections[end]]
    vi_fr_n = _PMD.var(pm, nw, :vi, f_bus)[f_connections[end]]
    vr_to_P = [_PMD.var(pm, nw, :vr, t_bus)[c] for c in t_connections[1:end-1]]
    vi_to_P = [_PMD.var(pm, nw, :vi, t_bus)[c] for c in t_connections[1:end-1]]
    vr_to_n = _PMD.var(pm, nw, :vr, t_bus)[t_connections[end]]
    vi_to_n = _PMD.var(pm, nw, :vi, t_bus)[t_connections[end]]
    
    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : _PMD.var(pm, nw, :tap, trans_id)[idx] for idx in 1:length(tm_fixed)]
    scale = (tm_scale*pol).*tm

    JuMP.@constraint(pm.model, (vr_fr_P.-vr_fr_n) .== scale.*(vr_to_P.-vr_to_n))
    JuMP.@constraint(pm.model, (vi_fr_P.-vi_fr_n) .== scale.*(vi_to_P.-vi_to_n))
end
function constraint_mc_transformer_voltage_dy(pm::_PMD.RectangularVoltageExplicitNeutralModels, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    vr_fr_P = [_PMD.var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    vi_fr_P = [_PMD.var(pm, nw, :vi, f_bus)[c] for c in f_connections]
    vr_to_P = [_PMD.var(pm, nw, :vr, t_bus)[c] for c in t_connections[1:end-1]]
    vi_to_P = [_PMD.var(pm, nw, :vi, t_bus)[c] for c in t_connections[1:end-1]]
    vr_to_n = _PMD.var(pm, nw, :vr, t_bus)[t_connections[end]]
    vi_to_n = _PMD.var(pm, nw, :vi, t_bus)[t_connections[end]]

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : _PMD.var(pm, nw, :tap, trans_id)[idx] for idx in 1:length(tm_fixed)]
    scale = (tm_scale*pol).*tm

    n_phases = length(tm)
    Md = _get_delta_transformation_matrix(n_phases)

    JuMP.@constraint(pm.model, Md*vr_fr_P .== scale.*(vr_to_P .- vr_to_n))
    JuMP.@constraint(pm.model, Md*vi_fr_P .== scale.*(vi_to_P .- vi_to_n))
end


function constraint_mc_transformer_current(pm::_PMD.AbstractExplicitNeutralIVRModel, i::Int; nw::Int=_IM.nw_id_default, fix_taps::Bool=true)
    # if ref(pm, nw_id_default, :conductors)!=3
    #     error("Transformers only work with networks with three conductors.")
    # end

    transformer = _PMD.ref(pm, nw, :transformer, i)
    f_bus = transformer["f_bus"]
    t_bus = transformer["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    configuration = transformer["configuration"]
    f_connections = transformer["f_connections"]
    t_connections = transformer["t_connections"]
    tm_set = transformer["tm_set"]
    tm_fixed = fix_taps ? ones(Bool, length(tm_set)) : transformer["tm_fix"]
    tm_scale = calculate_tm_scale(transformer, _PMD.ref(pm, nw, :bus, f_bus), _PMD.ref(pm, nw, :bus, t_bus))

    #TODO change data model
    # there is redundancy in specifying polarity seperately on from and to side
    #TODO change this once migrated to new data model
    pol = transformer["polarity"]
    if configuration == _PMD.WYE
        constraint_mc_transformer_current_yy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_connections, t_connections, pol, tm_set, tm_fixed, tm_scale)
    elseif configuration == _PMD.DELTA
        constraint_mc_transformer_current_dy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_connections, t_connections, pol, tm_set, tm_fixed, tm_scale)
    elseif configuration == "zig-zag"
        error("Zig-zag not yet supported.")
    end
end

function constraint_mc_transformer_current_yy(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    cr_fr_P = _PMD.var(pm, nw, :crt, f_idx)
    ci_fr_P = _PMD.var(pm, nw, :cit, f_idx)
    cr_to_P = _PMD.var(pm, nw, :crt, t_idx)
    ci_to_P = _PMD.var(pm, nw, :cit, t_idx)
    
    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : _PMD.var(pm, nw, :tap, trans_id)[idx] for idx in 1:length(tm_fixed)]
    scale = (tm_scale*pol).*tm

    JuMP.@constraint(pm.model, scale.*cr_fr_P .+ cr_to_P .== 0)
    JuMP.@constraint(pm.model, scale.*ci_fr_P .+ ci_to_P .== 0)

    _PMD.var(pm, nw, :crt_bus)[f_idx] = _merge_bus_flows(pm, [cr_fr_P..., -sum(cr_fr_P)], f_connections)
    _PMD.var(pm, nw, :cit_bus)[f_idx] = _merge_bus_flows(pm, [ci_fr_P..., -sum(ci_fr_P)], f_connections)
    _PMD.var(pm, nw, :crt_bus)[t_idx] = _merge_bus_flows(pm, [cr_to_P..., -sum(cr_to_P)], t_connections)
    _PMD.var(pm, nw, :cit_bus)[t_idx] = _merge_bus_flows(pm, [ci_to_P..., -sum(ci_to_P)], t_connections)
end

function constraint_mc_transformer_current_dy(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    cr_fr_P = _PMD.var(pm, nw, :crt, f_idx)
    ci_fr_P = _PMD.var(pm, nw, :cit, f_idx)
    cr_to_P = _PMD.var(pm, nw, :crt, t_idx)
    ci_to_P = _PMD.var(pm, nw, :cit, t_idx)
    
    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : _PMD.var(pm, nw, :tap, trans_id)[idx] for idx in 1:length(tm_fixed)]
    scale = (tm_scale*pol).*tm

    n_phases = length(tm)
    Md = _get_delta_transformation_matrix(n_phases)

    JuMP.@constraint(pm.model, scale.*cr_fr_P .+ cr_to_P .== 0)
    JuMP.@constraint(pm.model, scale.*ci_fr_P .+ ci_to_P .== 0)

    _PMD.var(pm, nw, :crt_bus)[f_idx] = _merge_bus_flows(pm, Md'*cr_fr_P, f_connections)
    _PMD.var(pm, nw, :cit_bus)[f_idx] = _merge_bus_flows(pm, Md'*ci_fr_P, f_connections)
    _PMD.var(pm, nw, :crt_bus)[t_idx] = _merge_bus_flows(pm, [cr_to_P..., -sum(cr_to_P)], t_connections)
    _PMD.var(pm, nw, :cit_bus)[t_idx] = _merge_bus_flows(pm, [ci_to_P..., -sum(ci_to_P)], t_connections)
end

"Merges flow variables that enter the same terminals, i.e. multiple neutrals of an underground cable connected to same neutral terminal"
function _merge_bus_flows(pm::_PMD.AbstractExplicitNeutralIVRModel, flows::Vector, connections::Vector)::JuMP.Containers.DenseAxisArray
    flows_merged = []
    conns_unique = unique(connections)
    for t in conns_unique
        idxs = findall(connections.==t)
        flows_t = flows[idxs]
        if length(flows_t)==1
            flows_merged_t = flows_t[1]
        else
            flows_merged_t = sum(flows_t)
        end
        push!(flows_merged, flows_merged_t)
    end
    JuMP.Containers.DenseAxisArray(flows_merged, conns_unique)
end

"creates a delta transformation matrix"
function _get_delta_transformation_matrix(n_phases::Int)::Matrix{Int}
    @assert(n_phases>2, "We only define delta transforms for three and more conductors.")
    Md = LinearAlgebra.diagm(0=>fill(1, n_phases), 1=>fill(-1, n_phases-1))
    Md[end,1] = -1
    return Md
end


"Probeert de transformer-id te vinden die bij deze virtual branch hoort."
function _find_transformer_for_virtual_branch(pm, nw::Int, branch::Dict{String,Any})
    bname = string(get(branch, "name", ""))

    # verwacht bv: _virtual_branch.transformer.tx3170_1
    m = match(r"_virtual_branch\.transformer\.([^_]+)_", bname)
    if m === nothing
        # fallback: probeer op source_id
        bsrc = string(get(branch, "source_id", ""))
        m = match(r"_virtual_branch\.transformer\.([^_]+)_", bsrc)
        m === nothing && return nothing
    end

    tag = m.captures[1]   # bv "tx3170"

    # zoek transformer met tag in name of source_id
    for (tr_id, tr) in _PMD.ref(pm, nw, :transformer)
        trname = string(get(tr, "name", ""))
        trsrc  = string(get(tr, "source_id", ""))
        if occursin(tag, trname) || occursin(tag, trsrc)
            return tr_id
        end
    end

    return nothing
end



function constraint_mc_bus_voltage_drop(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=_IM.nw_id_default)::Nothing
    branch = _PMD.ref(pm, nw, :branch, i)

    r0 = branch["br_r"]
    x0 = branch["br_x"]

    bname = string(get(branch, "name", ""))
    is_virtual_tr = startswith(bname, "_virtual_branch.transformer.") && endswith(bname, "_1")

    # default: constant r,x (zoals PMD)
    r_use = r0
    x_use = x0

    if is_virtual_tr
        tr_id = get(branch, "transformer_id", nothing)
        tr_id === nothing && (tr_id = _find_transformer_for_virtual_branch(pm, nw, branch))

        if tr_id !== nothing
            tr = _PMD.ref(pm, nw, :transformer, tr_id)

            tap_dict = get(_PMD.var(pm, nw), :tap, nothing)

            tap1 = (tap_dict !== nothing && haskey(tap_dict, tr_id)) ? tap_dict[tr_id][1] : tr["tm_set"][1]

            # Als tap1 een JuMP variabele/expression is -> NL tap-aware constraint
            if tap1 isa JuMP.AbstractJuMPScalar
                constraint_mc_bus_voltage_drop_tap(pm, nw, i,
                    branch["f_bus"], branch["t_bus"], (i, branch["f_bus"], branch["t_bus"]),
                    branch["f_connections"], branch["t_connections"],
                    r0, x0, tap1
                )
                return nothing
            end
        end
    end
    _PMD.constraint_mc_bus_voltage_drop(pm, nw, i,
        branch["f_bus"], branch["t_bus"], (i, branch["f_bus"], branch["t_bus"]),
        branch["f_connections"], branch["t_connections"],
        r_use, x_use
    )

    return nothing
end

function constraint_mc_bus_voltage_drop_tap(pm::_PMD.AbstractExplicitNeutralIVRModel,
    nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int},
    f_connections::Vector{Int}, t_connections::Vector{Int},
    r0::Matrix{<:Real}, x0::Matrix{<:Real}, tap1
)
    vr_fr = [_PMD.var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    vi_fr = [_PMD.var(pm, nw, :vi, f_bus)[c] for c in f_connections]
    vr_to = [_PMD.var(pm, nw, :vr, t_bus)[c] for c in t_connections]
    vi_to = [_PMD.var(pm, nw, :vi, t_bus)[c] for c in t_connections]

    csr_fr = _PMD.var(pm, nw, :csr, f_idx[1])  
    csi_fr = _PMD.var(pm, nw, :csi, f_idx[1])  

    n = length(f_connections)
    
    for p in 1:n
        r_eff = r0[p,p] / tap1

        JuMP.@constraint(pm.model,
            vr_to[p] == vr_fr[p] - r_eff*csr_fr[p] + x0[p,p]*csi_fr[p]
        )
        println("tap1: ", tap1, " r_eff: ", r_eff)
        JuMP.@constraint(pm.model,
            vi_to[p] == vi_fr[p] - r_eff*csi_fr[p] - x0[p,p]*csr_fr[p]
        )
    end
end