export solve_acp_mc_se_oltc, build_mc_se_oltc

"solves the AC state estimation with OLTC tap estimation"
function solve_acp_mc_se_oltc(data::Union{Dict{String,<:Any},String}, solver; kwargs...)
    return solve_mc_se_oltc(data, _PMD.ACPUPowerModel, solver; kwargs...)
end

function solve_mc_se_oltc(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    if haskey(data["se_settings"], "criterion")
        _PMDSE.assign_unique_individual_criterion!(data)
    end
    if !haskey(data["se_settings"], "rescaler")
        data["se_settings"]["rescaler"] = 1
    end
    if !haskey(data["se_settings"], "number_of_gaussian")
        data["se_settings"]["number_of_gaussian"] = 10
    end

    # --- KEY FIX: force multinetwork solving when data has "nw" ---
    is_mn = haskey(data, "nw") || get(data, "multinetwork", false) == true

    if is_mn
        return _PMD.solve_mc_model(data, model_type, solver, build_mc_se_oltc; multinetwork=true, kwargs...)
    else
        return _PMD.solve_mc_model(data, model_type, solver, build_mc_se_oltc; kwargs...)
    end
end

"Specification of the SE problem including Transformer Taps as variables"
function build_mc_se_oltc(pm::_PMD.IVRENPowerModel)
    
    # Time varying variables
    for (n, _) in _PMD.nws(pm)
        # Variables
        variable_mc_bus_voltage(pm, nw=n, bounded = true)
        _PMD.variable_mc_branch_current(pm, nw=n, bounded = true)
        _PMD.variable_mc_generator_current(pm, nw=n, bounded = true)
        _PMD.variable_mc_transformer_current(pm, nw=n, bounded = true)
        variable_mc_load_current(pm, nw=n, bounded = true)    
        variable_mc_residual(pm, nw=n, bounded = true)
        variable_mc_measurement(pm, nw=n, bounded = false)
        variable_mc_transformer_tap(pm, nw=n, bounded = true)    # --- ADDED: Tap Estimation Variable ---
        #_PMD.variable_mc_switch_current(pm, nw=n, bounded = true)  # --- ADDED: Switch Current Variable ---
    end
    
    

    # Constraints
    constraint_mc_transformer_tap_time_invariant(pm)   # --- ADDED: Tap Time-Invariant Constraint ---
    # Time varying constraints
    for (n, _) in _PMD.nws(pm)
        for i in _PMD.ids(pm, n, :bus)
            if i in _PMD.ids(pm, n, :ref_buses)
            _PMD.constraint_mc_voltage_reference(pm, i, nw = n)  # vm is not fixed
            end
        end
        
        for id in _PMD.ids(pm, n, :gen)
            constraint_mc_generator_current_se(pm, id, nw = n)
        end
#
        for i in _PMD.ids(pm, n, :transformer)
            constraint_mc_transformer_tap_equal_phase(pm, i, nw = n)
            #constraint_mc_transformer_tap_test(pm, i)   # --- ADDED: Tap Constraint ---
        end
        for i in _PMD.ids(pm, n, :transformer)
            constraint_mc_transformer_voltage(pm, i,fix_taps=false, nw = n)  # --- MODIFIED: fix_taps=false ---
            constraint_mc_transformer_current(pm, i,fix_taps=false, nw = n)# --- MODIFIED: fix_taps=false 
        end


        for i in _PMD.ids(pm, n, :branch)
            _PMD.constraint_mc_current_from(pm, i, nw = n)
            _PMD.constraint_mc_current_to(pm, i, nw = n)
            _PMD.constraint_mc_bus_voltage_drop(pm, i, nw = n)
        end
#
        for (i,bus) in _PMD.ref(pm, n, :bus)
            constraint_mc_current_balance_se(pm, i, nw = n)
            #constraint_mc_neutral_grounding(pm, i)  #TODO: make it only grounded if load is grounded
        end

        for (i,meas) in _PMD.ref(pm, n, :meas)
            constraint_mc_residual(pm, i, nw = n)
        end
    end

        
    # Variables
    #variable_mc_bus_voltage(pm, bounded = true)
    #_PMD.variable_mc_branch_current(pm, bounded = true)
    #_PMD.variable_mc_generator_current(pm, bounded = true)
    #_PMD.variable_mc_transformer_current(pm, bounded = true)
    #variable_mc_load_current(pm, bounded = true)    
    #variable_mc_residual(pm, bounded = true)
    #variable_mc_measurement(pm, bounded = false)
    #variable_mc_transformer_tap(pm,  bounded = true)    # --- ADDED: Tap Estimation Variable ---
    # Constraints
    
    
    #for i in _PMD.ids(pm, :bus)
    #    if i in _PMD.ids(pm, :ref_buses)
    #    _PMD.constraint_mc_voltage_reference(pm, i)  # vm is not fixed
    #    end
    #end
        
    #for id in _PMD.ids(pm, :gen)
     #       constraint_mc_generator_current_se(pm, id)
    #end

    #for i in _PMD.ids(pm, :transformer)
    #    constraint_mc_transformer_tap_equal_phase(pm, i)
    #    #constraint_mc_transformer_tap_test(pm, i)   # --- ADDED: Tap Constraint ---
    #end
    #for i in _PMD.ids(pm, :transformer)
    #    constraint_mc_transformer_voltage(pm, i,fix_taps=false)  # --- MODIFIED: fix_taps=false ---
    #    constraint_mc_transformer_current(pm, i,fix_taps=false) # --- MODIFIED: fix_taps=false 
    #end
    #for i in _PMD.ids(pm, :branch)
    #    _PMD.constraint_mc_current_from(pm, i)
     #   _PMD.constraint_mc_current_to(pm, i)
     #   _PMD.constraint_mc_bus_voltage_drop(pm, i)
    #end
    #for (i,bus) in _PMD.ref(pm, :bus)
    #    constraint_mc_current_balance_se(pm, i)
    #    #constraint_mc_neutral_grounding(pm, i)  #TODO: make it only grounded if load is grounded
    #end
    #for (i,meas) in _PMD.ref(pm, :meas)
    #    constraint_mc_residual(pm, i)
    #end

    # Objective
    objective_mc_se(pm)
end