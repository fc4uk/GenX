"""
GenX: An Configurable Capacity Expansion Model
Copyright (C) 2021,  Massachusetts Institute of Technology
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
A complete copy of the GNU General Public License v2 (GPLv2) is available
in LICENSE.txt.  Users uncompressing this from an archive may not have
received this license file.  If not, see <http://www.gnu.org/licenses/>.
"""

@doc raw"""
    function configure_ddp_dicts(setup::Dict, inputs::Dict)

This function instantiates Dictionary objects containing the names of linking expressions, constraints, and variables used in multi-stage modeling.

inputs:

* setup - Dictionary object containing GenX settings and key parameters.
* inputs – Dictionary of inputs for each model period, generated by the load\_inputs() method.

returns:

* start\_cap\_d – Dictionary which contains linking expression names as keys and linking constraint names as values, used for setting the end capacity in stage $p$ to the starting capacity in stage $p+1$.
* cap\_track\_d – Dictionary which contains linking variable names as keys and linking constraint names as values, used for enforcing endogenous retirements.
"""
function configure_ddp_dicts(setup::Dict, inputs::Dict)

    # start_cap_d dictionary contains key-value pairs of available capacity investment expressions
    # as keys and their corresponding linking constraints as values
    start_cap_d = Dict([(Symbol("eTotalCap"), Symbol("cExistingCap"))])

    if !isempty(inputs["STOR_ALL"])
        start_cap_d[Symbol("eTotalCapEnergy")] = Symbol("cExistingCapEnergy")
    end

    if !isempty(inputs["STOR_ASYMMETRIC"])
        start_cap_d[Symbol("eTotalCapCharge")] = Symbol("cExistingCapCharge")
    end

    if setup["NetworkExpansion"] == 1 && inputs["Z"] > 1
        start_cap_d[Symbol("eAvail_Trans_Cap")] = Symbol("cExistingTransCap")
    end

    # This dictionary contains the endogenous retirement constraint name as a key,
    # and a tuple consisting of the associated tracking array constraint and variable as the value
    cap_track_d = Dict([(Symbol("vCAPTRACK"), Symbol("cCapTrack"))])

    if !isempty(inputs["STOR_ALL"])
        cap_track_d[Symbol("vCAPTRACKENERGY")] = Symbol("cCapTrackEnergy")
    end

    if !isempty(inputs["STOR_ASYMMETRIC"])
        cap_track_d[Symbol("vCAPTRACKCHARGE")] = Symbol("cCapTrackCharge")
    end

    return start_cap_d, cap_track_d
end

@doc raw"""
	function run_ddp(models_d::Dict, setup::Dict, inputs_d::Dict)

This function run the dual dynamic programming (DDP) algorithm, as described in [Pereira and Pinto (1991)](https://doi.org/10.1007/BF01582895), and more recently, [Lara et al. (2018)](https://doi.org/10.1016/j.ejor.2018.05.039). Note that if the algorithm does not converge within 10,000 (currently hardcoded) iterations, this function will return models with sub-optimal solutions. However, results will still be printed as if the model is finished solving.

inputs:

  * models\_d – Dictionary which contains a JuMP model for each model period.
  * setup - Dictionary object containing GenX settings and key parameters.
  * inputs\_d – Dictionary of inputs for each model stage, generated by the load\_inputs() method.

returns:

  * models\_d – Dictionary which contains a JuMP model for each model stage, modified by this method.
  * stats\_d – Dictionary which contains the run time, upper bound, and lower bound of each DDP iteration.
  * inputs\_d – Dictionary of inputs for each model stage, generated by the load\_inputs() method, modified by this method.
"""
function run_ddp(models_d::Dict, setup::Dict, inputs_d::Dict)

    settings_d = setup["MultiStageSettingsDict"]
    num_stages = settings_d["NumStages"]  # Total number of investment planning stages
    EPSILON = settings_d["ConvergenceTolerance"] # Tolerance
    myopic = settings_d["Myopic"] == 1 # 1 if myopic (only one forward pass), 0 if full DDP

    start_cap_d, cap_track_d = configure_ddp_dicts(setup, inputs_d[1])

    ic = 0 # Iteration Counter

    results_d = Dict() # Dictionary to store the results to return
    stats_d = Dict() # Dictionary to store the statistics (total time, upper bound, and lower bound for each iteration)
    times_a = [] # Array to store the total time of each iteration
    upper_bounds_a = [] # Array to store the upper bound of each iteration
    lower_bounds_a = [] # Array to store the lower bound of each iteration

    # Step a.i) Initialize cost-to-go function for t = 1:num_stages
    for t in 1:num_stages
        models_d[t] = initialize_cost_to_go(settings_d, models_d[t], inputs_d[t])
    end

    # Step a.ii) Set objective upper bound
    global z_upper = Inf

    # Step b.i) Solve the approximate first-stage problem
    println("***********")
    println("Solving First Stage Problem")
    println("***********")


    t = 1 # Stage = 1
    solve_time_d = Dict()
    ddp_prev_time = time() # Begin tracking time of each iteration
    models_d[t], solve_time_d[t] = solve_model(models_d[t], setup)
    inputs_d[t]["solve_time"] = solve_time_d[t]

    # Step c.i) Initialize the lower bound, equal to the objective function value for the first period in the first iteration
    global z_lower = objective_value(models_d[t])

    # Step c.ii) If the relative difference between upper and lower bounds are small, break loop
    while ((z_upper - z_lower) / z_lower > EPSILON)

        ic = ic + 1 # Increase iteration counter by 1

        if (ic > 10000)
            println("***********")
            println("Exiting Without Covergence!")
            println(string("Upper Bound = ", z_upper))
            println(string("Lower Bound = ", z_lower))
            println("***********")

            stats_d["TIMES"] = times_a
            stats_d["UPPER_BOUNDS"] = upper_bounds_a
            statd_d["LOWER_BOUNDS"] = lower_bounds_a

            return models_d, stats_d, inputs_d
        end

        println("***********")
        println(string("Iteration Number: ", ic))
        println(string("Upper Bound = ", z_upper))
        println(string("Lower Bound = ", z_lower))
        println("***********")

        # Step d) Forward pass for t = 1:num_stages
        ## For first iteration we dont need to solve forward pass for first stage (we did that already above),
        ## but we need to update forward pass solution for the first stage for subsequent iterations
        if ic > 1
            t = 1 #  update forward pass solution for the first stage
            models_d[t], solve_time_d[t] = solve_model(models_d[t], setup)
            inputs_d[t]["solve_time"] = solve_time_d[t]
        end
        ## Forward pass for t=2:num_stages
        for t in 2:num_stages

            println("***********")
            println(string("Forward Pass t = ", t))
            println("***********")

            # Step d.i) Fix initial investments for model at time t given optimal solution for time t-1
            models_d[t] = fix_initial_investments(models_d[t-1], models_d[t], start_cap_d)

            # Step d.ii) Fix capacity tracking variables for endogenous retirements
            models_d[t] = fix_capacity_tracking(models_d[t-1], models_d[t], cap_track_d, t)

            # Step d.iii) Solve the model at time t
            models_d[t], solve_time_d[t] = solve_model(models_d[t], setup)
            inputs_d[t]["solve_time"] = solve_time_d[t]

        end

        ### For the myopic solution, algorithm should terminate here after the first forward pass calculation and then move to Outputs writing.
        if myopic
            println("***********")
            println("Exiting After First Forward Pass! (Myopic)")
            println(string("Upper Bound = ", z_upper))
            println(string("Lower Bound = ", z_lower))
            println("***********")

            stats_d["TIMES"] = times_a
            stats_d["UPPER_BOUNDS"] = upper_bounds_a
            stats_d["LOWER_BOUNDS"] = lower_bounds_a
            return models_d, stats_d, inputs_d
        end
        ###

        # Step e) Calculate the new upper bound
        z_upper_temp = 0
        for t in 1:num_stages
            z_upper_temp = z_upper_temp + (objective_value(models_d[t]) - value(models_d[t][:vALPHA]))
        end

        # If the upper bound decreased, set it as the new upper bound
        if z_upper_temp < z_upper
            z_upper = z_upper_temp
        end

        append!(upper_bounds_a, z_upper) # Store current iteration upper bound

        # Step f) Backward pass for t = num_stages:2
        for t in num_stages:-1:2

            println("***********")
            println(string("Backward Pass t = ", t))
            println("***********")

            # Step f.i) Add a cut to the previous time step using information from the current time step
            models_d[t-1] = add_cut(models_d[t-1], models_d[t], start_cap_d, cap_track_d)

            # Step f.ii) Solve the model with the additional cut at time t-1
            models_d[t-1], solve_time_d[t-1] = solve_model(models_d[t-1], setup)
            inputs_d[t-1]["solve_time"] = solve_time_d[t-1]
        end

        # Step g) Recalculate lower bound and go back to c)
        z_lower = objective_value(models_d[1])
        append!(lower_bounds_a, z_lower) # Store current iteration lower bound

        # Step h) Store the total time of the current iteration (in seconds)
        ddp_iteration_time = time() - ddp_prev_time
        append!(times_a, ddp_iteration_time)
        ddp_prev_time = time()
    end

    println("***********")
    println("Successful Convergence!")
    println(string("Upper Bound = ", z_upper))
    println(string("Lower Bound = ", z_lower))
    println("***********")


    ### STEP I) One final forward pass to guarantee convergence
    # Forward pass for t = 1:num_stages
    t = 1 #  update forward pass solution for the first stage
    models_d[t], solve_time_d[t] = solve_model(models_d[t], setup)
    inputs_d[t]["solve_time"] = solve_time_d[t]
    ## Forward pass for t=2:num_stages
    for t in 2:num_stages
        println("***********")
        println(string("Final Forward Pass t = ", t))
        println("***********")

        # Step d.i) Fix initial investments for model at time t given optimal solution for time t-1
        models_d[t] = fix_initial_investments(models_d[t-1], models_d[t], start_cap_d)

        # Step d.ii) Fix capacity tracking variables for endogenous retirements
        models_d[t] = fix_capacity_tracking(models_d[t-1], models_d[t], cap_track_d, t)

        # Step d.iii) Solve the model at time t
        models_d[t], solve_time_d[t] = solve_model(models_d[t], setup)
        inputs_d[t]["solve_time"] = solve_time_d[t]
    end
    ##### END of final forward pass

    stats_d["TIMES"] = times_a
    stats_d["UPPER_BOUNDS"] = upper_bounds_a
    stats_d["LOWER_BOUNDS"] = lower_bounds_a

    return models_d, stats_d, inputs_d
end

@doc raw"""
	function write_multi_stage_outputs(stats_d::Dict, outpath::String, settings_d::Dict)

This function calls various methods which write multi-stage modeling outputs as .csv files.

inputs:

  * stats\_d – Dictionary which contains the run time, upper bound, and lower bound of each DDP iteration.
  * outpath – String which represents the path to the Results directory.
  * settings\_d - Dictionary containing settings configured in the GenX settings genx\_settings.yml file as well as the multi-stage settings file multi\_stage\_settings.yml.
"""
function write_multi_stage_outputs(stats_d::Dict, outpath::String, settings_d::Dict, inputs_dict::Dict)

    multi_stage_settings_d = settings_d["MultiStageSettingsDict"]

    write_capacities_discharge(outpath, multi_stage_settings_d)
    #write_capacities_charge(outpath, multi_stage_settings_d)
    #write_capacities_energy(outpath, multi_stage_settings_d)
    #write_network_expansion(outpath, multi_stage_settings_d)
    write_costs(outpath, multi_stage_settings_d, inputs_dict)
    write_stats(outpath, stats_d)
    write_settings(outpath, settings_d)

end

@doc raw"""
	function write_capacities_discharge(outpath::String, settings_d::Dict)

This function writes the file capacities\_multi\_stage.csv to the Results directory. This file contains starting resource capcities from the first model stage and end resource capacities for the first and all subsequent model stages.

inputs:

  * outpath – String which represents the path to the Results directory.
  * settings\_d - Dictionary containing settings dictionary configured in the multi-stage settings file multi\_stage\_settings.yml.
"""
function write_capacities_discharge(outpath::String, settings_d::Dict)
    # TO DO - DO THIS FOR ENERGY CAPACITY AS WELL

    num_stages = settings_d["NumStages"] # Total number of investment planning stages
    capacities_d = Dict()

    for p in 1:num_stages
        inpath = joinpath(outpath, "Results_p$p")
        capacities_d[p] = DataFrame(CSV.File(joinpath(inpath, "capacity.csv"), header=true), copycols=true)
    end

    # Set first column of DataFrame as resource names from the first stage
    df_cap = DataFrame(Resource=capacities_d[1][!, :Resource], Zone=capacities_d[1][!, :Zone])

    # Store starting capacities from the first stage
    df_cap[!, Symbol("StartCap_p1")] = capacities_d[1][!, :StartCap]

    # Store end capacities for all stages
    for p in 1:num_stages
        df_cap[!, Symbol("EndCap_p$p")] = capacities_d[p][!, :EndCap]
    end

    CSV.write(joinpath(outpath, "capacities_multi_stage.csv"), df_cap)

end

function write_capacities_charge(outpath::String, settings_d::Dict) end

function write_capacities_energy(outpath::String, settings_d::Dict) end

function write_network_expansion(outpath::String, settings_d::Dict)
    # Include discounted NE costs and capacities for each model period
end

@doc raw"""
	function write_costs(outpath::String, settings_d::Dict)

This function writes the file costs\_multi\_stage.csv to the Results directory. This file contains variable, fixed, startup, network expansion, unmet reserve, and non-served energy costs discounted to year zero.

inputs:

  * outpath – String which represents the path to the Results directory.
  * settings\_d - Dictionary containing settings dictionary configured in the multi-stage settings file multi\_stage\_settings.yml.
"""
function write_costs(outpath::String, settings_d::Dict, inputs_dict::Dict)

    num_stages = settings_d["NumStages"] # Total number of DDP stages
    wacc = settings_d["WACC"] # Interest Rate and also the discount rate unless specified other wise
    stage_lens = settings_d["StageLengths"]
    myopic = settings_d["Myopic"] == 1 # 1 if myopic (only one forward pass), 0 if full DDP

    costs_d = Dict()
    for p in 1:num_stages
        cur_path = joinpath(outpath, "Results_p$p")
        costs_d[p] = DataFrame(CSV.File(joinpath(cur_path, "costs.csv"), header=true), copycols=true)
    end

    OPEXMULTS = [inputs_dict[j]["OPEXMULT"] for j in 1:num_stages] # Stage-wise OPEX multipliers to count multiple years between two model stages

    # Set first column of DataFrame as resource names from the first stage
    df_costs = DataFrame(Costs=costs_d[1][!, :Costs])

    # Store discounted total costs for each stage in a data frame
    for p in 1:num_stages
        if myopic
            DF = 1 # DF=1 because we do not apply discount factor in myopic case
        else
            DF = 1 / (1 + wacc)^(stage_lens[p] * (p - 1))  # Discount factor applied to ALL costs in each stage
        end
        df_costs[!, Symbol("TotalCosts_p$p")] = DF .* costs_d[p][!, Symbol("Total")]
    end

    # For OPEX costs, apply additional discounting
    for cost in ["cFuel", "cVOM", "cNSE", "cStart", "cUnmetRsv"]
        if cost in df_costs[!, :Costs]
            df_costs[df_costs[!, :Costs].==cost, 2:end] = transpose(OPEXMULTS) .* df_costs[df_costs[!, :Costs].==cost, 2:end]
        end
    end

    # Remove "cTotal" from results (as this includes Cost-to-Go)
    df_costs = df_costs[df_costs[!, :Costs].!="cTotal", :]

    CSV.write(joinpath(outpath, "costs_multi_stage.csv"), df_costs)

end

@doc raw"""
	function write_stats(outpath::String, stats_d::Dict)

This function writes the file stats\_multi\_stage.csv. to the Results directory. This file contains the runtime, upper bound, lower bound, and relative optimality gap for each iteration of the DDP algorithm.

inputs:

  * outpath – String which represents the path to the Results directory.
  * stats\_d – Dictionary which contains the run time, upper bound, and lower bound of each DDP iteration.
"""
function write_stats(outpath::String, stats_d::Dict)

    times_a = stats_d["TIMES"] # Time (seconds) of each iteration
    upper_bounds_a = stats_d["UPPER_BOUNDS"] # Upper bound of each iteration
    lower_bounds_a = stats_d["LOWER_BOUNDS"] # Lower bound of each iteration

    # Create an array of numbers 1 through total number of iterations
    iteration_count_a = collect(1:length(times_a))

    realtive_gap_a = (upper_bounds_a .- lower_bounds_a) ./ lower_bounds_a

    # Construct dataframe where first column is iteration number, second is iteration time
    df_stats = DataFrame(Iteration_Number=iteration_count_a,
        Seconds=times_a,
        Upper_Bound=upper_bounds_a,
        Lower_Bound=lower_bounds_a,
        Relative_Gap=realtive_gap_a)

    CSV.write(joinpath(outpath, "stats_multi_stage.csv"), df_stats)

end

@doc raw"""
	function fix_initial_investments(EP_prev::Model, EP_cur::Model, start_cap_d::Dict)

This function sets the right hand side values of the existing capacity linking constraints in the current stage $p$ to the realized values of the total available end capacity linking variable expressions from the previous stage $p-1$ as part of the forward pass.

inputs:

  * EP\_prev - JuMP model from the previous model stage $p-1$.
  * EP\_cur - JuMP model from the current model stage $p$.
  * start\_cap\_d – Dictionary which contains key-value pairs of available capacity investment expression names (as Symbols) as keys and their corresponding linking constraint names (as Symbols) as values.

returns: JuMP model with updated linking constraints.
"""
function fix_initial_investments(EP_prev::Model, EP_cur::Model, start_cap_d::Dict)

    # start_cap_d dictionary contains the starting capacity expression name (e) as a key,
    # and the associated linking constraint name (c) as a value
    for (e, c) in start_cap_d
        for y in keys(EP_cur[c])

            # Set the right hand side value of the linking initial capacity constraint in the current
            # stage to the value of the available capacity variable solved for in the previous stages
            set_normalized_rhs(EP_cur[c][y], value(EP_prev[e][y]))
        end
    end
    return EP_cur
end

@doc raw"""
	function fix_capacity_tracking(EP_prev::Model, EP_cur::Model, cap_track_d::Dict, cur_stage::Int)

This function sets the right hand side values of the new and retired capacity tracking linking constraints in the current stage $p$ to the realized values of the new and retired capacity tracking linking variables from the previous stage $p-1$ as part of the forward pass.
where tracking linking variables are defined variables for tracking, linking and passing realized expansion and retirement of capacities of each stage to the next stage.
Tracking linking variables are each defined in endogenous\_retirement\_discharge, endogenous\_retirement\_energy, and endogenous\_retirement\_charge functions. Three examples are "vCAPTRACK", "vCAPTRACKCHARGE", and ""vCAPTRACKENERGY"

inputs:

  * EP\_prev - JuMP model from the previous model stage $p-1$.
  * EP\_cur - JuMP model from the current model stage $p$.
  * cap\_track\_d – Dictionary which contains key-value pairs of capacity addition and retirement tracking variable names (as Symbols) as keys and their corresponding linking constraint names (as Symbols) as values.
  * cur\_period – Int representing the current model stage $p$.

returns: JuMP model with updated linking constraints.
"""
function fix_capacity_tracking(EP_prev::Model, EP_cur::Model, cap_track_d::Dict, cur_stage::Int)

    # cap_track_d dictionary contains the endogenous retirement tracking array variable name (v) as a key,
    # and the associated linking constraint name (c) as a value
    for (v, c) in cap_track_d

        # Tracking variables and constraints for retired capacity are named identicaly to those for newly
        # built capacity, except have the prefex "vRET" and "cRet", accordingly
        rv = Symbol("vRET", string(v)[2:end]) # Retired capacity tracking variable name (rv)
        rc = Symbol("cRet", string(c)[2:end]) # Retired capacity tracking constraint name (rc)

        for i in keys(EP_cur[c])
            i = i[1] # Extract integer index value from keys tuple - corresponding to generator index

            # For all previous stages, set the right hand side value of the tracking constraint in the current
            # stage to the value of the tracking constraint observed in the previous stage
            for p in 1:(cur_stage-1)
                # Tracking newly buily capacity over all previous stages
                JuMP.set_normalized_rhs(EP_cur[c][i, p], value(EP_prev[v][i, p]))
                # Tracking retired capacity over all previous stages
                JuMP.set_normalized_rhs(EP_cur[rc][i, p], value(EP_prev[rv][i, p]))
            end
        end
    end

    return EP_cur
end

@doc raw"""
	function add_cut(EP_cur::Model, EP_next::Model, start_cap_d::Dict, cap_track_d::Dict)

inputs:

  * EP\_cur - JuMP model from the current model stage $p$.
  * EP\_next - JuMP model from the next model stage $p+1$..
  * start\_cap\_d – Dictionary which contains key-value pairs of available capacity investment expression names (as Symbols) as keys and their corresponding linking constraint names (as Symbols) as values.
  * cap\_track\_d – Dictionary which contains key-value pairs of capacity addition and retirement tracking variable names (as Symbols) as keys and their corresponding linking constraint names (as Symbols) as values.

returns: JuMP expression representing a sum of Benders cuts for linking capacity investment variables to be added to the cost-to-go function.
"""
function add_cut(EP_cur::Model, EP_next::Model, start_cap_d::Dict, cap_track_d::Dict)

    next_obj_value = objective_value(EP_next) # Get the objective function value for the next investment planning stage

    eRHS = @expression(EP_cur, 0) # Initialize RHS of cut to 0

    # Generate cut components for investment decisions

    # start_cap_d dictionary contains the starting capacity expression name (e) as a key,
    # and the associated linking constraint name (c) as a value
    for (e, c) in start_cap_d

        # Continue if nothing to add to the cut
        if isempty(EP_next[e])
            continue
        end

        # Generate the cut component
        eCurRHS = generate_cut_component_inv(EP_cur, EP_next, e, c)

        # Add the cut component to the RHS
        eRHS = eRHS + eCurRHS
    end

    # Generate cut components for endogenous retirements.

    # cap_track_d dictionary contains the endogenous retirement tracking array variable name (v) as a key,
    # and the associated linking constraint name (c) as a value
    for (v, c) in cap_track_d

        # Continue if nothing to add to the cut
        if isempty(EP_next[c])
            continue
        end

        # Generate the cut component for new capacity
        eCurRHS_cap = generate_cut_component_track(EP_cur, EP_next, v, c)

        rv = Symbol("vRET", string(v)[2:end]) # Retired capacity tracking variable (rv)
        rc = Symbol("cRet", string(c)[2:end]) # Retired capacity tracking constraint (rc)

        # Generate the cut component for retired capacity
        eCurRHS_ret = generate_cut_component_track(EP_cur, EP_next, rv, rc)

        # Add the cut component to the RHS
        eRHS = eRHS + eCurRHS_cap + eCurRHS_ret
    end

    # Add the cut to the model
    @constraint(EP_cur, EP_cur[:vALPHA] >= next_obj_value - eRHS)

    return EP_cur
end

@doc raw"""
	function generate_cut_component_inv(EP_cur::Model, EP_next::Model, expr_name::Symbol, constr_name::Symbol)

This function generates Bender's cut expressions for total new or retired capacity tracking linking variables in the form:
```math
\begin{aligned}
        \mu_{next}^{\top}(\hat{x}_{cur} - x_{cur})
\end{aligned}
```
where $\mu_{next}$ is a vector of dual values of the linking constraints defined by constr\_name in EP\_next, $\hat{x}_{cur}$ is a vector of realized values from the forward pass of the new or retired capacity tracking linking variables var\_name from EP\_cur, and $x_{cur}$ is a vector of unrealized new or retired capacity tracking linking variables from EP\_cur.

inputs:

  * EP\_cur - JuMP model from the current model stage $p$.
  * EP\_next - JuMP model from the next model stage $p+1$.
  * var\_name – Symbol representing the name of a JuMP variable array which contains total new or retired capacity tracking linking variables.
  * constr\_name – Symbol representing the name of the array of linking JuMP constraints which contain total new or retired capacity tracking linking variables.

returns: JuMP expression representing a sum of Benders cuts for linking capacity investment variables to be added to the cost-to-go function.
"""
function generate_cut_component_track(EP_cur::Model, EP_next::Model, var_name::Symbol, constr_name::Symbol)

    next_dual_value = Float64[]
    cur_inv_value = Float64[]
    cur_inv_var = []

    for k in keys(EP_next[constr_name])
        y = k[1] # Index representing resource
        p = k[2] # Index representing stage

        push!(next_dual_value, getdual(EP_next[constr_name][y, p]))
        push!(cur_inv_value, getvalue(EP_cur[var_name][y, p]))
        push!(cur_inv_var, EP_cur[var_name][y, p])
    end

    eCutComponent = @expression(EP_cur, dot(next_dual_value, (cur_inv_value .- cur_inv_var)))

    return eCutComponent
end

@doc raw"""
	function generate_cut_component_inv(EP_cur::Model, EP_next::Model, expr_name::Symbol, constr_name::Symbol)

This function generates Bender's cut expressions for linking capacity investment variable expression in the form:
```math
\begin{aligned}
        \mu_{next}^{\top}(\hat{x}_{cur} - x_{cur})
\end{aligned}
```
where $\mu_{next}$ is a vector of dual values of the linking constraints defined by constr\_name in EP\_next, $\hat{x}_{cur}$ is a vector of realized values from the forward pass of the linking capacity investment variable expressions expr\_name from EP\_cur, and $x_{cur}$ is a vector of unrealized linking capacity investment variable expressions from EP\_cur. inputs:

inputs:

  * EP\_cur - JuMP model from the current model stage $p$, solved in the forward pass.
  * EP\_next - JuMP model from the next model stage $p+1$, solved in the forward pass.
  * expr\_name – Symbol representing the name of a JuMP expression array which contains linking capacity investment variables.
  * constr\_name – Symbol representing the name of the array of linking JuMP constraints which contain the linking capacity investment variables.

returns: JuMP expression representing a sum of Benders cuts for linking capacity investment variables to be added to the cost-to-go function.
"""
function generate_cut_component_inv(EP_cur::Model, EP_next::Model, expr_name::Symbol, constr_name::Symbol)

    next_dual_value = Float64[]
    cur_inv_value = Float64[]
    cur_inv_var = []

    for y in keys(EP_next[constr_name])

        push!(next_dual_value, dual(EP_next[constr_name][y]))
        push!(cur_inv_value, value(EP_cur[expr_name][y]))
        push!(cur_inv_var, EP_cur[expr_name][y])
    end

    eCutComponent = @expression(EP_cur, dot(next_dual_value, (cur_inv_value .- cur_inv_var)))

    return eCutComponent
end

@doc raw"""
	function initialize_cost_to_go(settings_d::Dict, EP::Model)

This function scales the model objective function so that costs are consistent with multi-stage modeling and introduces a cost-to-go function variable to the objective function.

The updated objective function $OBJ^{*}$ returned by this method takes the form:
```math
\begin{aligned}
    OBJ^{*} = DF * OPEXMULT * OBJ + \alpha
\end{aligned}
```
where $OBJ$ is the original objective function. $OBJ$ is scaled by two terms. The first is a discount factor (applied only in the non-myopic case), which discounts costs associated with the model stage $p$ to year-0 dollars:
```math
\begin{aligned}
    DF = \frac{1}{(1+WACC)^{L*(p-1)}}
\end{aligned}
```
where $WACC$ is the weighted average cost of capital, and $L$ is the length of each stage in years (both set in multi\_stage\_settings.yml)

The second term is a discounted sum of annual operational expenses incurred each year of a multi-year model stage:
```math
\begin{aligned}
    & OPEXMULT = \sum^{L}_{l=1}\frac{1}{(1+WACC)^{l-1}}
\end{aligned}
```
Note that although the objective function contains investment costs, which occur only once and thus do not need to be scaled by OPEXMULT, these costs are multiplied by a factor of $\frac{1}{WACC}$ before being added to the objective function in investment\_discharge\_multi\_stage(), investment\_charge\_multi\_stage(), investment\_energy\_multi\_stage(), and transmission\_multi\_stage(). Thus, this step scales these costs back to their correct value.

The cost-to-go function $\alpha$ represents an approximation of future costs given the investment and retirement decisions in the current stage. It is constructed through the addition of cuts to the cost-to-go function $\alpha$ during the backwards pass.

inputs:

  * settings\_d - Dictionary containing settings dictionary configured in the multi-stage settings file multi\_stage\_settings.yml.
  * EP – JuMP model.

returns: JuMP model with updated objective function.
"""
function initialize_cost_to_go(settings_d::Dict, EP::Model, inputs::Dict)

    cur_stage = settings_d["CurStage"] # Current DDP Investment Planning Stage
    stage_len = settings_d["StageLengths"][cur_stage]
    wacc = settings_d["WACC"] # Interest Rate  and also the discount rate unless specified other wise
    myopic = settings_d["Myopic"] == 1 # 1 if myopic (only one forward pass), 0 if full DDP
    OPEXMULT = inputs["OPEXMULT"] # OPEX multiplier to count multiple years between two model stages, set in configure_multi_stage_inputs.jl

    # Overwrite the objective function to include the cost-to-go variable (not in myopic case)
    # Multiply discount factor to all terms except the alpha term or the cost-to-go function
    # All OPEX terms get an additional adjustment factor
    if myopic
        ### No discount factor or OPEX multiplier applied in myopic case as costs are left annualized.
        @objective(EP, Min, EP[:eObj])
    else
        DF = 1 / (1 + wacc)^(stage_len * (cur_stage - 1))  # Discount factor applied all to costs in each stage ###
        # Initialize the cost-to-go variable
        @variable(EP, vALPHA >= 0)
        @objective(EP, Min, DF * OPEXMULT * EP[:eObj] + vALPHA)
    end

    return EP

end
