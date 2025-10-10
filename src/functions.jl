using CommonSolve

"""
    init(prob::MermaidProblem, alg::MermaidSolver; kwargs...)

Creates a [MermaidIntegrator](@ref) from a [MermaidProblem](@ref).

# Arguments
- `prob::MermaidProblem`: The hybrid problem to be solved.
- `alg::AbstractMermaidSolver`: The [AbstractMermaidSolver](@ref) algorithm to be used for solving the problem.
- `kwargs...`: Additional keyword arguments for the solver.

# Returns
- `MermaidIntegrator`: The initialized integrator for the problem.
"""
function CommonSolve.init(prob::MermaidProblem, alg::AbstractMermaidSolver; save_vars=[])
    # Sort connectors so that the ones with # in outputs are first
    connectors = sort(prob.connectors; by=x -> any([contains(conn.variable, "#") for conn in x.outputs]), rev=true)
    # Initialize the solver
    integrators = Vector{Any}()
    for c in prob.components
        integrator = CommonSolve.init(c, connectors)
        push!(integrators, integrator)
    end
    return MermaidIntegrator(integrators, connectors, prob.max_t, 0.0, alg, save_vars)
end

"""
    step!(int::MermaidIntegrator, dt::Float64)

Steps the integrator for the given time step. The method is defined by the `alg` field of the [MermaidIntegrator](@ref).

# Arguments
- `int::MermaidIntegrator`: The integrator to be stepped.
- `dt::Float64`: The time step for the integrator.
"""
function CommonSolve.step!(merInt::MermaidIntegrator, dt)
    merInt.alg(merInt, dt)
end

"""
    solve!(int::MermaidIntegrator)

Solves the problem using the [MermaidIntegrator](@ref).
This handles all the message passing and calls step! on the [MermaidIntegrator](@ref).

# Arguments
- `int::MermaidIntegrator`: The integrator to be solved.

# Returns
- `MermaidSolution`: The [solution](@ref MermaidSolution) of the problem.
"""
function CommonSolve.solve!(int::MermaidIntegrator)
    t = [0.0]
    dt = minimum([i.component.time_step for i in int.integrators]) # Minimum isnt sufficient to guarantee we don't jump over anything
    sol = MermaidSolution(int)
    update_solution!(sol, int)
    while int.currtime < int.maxt
        CommonSolve.step!(int, dt)
        update_solution!(sol, int)
    end
    return sol
end

"""
    getstate(merInt::MermaidIntegrator, key::ConnectedVariable)

Retrieve the state of a component within a [MermaidIntegrator](@ref) based on the provided [ConnectedVariable](@ref).

# Arguments
- `merInt::MermaidIntegrator`: The integrator containing multiple component integrators.
- `key::ConnectedVariable`: The key specifying which component's state to retrieve.

# Returns
- The state associated with the specified key, or `nothing` if not found.
"""
function getstate(merInt::MermaidIntegrator, key::ConnectedVariable)
    # Get the state of the component based on the key
    for integrator in merInt.integrators
        if integrator.component.name == key.component
            return getstate(integrator, key)
        end
    end
end

"""
    setstate!(merInt::MermaidIntegrator, key::ConnectedVariable, value)

Set the state of a component integrator inside a [MermaidIntegrator](@ref) based on the provided key and value.

# Arguments
- `merInt::MermaidIntegrator`: The integrator containing multiple component integrators.
- `key::ConnectedVariable`: The key specifying which component's state to set.
- `value`: The value to assign to the specified state.
"""
function setstate!(merInt::MermaidIntegrator, key::ConnectedVariable, value)
    # Set the state of the component based on the key
    for integrator in merInt.integrators
        if integrator.component.name == key.component
            setstate!(integrator, key, value)
            return nothing
        end
    end
end

"""
    update_outputs!(compInt::ComponentIntegrator)

Update the outputs field of a [ComponentIntegrator](@ref) based on its current state.

# Arguments
- `compInt::ComponentIntegrator`: The component integrator whose outputs are to be updated.
"""
function update_outputs!(compInt::ComponentIntegrator)
    # Update the outputs of the component based on the current state
    for output_key in keys(compInt.outputs)
        compInt.outputs[output_key] = getstate(compInt, output_key)
    end
end

"""
    update_inputs!(mermaidInt::MermaidIntegrator)

Updates the input values of each [ComponentIntegrator](@ref) within the given [MermaidIntegrator](@ref) instance based on the outputs of other components and the defined [Connectors](@ref Connector).

# Arguments
- `mermaidInt::MermaidIntegrator`: The integrator containing components, connectors, and their current states.

# Description
For each connector in `mermaidInt.connectors`, this function:
- Collects the outputs from the source components specified in the connector's inputs.
- If a function (`conn.func`) is defined for the connector, applies it to the collected inputs to compute the outputs; otherwise, passes the inputs directly.
- Assigns the resulting outputs to the appropriate input fields of the destination components specified in the connector's outputs.
"""
function update_inputs!(mermaidInt::MermaidIntegrator)
    # Update the inputs of the ODE component based on the outputs of other components
    for conn in mermaidInt.connectors
        # Get the values of the connectors inputs
        inputs = []
        for input in conn.inputs
            # Find the corresponding integrator
            index = findfirst(i -> i.component.name == input.component, mermaidInt.integrators)
            if index !== nothing
                integrator = mermaidInt.integrators[index]
                # Get the value of the input from the integrator
                push!(inputs, integrator.outputs[input])
            end
        end
        if isnothing(conn.func)
            if length(inputs) == 1
                outputs = inputs[1]
            else
                outputs = inputs
            end
        else
            outputs = conn.func(inputs...)
        end
        # Set the inputs of the corresponding integrators
        for output in conn.outputs
            # Find the corresponding integrator
            index = findfirst(i -> i.component.name == output.component, mermaidInt.integrators)
            if index !== nothing
                integrator = mermaidInt.integrators[index]
                # Set the input value for the integrator
                integrator.inputs[output] = outputs
            end
        end
    end
end

"""
    update_solution!(sol::MermaidSolution, merInt::MermaidIntegrator)

Update the [MermaidSolution](@ref) `sol` with the current time and state from the [MermaidIntegrator](@ref).

# Arguments
- `sol::MermaidSolution`: The [MermaidSolution](@ref) to be updated. It contains time points (`t`) and a dictionary of state histories (`u`).
- `merInt::MermaidIntegrator`: The integrator object providing the current time (`currtime`) and state access via `getstate`.

# Description
Appends the current time from `merInt` to `sol.t`. For each key in `sol.u`, retrieves the corresponding state from `merInt` using `getstate` and appends it to the respective vector in `sol.u`.
"""
function update_solution!(sol::MermaidSolution, merInt::MermaidIntegrator)
    # Update the solution with the current time and state
    push!(sol.t, merInt.currtime)
    for key in keys(sol.u)
        push!(sol.u[key], getstate(merInt, key))
    end
end

"""
    parsevariable(name::AbstractString) -> ConnectedVariable

Parses a variable name as a foramtted string to a [ConnectedVariable](@ref).

# Arguments
- `name::AbstractString`: The variable name to parse. It can include an optional index, which may be a single integer (e.g., `"comp.var[3]"`) or a range (e.g., `"comp.var[1:5]"`).

# Returns
- `ConnectedVariable`: The corresponding [ConnectedVariable](@ref), containing a `component`, `variable`, `index` (which can be `nothing`, an `Int`, or a `UnitRange{Int}`), and the original `name`.
"""
function parsevariable(name)
    # Parse the variable name to extract its parts
    component, variable = split(name, ".")
    # Is there a variable index
    if contains(variable, "[")
        variable, index = split(variable, "[")
        # Strip the final "]"
        index = strip(index, ']')
        # This will parse "1:5" to a UnitRange{Int} or "3" to an Int
        index = eval(Meta.parse(index))
    else
        # No index
        index = nothing
    end
    # Is there a duplicated index
    if contains(component, "[")
        component, dupindex = split(component, "[")
        # Strip the final "]"
        dupindex = strip(dupindex, ']')
        # This will parse "1:5" to a UnitRange{Int} or "3" to an Int
        dupindex = eval(Meta.parse(dupindex))
    else
        dupindex = nothing
    end
    return ConnectedVariable(component, variable, index, dupindex)
end

function fullname(var::ConnectedVariable)
    # Construct the full name of the variable
    comp = var.component
    dupindex = isnothing(var.duplicatedindex) ? "" : "[" * string(var.duplicatedindex) * "]"
    variable = var.variable
    index = isnothing(var.variableindex) ? "" : "[" * string(var.variableindex) * "]"
    return comp * dupindex * "." * variable * index
end

"""
    Base.getindex(sol::MermaidSolution, var::AbstractString)

Get the solution array for a variable from a [MermaidSolution](@ref).

# Arguments
- `sol::MermaidSolution`: The solution object.
- `var::AbstractString`: The variable name.

# Returns
- The solution array for the specified variable.
"""
function Base.getindex(sol::MermaidSolution, var::AbstractString)
    var = parsevariable(var)
    return Base.getindex(sol, var)
end

"""
    Base.getindex(sol::MermaidSolution, var::ConnectedVariable)

Get the solution array for a variable from a [MermaidSolution](@ref).

# Arguments
- `sol::MermaidSolution`: The solution object.
- `var::ConnectedVariable`: The variable name.

# Returns
- The solution array for the specified variable.
"""
function Base.getindex(sol::MermaidSolution, var::ConnectedVariable)
    if haskey(sol.u, var)
        return sol.u[var]
    else
        # See if we have a key without an index
        var_no_index = ConnectedVariable(var.component, var.variable, nothing, nothing) # TODO I'm not sure how the duplicatedindex data is stored in the solution
        if haskey(sol.u, var_no_index)
            return [i[var.variableindex] for i in sol.u[var_no_index]]
        end
        for key in keys(sol.u)
            if !isnothing(key.variableindex) && !isnothing(var.variableindex)
                if key.variable == var.variable && key.component == var.component && issubset(var.variableindex, key.variableindex)
                    if length(var.variableindex) == 1
                        # If the variableindex is a single value, return the corresponding value
                        return [i[findfirst(x->x==var.variableindex[1], key.variableindex)] for i in sol.u[key]]
                    else
                        return [[i[findfirst(x->x==v, key.variableindex)] for v in var.variableindex] for i in sol.u[key]]
                    end
                end
            end
        end
    end
    throw(KeyError(var))
end

"""
    Base.getindex(sol::MermaidSolution, index::Int)

Get the solution array for a variable from a [MermaidSolution](@ref) at time sol.t[index].

# Arguments
- `sol::MermaidSolution`: The solution object.
- `index::Int`: The time index of the variable.

# Returns
- The solution array for the specified index.
"""
function Base.getindex(sol::MermaidSolution, index::Int)
    if index < 1 || index > length(sol.t)
        throw(BoundsError(sol.t, index))
    end
    return MermaidSolution([sol.t[index]], Dict([var => sol.u[var][index] for var in keys(sol.u)]))
end

"""
(sol::MermaidSolution)(t::Real)

Interpolates the solution at a given time `t` using linear interpolation.

# Arguments
- `sol::MermaidSolution`: The solution object containing time points and state histories.
- `t::Real`: The time at which to interpolate the solution.

# Returns
- `MermaidSolution`: A new [MermaidSolution](@ref) object containing the interpolated time and state.
"""
function (sol::MermaidSolution)(t::Real)
    if t < sol.t[1] || t > sol.t[end]
        throw("Time $t is out of bounds for the solution range [$(sol.t[1]), $(sol.t[end])].")
    end
    lb = findlast(x -> x <= t, sol.t)
    ub = findfirst(x -> x >= t, sol.t)
    if lb == ub
        return sol[lb]
    end
    change = (t-sol.t[lb])/(sol.t[ub]-sol.t[lb])
    return MermaidSolution(
        [sol.t[lb] + change * (sol.t[ub] - sol.t[lb])],
        Dict([var => sol.u[var][lb] .+ change * (sol.u[var][ub] .- sol.u[var][lb]) for var in keys(sol.u)])
    )
end

"""
    inputsandoutputs(integrator::ComponentIntegrator, conns::Vector{Connector}, compName::AbstractString)

Generates the inputs and outputs of a component integrator based on its connections.

# Arguments
- `integrator::ComponentIntegrator`: The component integrator whose inputs and outputs are to be extracted.
- `conns::Vector{Connector}`: The connectors that define the inputs and outputs of the component.

# Returns
- `outputs::OrderedDict{ConnectedVariable,Any}`: An ordered dictionary mapping [ConnectedVariable](@ref) names to their initial values from the component.
- `inputs::OrderedDict{ConnectedVariable,Any}`: An ordered dictionary mapping [ConnectedVariable](@ref) names to their current values (initially 0).
"""
function inputsandoutputs(integrator::ComponentIntegrator, conns::Vector{Connector})
    outputs = OrderedDict{ConnectedVariable,Any}() # Full variable name => Initial value from component
    inputs = OrderedDict{ConnectedVariable,Any}() # Full variable name => Value (initially 0)
    for conn in conns
        # If connection has an input from this component, store its index and function as a ComponentIntegrator.output
        for input in conn.inputs
            if input.component == integrator.component.name
                outputs[input] = getstate(integrator, input)
            end
        end
        for output in conn.outputs
            if output.component == integrator.component.name
                inputs[output] = isnothing(output.variableindex) ? 0 : [0 for _ in output.variableindex]
            end
        end
    end
    return inputs, outputs
end

gettime(integrator::ComponentIntegrator) = getstate(integrator, ConnectedVariable(integrator.component.name, "#time", nothing, nothing))
settime!(integrator::ComponentIntegrator, t) = setstate!(integrator, ConnectedVariable(integrator.component.name, "#time", nothing, nothing), t)

variables(integrator::ComponentIntegrator) = variables(integrator.component)
