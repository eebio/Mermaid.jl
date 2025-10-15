
"""
    MermaidSolution

Struct for storing the solution of a hybrid Mermaid simulation.

# Fields
- `t::Vector`: Time points.
- `u::Dict`: Dictionary mapping variables to their solution arrays.
"""
struct MermaidSolution <: AbstractMermaidSolution
    t::Vector
    u::Dict
end

"""
    MermaidSolution(int::MermaidIntegrator)

Construct a [MermaidSolution](@ref) from a [MermaidIntegrator](@ref).

# Arguments
- `int::MermaidIntegrator`: The integrator to extract solution structure from.

# Returns
- `MermaidSolution`: The initialized solution object.
"""
function MermaidSolution(int::AbstractMermaidIntegrator)
    u = Dict()
    if length(int.save_vars) == 0
        for i in int.integrators
            for key in variables(i)
                if contains(key, "#")
                    continue
                end
                fullname = join([i.component.name, key], ".")
                u[parsevariable(fullname)] = []
            end
        end
    else
        for var in int.save_vars
            u[parsevariable(var)] = []
        end
    end
    return MermaidSolution([], u)
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
function update_solution!(sol::AbstractMermaidSolution, merInt::AbstractMermaidIntegrator)
    # Update the solution with the current time and state
    push!(sol.t, merInt.currtime)
    for key in keys(sol.u)
        push!(sol.u[key], getstate(merInt, key))
    end
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
function Base.getindex(sol::AbstractMermaidSolution, var::AbstractString)
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
function Base.getindex(sol::AbstractMermaidSolution, var::AbstractConnectedVariable)
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
                if key.variable == var.variable && key.component == var.component &&
                   issubset(var.variableindex, key.variableindex)
                    if length(var.variableindex) == 1
                        # If the variableindex is a single value, return the corresponding value
                        return [i[findfirst(
                                    x -> x == var.variableindex[1], key.variableindex)]
                                for i in sol.u[key]]
                    else
                        return [[i[findfirst(x -> x == v, key.variableindex)]
                                 for v in var.variableindex] for i in sol.u[key]]
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
function Base.getindex(sol::AbstractMermaidSolution, index::Int)
    if index < 1 || index > length(sol.t)
        throw(BoundsError(sol.t, index))
    end
    return MermaidSolution(
        [sol.t[index]], Dict([var => sol.u[var][index] for var in keys(sol.u)]))
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
function (sol::AbstractMermaidSolution)(t::Real)
    if t < sol.t[1] || t > sol.t[end]
        throw("Time $t is out of bounds for the solution range [$(sol.t[1]), $(sol.t[end])].")
    end
    lb = findlast(x -> x <= t, sol.t)
    ub = findfirst(x -> x >= t, sol.t)
    if lb == ub
        return sol[lb]
    end
    change = (t - sol.t[lb]) / (sol.t[ub] - sol.t[lb])
    return MermaidSolution(
        [sol.t[lb] + change * (sol.t[ub] - sol.t[lb])],
        Dict([var => sol.u[var][lb] .+ change * (sol.u[var][ub] .- sol.u[var][lb])
              for var in keys(sol.u)])
    )
end
