using CommonSolve

"""
    MermaidProblem <: AbstractMermaidProblem
    MermaidProblem(;
        components::Vector{AbstractComponent},
        connectors::Vector{AbstractConnector},
        tspan::Tuple{Float64, Float64},
        timescales::Vector{Float64}=ones(length(components)))

Defines a Mermaid hybrid simulation problem.

# Keyword Arguments
- `components::Vector{<:AbstractComponent}`: Vector of Components.
- `connectors::Vector{<:AbstractConnector}`: Vector of [Connector](@ref).
- `tspan::Tuple{Float64, Float64}`: The time span of the simulation.
- `timescales::Vector{Float64}`=ones(length(components)): Timescales for each component.
    Each component's timescale will be multiplied by the component's time to convert it to
    the universal simulation time.
"""
@kwdef struct MermaidProblem <: AbstractMermaidProblem
    components::Vector{AbstractComponent}
    connectors::Vector{AbstractConnector}
    tspan::Tuple{Float64, Float64}
    timescales::Vector{Float64} = ones(length(components))
end

mutable struct MermaidIntegrator{X <: AbstractMermaidSolver, S <: Union{Function, AbstractVector}} <: AbstractMermaidIntegrator
    integrators::Vector{<:AbstractComponentIntegrator}
    connectors::Vector{<:AbstractConnector}
    tspan::Tuple{Float64, Float64}
    currtime::Float64
    alg::X
    save_vars::Vector{<:AbstractString}
    saveat::S
    timescales::Vector{Float64}
end

"""
    init(prob::AbstractMermaidProblem, alg::AbstractMermaidSolver;
    save_vars = nothing, saveat = nothing)

Defines the integrator for a Mermaid hybrid simulation.

# Arguments
- `prob::AbstractMermaidProblem`: The problem to be solved.
- `alg::AbstractMermaidSolver`: The Mermaid solver algorithm to be used.
- `save_vars`: Variables to be saved during the simulation. Options include:
    - `:all`: Save all variables.
    - `:none`: Save no variables. Equivalent to an empty vector.
    - `Vector{String}`: A vector of connected variable names to save.
    - `nothing`: Default. Save all non-special variables (i.e., variables that do not start with '#').
- `saveat::Function`: When to save the variables during the simulation. Can be a function
    that takes the integrator and current time and returns a boolean, or a vector of time
    points at which to save, or a number representing the time interval (save at
    `tspan[1]:saveat:tspan[2]`).
"""
function CommonSolve.init(prob::AbstractMermaidProblem, alg::AbstractMermaidSolver;
        save_vars = nothing, saveat = nothing)
    # Initialize the solver
    integrators = [init(c) for c in prob.components]

    # Process save_vars
    if isnothing(save_vars) || save_vars == :all
        tmp = String[]
        for int in integrators
            for var in variables(int)
                if var[1] != '#' || save_vars == :all
                    push!(tmp, string(name(int), ".", var))
                end
            end
        end
        save_vars = tmp
    end
    if (save_vars isa AbstractVector && length(save_vars) == 0) || save_vars == :none
        save_vars = String[]
    end

    # Process saveat
    if isnothing(saveat)
        saveat = (integrator, t) -> true
    end
    if saveat isa Number
        saveat = prob.tspan[1]:saveat:prob.tspan[2]
    end
    return MermaidIntegrator(
        integrators, prob.connectors, prob.tspan, 0.0, alg, save_vars, saveat, prob.timescales)
end

function CommonSolve.step!(merInt::AbstractMermaidIntegrator)
    step!(merInt, merInt.alg)
end

"""
    solve!(merInt::AbstractMermaidIntegrator)

Solves the problem using the [MermaidIntegrator](@ref). This handles all the message passing
    and calls step! on the [MermaidIntegrator](@ref).

# Arguments
- `merInt::AbstractMermaidIntegrator`: The integrator to be solved.

# Returns
- `MermaidSolution`: The [solution](@ref MermaidSolution) of the problem.
"""
function CommonSolve.solve!(merInt::AbstractMermaidIntegrator)
    sol = MermaidSolution(merInt)
    if should_save(merInt, merInt.saveat)
        update_solution!(sol, merInt)
    end
    while merInt.currtime < merInt.tspan[2]
        step!(merInt)
        if should_save(merInt, merInt.saveat)
            update_solution!(sol, merInt)
        end
    end
    return sol
end

function should_save(merInt::AbstractMermaidIntegrator, saveat::AbstractVector)
    if merInt.currtime in saveat
        return true
    else
        return false
    end
end

function should_save(merInt::AbstractMermaidIntegrator, saveat::Function)
    return saveat(merInt, merInt.currtime)
end

function getstate(merInt::AbstractMermaidIntegrator, key::AbstractConnectedVariable)
    # Get the state of the component based on the key
    for integrator in merInt.integrators
        if name(integrator) == key.component
            return getstate(integrator, key)
        end
    end
end

function gettime(merInt::AbstractMermaidIntegrator)
    # Get the current time of the integrator
    return merInt.currtime
end

function setstate!(merInt::AbstractMermaidIntegrator, key::AbstractConnectedVariable, value)
    # Set the state of the component based on the key
    for integrator in merInt.integrators
        if name(integrator) == key.component
            setstate!(integrator, key, value)
            return nothing
        end
    end
end

"""
    gettime(merInt::AbstractComponentIntegrator)

Get the current time of the integrator.

# Arguments
- `int::AbstractComponentIntegrator`: The integrator whose time is to be retrieved.

# Returns
- The current time of the integrator.
"""
function gettime(int::AbstractComponentIntegrator)
    getstate(
        int, ConnectedVariable(name(int), "#time", nothing, nothing))
end

"""
    settime!(merInt::AbstractComponentIntegrator, t)

Set the current time of the integrator.

# Arguments
- `int::AbstractComponentIntegrator`: The integrator whose time is to be set.
- `t`: The time to set.
"""
function settime!(int::AbstractComponentIntegrator, t)
    setstate!(int,
        ConnectedVariable(name(int), "#time", nothing, nothing), t)
end

"""
    variables(integrator::AbstractComponent)
    variables(integrator::AbstractComponentIntegrator)

Get the variables names of the component or integrator that can be accessed through getstate
    and setstate!. This includes special variables like `#time`.
"""
variables(integrator::AbstractComponentIntegrator) = variables(integrator.component)

function getstate(int::AbstractComponentIntegrator; copy = false)
    if copy
        return deepcopy(getstate(int))
    else
        return getstate(int)
    end
end

"""
    timestep(int::AbstractComponent)
    timestep(comp::AbstractComponentIntegrator)

Get the proposed time step of the integrator or component. It can depend on the current
    state.
"""
timestep(comp::AbstractComponent) = comp.timestep
timestep(int::AbstractComponentIntegrator) = timestep(int.component)

"""
    name(int::AbstractComponentIntegrator)
    name(comp::AbstractComponent)

Get the name of the integrator or component.
"""
name(comp::AbstractComponent) = comp.name
name(int::AbstractComponentIntegrator) = name(int.component)
