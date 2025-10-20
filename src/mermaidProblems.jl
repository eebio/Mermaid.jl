using CommonSolve

"""
    MermaidProblem <: AbstractMermaidProblem
    MermaidProblem(;
        components::Vector{AbstractComponent},
        connectors::Vector{AbstractConnector},
        max_t::Float64=1.0)

Defines a Mermaid hybrid simulation problem.

# Keyword Arguments
- `components::Vector{<:AbstractComponent}`: Vector of Components.
- `connectors::Vector{<:AbstractConnector}`: Vector of [Connector](@ref).
- `max_t::Float64`=1.0: Maximum simulation time. Defaults to 1.0.
"""
@kwdef struct MermaidProblem <: AbstractMermaidProblem
    components::Vector{AbstractComponent}
    connectors::Vector{AbstractConnector}
    max_t::Float64 = 1.0
end

"""
    MermaidIntegrator <: AbstractMermaidIntegrator
    MermaidIntegrator(
        integrators::Vector{AbstractComponentIntegrator},
        connectors::Vector{AbstractConnector},
        maxt::Float64,
        currtime::Float64,
        alg::AbstractMermaidSolver,
        save_vars::Vector{String})

Defines the integrator for a Mermaid hybrid simulation.

# Fields
- `integrators::Vector{<:AbstractComponentIntegrator}`: Vector of component integrators.
- `connectors::Vector{<:AbstractConnector}`: Vector of [Connector](@ref) between components.
- `maxt::Float64`: Maximum simulation time.
- `currtime::Float64`: Current simulation time.
- `alg::AbstractMermaidSolver`: The Mermaid solver algorithm to be used.
- `save_vars::Vector{String}`: Variables to be saved during the simulation.
"""
mutable struct MermaidIntegrator <: AbstractMermaidIntegrator
    integrators::Vector{AbstractComponentIntegrator}
    connectors::Vector{AbstractConnector}
    maxt::Float64
    currtime::Float64
    alg::AbstractMermaidSolver
    save_vars::Vector{String}
end

function CommonSolve.init(
        prob::AbstractMermaidProblem, alg::AbstractMermaidSolver; save_vars = [])
    # Sort connectors so that the ones with # in outputs are first
    connectors = sort(prob.connectors;
        by = x -> any([contains(conn.variable, "#") for conn in x.outputs]), rev = true)
    # Initialize the solver
    integrators = [CommonSolve.init(c, connectors) for c in prob.components]
    return MermaidIntegrator(integrators, connectors, prob.max_t, 0.0, alg, save_vars)
end

function CommonSolve.step!(merInt::AbstractMermaidIntegrator)
    CommonSolve.step!(merInt, merInt.alg)
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
    update_solution!(sol, merInt)
    while merInt.currtime < merInt.maxt
        CommonSolve.step!(merInt)
        update_solution!(sol, merInt)
    end
    return sol
end

function getstate(merInt::AbstractMermaidIntegrator, key::AbstractConnectedVariable)
    # Get the state of the component based on the key
    for integrator in merInt.integrators
        if integrator.component.name == key.component
            return getstate(integrator, key)
        end
    end
end

function setstate!(merInt::AbstractMermaidIntegrator, key::AbstractConnectedVariable, value)
    # Set the state of the component based on the key
    for integrator in merInt.integrators
        if integrator.component.name == key.component
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
        int, ConnectedVariable(int.component.name, "#time", nothing, nothing))
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
        ConnectedVariable(int.component.name, "#time", nothing, nothing), t)
end

variables(integrator::AbstractComponentIntegrator) = variables(integrator.component)
