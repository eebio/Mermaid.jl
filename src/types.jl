using OrdinaryDiffEqCore: ODEIntegrator

@kwdef struct Connector
    inputs::Vector{String}
    outputs::Vector{String}
    func = nothing
end

abstract type AbstractComponent end
# Required fields:
# - model: The model of the component
# - name: The name of the component
# - outputs: The outputs of the component
# - inputs: The inputs of the component
# - state: The state of the component

abstract type AbstractTimeIndependentComponent <: AbstractComponent end
# Required fields:
# As above, plus:
# - model: should be a function which takes inputs and state and returns outputs and edits state in place
# - precompute: Bool indicating whether to compute the solution prior to everything else

abstract type AbstractTimeDependentComponent <: AbstractComponent end

abstract type AbstractMermaidSolver end

abstract type ComponentIntegrator end

mutable struct MermaidIntegrator
    integrators::Vector
    connectors::Vector{Connector} # Change name to connectors
    maxt::Float64
    currtime::Float64
    alg::AbstractMermaidSolver
    save_vars::Vector{String}
end

@kwdef struct MermaidProblem
    components::Vector
    connectors::Vector{Connector}
    max_t::Float64 = 1.0
end

struct MermaidSolution
    t::Vector
    u::Dict
end

function MermaidSolution(int::MermaidIntegrator)
    u = Dict()
    for i in int.integrators
        for key in keys(i.component.state_names)
            fullname = join([i.component.name, key], ".")
            if length(int.save_vars) == 0 || fullname in int.save_vars
                u[fullname] = []
            end
        end
    end
    return MermaidSolution([], u)
end
