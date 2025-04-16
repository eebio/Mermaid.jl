using OrdinaryDiffEqCore: ODEIntegrator

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


# Predefined concrete types
@kwdef mutable struct ODEComponent <: AbstractComponent # No longer needs to be mutable
    model::ODEProblem
    name::String = "ODE Component"
    outputs::Dict{String,Any} = NamedTuple{String,Any}()
    inputs::Dict{String,Any} = NamedTuple{String,Any}()
    output_indices::Dict{String,Any} = NamedTuple{String,Any}()
    state = Vector{Float64}() # Can be removed
    time_step::Float64 = 1.0
    time::Float64 = 0.0 # Can be removed
end

abstract type ComponentIntegrator end

mutable struct ODEComponentIntegrator <: ComponentIntegrator
    integrator::OrdinaryDiffEqCore.ODEIntegrator
    component::ODEComponent
end

mutable struct MermaidIntegrator
    integrators::Vector
    maxt::Float64
    currtime::Float64
end

@kwdef struct MermaidProblem
    components::Vector
    max_t::Float64 = 1.0
end

struct MermaidSolver
end

struct MermaidSolution
    t::Vector
    u::Dict
end
