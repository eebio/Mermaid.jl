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

abstract type AbstractMermaidSolver end

# Predefined concrete types
@kwdef struct ODEComponent <: AbstractTimeDependentComponent
    model::ODEProblem
    name::String = "ODE Component"
    output_indices::Dict{String,Any} = NamedTuple{String,Any}()
    input_names::Vector{String} = []
    time_step::Float64 = 1.0
    alg = Rodas5()
    intkwargs::Tuple{Pair{Symbol,Any}} = ()
end

abstract type ComponentIntegrator end

mutable struct ODEComponentIntegrator <: ComponentIntegrator
    integrator::OrdinaryDiffEqCore.ODEIntegrator
    component::ODEComponent
    outputs::Dict{String,Any}
    inputs::Dict{String,Any}
end

mutable struct MermaidIntegrator
    integrators::Vector
    maxt::Float64
    currtime::Float64
    alg::AbstractMermaidSolver
end

@kwdef struct MermaidProblem
    components::Vector
    max_t::Float64 = 1.0
end

struct MermaidSolution
    t::Vector
    u::Dict
end
