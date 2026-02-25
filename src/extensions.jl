export DEComponent, DEComponentIntegrator
export AgentsComponent, AgentsComponentIntegrator
export MOLComponent, MOLComponentIntegrator
export SurrogateComponent, SurrogateComponentIntegrator

struct DEComponent{A, B, C, D, E, F} <: AbstractTimeDependentComponent
    model::A
    name::B
    state_names::C
    timestep::D
    alg::E
    intkwargs::F
end

mutable struct DEComponentIntegrator{A, B} <: AbstractComponentIntegrator
    integrator::A
    component::B
end

struct AgentsComponent{A, B, C, D} <: AbstractTimeDependentComponent
    model::A
    name::B
    state_names::C
    timestep::D
end

mutable struct AgentsComponentIntegrator{A, B} <: AbstractComponentIntegrator
    integrator::A
    component::B
end

struct MOLComponent{A, B, C, D, E, F} <: AbstractTimeDependentComponent
    model::A
    name::B
    state_names::C
    timestep::D
    alg::E
    intkwargs::F
end

mutable struct MOLComponentIntegrator{A, B} <: AbstractComponentIntegrator
    integrator::A
    component::B
end

struct SurrogateComponent{A, B, C, D, E, F, G, H, I} <: AbstractTimeDependentComponent
    component::A
    name::B
    surrogate::C
    timestep::D
    state_names::E
    lower_bound::F
    upper_bound::G
    n_samples::H
    kwargs::I
end

mutable struct SurrogateComponentIntegrator{A, B, C, D, E} <: AbstractComponentIntegrator
    integrator::A
    component::B
    state::C
    time::D
    surrogate::E
end
