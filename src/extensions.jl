export DEComponent, DEComponentIntegrator
export AgentsComponent, AgentsComponentIntegrator
export MOLComponent, MOLComponentIntegrator
export SurrogateComponent, SurrogateComponentIntegrator

struct DEComponent{A, B, C, D, E, F} <: AbstractTimeDependentComponent
    model::A
    name::B
    state_names::C
    time_step::D
    alg::E
    intkwargs::F
end

mutable struct DEComponentIntegrator{A, B, C, D} <: AbstractComponentIntegrator
    integrator::A
    component::B
    outputs::C
    inputs::D
end

struct AgentsComponent{A, B, C, D} <: AbstractTimeDependentComponent
    model::A
    name::B
    state_names::C
    time_step::D
end

mutable struct AgentsComponentIntegrator{A, B, C, D} <: AbstractComponentIntegrator
    integrator::A
    component::B
    outputs::C
    inputs::D
end

struct MOLComponent{A, B, C, D, E, F} <: AbstractTimeDependentComponent
    model::A
    name::B
    state_names::C
    time_step::D
    alg::E
    intkwargs::F
end

mutable struct MOLComponentIntegrator{A, B, C, D} <: AbstractComponentIntegrator
    integrator::A
    component::B
    outputs::C
    inputs::D
end

struct SurrogateComponent{A, B, C, D, E, F, G, H, I} <: AbstractTimeDependentComponent
    component::A
    name::B
    surrogate::C
    time_step::D
    state_names::E
    lower_bound::F
    upper_bound::G
    n_samples::H
    kwargs::I
end

mutable struct SurrogateComponentIntegrator{A, B, C, D, E, F, G} <:
               AbstractComponentIntegrator
    integrator::A
    component::B
    outputs::C
    inputs::D
    state::E
    time::F
    surrogate::G
end
