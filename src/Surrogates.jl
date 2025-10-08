using Surrogates
using Flux
# TODO time dependence: what if an ODE non-autonomous? want an option to include time in the surrogate?

"""
    SurrogateComponent <: AbstractComponent

Represents a component that is replaced with a surrogate in the simulation, speeding up computation of a complex step! function.

# Fields
- `component::AbstractTimeDependentComponent`: The original component to be replaced with a surrogate.
- `name::String`: Name of the surrogate component.
"""
@kwdef struct SurrogateComponent <: AbstractTimeDependentComponent
    component::AbstractTimeDependentComponent
    name::String = component.name
    time_step::Float64 = component.time_step
    state_names = component.state_names
    lower_bound
    upper_bound
    model = nothing
    n_samples = 100
    n_epochs = 100
end

mutable struct SurrogateComponentIntegrator <: ComponentIntegrator
    integrator::ComponentIntegrator
    component::SurrogateComponent
    outputs::OrderedDict{ConnectedVariable,Any}
    inputs::OrderedDict{ConnectedVariable,Any}
    state
    time::Float64
    surrogate
end

function CommonSolve.init(c::SurrogateComponent, conns::Vector{Connector})
    lower_bound = c.lower_bound
    upper_bound = c.upper_bound
    n_samples = c.n_samples
    integrator = init(c.component, conns)
    inputs = integrator.inputs
    outputs = integrator.outputs
    initial_state = getstate(integrator)
    function step(x)
        setstate!(integrator, x) # TODO x will be a tuple, will that be okay?
        step!(integrator)
        return getstate(integrator)
    end
    if isnothing(c.model)
        model1 = f64(Chain(
            Dense(length(lower_bound), 32, relu),
            Dense(32, 32, relu),
            Dense(32, length(lower_bound))
        ))
    else
        model1 = c.model
    end
    xys = sample(n_samples, lower_bound, upper_bound, SobolSample())
    zs = step.(xys)
    if length(lower_bound) == 1
        xys = [[x] for x in xys]
    end

    learning_rate = 0.1
    optimizer = Descent(learning_rate)  # Simple gradient descent. See Flux documentation for other options.

    sgt = NeuralSurrogate(xys, zs, lower_bound, upper_bound; model=model1,
        opt=optimizer, n_epochs=c.n_epochs)

    return SurrogateComponentIntegrator(integrator, c, inputs, outputs, initial_state, 0.0, sgt)
end

function CommonSolve.step!(compInt::SurrogateComponentIntegrator)
    surr = compInt.surrogate(compInt.state)
    if length(surr) == 1
        compInt.state = surr[1]  # Convert 1-element array to scalar
    else
        compInt.state = vec(surr)
    end
    compInt.time += compInt.component.time_step
end

function getstate(compInt::SurrogateComponentIntegrator, key)
    setstate!(compInt.integrator, compInt.state)
    return getstate(compInt.integrator, key)
end

function getstate(compInt::SurrogateComponentIntegrator)
    return compInt.state
end

function setstate!(compInt::SurrogateComponentIntegrator, state)
    compInt.state = state
end

function setstate!(compInt::SurrogateComponentIntegrator, key, value)
    setstate!(compInt.integrator, compInt.state)
    setstate!(compInt.integrator, key, value)
    compInt.state = getstate(compInt.integrator)
end

function gettime(compInt::SurrogateComponentIntegrator)
    return compInt.time
end

function settime!(compInt::SurrogateComponentIntegrator, t)
    compInt.time = t
end

function variables(component::SurrogateComponent)
    return variables(component.component)
end
