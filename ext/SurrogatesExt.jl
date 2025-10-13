module SurrogatesExt

using Mermaid
using Surrogates
using Flux
using CommonSolve
using OrderedCollections: OrderedDict
# TODO time dependence: what if an ODE non-autonomous? want an option to include time in the surrogate?

"""
    SurrogateComponent(component::AbstractTimeDependentComponent, lower_bound, upper_bound; name=component.name,
    time_step::Float64=component.time_step,
    state_names::Dict{String,Any}=component.state_names,
    model=nothing,
    n_samples=1000,
    n_epochs=1000)

Represents a component that is replaced with a surrogate in the simulation, speeding up computation of a complex step! function.

# Arguments
- `component`: The original component to be replaced with a surrogate.
- `lower_bound`: Lower bounds for each state variable for surrogate sampling.
- `upper_bound`: Upper bounds for each state variable for surrogate sampling.

# Keyword Arguments
- `name`: Name of the component (default: same as original component)
- `time_step`: Time step for the component (default: same as original component)
- `state_names`: Dictionary mapping variable names (as strings) to their corresponding indices in the
    state vector or symbols from ModelingToolkit/Symbolics (default: same as original component)
- `model`: A Flux.jl model to use as the surrogate. If `nothing`, a default feedforward neural network is created.
- `n_samples`: Number of samples to use for training the surrogate (default: 1000)
- `n_epochs`: Number of training epochs for the surrogate (default: 1000)
"""
function Mermaid.SurrogateComponent(component::AbstractTimeDependentComponent, lower_bound, upper_bound;
    name=component.name,
    time_step::Real=component.time_step,
    state_names=component.state_names,
    model=nothing,
    n_samples=1000,
    n_epochs=1000)
    return Mermaid.SurrogateComponent(component, name, time_step, state_names, lower_bound, upper_bound, model, n_samples, n_epochs)
end

function CommonSolve.init(c::SurrogateComponent, conns::Vector{Connector})
    lower_bound = c.lower_bound
    upper_bound = c.upper_bound
    n_samples = c.n_samples
    integrator = CommonSolve.init(c.component, conns)
    inputs = integrator.inputs
    outputs = integrator.outputs
    initial_state = getstate(integrator)
    function step(x)
        if x isa Tuple
            x = collect(x)
        end
        setstate!(integrator, x)
        CommonSolve.step!(integrator)
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

function Mermaid.getstate(compInt::SurrogateComponentIntegrator, key)
    if first(key.variable) == '#'
        if key.variable == "#time"
            return compInt.time
        end
    end
    setstate!(compInt.integrator, compInt.state)
    return getstate(compInt.integrator, key)
end

function Mermaid.getstate(compInt::SurrogateComponentIntegrator)
    return compInt.state
end

function Mermaid.setstate!(compInt::SurrogateComponentIntegrator, state)
    compInt.state = state
end

function Mermaid.setstate!(compInt::SurrogateComponentIntegrator, key, value)
    if first(key.variable) == '#'
        if key.variable == "#time"
            compInt.time = value
            return nothing
        end
    end
    setstate!(compInt.integrator, compInt.state)
    setstate!(compInt.integrator, key, value)
    compInt.state = getstate(compInt.integrator)
end

function Mermaid.variables(component::SurrogateComponent)
    return variables(component.component)
end

end
