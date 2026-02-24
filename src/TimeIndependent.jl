"""
    TimeIndependentComponent(name::String, func::Function, initial_state)

TimeIndependentComponent represents a component that does not evolve in time, but instead
    computes its state based on its inputs, irregardless of the change in time.

# Arguments
- `name::String`: The name of the component.
- `func::Function`: The function that computes the component's state based on its inputs. It
    should take a single input of the same type as initial_state.
- `initial_state`: The initial state of the component. This should be a valid input to
    `func`.
"""
struct TimeIndependentComponent{F, S} <: AbstractTimeIndependentComponent where {F <: Function, S}
    name::String
    func::F
    initial_state::S
end


mutable struct TimeIndependentComponentIntegrator{I, O, F <: Function} <:
               AbstractComponentIntegrator
    component::TimeIndependentComponent
    input_state::I
    output_state::O
    func::F
    mode::Symbol
end

function CommonSolve.init(c::TimeIndependentComponent)
    return TimeIndependentComponentIntegrator(c, c.initial_state, c.func(c.initial_state), c.func, :output)
end

function getstate(compInt::TimeIndependentComponentIntegrator, key)
    if key.variable == "#state"
        return getstate(compInt)
    end
    return nothing
end

function getstate(compInt::TimeIndependentComponentIntegrator)
    if compInt.mode == :input
        compInt.output_state = compInt.func(compInt.input_state)
        compInt.mode = :output
    end
    return compInt.output_state
end

function setstate!(compInt::TimeIndependentComponentIntegrator, key, value)
    if key.variable == "#state"
        setstate!(compInt, value)
    end
    return nothing
end

function setstate!(compInt::TimeIndependentComponentIntegrator, state)
    compInt.input_state = state
    compInt.mode = :input
    return nothing
end

function variables(compInt::TimeIndependentComponent)
    return ["#state"]
end
