using CommonSolve

"""
    MinimumTimeStepper()

A solver that steps the integrators in a `MermaidIntegrator` advancing to the next event.
"""
struct MinimumTimeStepper <: AbstractMermaidSolver
end

function CommonSolve.step!(merInt::MermaidIntegrator, ::MinimumTimeStepper)
    # Update the current time
    min_t = Inf
    for int in merInt.integrators
        next_t = gettime(int) + int.component.time_step
        if next_t < min_t
            min_t = next_t
        end
    end
    merInt.currtime = min_t
    # Update the inputs of all components
    update_inputs!(merInt)
    # Step the integrator
    for int in merInt.integrators
        if gettime(int) + int.component.time_step <= merInt.currtime
            # Step the integrator
            CommonSolve.step!(int)
        end
    end
    for int in merInt.integrators
        # Update the outputs of the component
        update_outputs!(int)
    end
end
