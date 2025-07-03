"""
    MinimumTimeStepper()

A solver that steps the integrators in a `MermaidIntegrator` at the minimum time step of all components.
"""
struct MinimumTimeStepper <: AbstractMermaidSolver
end

function (m::MinimumTimeStepper)(merInt::MermaidIntegrator, dt::Float64)
    # Update the current time
    merInt.currtime += dt
    # Update the inputs of all components
    update_inputs!(merInt)
    # Step the integrator
    for int in merInt.integrators
        while gettime(int) + int.component.time_step <= merInt.currtime
            # Step the integrator
            CommonSolve.step!(int)
        end
    end
    for int in merInt.integrators
        # Update the outputs of the component
        update_outputs!(int)
    end
end
