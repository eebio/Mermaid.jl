struct MinimumTimeStepper <: AbstractMermaidSolver
end

function (m::MinimumTimeStepper)(merInt::MermaidIntegrator, dt::Float64)
    # Update the current time
    merInt.currtime += dt
    # Step the integrator
    for int in merInt.integrators
        # Update the inputs of the component
        update_inputs!(int, merInt)
        while int.integrator.t + int.component.time_step <= merInt.currtime
            # Step the integrator
            CommonSolve.step!(int)
        end
    end
    for int in merInt.integrators
        # Update the outputs of the component
        update_outputs!(int)
    end
end
