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
    # Apply connections
    for conn in merInt.connectors
        runconnection!(merInt, conn)
    end
    # Step the integrator
    for int in merInt.integrators
        if gettime(int) + int.component.time_step <= merInt.currtime
            # Step the integrator
            CommonSolve.step!(int)
        end
    end
end
