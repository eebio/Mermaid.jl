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
    for (int, timescale) in zip(merInt.integrators, merInt.timescales)
        next_t = (gettime(int) + timestep(int)) * timescale
        if next_t < min_t
            min_t = next_t
        end
    end
    # Stop early if the user requested to save somewhere
    if merInt.saveat isa AbstractVector && any(merInt.currtime .< merInt.saveat .< min_t)
        min_t = first(merInt.saveat[merInt.saveat .> merInt.currtime])
    end
    merInt.currtime = min_t
    # Apply connections
    for conn in merInt.connectors
        runconnection!(merInt, conn)
    end
    # Step the integrator
    for (int, timescale) in zip(merInt.integrators, merInt.timescales)
        if (gettime(int) + timestep(int)) * timescale <= nextfloat(merInt.currtime)
            # Step the integrator
            step!(int)
            # Force time synchronization after stepping to avoid floating point issues.
            # Especially important for handling multiple timescales.
            if gettime(int) * timescale == nextfloat(merInt.currtime)
                settime!(int, merInt.currtime / timescale)
            end
        end
    end
end
