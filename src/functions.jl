using CommonSolve

# Common solve interfaces
"""
    init(prob::MermaidProblem, alg::MermaidSolver; kwargs...)

Converts a MermaidProblem to a MermaidIntegrator.

# Arguments
- `prob::MermaidProblem`: The problem to be solved.
- `alg::MermaidSolver`: The algorithm to be used for solving the problem.
- `kwargs...`: Additional keyword arguments for the solver.

# Returns
- `MermaidIntegrator`: The initialized integrator for the problem.
"""
function CommonSolve.init(prob::MermaidProblem, alg::MermaidSolver)
    # Initialize the solver
    integrators = Vector{ODEComponentIntegrator}()
    for c in prob.components
        integrator = CommonSolve.init(c)
        push!(integrators, integrator)
    end
    return MermaidIntegrator(integrators, prob.max_t, 0.0)
end


"""
    init(c::ODEComponent)
"""
function CommonSolve.init(c::ODEComponent)
    # Initialize the integrator for the component
    outputs = Dict{String,Any}([key => c.model.u0[value] for (key, value) in c.output_indices])
    inputs = Dict{String,Any}([key => 0.0 for key in c.input_names])
    integrator = ODEComponentIntegrator(init(c.model, c.alg; dt=c.time_step, c.intkwargs...), c, outputs, inputs)
    return integrator
end

"""
    step!(int::MermaidIntegrator, dt::Float64)

Steps the integrator for the given time step.

# Arguments
- `int::MermaidIntegrator`: The integrator to be stepped.
- `dt::Float64`: The time step for the integrator.
"""
function CommonSolve.step!(merInt::MermaidIntegrator, dt)
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

"""
    solve!(int::MermaidIntegrator)

Solves the problem using the integrator.

This handles all the message passing and calls step! on the MermaidIntegrator.

# Arguments
- `int::MermaidIntegrator`: The integrator to be solved.

# Returns
- `MermaidSolution`: The solution of the problem.
"""
function CommonSolve.solve!(int::MermaidIntegrator)
    t = [0.0]
    dt = minimum([i.component.time_step for i in int.integrators]) # Minimum isnt sufficient to guarantee we don't jump over anything
    u = Dict()
    for i in int.integrators
        for key in keys(i.outputs)
            u[i.component.name * "." * key] = []
        end
    end
    for i in int.integrators
        for key in keys(i.outputs)
            push!(u[i.component.name*"."*key], i.outputs[key])
        end
    end
    while int.currtime < int.maxt
        CommonSolve.step!(int, dt)
        push!(t, int.currtime)
        for i in int.integrators
            for key in keys(i.outputs)
                push!(u[i.component.name * "." * key], i.outputs[key])
            end
        end
    end
    return MermaidSolution(t, u)
end

"""
    step!(compInt::ODEComponentIntegrator)

Steps the ODE component integrator.
# Arguments
- `compInt::ODEComponentIntegrator`: The ODE component integrator to be stepped.
"""
function CommonSolve.step!(compInt::ODEComponentIntegrator)
    # Update compInt.integrator with the parameters from inputs
    p = collect(values(compInt.inputs)) # Need to be smarter about this
    compInt.integrator.p = p
    u_modified!(compInt.integrator, true)
    CommonSolve.step!(compInt.integrator)
end

function update_outputs!(compInt::ODEComponentIntegrator)
    # Update the outputs of the ODE component based on the current state
    s = compInt.integrator.u
    for output_key in keys(compInt.outputs)
        # Update the output data for the component
        compInt.outputs[output_key] = s[compInt.component.output_indices[output_key]]
    end
end

function update_inputs!(compInt::ComponentIntegrator, mermaidInt::MermaidIntegrator)
    # Update the inputs of the ODE component based on the outputs of other components
    for key in keys(compInt.inputs)
        comp_name, variable_name = split(key, ".")
        for int in mermaidInt.integrators
            if int.component.name == comp_name
                # Update the input data for the component
                o = int.outputs
                compInt.inputs[key] = o[variable_name]
            end
        end
    end
end
