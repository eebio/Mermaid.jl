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
function CommonSolve.init(prob::MermaidProblem, alg::AbstractMermaidSolver)
    # Initialize the solver
    integrators = Vector{ODEComponentIntegrator}()
    for c in prob.components
        integrator = CommonSolve.init(c, prob.connectors)
        push!(integrators, integrator)
    end
    return MermaidIntegrator(integrators, prob.connectors, prob.max_t, 0.0, alg)
end

"""
    step!(int::MermaidIntegrator, dt::Float64)

Steps the integrator for the given time step.

# Arguments
- `int::MermaidIntegrator`: The integrator to be stepped.
- `dt::Float64`: The time step for the integrator.
"""
function CommonSolve.step!(merInt::MermaidIntegrator, dt)
    merInt.alg(merInt, dt)
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
            u[key] = []
        end
    end
    for i in int.integrators
        for key in keys(i.outputs)
            push!(u[key], i.outputs[key])
        end
    end
    while int.currtime < int.maxt
        CommonSolve.step!(int, dt)
        push!(t, int.currtime)
        for i in int.integrators
            for key in keys(i.outputs)
                push!(u[key], i.outputs[key])
            end
        end
    end
    return MermaidSolution(t, u)
end

function update_outputs!(compInt::ODEComponentIntegrator)
    # Update the outputs of the ODE component based on the current state
    s = compInt.integrator.u
    for output_key in keys(compInt.outputs)
        # Update the output data for the component
        _, var_name = split(output_key, ".")
        index = compInt.component.state_names[var_name]
        # If the index is a MTK symbol then get the variable index
        if !isa(index, Integer)
            index = variable_index(compInt.component.model.f.sys, index)
        end
        compInt.outputs[output_key] = s[index]
    end
end

function update_inputs!(mermaidInt::MermaidIntegrator)
    # Update the inputs of the ODE component based on the outputs of other components
    for conn in mermaidInt.connectors
        # Get the values of the connectors inputs
        inputs = []
        for input in conn.inputs
            # Get the component name and variable name
            comp_name, var_name = split(input, ".")
            # Find the corresponding integrator
            index = findfirst(i -> i.component.name == comp_name, mermaidInt.integrators)
            if index !== nothing
                integrator = mermaidInt.integrators[index]
                # Get the value of the input from the integrator
                push!(inputs, integrator.outputs[input])
            end
        end
        outputs = conn.func(inputs...)
        # Set the inputs of the corresponding integrators
        for output in conn.outputs
            # Get the component name and variable name
            comp_name, var_name = split(output, ".")
            # Find the corresponding integrator
            index = findfirst(i -> i.component.name == comp_name, mermaidInt.integrators)
            if index !== nothing
                integrator = mermaidInt.integrators[index]
                # Set the input value for the integrator
                integrator.inputs[integrator.component.name * "." * var_name] = outputs
            end
        end
    end
end
