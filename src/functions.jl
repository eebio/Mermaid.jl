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
function CommonSolve.init(prob::MermaidProblem, alg::AbstractMermaidSolver; save_vars=[])
    # Initialize the solver
    integrators = Vector{Any}()
    for c in prob.components
        integrator = CommonSolve.init(c, prob.connectors)
        push!(integrators, integrator)
    end
    return MermaidIntegrator(integrators, prob.connectors, prob.max_t, 0.0, alg, save_vars)
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
    sol = MermaidSolution(int)
    update_solution!(sol, int)
    while int.currtime < int.maxt
        CommonSolve.step!(int, dt)
        update_solution!(sol, int)
    end
    return sol
end

function update_outputs!(compInt::ComponentIntegrator)
    # Update the outputs of the component based on the current state
    for output_key in keys(compInt.outputs)
        compInt.outputs[output_key] = getstate(compInt, output_key)
    end
end

function update_inputs!(mermaidInt::MermaidIntegrator)
    # Update the inputs of the ODE component based on the outputs of other components
    for conn in mermaidInt.connectors
        # Get the values of the connectors inputs
        inputs = []
        for input in conn.inputs
            # Find the corresponding integrator
            index = findfirst(i -> i.component.name == input.component, mermaidInt.integrators)
            if index !== nothing
                integrator = mermaidInt.integrators[index]
                # Get the value of the input from the integrator
                push!(inputs, integrator.outputs[input])
            end
        end
        if isnothing(conn.func)
            if length(inputs) == 1
                outputs = inputs[1]
            else
                outputs = inputs
            end
        else
            outputs = conn.func(inputs...)
        end
        # Set the inputs of the corresponding integrators
        for output in conn.outputs
            # Find the corresponding integrator
            index = findfirst(i -> i.component.name == output.component, mermaidInt.integrators)
            if index !== nothing
                integrator = mermaidInt.integrators[index]
                # Set the input value for the integrator
                integrator.inputs[output] = outputs
            end
        end
    end
end

function update_solution!(sol::MermaidSolution, merInt::MermaidIntegrator)
    # Update the solution with the current time and state
    # TODO this still uses the old style of variable names
    push!(sol.t, merInt.currtime)
    for i in merInt.integrators
        for key in keys(i.component.state_names)
            fullname = join([i.component.name, key], ".")
            if fullname in keys(sol.u)
                push!(sol.u[join([i.component.name, key], ".")], getstate(i, parsevariable(join([i.component.name, key], "."))))
            end
        end
    end
end

function parsevariable(name)
    # Parse the variable name to extract its parts
    component, variable = split(name, ".")
    # Is there a variable index
    if contains(variable, "[")
        variable, index = split(variable, "[")
        # Strip the final "]"
        index = strip(index, ']')
        # Is the index a range
        if contains(index, ":")
            # Extract the range
            start, stop = split(range, ":")
            start = parse(Int, start)
            stop = parse(Int, stop)
            index = start:stop
        else
            # Index is just an integer
            index = parse(Int, index)
        end
    else
        # No index
        index = nothing
    end
    return ConnectedVariable(component, variable, index)
end
