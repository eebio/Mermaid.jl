using CommonSolve

function step!(c::AbstractTimeIndependentComponent)
    i = c.inputs
    s = c.state
    m = c.model
    o, s = m(s, i)
    c.outputs = o
    c.state = s
    return nothing
end

function update_inputs!(c::AbstractComponent, component_list::Vector{T}) where T <: AbstractComponent
    # Update the inputs of the ODE component based on the outputs of other components
    i = c.inputs
    for key in keys(i)
        comp_name, variable_name = split(key, ".")
        for comp in component_list
            if comp.name == comp_name
                # Update the input data for the component
                o = comp.outputs
                i[key] = o[variable_name]
            end
        end
    end
end

function update_outputs!(c::ODEComponent)
    # Update the outputs of the ODE component based on the current state
    s = c.state
    for output_key in keys(c.outputs)
        # Update the output data for the component
        c.outputs[output_key] = s[c.output_indices[output_key]]
    end
end

function step!(c::ODEComponent)
    sol = solve(c.model, Euler(), u0=c.state, p=collect(values(c.inputs)), tspan=(0.0, c.time_step), adaptive=false, dt=c.time_step)
    c.state = sol.u[end]
    return nothing
end

function solve!(v::Vector{ODEComponent}, max_t)
    # Solve the system of components
    curr_t = 0.0
    dt = minimum([c.time_step for c in v])
    # Create an array to store the outputs
    outputs = zeros(Int(round(max_t / dt))+1,3)
    for _ in eachrow(outputs)
        for c in v
            if c.time < curr_t
                update_inputs!(c, v)
                step!(c)
                c.time += c.time_step
            end
        end
        for c in v
            update_outputs!(c)
        end
        outputs[Int(round(curr_t / dt)) + 1, :] = pushfirst!([c.outputs["pop"] for c in v], curr_t)
        curr_t += dt
    end
    return v, outputs
end




# Common solve intefaces
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
    integrator = ODEComponentIntegrator(init(c.model, Euler(), adaptive=false, dt=c.time_step), c)
    return integrator
end

"""
    step!(int::MermaidIntegrator, dt::Float64)

Steps the integrator for the given time step.

# Arguments
- `int::MermaidIntegrator`: The integrator to be stepped.
- `dt::Float64`: The time step for the integrator.
"""
function CommonSolve.step!(merInt::MermaidIntegrator, dt::Float64)
    # Step the integrator
    for int in merInt.integrators
        # Update the inputs of the component
        update_inputs!(int, merInt)
        while int.integrator.t <= merInt.currtime
            # Step the integrator
            CommonSolve.step!(int)
        end
    end
    for int in merInt.integrators
        # Update the outputs of the component
        update_outputs!(int)
    end
    # Update the current time
    merInt.currtime += dt
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
    t = []
    u = Dict()
    for i in int.integrators
        for key in keys(i.component.outputs)
            u[i.component.name * "." * key] = []
        end
    end
    while int.currtime <= int.maxt
        push!(t, int.currtime)
        for i in int.integrators
            for key in keys(i.component.outputs)
                push!(u[i.component.name * "." * key], i.component.outputs[key])
            end
        end
        CommonSolve.step!(int, int.integrators[1].component.time_step)
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
    p = collect(values(compInt.component.inputs)) # Need to be smarter about this
    compInt.integrator.p = p
    u_modified!(compInt.integrator, true)
    CommonSolve.step!(compInt.integrator, compInt.component.time_step)
end

function update_outputs!(compInt::ODEComponentIntegrator)
    # Update the outputs of the ODE component based on the current state
    s = compInt.integrator.u
    for output_key in keys(compInt.component.outputs)
        # Update the output data for the component
        compInt.component.outputs[output_key] = s[compInt.component.output_indices[output_key]]
    end
end

function update_inputs!(compInt::ComponentIntegrator, mermaidInt::MermaidIntegrator)
    # Update the inputs of the ODE component based on the outputs of other components
    for key in keys(compInt.component.inputs)
        comp_name, variable_name = split(key, ".")
        for int in mermaidInt.integrators
            if int.component.name == comp_name
                # Update the input data for the component
                o = int.component.outputs
                compInt.component.inputs[key] = o[variable_name]
            end
        end
    end
end

# TODO I think we might still be preturbing the components because the integrators contains references to the component and the component is changed, I probably need a copy for the components during compInt init
