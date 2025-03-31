function step!(c::AbstractTimeIndependentComponent)
    i = c.inputs
    s = c.state
    m = c.model
    o, s = m(s, i)
    c.outputs = o
    c.state = s
    return nothing
end

function update_inputs!(c::AbstractComponent, component_list::Vector{AbstractComponent})
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
    sol = solve(c.model, Tsit5(), u0=c.state, p=id, tspan=(0.0, c.time_step), adaptive=false, dt=c.timestep)
    c.state = sol.u[end]
    return nothing
end

function solve!(v::Vector{ODEComponent}, max_t)
    # Solve the system of components
    curr_t = 0.0
    while curr_t < max_t
        for c in v
            if c.time < curr_t
                update_inputs!(c, v)
                step!(c)
            end
        end
    end
    return v
end
