using Mermaid
using StochasticDiffEq

includet("gfp.jl")
includet("cells.jl")

maxt = 600.0

repressilator = get_repressilator()
sde = SDEProblem(repressilator, u0, tspan, ps)

agents = initialize_cell_model()

rep = DEComponent(model=sde, name="repressilator", state_names=Dict("gfp" => repressilator.gfp),
    alg=EM(), time_step=agents.dt)

abm = AgentsComponent(
    model=agents,
    name="cells",
    state_names=Dict("gfp" => :gfp),
    time_step=agents.dt
)

dup = DuplicatedComponent(
    component=rep,
    default_state=sde.u0,
    initial_states=[],
)

conn1 = Connector(
    inputs=["cells.#ids"],
    outputs=["repressilator.#ids"],
)

conn2 = Connector(
    inputs=["repressilator.gfp"],
    outputs=["cells.gfp"],
)

voronoi_marker = (model, cell) -> begin
    verts = get_polygon_coordinates(model.tessellation, cell.index)
    return Makie.Polygon([Point2f(getxy(q) .- cell.pos) for q in verts])
end
voronoi_color(cell) = get(cgrad([:black, :green]), cell.gfp / 1000.0)
fig, ax = abmplot(agents, agent_marker=cell -> voronoi_marker(agents, cell), agent_color=voronoi_color,
    agentsplotkwargs=(strokewidth=1,), figure=(; size=(1600, 800), fontsize=34),
    axis=(; width=800, height=800))
t = [0.0]
gfp1 = [agents[1].gfp]
gfp2 = [agents[2].gfp]
plot_layout = fig[:, end+1] = GridLayout()
gfp1_layout = plot_layout[1, 1] = GridLayout()
ax_1 = Axis(gfp1_layout[1, 1], xlabel="Time", ylabel="GFP 1", width=600, height=400)
lines!(ax_1, t, gfp1, color=:black, label="Total", linewidth=3)
vlines!(ax_1, 0.0, color=:grey, linestyle=:dash, linewidth=3)
Makie.xlims!(ax_1, 0, maxt)
Makie.ylims!(ax_1, 0, 7000)
gfp2_layout = plot_layout[2, 1] = GridLayout()
ax_2 = Axis(gfp2_layout[1, 1], xlabel="Time", ylabel="GFP 2", width=600, height=400)
lines!(ax_2, t, gfp2, color=:black, label="Total", linewidth=3)
vlines!(ax_2, 0.0, color=:grey, linestyle=:dash, linewidth=3)
Makie.xlims!(ax_2, 0, maxt)
Makie.ylims!(ax_2, 0, 7000)
resize_to_layout!(fig)
io = VideoStream(fig; framerate=20)
function plot_input(model)
    if abmtime(model) % 100 == 0
        empty!(ax)
        empty!(ax_1)
        empty!(ax_2)
        abmplot!(ax, model; agent_marker=cell -> voronoi_marker(model, cell), agent_color=voronoi_color,
            agentsplotkwargs=(strokewidth=1,))
        push!(t, abmtime(model) * model.dt)
        push!(gfp1, model[1].gfp)
        push!(gfp2, model[2].gfp)
        lines!(ax_1, t, gfp1, color=:black, label="Total", linewidth=3)
        vlines!(ax_1, t[end], color=:grey, linestyle=:dash, linewidth=3)
        lines!(ax_2, t, gfp2, color=:black, label="Total", linewidth=3)
        vlines!(ax_2, t[end]    , color=:grey, linestyle=:dash, linewidth=3)
        recordframe!(io)
        @show abmtime(model) * model.dt
        @show nagents(model)
    end
end

conn3 = Connector(
    inputs=["cells.#model"],
    outputs=Vector{String}(),
    func=(model) -> plot_input(model)
)

function set_initial_states!(states, ids, model) # Do mutating functions work in Mermaid connectors?
    # init_states is returned, states is mutated
    init_states = Dict{Int,Vector{Float64}}()
    for id in allids(model)
        if id ∉ ids
            # New cell, so either divided or freshing created at start of simulation
            parent = model[id].parent
            if !isnothing(parent)
                states[ids .== parent.id] = states[ids .== parent.id] ./ 2
                init_states[id] = first(deepcopy(states[ids .== parent.id]))
            end
        end
    end
    return init_states
end

conn4 = Connector(
    inputs=["repressilator.#states", "repressilator.#ids", "cells.#model"],
    outputs=["repressilator.#init_states"],
    func = set_initial_states!
)

mp = MermaidProblem(
    components=[dup, abm],
    connectors=[conn4, conn1, conn2, conn3],
    max_t=maxt,
)

alg = MinimumTimeStepper()
sol = solve(mp, alg)

save("examples/outputs/repressilator.mp4", io)

# TODO There might be an assumption in Duplicated that the ids are consecutive, I have an error that tried to set the state of agents at 2365 - BoundsError: attempt to access 2364-element Vector{Any} at index [2365]
