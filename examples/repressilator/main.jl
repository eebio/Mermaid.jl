using Mermaid
using StochasticDiffEq

includet("gfp.jl")
includet("cells.jl")
includet("growth.jl")
includet("improved.jl")

using Random
Random.seed!(1234)

using .Repressilator

maxt = 5.0

repressilator = Repressilator.repressilator
sde = SDEProblem(repressilator, Repressilator.u0, Repressilator.tspan, Repressilator.ps)

sde_improved = sde_improved_repressilator()
improved = sde_improved.f.sys

agents = initialize_cell_model()

growth = ode_growth()
g_model = get_growth_model()

rep = DEComponent(model=sde, name="repressilator", state_names=Dict("gfp" => repressilator.gfp,
    "growth_rate" => repressilator.gr),
    alg=EM(), time_step=agents.dt, intkwargs=(:maxiters => Inf, :save_everystep => false))

rep_imp = DEComponent(model=sde_improved, name="repressilator", state_names=Dict("gfp" => improved.gfp,
        "growth_rate" => improved.gr, "volume" => improved.V),
    alg=Tsit5(), time_step=agents.dt, intkwargs=(:maxiters => Inf, :save_everystep => false))

abm = AgentsComponent(
    model=agents,
    name="cells",
    state_names=Dict("gfp" => :gfp, "size" => :size, "nutrients" => :nuts, "nut_import_rate" => :nut_import_rate),
    time_step=agents.dt
)

gro = DEComponent(model=growth, name="growth", state_names=Dict("s" => g_model.s, "λ" => g_model.λ, "mass" => g_model.M, "import" => g_model.ν_imp),
    alg=Rosenbrock23(), time_step=agents.dt, intkwargs=(:maxiters => Inf, :isoutofdomain => (u, p, t) -> any(x -> x < 0, u), :save_everystep => false))

dup_r = DuplicatedComponent(
    component=rep,
    default_state=sde.u0,
    initial_states=[],
)

dup_i = DuplicatedComponent(
    component=rep_imp,
    default_state=sde_improved.u0,
    initial_states=[],
)

dup_g = DuplicatedComponent(
    component=gro,
    default_state=growth.u0,
    initial_states=[],
)

conn_ids_1 = Connector(
    inputs=["cells.#ids"],
    outputs=["repressilator.#ids"],
)

conn_ids_2 = Connector(
    inputs=["cells.#ids"],
    outputs=["growth.#ids"]
)

conn_gfp = Connector(
    inputs=["repressilator.gfp"],
    outputs=["cells.gfp"],
)

voronoi_marker = (model, cell) -> begin
    #return :+
    verts = get_polygon_coordinates(model.tessellation, cell.index)
    return Makie.Polygon([Point2f(getxy(q) .- cell.pos) for q in verts])
end
voronoi_color(cell) = get(cgrad([:black, :green]), cell.gfp/cell.size^3 / 200.0)
fig, ax = abmplot(agents, agent_marker=cell -> voronoi_marker(agents, cell), agent_color=voronoi_color,
    agentsplotkwargs=(strokewidth=1,), figure=(; size=(1600, 800), fontsize=34),
    axis=(; width=800, height=800), heatarray=:nutrients, heatkwargs=(colorrange=(0.0, 1.0),))
abmplot!(ax, agents; agent_marker=:xcross, agent_color=:red, agent_size=cell -> cell.id ∈ [1, 2] ? 10 : 0)
t = [0.0]
gfp1 = [agents[1].gfp/agents[1].size^3]
gfp2 = [agents[2].gfp/agents[2].size^3]
nut1 = [agents[1].nuts]
nut2 = [agents[2].nuts]
size1 = [agents[1].size^3]
size2 = [agents[2].size^3]
plot_layout = fig[:, end+1] = GridLayout()
gfp1_layout = plot_layout[1, 1] = GridLayout()
ax_1 = Axis(gfp1_layout[1, 1], xlabel="Time", ylabel="GFP 1", width=600, height=400)
ax_1_2 = Axis(gfp1_layout[1, 1], ylabel="size 1", yaxisposition=:right, yticklabelcolor=:blue)
lines!(ax_1, t, gfp1, color=:black, label="Total", linewidth=3)
lines!(ax_1_2, t, size1, color=:blue, label="size", linewidth=3)
lines!(ax_1_2, t, nut1, color=:green, label="nutrients", linewidth=3)
vlines!(ax_1, 0.0, color=:grey, linestyle=:dash, linewidth=3)
Makie.xlims!(ax_1, 0, maxt)
Makie.ylims!(ax_1, 0, 2000)
Makie.xlims!(ax_1_2, 0, maxt)
Makie.ylims!(ax_1_2, 0, 5.0)
gfp2_layout = plot_layout[2, 1] = GridLayout()
ax_2 = Axis(gfp2_layout[1, 1], xlabel="Time", ylabel="GFP 2", width=600, height=400)
ax_2_2 = Axis(gfp2_layout[1, 1], ylabel="size 2", yaxisposition=:right, yticklabelcolor=:blue)
lines!(ax_2, t, gfp2, color=:black, label="Total", linewidth=3)
lines!(ax_2_2, t, size2, color=:blue, label="size", linewidth=3)
lines!(ax_2_2, t, nut2, color=:green, label="nutrients", linewidth=3)
vlines!(ax_2, 0.0, color=:grey, linestyle=:dash, linewidth=3)
Makie.xlims!(ax_2, 0, maxt)
Makie.ylims!(ax_2, 0, 2000)
Makie.xlims!(ax_2_2, 0, maxt)
Makie.ylims!(ax_2_2, 0, 5.0)
resize_to_layout!(fig)
io = VideoStream(fig; framerate=20)
function plot_input(model)
    push!(t, abmtime(model) * model.dt)
    push!(gfp1, model[1].gfp/model[1].size^3)
    push!(gfp2, model[2].gfp/model[2].size^3)
    push!(nut1, model[1].nuts)
    push!(nut2, model[2].nuts)
    push!(size1, model[1].size^3)
    push!(size2, model[2].size^3)
    if abmtime(model) % 10 == 0
        empty!(ax)
        empty!(ax_1)
        empty!(ax_2)
        empty!(ax_1_2)
        empty!(ax_2_2)
        abmplot!(ax, model; agent_marker=cell -> voronoi_marker(model, cell), agent_color=voronoi_color,
            agentsplotkwargs=(strokewidth=1,), heatarray=:nutrients, heatkwargs=(colorrange=(0.0, 1.0),))
        abmplot!(ax, model; agent_marker=:xcross, agent_color=:red, agent_size= cell->cell.id ∈ [1,2] ? 10 : 0)
        lines!(ax_1, t, gfp1, color=:black, label="Total", linewidth=3)
        lines!(ax_1_2, t, size1, color=:blue, label="size", linewidth=3)
        lines!(ax_1_2, t, nut1, color=:green, label="nutrients", linewidth=3)
        vlines!(ax_1, t[end], color=:grey, linestyle=:dash, linewidth=3)
        lines!(ax_2, t, gfp2, color=:black, label="Total", linewidth=3)
        lines!(ax_2_2, t, size2, color=:blue, label="size", linewidth=3)
        lines!(ax_2_2, t, nut2, color=:green, label="nutrients", linewidth=3)
        vlines!(ax_2, t[end], color=:grey, linestyle=:dash, linewidth=3)
        recordframe!(io)
        @show abmtime(model) * model.dt
        @show nagents(model)
    end
end

conn_plot = Connector(
    inputs=["cells.#model"],
    outputs=Vector{String}(),
    func=(model) -> plot_input(model)
)

function set_initial_states!(states, ids, model) # Do mutating functions work in Mermaid connectors?
    # init_states is returned, states is mutated
    init_states = Dict{Int,Vector{Float64}}()
    for id in allids(model)
        if id ∉ ids
            # New cell, so either divided or freshly created at start of simulation
            parent = model[id].parent
            if !isnothing(parent)
                states[ids .== parent.id] = states[ids .== parent.id] ./ 2
                init_states[id] = first(deepcopy(states[ids .== parent.id]))
            end
        end
    end
    return init_states
end

conn_init_states_rep = Connector(
    inputs=["repressilator.#states", "repressilator.#ids", "cells.#model"],
    outputs=["repressilator.#init_states"],
    func = set_initial_states!
)

conn_init_states_growth = Connector(
    inputs=["growth.#states", "growth.#ids", "cells.#model"],
    outputs=["growth.#init_states"],
    func=set_initial_states!
)

conn_gr = Connector(
    inputs=["growth.λ"],
    outputs=["repressilator.growth_rate"],
    #func = (λ) -> ones(length(λ))*0.6,
)

conn_size = Connector(
    inputs=["growth.mass"],
    outputs=["cells.size"],
    func = (M) -> (M ./ 1e8).^(1/3)
)

conn_volume = Connector(
    inputs=["growth.mass"],
    outputs=["repressilator.volume"],
    func = (V) -> (V ./ 1e8)
)

conn_nuts = Connector(
    inputs=["cells.nutrients"],
    outputs=["growth.s"],
    func = (nutrients) -> nutrients .* 1e4,
)

conn_nuts_imp = Connector(
    inputs=["growth.import"],
    outputs=["cells.nut_import_rate"],
    func = (imp) -> imp ./ 1e8*5*0, # Scaling is arbitrary
)

mp = MermaidProblem(
    components=[dup_g, dup_i, abm], # TODO get awkward error when repressilator is run before growth
    connectors=[conn_init_states_rep, conn_init_states_growth, conn_ids_1, conn_ids_2, conn_gfp, conn_gr, conn_size, conn_nuts, conn_volume, conn_nuts_imp, conn_plot],
    max_t=maxt,
)

alg = MinimumTimeStepper()
sol = solve(mp, alg)

save("examples/outputs/repressilator.mp4", io)

# TODO There might be an assumption in Duplicated that the ids are consecutive, I have an error that tried to set the state of agents at 2365 - BoundsError: attempt to access 2364-element Vector{Any} at index [2365]
