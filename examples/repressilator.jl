# using Catalyst

# # Constants
# eff = 20
# n = 2
# KM = 40
# tau_mRNA = 2
# tau_prot = 10
# ps_a = 0.5
# ps_0 = 5e-4

# # Derived
# t_ave = tau_mRNA / log(2)
# beta = tau_mRNA / tau_prot
# k_tl = eff / t_ave
# a_tr = (ps_a-ps_0)*60
# a0_tr = ps_0*60
# kd_prot = log(2) / tau_prot
# kd_mRNA = log(2) / tau_mRNA
# alpha = a_tr * eff * tau_prot / (KM * log(2))

# repressilator = @reaction_network repressilator begin
#     hillr(P₃, a_tr, KM, n) + a0_tr, ∅ --> m₁
#     hillr(P₁, a_tr, KM, n) + a0_tr, ∅ --> m₂
#     hillr(P₂, a_tr, KM, n) + a0_tr, ∅ --> m₃
#     kd_mRNA, m₁ --> ∅
#     kd_mRNA, m₂ --> ∅
#     kd_mRNA, m₃ --> ∅
#     k_tl, m₁ --> m₁ + P₁
#     k_tl, m₂ --> m₂ + P₂
#     k_tl, m₃ --> m₃ + P₃
#     kd_prot, P₁ --> ∅
#     kd_prot, P₂ --> ∅
#     kd_prot, P₃ --> ∅
#     hillr(P₁, a_tr, KM, n) + a0_tr, ∅ --> m_gfp
#     kd_mRNA, m_gfp --> ∅
#     k_tl, m_gfp --> m_gfp + gfp
#     kd_prot, gfp --> ∅
# end

# using DifferentialEquations, Plots

# u0 = [:P₁ => 150.0, :P₂ => 140.0, :P₃ => 130.0, :m₁ => 10.0, :m₂ => 12.0, :m₃ => 11.0, :m_gfp => 0.0, :gfp => 0.0]
# tspan = (0., 1500.)
# ps = [:k_tl => k_tl, :KM => KM, :a0_tr => a0_tr, :a_tr => a_tr, :kd_prot => kd_prot, :n => n, :kd_mRNA => kd_mRNA]

# prob = ODEProblem(repressilator, u0, tspan, ps)
# sol = solve(prob, Tsit5())

# plot(sol, vars=[:gfp], xlabel="Time", ylabel="Concentration", label="gfp - Cell 1", title="Repressilator Protein Dynamics", linewidth=2)
# sol = solve(prob, Tsit5())
# display(plot!(sol, vars=[:gfp], label="gfp - Cell 2", linestyle=:dash, linewidth=2))

# using StochasticDiffEq
# sde = SDEProblem(repressilator, u0, tspan, ps)
# sol_sde = solve(sde, EM(), dt=0.1)
# plot(sol_sde, vars=[:gfp], xlabel="Time", ylabel="Concentration", label="gfp - Cell 1", title="Repressilator Protein Dynamics (Stochastic)", linewidth=2)
# sol_sde = solve(sde, EM(), dt=0.1)
# display(plot!(sol_sde, vars=[:gfp], label="gfp - Cell 2", linestyle=:dash, linewidth=2))

# using Mermaid

# comp1 = DEComponent(model=sde, name="Repressilator", state_names=Dict{String,Any}("gfp" => repressilator.gfp))

using Agents, StaticArrays

@agent struct Cell(ContinuousAgent{2,Float64})
    gfp::Float64
    const birth::Float64
    death::Float64 = Inf
end

using DelaunayTriangulation
const DT = DelaunayTriangulation
DT.getx(cell::Cell) = cell.pos[1]
DT.gety(cell::Cell) = cell.pos[2]
DT.number_type(::Type{Cell}) = Float64
DT.number_type(::Type{Vector{Cell}}) = Float64
DT.is_point2(::Cell) = true

spring_constant(p, q) = 20.0 # μ
heterotypic_spring_constant(p, q) = 1.0 # μₕₑₜ
drag_coefficient(p) = 1 / 2 # η
mature_cell_spring_rest_length(p, q) = 0.4 # s
expansion_rate(p, q) = 0.05 * mature_cell_spring_rest_length(p, q) # ε
cutoff_distance(p, q) = 1.5 # ℓₘₐₓ
intrinsic_proliferation_rate(p) = 0.4 # β
carrying_capacity_density(p) = 100.0^2 # K
min_division_age(p) = 1.0 # tₘᵢₙ
max_division_age(p) = 15.0 # tₘₐₓ
max_age(p) = 10.0 # dₘₐₓ
death_rate(p) = 0.001 # psick
mutation_probability(p) = 0.3 # pₘᵤₜ
min_area(p) = 1e-2 # Aₘᵢₙ

using LinearAlgebra
spring_constant(model, i::Int, j::Int, t) = spring_constant(model, model[i], model[j], t)
function spring_constant(model, p, q, t)
    μ = spring_constant(p, q)
    return μ # no adhesion for the initial population
end

rest_length(model, i::Int, j::Int, t) = rest_length(model, model[i], model[j]..., t)
function rest_length(model, p, q, t)
    s = mature_cell_spring_rest_length(p, q)
    ε = expansion_rate(p, q)
    return min(s, (s - ε) * t + ε)
end

function proliferation_rate(model, i::Int, t)
    p = model[i]
    age = t - p.birth
    tₘᵢₙ = min_division_age(p)
    tₘₐₓ = max_division_age(p)
    A = get_area(model.tessellation, i)
    if age ≤ tₘᵢₙ || age ≥ tₘₐₓ || A < min_area(p)
        return 0.0
    end
    vorn = model.tessellation
    Aᵢ = get_area(vorn, i)
    β = intrinsic_proliferation_rate(p)
    K = carrying_capacity_density(p)
    return max(0.0, β * (1 - 1 / (K * Aᵢ)))
end

force(model, i::Int, j::Int, t) = force(model, model[i], model[j], t)
function force(model, p, q, t)
    δ = norm(p.pos - q.pos)
    if δ > cutoff_distance(p, q)
        return SVector(0.0, 0.0)
    end
    μ = spring_constant(model, p, q, t)
    s = rest_length(model, p, q, t)
    rᵢⱼ = q.pos - p.pos
    return μ * (norm(rᵢⱼ) - s) * rᵢⱼ / norm(rᵢⱼ)
end
function force(model, i::Int, t)
    F = SVector(0.0, 0.0)
    for j in get_neighbours(model.triangulation, i) # TODO Does this include dead cells?
        DT.is_ghost_vertex(j) && continue
        F = F + force(model, i, j, t)
    end
    return F
end

velocity(model, i, t) = force(model, i, t) / drag_coefficient(model[i])

function update_velocities!(model, t)
    for i in each_solid_vertex(model.triangulation)
        model[i].vel = velocity(model, i, t)
    end
    return model
end
function update_positions!(model, t)
    update_velocities!(model, t)
    for i in each_solid_vertex(model.triangulation)
        move_agent!(model[i], model, model.dt)
    end
    return model
end

function proliferation_probability(model, t)
    Δt = model.dt
    probs = zeros(nagents(model)) # Technically nagents is not the number of alive agents, but with the way we are handling agents this is correct
    for i in allids(model)
        if !DT.has_vertex(model.triangulation, i) || i in model.dead_cells
            i > 1 && (probs[i] = probs[i-1])
            continue
        end
        Gᵢ = proliferation_rate(model, i, t)
        if i > 1
            probs[i] = probs[i-1] + Gᵢ * Δt
        else
            probs[i] = Gᵢ * Δt
        end
    end
    return probs
end
function select_proliferative_cell(model, probs)
    E = probs[end]
    u = rand() * E
    i = searchsortedlast(probs, u) + 1 # searchsortedlast instead of searchsortedfirst since we skip over some agents in probs
    return i
end

function sample_triangle(tri::Triangulation, T)
    i, j, k = triangle_vertices(T)
    p, q, r = get_point(tri, i, j, k)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    a = (qx - px, qy - py)
    b = (rx - px, ry - py)
    u₁, u₂ = rand(), rand()
    if u₁ + u₂ > 1
        u₁, u₂ = 1 - u₁, 1 - u₂
    end
    ax, ay = getxy(a)
    bx, by = getxy(b)
    wx, wy = u₁ * ax + u₂ * bx, u₁ * ay + u₂ * by
    return SVector(px + wx, py + wy)
end

using StreamSampling
function random_triangle(tri::Triangulation)
    triangles = DT.each_solid_triangle(tri)
    area(T) = DT.triangle_area(get_point(tri, triangle_vertices(T)...)...)
    T = itsample(triangles, area)
    return T
end

function triangulate_voronoi_cell(vorn::VoronoiTessellation, i)
    S = @view get_polygon(vorn, i)[1:end-1]
    points = DT.get_polygon_points(vorn)
    return triangulate_convex(points, S)
end
function sample_voronoi_cell(vorn::VoronoiTessellation, i)
    tri = triangulate_voronoi_cell(vorn, i)
    T = random_triangle(tri)
    return sample_triangle(tri, T)
end

function place_daughter_cell!(model, i, t)
    parent = model[i]
    daughter = sample_voronoi_cell(model.tessellation, i) # this is an SVector, not a Cell
    add_agent!(daughter, model; birth=t, gfp = parent.gfp, vel=SVector(0.0, 0.0)) # HERE we want the cell to be nearby, not sure this does it? maybe daughter is a nearby position
    return daughter
end
function proliferate_cells!(model, t)
    probs = proliferation_probability(model, t)
    u = rand()
    event = u < probs[end]
    !event && return false
    i = select_proliferative_cell(model, probs)
    daughter = place_daughter_cell!(model, i, t)
    return true
end

function cull_cell!(model, i, t)
    p = model[i]
    elder = t - p.birth > max_age(p)
    sick = rand() < model.dt * death_rate(p)
    xmax, ymax = spacesize(model)
    x, y = p.pos
    distance_to_origin = norm(p.pos-SVector(xmax / 2, ymax / 2)) # distance to the center of the space
    outside = distance_to_origin > 0.5 * min(xmax, ymax)
    if elder || sick || outside
        push!(model.dead_cells, i)
        p.death = t
    end
    return model
end
function cull_cells!(model, t)
    for i in each_solid_vertex(model.triangulation)
        cull_cell!(model, i, t)
    end
    return model
end

function model_step!(model)
    stepn = abmtime(model)
    t = stepn * model.dt
    cull_cells!(model, t)
    proliferate_cells!(model, t)
    update_positions!(model, t)
    model.triangulation = retriangulate(model.triangulation, collect(allagents(model)))
    model.tessellation = voronoi(model.triangulation, clip=true)
    return model
end

using Random
function initialize_cell_model(;
    ninit=10,
    radius=2.0,
    dt=0.01,
    sides=SVector(20.0, 20.0))
    Random.seed!(0)
    # Generate the initial random positions
    cent = SVector(sides[1] / 2, sides[2] / 2)
    cells = map(1:ninit) do i
        θ = 2π * rand()
        r = radius * sqrt(rand())
        pos = cent + SVector(r * cos(θ), r * sin(θ))
        cell = Cell(; id=i, pos=pos,
            birth=0.0, gfp = rand(), vel=SVector(0.0, 0.0))
    end
    positions = [cell.pos for cell in cells]

    # Compute the triangulation and the tessellation
    triangulation = triangulate(cells)
    tessellation = voronoi(triangulation, clip=true)

    # Define the model parameters
    properties = Dict(
        :triangulation => triangulation,
        :tessellation => tessellation,
        :dt => dt,
        :dead_cells => Set{Int}()
    )

    # Define the space
    space = ContinuousSpace(sides; periodic=false)

    # Define the model
    model = StandardABM(Cell, space; model_step!, properties)

    # Add the agents
    for (id, pos) in pairs(positions)
        add_agent!(pos, model; birth=0.0, gfp=rand(), vel=SVector(0.0, 0.0))
    end

    return model
end

count_total(model) = num_solid_vertices(model.triangulation)

using StatsBase
function average_cell_area(model)
    area_itr = (get_area(model.tessellation, i) for i in each_solid_vertex(model.triangulation))
    mean_area = mean(area_itr)
    return mean_area
end
function cell_diameter(vorn, i)
    S = get_polygon(vorn, i)
    # This is an O(|S|^2) method, but |S| is small so it is fine
    max_d = 0.0
    for i in S
        p = get_polygon_point(vorn, i)
        for j in S
            i == j && continue
            q = get_polygon_point(vorn, j)
            d = norm(getxy(p) .- getxy(q))
            max_d = max(max_d, d)
        end
    end
    return max_d
end
function average_cell_diameter(model)
    diam_itr = (cell_diameter(model.tessellation, i) for i in each_solid_vertex(model.triangulation))
    mean_diam = mean(diam_itr)
    return mean_diam
end
function average_spring_length(model)
    spring_itr = (norm(model[i].pos - model[j].pos) for (i, j) in each_solid_edge(model.triangulation))
    mean_spring = mean(spring_itr)
    return mean_spring
end

finalT = 50.0
model = initialize_cell_model()
nsteps = Int(finalT / model.dt)
mdata = [count_total,
    average_cell_area, average_cell_diameter, average_spring_length]
agent_df, model_df = run!(model, nsteps; mdata);

using CairoMakie
time = 0:model.dt:finalT
fig = Figure(fontsize=24)
ax = Axis(fig[1, 1], xlabel="Time", ylabel="Count", width=600, height=400)
lines!(ax, time, model_df[!, :count_total], color=:black, label="Total", linewidth=3)
axislegend(ax, position=:lt)
ax = Axis(fig[1, 2], xlabel="Time", ylabel="Average", width=600, height=400)
lines!(ax, time, model_df[!, :average_cell_area], color=:black, label="Cell area", linewidth=3)
lines!(ax, time, model_df[!, :average_cell_diameter], color=:magenta, label="Cell diameter", linewidth=3)
lines!(ax, time, model_df[!, :average_spring_length], color=:red, label="Spring length", linewidth=3)
axislegend(ax, position=:rb)
resize_to_layout!(fig)
fig

voronoi_marker = (model, cell) -> begin
    id = cell.id
    verts = get_polygon_coordinates(model.tessellation, id)
    return Makie.Polygon([Point2f(getxy(q) .- cell.pos) for q in verts])
end
voronoi_color(cell) = get(cgrad([:black, :green]), cell.gfp)
model = initialize_cell_model() # reinitialise the model for the animation
fig, ax, amobs = abmplot(model, agent_marker=cell -> voronoi_marker(model, cell), agent_color=voronoi_color,
    agentsplotkwargs=(strokewidth=1,), figure=(; size=(1600, 800), fontsize=34), mdata=mdata,
    axis=(; width=800, height=800), when=10)
current_time = Observable(0.0)
t = Observable([0.0])
ntotal = Observable(amobs.mdf[][!, :count_total])
avg_area = Observable(amobs.mdf[][!, :average_cell_area])
avg_diam = Observable(amobs.mdf[][!, :average_cell_diameter])
avg_spring = Observable(amobs.mdf[][!, :average_spring_length])
plot_layout = fig[:, end+1] = GridLayout()
count_layout = plot_layout[1, 1] = GridLayout()
ax_count = Axis(count_layout[1, 1], xlabel="Time", ylabel="Count", width=600, height=400)
lines!(ax_count, t, ntotal, color=:black, label="Total", linewidth=3)
vlines!(ax_count, current_time, color=:grey, linestyle=:dash, linewidth=3)
xlims!(ax_count, 0, finalT)
ylims!(ax_count, 0, 800)
avg_layout = plot_layout[2, 1] = GridLayout()
ax_avg = Axis(avg_layout[1, 1], xlabel="Time", ylabel="Average", width=600, height=400)
lines!(ax_avg, t, avg_area, color=:black, label="Cell area", linewidth=3)
lines!(ax_avg, t, avg_diam, color=:magenta, label="Cell diameter", linewidth=3)
lines!(ax_avg, t, avg_spring, color=:red, label="Spring length", linewidth=3)
vlines!(ax_avg, current_time, color=:grey, linestyle=:dash, linewidth=3)
axislegend(ax_avg, position=:rt)
xlims!(ax_avg, 0, finalT)
ylims!(ax_avg, 0, 2)
resize_to_layout!(fig)
on(amobs.mdf) do mdf
    current_time[] = abmtime(amobs.model[]) * model.dt
    t.val = mdf[!, :time] .* model.dt
    ntotal[] = mdf[!, :count_total]
    avg_area[] = mdf[!, :average_cell_area]
    avg_diam[] = mdf[!, :average_cell_diameter]
    avg_spring[] = mdf[!, :average_spring_length]
end

record(fig, "delaunay_model.mp4", 1:(nsteps÷10), framerate=24) do i
    step!(amobs, 10)
end
