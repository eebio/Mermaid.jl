using Agents, StaticArrays
using Random
using LinearAlgebra
using StatsBase
using CairoMakie
using StreamSampling
using DelaunayTriangulation
const DT = DelaunayTriangulation

@agent struct Cell(ContinuousAgent{2,Float64})
    gfp::Float64
    const birth::Float64
    death::Float64 = Inf
    index::Int = typemin(Int)
    const parent::Union{Cell, Nothing} = nothing
    nut_import_rate::Float64 = 0.0
    nuts::Float64
    size::Float64
end

DT.getx(cell::Cell) = cell.pos[1]
DT.gety(cell::Cell) = cell.pos[2]
DT.number_type(::Type{Cell}) = Float64
DT.number_type(::Type{Vector{Cell}}) = Float64
DT.is_point2(::Cell) = true

spring_constant(p, q) = 20.0 # μ
heterotypic_spring_constant(p, q) = 1.0 # μₕₑₜ
drag_coefficient(p) = 1.0 # η
mature_cell_spring_rest_length(p, q) = p.size + q.size # s
expansion_rate(p, q) = 0.05 * mature_cell_spring_rest_length(p, q) # ε
cutoff_distance(p, q) = 1.5 # ℓₘₐₓ
intrinsic_proliferation_rate(p) = 0.05 # β
carrying_capacity_density(p) = 100.0^2 # K
max_age(p) = 10.0 # dₘₐₓ
min_division_age(p) = 15.0 # tₘᵢₙ
max_division_age(p) = 10000.0 # tₘₐₓ
death_rate(p) = 0.001 # psick
mutation_probability(p) = 0.3 # pₘᵤₜ
min_area(p) = 1e-2 # Aₘᵢₙ

spring_constant(model, i::Int, j::Int, t) = spring_constant(model, model[i], model[j], t)
function spring_constant(model, p, q, t)
    μ = spring_constant(p, q)
    return μ # no adhesion for the initial population
end

rest_length(model, i::Int, j::Int, t) = rest_length(model, model[i], model[j]..., t)
function rest_length(model, p, q, t)
    return mature_cell_spring_rest_length(p, q) * 0.2
end

function proliferation_rate(model, i::Int, t)
    return model[i].size^3 > 1.0 ? 3.0 : 0.0
end

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
function force(model, p::Cell, t)
    F = SVector(0.0, 0.0)
    for j in get_neighbours(model.triangulation, p.index)
        DT.is_ghost_vertex(j) && continue
        # Find the agent with index j
        F = F + force(model, p, tri2agent(model, j), t)
    end
    return F
end

tri2agent(model, p::Int) = model.triangulation.points[p]

velocity(model, p::Cell, t) = force(model, p, t) / drag_coefficient(p)

function update_velocities!(model, t)
    for p in allagents(model)
        p.vel = velocity(model, p, t)
    end
    return model
end
function update_positions!(model, t)
    update_velocities!(model, t)
    for p in allagents(model)
        move_agent!(p, model, model.dt)
    end
    return model
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
    daughter = SVector(rand()-0.5, rand()-0.5) * 2 * 0.2 * parent.size
    while norm(daughter) > 0.2 * parent.size
        daughter = SVector(rand()-0.5, rand()-0.5) * 2 * 0.2 * parent.size
    end
    daughter = parent.pos + daughter
    add_agent!(daughter, model; birth=t, gfp=parent.gfp, vel=SVector(0.0, 0.0), parent=parent, size=parent.size*0.5^(1/3), nuts=parent.nuts)
    return daughter
end
function proliferate_cells!(model, t)
    for p in collect(allagents(model))
        Gᵢ = proliferation_rate(model, p.id, t)
        u = rand()
        if u < Gᵢ * model.dt
            place_daughter_cell!(model, p.id, t)
            p.size *= 0.5^(1/3)
        end
    end
    return true
end

function cull_cells!(model, t)
    for p in allagents(model)
        elder = t - p.birth > max_age(p)
        sick = rand() < model.dt * death_rate(p)
        xmax, ymax = spacesize(model)
        x, y = p.pos
        distance_to_origin = norm(p.pos - SVector(xmax / 2, ymax / 2))
        outside = distance_to_origin > 0.5 * min(xmax, ymax)
        if (elder || sick || outside)
            remove_agent!(p.id, model)
            p.death = t
        end
    end
    return model
end

function update_nutrients!(model, t)
    # Clear previous growth
    for p in allagents(model)
        # Get the triangle of that voronoi cell
        num_sample_points = 10
        inds = [get_spatial_index(p.pos + 0.2 * 2 * p.size * SVector(rand() - 0.5, rand() - 0.5),
                    model.nutrients, model) for _ in 1:num_sample_points]
        # Points will be duplicated if they are closer to the centre of the cell
        p.nuts = mean(model.nutrients[inds])
        for ind in inds
            #model.nutrients[ind] -= p.nut_import_rate * model.dt/num_sample_points
            #model.nutrients[ind] = max(0, model.nutrients[ind])
        end
    end
    # Diffuse nutrients
    #tmp = copy(model.nutrients)
    #for x in 2:size(model.nutrients, 1)-1, y in 2:size(model.nutrients, 2)-1
        #tmp[x,y] = max(0, 0.99*model.nutrients[x,y] + 0.01*(sum(model.nutrients[x-1:x+1,y-1:y+1]) - model.nutrients[x,y]) / 8)
    #end
    #model.nutrients = tmp
    return model
end

function model_step!(model)
    stepn = abmtime(model)
    t = stepn * model.dt
    update_positions!(model, t)
    update_nutrients!(model, t)
    #cull_cells!(model, t)
    proliferate_cells!(model, t)
    model.triangulation = retriangulate(model.triangulation, collect(allagents(model)))
    set_indexes(model)
    return model
end

function set_indexes(model)
    for (i, p) in enumerate(model.triangulation.points)
        p.index = i
    end
end

function initialize_cell_model(;
    ninit=3,
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
            birth=0.0, gfp=rand(), vel=SVector(0.0, 0.0), nuts=1.0, size=1.0)
    end
    positions = [cell.pos for cell in cells]

    # Compute the triangulation and the tessellation
    triangulation = triangulate(cells)

    # Define the model parameters
    properties = Dict(
        :triangulation => triangulation,
        :dt => dt,
        :nutrients => ones(100,100),
    )

    # Define the space
    space = ContinuousSpace(sides; periodic=false)

    # Define the model
    model = StandardABM(Cell, space; model_step!, properties, container=Vector)

    # Add the agents
    for (id, pos) in pairs(positions)
        add_agent!(pos, model; birth=0.0, gfp=rand(), vel=SVector(0.0, 0.0), index=id, nuts=1.0, size=1.0)
    end

    return model
end
