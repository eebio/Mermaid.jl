using Agents, StaticArrays
using Random
using LinearAlgebra
using StatsBase
using DelaunayTriangulation
const DT = DelaunayTriangulation

@agent struct Cell(ContinuousAgent{2,Float64})
    gfp::Float64
    index::Int = typemin(Int)
    const parent::Union{Cell, Nothing} = nothing
    nut_import_rate::Float64 = 0.0
    nuts::Float64
    size::Float64
end

mutable struct ModelProperties{X <: DT.Triangulation}
    triangulation::X
    dt::Float64
    nutrients::Matrix{Float64}
    nutrient_import_rate::Matrix{Float64}
end

DT.getx(cell::Cell) = cell.pos[1]
DT.gety(cell::Cell) = cell.pos[2]
DT.number_type(::Type{Cell}) = Float64
DT.number_type(::Type{Vector{Cell}}) = Float64
DT.is_point2(::Cell) = true

spring_constant() = 20.0 # μ
drag_coefficient() = 1.0 # η
cutoff_distance(p, q) = 1.5 # ℓₘₐₓ
rest_length(p, q) = (p.size + q.size) * 0.2 # ℓ₀
proliferation_rate(model, i::Int) = model[i].size^3 > 1.0 ? 3.0 : 0.0

function force(p::Cell, q::Cell)
    δ = norm(p.pos - q.pos)
    if δ > cutoff_distance(p, q)
        return SVector(0.0, 0.0)
    end
    μ = spring_constant()
    s = rest_length(p, q)
    rᵢⱼ = q.pos - p.pos
    return μ * (norm(rᵢⱼ) - s) * rᵢⱼ / norm(rᵢⱼ)
end
function force(model, p::Cell)
    F = SizedVector(0.0, 0.0)
    for j in get_neighbours(model.triangulation, p.index)
        DT.is_ghost_vertex(j) && continue
        # Find the agent with index j
        F .+= force(p, tri2agent(model.triangulation, j))
    end
    return F
end

tri2agent(triangulation, p::Int) = triangulation.points[p]

velocity(model, p::Cell) = force(model, p) / drag_coefficient()

function update_velocities!(model)
    for p in allagents(model)
        p.vel = velocity(model, p)
    end
    return model
end
function update_positions!(model)
    update_velocities!(model)
    for p in allagents(model)
        move_agent!(p, model, model.dt)
    end
    return model
end

function place_daughter_cell!(model, i)
    parent = model[i]
    daughter = SVector(rand()-0.5, rand()-0.5) * 2 * 0.2 * parent.size
    while norm(daughter) > 0.2 * parent.size
        daughter = SVector(rand()-0.5, rand()-0.5) * 2 * 0.2 * parent.size
    end
    daughter = parent.pos + daughter
    add_agent!(daughter, model; gfp=parent.gfp, vel=SVector(0.0, 0.0), parent=parent, size=parent.size*0.5^(1/3), nuts=parent.nuts)
    return daughter
end
function proliferate_cells!(model)
    for p in collect(allagents(model))
        Gᵢ = proliferation_rate(model, p.id)
        u = rand()
        if u < Gᵢ * model.dt
            place_daughter_cell!(model, p.id)
            p.size *= 0.5^(1/3)
        end
    end
    return true
end

function get_spatial_index_type_stable(pos, property::AbstractArray{T, D}, model::ABM) where {T, D}
    ssize = spacesize(model)
    propertysize = size(property)
    εs = ssize ./ propertysize
    idxs = floor.(Int, pos ./ εs) .+ 1
    return CartesianIndex(Tuple{Vararg{Int, D}}(idxs))
    #return CartesianIndex(Tuple(idxs))
end

function update_nutrients!(model)
    # Clear previous growth
    model.nutrient_import_rate .= 0.0
    tmp_pos = MVector{2, Float64}(0.0, 0.0)
    for p in allagents(model)
        # Get the triangle of that voronoi cell
        num_sample_points = 10
        inds = Vector{CartesianIndex{2}}(undef, num_sample_points)
        for i in 1:num_sample_points
            # Generate random point inside the cell's bounding box
            tmp_pos .= p.pos .+ 0.2 * 2 * p.size .* (rand(MVector{2}) .- 0.5)
            # Clamp manually
            tmp_pos[1] = clamp(tmp_pos[1], 0.0, abmspace(model).extent[1])
            tmp_pos[2] = clamp(tmp_pos[2], 0.0, abmspace(model).extent[2])
            # Get index
            inds[i] = get_spatial_index_type_stable(tmp_pos, model.nutrients, model)
        end
        # Points will be duplicated if they are closer to the centre of the cell
        p.nuts = mean(model.nutrients[inds])
        for ind in inds
            model.nutrient_import_rate[ind] += p.nut_import_rate * model.dt/num_sample_points
        end
    end
    return model
end

function model_step!(model)
    update_positions!(model)
    update_nutrients!(model)
    proliferate_cells!(model)
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
    ninit=5,
    radius=0.5,
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
            gfp=rand(), vel=SVector(0.0, 0.0), nuts=1.0, size=1.0)
    end
    positions = [cell.pos for cell in cells]

    # Compute the triangulation and the tessellation
    triangulation = triangulate(cells)

    # Define the model parameters
    properties = ModelProperties(
        triangulation,
        dt,
        ones(100,100),
        zeros(100,100),
    )

    # Define the space
    space = ContinuousSpace(sides; periodic=false)

    # Define the model
    model = StandardABM(Cell, space; model_step!, properties, container=Vector)

    # Add the agents
    for (id, pos) in pairs(positions)
        add_agent!(pos, model; gfp=rand(), vel=SVector(0.0, 0.0), index=id, nuts=1.0, size=1.0)
    end

    return model
end
