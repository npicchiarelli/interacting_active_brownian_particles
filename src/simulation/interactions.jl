#Functions to apply interactions in ABPs systems
using Distances
using GeometryBasics
using VoronoiCells
include("boundary_conditions.jl")
include("force_functions.jl")
include("initialization.jl")
include(joinpath("..","utilities","geo_toolbox.jl"))

#Function to calculate force vectors

# Buffers of the range forces for one chunk of particles, allocated once per simulation and reused at every step
struct RangeBuffers
    Δx::Vector{Float64}               # displacements and distances of the particles within range of one particle
    Δy::Vector{Float64}
    d::Vector{Float64}
    forces::Vector{Matrix{Float64}}   # forces[n]: n×2 force terms of a particle with n neighbours, created when first needed
    ΣF::Matrix{Float64}               # 1×2 sum of the force terms
end
RangeBuffers(Np::Int) = RangeBuffers(zeros(Np), zeros(Np), zeros(Np), Vector{Matrix{Float64}}(undef, Np), zeros(1,2))

# One set of buffers per thread
RangeWorkspace(Np::Int) = [RangeBuffers(Np) for _ in 1:min(nthreads(), Np)]

function interactions_range(xy::Array{Float64, 2}, L::Float64, l::Float64, Np::Int, int_func::Function, int_params...; workspace::Vector{RangeBuffers}=RangeWorkspace(Np))
    # Preallocate the result array for efficiency
    ΣFtot = zeros(Float64, Np, 2)
    nchunks = length(workspace)

    @threads for c in 1:nchunks
        range_forces!(ΣFtot, c:nchunks:Np, xy, L, l, int_func, int_params, workspace[c])
    end
    return ΣFtot
end

# Total force on each particle in `particles` from all the particles within range l, written in ΣFtot
function range_forces!(ΣFtot::Array{Float64, 2}, particles, xy::Array{Float64, 2}, L::Float64, l::Float64, int_func::Function, int_params::Tuple, buffers::RangeBuffers)
    Np = size(xy,1)
    # Displacements and distances of the particles inside the interaction range
    Δx_in, Δy_in, d_in = buffers.Δx, buffers.Δy, buffers.d

    for i in particles
        n = 0
        for j in 1:Np
            j == i && continue
            # Compute distances and apply periodic boundary conditions, as in periodic_BC_array!
            Δx = xy[j,1] - xy[i,1]
            Δy = xy[j,2] - xy[i,2]
            if abs(Δx) > L/2
                Δx = Δx - sign(Δx) * L
            end
            if abs(Δy) > L/2
                Δy = Δy - sign(Δy) * L
            end
            d = sqrt(abs2(Δx) + abs2(Δy))
            if d <= l
                n += 1
                Δx_in[n], Δy_in[n], d_in[n] = Δx, Δy, d
            end
        end
        # If no particles are inside the interaction range, force remains zero
        n == 0 && continue

        isassigned(buffers.forces, n) || (buffers.forces[n] = Matrix{Float64}(undef, n, 2))
        forces = buffers.forces[n]
        for k in 1:n
            f = int_func(d_in[k], int_params...)
            forces[k,1] = f * (Δx_in[k] / d_in[k])
            forces[k,2] = f * (Δy_in[k] / d_in[k])
        end
        # Summing an n×2 Matrix over dims=1 adds the terms in the same order as the original sum(forces .* dirs, dims=1)
        sum!(buffers.ΣF, forces)
        ΣFtot[i,1] = -buffers.ΣF[1]
        ΣFtot[i,2] = -buffers.ΣF[2]
    end
    return nothing
end

#Functions to calculate force vectors with voronoi rationale.
function get_force_from_adjacency(adjacency::Set{Tuple{Int64, Int64}}, xy::Array{Float64, 2}, Np::Int, int_func::Function, int_params...)
    ΣFtot = zeros(Float64, Np, 2)
    for a in adjacency
        # Since tuples are ordered and the set is composed starting from points inside the simulation box, we are sure that a[1] <= Np
        dxy = xy[a[1],:] - xy[a[2],:]
        dist = sqrt(sum(abs2, dxy))
        dir = dxy / dist
        force = int_func(dist, int_params...)
        ΣFtot[a[1],:] .+= force * dir
        # add the opposite force to the second point if it is within the simulation box
        if a[2]<=Np
            ΣFtot[a[2],:] .-= force * dir
        end
    end
    return ΣFtot
end

# Periodic images of the particles that fall less than m outside the box [-L/2, L/2]², corners included.
# Returns their positions and, for each image, its particle and shift in box lengths.
function periodic_images(xy::Array{Float64, 2}, L::Real, m::Real)
    n = ceil(Int, m/L)   # box lengths spanned by the margin
    W = L/2 + m
    xs, ys, images = Float64[], Float64[], NTuple{3,Int}[]
    for k in axes(xy, 1), sx in -n:n, sy in -n:n
        (sx == 0 && sy == 0) && continue
        x = xy[k,1] + sx*L
        y = xy[k,2] + sy*L
        if abs(x) < W && abs(y) < W
            push!(xs, x); push!(ys, y); push!(images, (k, sx, sy))
        end
    end
    return xs, ys, images
end

# Voronoi tessellation of the particles in the periodic box, with the images needed for the cells of the particles
# to be exactly the periodic ones. A point cuts the cell of particle i only if it lies within |v - xᵢ| of one of its
# vertices v, so the images are added until none is missing from these disks. Starts from the images within `margin`
# of the box. Returns the tessellation, the positions of all the generators (the Np particles first) and the
# particle of each generator.
function periodic_voronoi(xy::Array{Float64, 2}, L::Real; margin::Real = 5L/sqrt(size(xy,1)))
    Np = size(xy, 1)
    m = min(margin, 1.5L)   # wider margins are never needed and would not fit in the rectangle
    W0 = L/2 + m        # all the images inside this square are present
    Wr = L/2 + 2L       # tessellation rectangle; a cell bounded by the images of its own particle stays within L/2 + 1.71L
    rect = Rectangle(Point2(-Wr, -Wr), Point2(Wr, Wr))
    xs, ys, images = periodic_images(xy, L, m)
    xs = vcat(xy[:,1], xs); ys = vcat(xy[:,2], ys)
    sources = vcat(1:Np, first.(images))
    added = Set{NTuple{3,Int}}()
    added_now = NTuple{3,Int}[]   # images added since the last tessellation
    present(k, sx, sy) = (sx == 0 && sy == 0) || (abs(xy[k,1] + sx*L) < W0 && abs(xy[k,2] + sy*L) < W0) || (k, sx, sy) in added
    function add!(k, sx, sy)
        push!(xs, xy[k,1] + sx*L); push!(ys, xy[k,2] + sy*L); push!(sources, k); push!(added, (k, sx, sy)); push!(added_now, (k, sx, sy))
    end
    while true
        tess = voronoicells(xy_to_points([xs ys]), rect)
        empty!(added_now)
        for i in 1:Np
            cell = tess.Cells[i]
            # A cell clipped by the rectangle, or with disks leaving it, is bounded by the images of its own particle
            if any(v -> max(abs(v[1]), abs(v[2])) + sqrt((v[1]-xy[i,1])^2 + (v[2]-xy[i,2])^2) >= Wr, cell)
                for sx in -1:1, sy in -1:1
                    present(i, sx, sy) || add!(i, sx, sy)
                end
                # with its own images in the tessellation the cell lies within L/2 of the particle, so this cannot happen
                any(k -> k[1] == i, added_now) || error("periodic_voronoi: cell $i is not bounded by the images of its particle")
                continue
            end
            for v in cell
                r2 = (v[1]-xy[i,1])^2 + (v[2]-xy[i,2])^2
                r = sqrt(r2)
                # disks inside W0 hold no missing image
                max(abs(v[1]), abs(v[2])) + r < W0 - 1e-9*L && continue
                for k in 1:Np
                    for sx in ceil(Int, (v[1] - r - xy[k,1])/L):floor(Int, (v[1] + r - xy[k,1])/L),
                        sy in ceil(Int, (v[2] - r - xy[k,2])/L):floor(Int, (v[2] + r - xy[k,2])/L)
                        present(k, sx, sy) && continue
                        # slightly widened disk: an extra image does not change the cell, a missing one would
                        if (xy[k,1] + sx*L - v[1])^2 + (xy[k,2] + sy*L - v[2])^2 < r2*(1 + 1e-9)
                            add!(k, sx, sy)
                        end
                    end
                end
            end
        end
        isempty(added_now) && return tess, [xs ys], sources
    end
end

function interaction_voronoi(xy::Array{Float64, 2}, θ::Array{Float64,1}, oc_length::Float64, L::Real, Np::Int, int_func::Function, int_params...)
    # Step 1: Get the periodic tessellation
    tess, xy_pp, sources = periodic_voronoi(xy, L)
    θ_pp = θ[sources]

    # Step 2: Get the adjacency from the periodic tessellation
    adjacency = get_adjacency_from_points(tess.Cells, Np)

    # Step 3: Compute the forces
    xy_pp_oc = xy_pp .+ oc_length* [cos.(θ_pp) sin.(θ_pp)]
    return get_force_from_adjacency(adjacency, xy_pp_oc, Np, int_func, int_params...)
end

#Function used to compute torques in aligning interactions
function force_torque(xy, θ, L, oc_length, int_func, int_params...; method::Symbol=:range, range::Real=0.0, workspace::Vector{RangeBuffers}=RangeWorkspace(size(xy,1)))
    if method == :voronoi
        forces = interaction_voronoi(xy, θ, oc_length, L, size(xy,1), int_func, int_params...)
    elseif method == :range
        xy_chgcen = xy .+ oc_length * [cos.(θ) sin.(θ)]
        forces = interactions_range(xy_chgcen, L, range, size(xy,1), int_func, int_params...; workspace)
    else
        throw(ArgumentError("Unknown method: $method. Use :range or :voronoi."))
    end
    torques = oc_length * (forces[:,2] .* cos.(θ) .- forces[:,1] .* sin.(θ))
    return forces, torques
end

# Correction for the first step when offcenter interactions are present, and a diverging force field is used (e.g. WCA). This function randomly reorients particles until no superposition is present, which is the only way to avoid divergences in the first step.
function offcenter_nosuperpose!(abpe::ABPE2, δt::Float64, offcenter::Float64, range::Float64, int_func::Function, int_params...)
    xy = position(abpe)
    θ = orientation(abpe)
    R = abpe.R
    L = abpe.L
    γₜ = diffusion_coeff(1e-6*R, T)[3] #Output in international system units kg/s
    γᵣ = (8e-12γₜ / 6) * R^2                #Output in international system units
    oc = offcenter*R

    xy_chgcen = xy .+ oc* [cos.(θ) sin.(θ)]

    ocdists = pairwise(Euclidean(), xy_chgcen, xy_chgcen, dims=1)
    superpositions = sum(0.0 .< ocdists .< 2R)/2
    j = 0
    while superpositions > 0

        ocdists = pairwise(Euclidean(), xy_chgcen, xy_chgcen, dims=1)
        sup_per_particle = sum(0.0 .< ocdists .< 2R, dims=1)
        id = vec(sup_per_particle .> 0)
        abpe.θ[id] .+= rand(abpe.θ[id.>0]).*2π
        xy_chgcen = xy .+ oc* [cos.(θ) sin.(θ)]
        superpositions = sum(0.0 .< pairwise(Euclidean(), xy_chgcen, xy_chgcen, dims=1) .< 2R)/2
        j +=1
    end
    return nothing
end
