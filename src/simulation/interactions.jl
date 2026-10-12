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

function interaction_voronoi(xy::Array{Float64, 2}, θ::Array{Float64,1}, oc_length::Float64, L::Real, Np::Int, int_func::Function, int_params...)
    #Step 1: Get the tessellation

    rect = Rectangle(Point2(-L/2, -L/2), Point2(L/2, L/2))

    tess = voronoicells(xy_to_points(xy), rect)

    # Step 2: Get the periodic projection of the tessellation

    xy_periodic_projection, proj_inds = find_boundary_points(xy, tess.Cells)
    θ_pp = vcat(θ, θ[proj_inds])
    xy_periodic_projection = vcat(xy, xy_periodic_projection)

    # Get the rectangle of the periodic projection
    minpoint = Point2(minimum(xy_periodic_projection[:,1]), minimum(xy_periodic_projection[:,2]))
    maxpoint = Point2(maximum(xy_periodic_projection[:,1]), maximum(xy_periodic_projection[:,2]))
    rect_periodic = Rectangle(minpoint, maxpoint)
    tess_periodic = voronoicells(xy_to_points(xy_periodic_projection), rect_periodic)

    # Step 3: Get the adjacency from the periodic tessellation
    adjacency = get_adjacency_from_points(tess_periodic.Cells, Np)
    
    #Step 4: Compute the forces
    xy_pp_oc = xy_periodic_projection .+ oc_length* [cos.(θ_pp) sin.(θ_pp)]
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
