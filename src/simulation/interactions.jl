#Functions to apply interactions in ABPs systems
using Distances
using GeometryBasics
using VoronoiCells
include("boundary_conditions.jl")
include("force_functions.jl")
include("initialization.jl")
include(joinpath("..","utilities","geo_toolbox.jl"))

#Function to calculate force vectors

function interactions_range(xy::Array{Float64, 2}, L::Float64, l::Float64, Np::Int, int_func::Function, int_params...)
    # Preallocate the result array for efficiency
    ΣFtot = zeros(Float64, Np, 2)
    
    @threads for i in axes(xy,1)
        # Compute distances and apply periodic boundary conditions
        xy_shifted = xy .- [xy[i,1] xy[i,2]]
        periodic_BC_array!(xy_shifted, L, Np)
        xy_nonzero = @view xy_shifted[1:end .!=i, :]

        dists = d2(xy_nonzero)
        inside = vec(dists .<= l)
        # If no particles are inside the interaction range, force remains zero
        !any(inside) && continue
        xy_inside = @view xy_nonzero[inside, :]
        dirs = zeros(Float64, sum(inside), 2)

        # Compute distances and interaction forces
        dists_nonzero = @view dists[inside]
        forces = int_func.(dists_nonzero, int_params...)
        
        # Compute normalized direction vectors
        dirs .= xy_inside ./ dists_nonzero

        # Sum forces for each direction and assign to ΣFtot
        ΣFtot[i, :] = -sum(forces .* dirs, dims=1)'
    end
    return ΣFtot
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

    rect = Rect(Point2(-L/2, -L/2), Point2(L/2, L/2))

    tess = voronoicells(xy_to_points(xy), rect)

    # Step 2: Get the periodic projection of the tessellation

    xy_periodic_projection, proj_inds = find_boundary_points(xy, tess.Cells)
    θ_pp = vcat(θ, θ[proj_inds])
    xy_periodic_projection = vcat(xy, xy_periodic_projection)

    # Get the rectangle of the periodic projection
    minpoint = Point2(minimum(xy_periodic_projection[:,1]), minimum(xy_periodic_projection[:,2]))
    maxpoint = Point2(maximum(xy_periodic_projection[:,1]), maximum(xy_periodic_projection[:,2]))
    rect_periodic = Rect(minpoint, maxpoint)
    tess_periodic = voronoicells(xy_to_points(xy_periodic_projection), rect_periodic)

    # Step 3: Get the adjacency from the periodic tessellation
    adjacency = get_adjacency_from_points(tess_periodic.Cells, Np)
    
    #Step 4: Compute the forces
    xy_pp_oc = xy_periodic_projection .+ oc_length* [cos.(θ_pp) sin.(θ_pp)]
    return get_force_from_adjacency(adjacency, xy_pp_oc, Np, int_func, int_params...)
end

#Function used to compute torques in aligning interactions
function force_torque(xy, θ, L, oc_length, int_func, int_params...; method::Symbol=:range, range::Real=0.0)
    if method == :voronoi
        forces = interaction_voronoi(xy, θ, oc_length, L, size(xy,1), int_func, int_params...)
    elseif method == :range
        xy_chgcen = xy .+ oc_length * [cos.(θ) sin.(θ)]
        forces = interactions_range(xy_chgcen, L, range, size(xy,1), int_func, int_params...)
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
