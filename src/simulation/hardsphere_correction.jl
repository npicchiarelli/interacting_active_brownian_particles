using .Threads
using Distances
using LinearAlgebra
include("boundary_conditions.jl")

# Functions for the hard sphere corrections

# Buffers of the hard sphere correction, allocated once per simulation and reused at every step
struct HardSphereWorkspace
    Δxy::Matrix{Float64}      # displacement that separates each particle from the ones it overlaps
    ncells::Int               # cells per side of the cell list
    head::Matrix{Int}         # first particle of each cell (0 if the cell is empty)
    next::Vector{Int}         # next particle in the same cell (0 at the end of the list)
    cell::Matrix{Int}         # cell coordinates of each particle
    partners::Vector{Int}     # overlapping partners of one particle
end

function HardSphereWorkspace(Np::Int, R::Float64, L::Float64)
    # Cells at least 2R wide, so overlapping particles always sit in neighbouring cells; at most about 4Np cells
    ncells = min(floor(Int, L / (2R*(1+1e-6))), ceil(Int, 2sqrt(Np)))
    return HardSphereWorkspace(zeros(Np,2), ncells, zeros(Int,ncells,ncells), zeros(Int,Np), zeros(Int,Np,2), zeros(Int,Np))
end

function hardsphere_periodic!(xy::Array{Float64,2}, workspace::HardSphereWorkspace, R::Float64,L::Float64; tol::Float64=0.)
    superpositions = 1
    counter = 0
    while superpositions > 0
        superpositions = hardsphere_correction_periodic!(workspace, xy, R, L, tol=1e-3)
        if superpositions > 0
            xy .+= workspace.Δxy
            periodic_BC_array!(xy,L,size(xy,1))
        end
        counter += 1
        if counter >= 100
            println("$superpositions superpositions remaining after 100 cycles")
            break
        end
    end
    return nothing
end

# Minimum image convention for a displacement along one axis
minimum_image(Δ::Float64, L::Float64) = abs(Δ) > abs(abs(Δ)-L) ? -sign(Δ)*abs(abs(Δ)-L) : Δ

# Periodic displacement and distance between particles a and b
@inline function pair_displacement(xy::Array{Float64,2}, a::Int, b::Int, L::Float64)
    Δx = minimum_image(xy[a,1] - xy[b,1], L)
    Δy = minimum_image(xy[a,2] - xy[b,2], L)
    return Δx, Δy, sqrt(Δx^2 + Δy^2)
end

# Adds to Δxy the displacements that separate the overlapping particles a and b
@inline function separate!(Δxy::Array{Float64,2}, xy::Array{Float64,2}, a::Int, b::Int, R::Float64, L::Float64, tol::Float64)
    Δx, Δy, d = pair_displacement(xy, a, b, L)
    s = ((1+tol)*2R / d - 1) / 2
    Δxy[a,1] += Δx*s
    Δxy[a,2] += Δy*s
    Δxy[b,1] += -Δx*s
    Δxy[b,2] += -Δy*s
    return nothing
end

# Fills Δxy with the displacements that separate overlapping particles and returns the number of overlapping pairs.
# Pairs (a, b) are taken in increasing a and then b, so each particle accumulates its displacement in increasing
# partner index, as in the original matrix row sums.
function hardsphere_correction_periodic!(Δxy::Array{Float64,2}, xy::Array{Float64,2}, R::Float64, L::Float64; tol::Float64=1e-3)
    Np = size(xy,1)
    fill!(Δxy, 0.0)
    superpositions = 0
    @inbounds for a in 1:Np-1, b in a+1:Np
        if 0. < pair_displacement(xy, a, b, L)[3] < 2R
            superpositions += 1
            separate!(Δxy, xy, a, b, R, L, tol)
        end
    end
    return superpositions
end

# Same result as above, looking for overlaps only among particles in neighbouring cells
function hardsphere_correction_periodic!(workspace::HardSphereWorkspace, xy::Array{Float64,2}, R::Float64, L::Float64; tol::Float64=1e-3)
    # The cell list needs at least 3 cells per side and particles within one box length from the centre
    if workspace.ncells < 3 || !all(x -> abs(x) <= L, xy)
        return hardsphere_correction_periodic!(workspace.Δxy, xy, R, L; tol=tol)
    end
    (; Δxy, ncells, head, next, cell, partners) = workspace
    Np = size(xy,1)
    w = L / ncells
    fill!(Δxy, 0.0)
    fill!(head, 0)
    @inbounds for i in Np:-1:1
        cx = mod(floor(Int, (xy[i,1] + L/2) / w), ncells) + 1
        cy = mod(floor(Int, (xy[i,2] + L/2) / w), ncells) + 1
        cell[i,1], cell[i,2] = cx, cy
        next[i] = head[cx,cy]
        head[cx,cy] = i
    end
    superpositions = 0
    @inbounds for a in 1:Np
        k = 0
        for i in -1:1, j in -1:1
            b = head[mod1(cell[a,1]+i, ncells), mod1(cell[a,2]+j, ncells)]
            while b != 0
                if b > a && 0. < pair_displacement(xy, a, b, L)[3] < 2R
                    k += 1
                    partners[k] = b
                end
                b = next[b]
            end
        end
        # Partners in increasing index, the order of the full scan
        sort!(view(partners, 1:k))
        for p in 1:k
            separate!(Δxy, xy, a, partners[p], R, L, tol)
        end
        superpositions += k
    end
    return superpositions
end

function hardsphere_periodic(xy::Array{Float64,2}, R::Float64, L::Float64; tol::Float64=1e-3) # called in initABPE
    workspace = HardSphereWorkspace(size(xy,1), R, L)
    hardsphere_periodic!(xy, workspace, R,L; tol=tol)
    return xy, workspace
end