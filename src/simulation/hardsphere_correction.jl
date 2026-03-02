using .Threads
using Distances
using LinearAlgebra
include("boundary_conditions.jl")

# Functions for the hard sphere corrections

function hardsphere_periodic!(xy::Array{Float64,2}, periodicdists::Array{Float64,2}, superpose::BitArray{2}, R::Float64,L::Float64; tol::Float64=0.)
    superpositions = 1
    counter = 0
    Δxy = zeros(size(periodicdists)...,2)
    while superpositions > 0
        Threads.@threads for i in 1:2
            Δxy[:,:,i] .= pairwise(-,xy[:,i],xy[:,i])
            ix = abs.(Δxy[:,:,i]) .> abs.(abs.(Δxy[:,:,i]).-L)
            Δxy[ix,i] .= -sign.(Δxy[ix,i]).*abs.(abs.(Δxy[ix,i]).-L)
        end
        periodicdists .= sqrt.(Δxy[:,:,1].^2 .+ Δxy[:,:,2].^2)
        superpose .= 0. .< periodicdists.<2R
        superpositions = sum(superpose)÷2
        if superpositions > 0
            hardsphere_correction_periodic!(xy,Δxy,periodicdists,superpose,R,tol=1e-3)
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

function hardsphere_correction_periodic!(xy::Array{Float64,2}, Δxy::Array{Float64,3}, dists::Array{Float64,2}, superpose::BitArray{2}, R::Float64; tol::Float64=1e-3)
    Δpi(Δxi,d,s) = s==0 ? 0.0 : Δxi*(((1+tol)*2R / d - 1) / 2 )
    Threads.@threads for i in 1:2
        xy[:,i] .+= sum(Δpi.(Δxy[:,:,i],dists,superpose),dims=2)
    end
    return nothing
end

function hardsphere_periodic(xy::Array{Float64,2}, R::Float64, L::Float64; tol::Float64=1e-3) # called in initABPE
    Np = size(xy,1)
    periodicdists = zeros(Np,Np)
    superpose = falses(Np,Np)
    uptriang = triu(trues(Np,Np),1)
    hardsphere_periodic!(xy, periodicdists, superpose, R,L; tol=tol)
    return xy, periodicdists, superpose, uptriang
end