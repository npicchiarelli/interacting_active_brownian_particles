using  Distributions, Random

include("hardsphere_correction.jl")
# Define an "ABPsEnsemble" Type
abstract type ABPsEnsemble end

# Define a specific type for 2D ABPsEnsembles (CURRENTLY ASSUMING ALL PARTICLES ARE EQUAL)
struct ABPE2 <: ABPsEnsemble
    Np::Int64                      # number of particles®®
    L::Float64                      # size of observation space (μm)
	R::Float64  # Radius (μm)                                   --> Vector{Float64}(undef,Np)
    T::Float64 #temperature (K)
	v::Vector{Float64}  	# velocity (μm/s)                   --> Vector{Float64}(undef,Np)
    ω::Vector{Float64} #Angular velocity (rad/s)                --> Vector{Float64}(undef,Np)            
	DT::Float64 # translational diffusion coefficient (μm^2/s)  --> Vector{Float64}(undef,Np)
	DR::Float64 # rotational diffusion coefficient (rad^2/s)    --> Vector{Float64}(undef,Np)
	x::Vector{Float64}    # x position (μm)
	y::Vector{Float64}    # y position (μm)
	θ::Vector{Float64}    # orientation (rad)
end

## Initialize ABP ensemble (CURRENTLY ONLY 2D) 
function initABPE(Np::Int64, L::Float64, R::Float64, T::Float64, vd::Union{Float64,Array{Float64,1},Distribution}, ωd::Union{Float64,Array{Float64,1},Distribution}, int_func::Function, offcenter::Float64, range::Float64, int_params...; η::Float64=1e-3)
    # translational diffusion coefficient [m^2/s] & rotational diffusion coefficient [rad^2/s] - R [m]
    # Intial condition will be choosen as per the geometry under study
    (vd isa Float64) ? vd = [vd] : Nt
    (ωd isa Float64) ? ωd = [ωd] : Nt
    DT, DR = diffusion_coeff(1e-6R, T)
    xyθ = (rand(Np,3).-0.5).*repeat([L L 2π],Np)
    xyθ[:,1:2], dists, superpose, uptriang = hardsphere_periodic(xyθ[:,1:2], R, L) #xyθ[:,1:2] gives x and y positions of intitial particles
    v = rand(vd, Np)
    ω = rand(ωd,Np)
    abpe = ABPE2( Np, L, R, T, v, ω, 1e12DT, DR, xyθ[:,1], xyθ[:,2], xyθ[:,3])
    return abpe, (dists, superpose, uptriang)
end

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Calculate diffusion coefficient and friction coefficient
function diffusion_coeff(R::Float64, T::Float64=300.0, η::Float64=1e-3)
    # Boltzmann constant [J/K]
    kB = 1.38e-23
    # friction coefficient [Ns/m]
    γ = 6*pi*R*η
    # translational diffusion coefficient [m^2/s]
    DT = kB*T/γ
    # rotational diffusion coefficient [rad^2/s]
    DR = 6*DT/(8*R^2)
    return DT, DR, γ
end
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

position(abpe::ABPE2) = [ abpe.x abpe.y ]
orientation(abpe::ABPE2) = abpe.θ
force(abpe::ABPE2) = [abpe.fx abpe.fy]
torque(abpe::ABPE2) = abpe.torque

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

