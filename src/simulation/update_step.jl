using .Threads
include("boundary_conditions.jl")
include("hardsphere_correction.jl")
include("initialization.jl")
include("interactions.jl")

# Functions to update particles for the next step

function update_heun(abpe::ABPE, matrices::Tuple{Matrix{Float64}, BitMatrix, BitMatrix}, δt::Float64, offcenter::Float64, range::Float64, int_func::Function, int_params...) where {ABPE <: ABPsEnsemble}

    δp_i = Array{Float64,2}(undef,abpe.Np,2)
    δθ_i = Array{Float64,1}(undef,abpe.Np)
    δp_f = Array{Float64,2}(undef,abpe.Np,2)
    δθ_f = Array{Float64,1}(undef,abpe.Np)
    γₜ = diffusion_coeff(1e-6*abpe.R, abpe.T)[3] #Output in international system units kg/s
    γᵣ = (8e-12γₜ / 6) * abpe.R^2                #Output in international system units
    oc_length = offcenter*abpe.R
    
    #intermediate step
    if (!isapprox(offcenter,0.0))
        f_i, t_i = force_torque(position(abpe), orientation(abpe), abpe.L, oc_length, range, int_func, int_params...)
    else
        f_i = interactions_range(position(abpe), abpe.L, range, abpe.Np, int_func, int_params...)
        t_i = zeros(abpe.Np)
    end 

    det_part_i = (abpe.v.*[cos.(abpe.θ) sin.(abpe.θ)] .+ 1e-6f_i/γₜ, abpe.ω .+ 1e-18t_i/γᵣ)
    noise_part =  (sqrt.(2*δt*abpe.DT)*randn(abpe.Np,2), sqrt(2*abpe.DR*δt)*randn(abpe.Np))
    δp_i .= noise_part[1] .+ δt.*det_part_i[1]
    δθ_i .= noise_part[2] .+ δt.*det_part_i[2]

    pθ_i = (position(abpe), orientation(abpe)) .+ (δp_i, δθ_i)

    periodic_BC_array!(pθ_i[1], abpe.L, abpe.Np)

    #final step
    if (!isapprox(offcenter,0.0))
        f_f, t_f = force_torque(pθ_i..., abpe.L, oc_length, range, int_func, int_params...)
    else
        f_f = interactions_range(pθ_i[1], abpe.L, range, abpe.Np, int_func, int_params...)
        t_f = zeros(abpe.Np)
    end 

    det_part_f = (abpe.v.*[cos.(pθ_i[2]) sin.(pθ_i[2])] .+ 1e-6f_f/γₜ, abpe.ω .+ 1e-18t_f/γᵣ)
    δp_f .= noise_part[1] .+ δt.*(det_part_f[1] .+ det_part_i[1])/2
    δθ_f .= noise_part[2] .+ δt.*(det_part_f[2] + det_part_i[2])/2

    pθ = (position(abpe), orientation(abpe)) .+ (δp_f, δθ_f)

    periodic_BC_array!(pθ[1],abpe.L, abpe.Np)
    hardsphere_periodic!(pθ[1], matrices[1], matrices[2],abpe.R, abpe.L)

    new_abpe = ABPE2( abpe.Np, abpe.L, abpe.R, abpe.T, abpe.v, abpe.ω, abpe.DT, abpe.DR, pθ[1][:,1], pθ[1][:,2], pθ[2] )

    return new_abpe
end