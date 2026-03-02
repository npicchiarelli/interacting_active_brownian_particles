using ArgParse
using CSV
using DataFrames
using Dates
using DelimitedFiles
using Distributions
using Logging
using Printf
using Random

# Simulation files
include("simulation/boundary_conditions.jl")
include("simulation/force_functions.jl")
include("simulation/hardsphere_correction.jl")
include("simulation/initialization.jl")
include("simulation/interactions.jl")
include("simulation/update_step.jl")

# Utilities
include("utilities/geo_toolbox.jl")

function parse_cl()
    s = ArgParseSettings()
    @add_arg_table s begin
        "simfolder"
            help = "the path to the simulation folder"
            required = true
    end
    return parse_args(s)
end

parsed_args = parse_cl()
simulation_folder = parsed_args["simfolder"]
settings_file = joinpath(simulation_folder, "settings.jl")

if isfile(settings_file)
    include(settings_file)
else
    error("Settings file not found at $settings_file")
end
println("Loaded settings from $settings_file")

# Preliminary calculations

DT, DR, γ = diffusion_coeff(1e-6R,T).*[1e12, 1, 1] # Translational and Rotational diffusion coefficients, drag coefficient
T_tot = δt*Nt
actual_steps = (Nt÷measevery)+1

if box_shape == :square
    packing_fraction = (pi*R^2/L^2) # Largest pf for spherical beads π/4 = 0.7853981633974483
    density = Np/L^2
    a,b = L,L
    packing_fraction > π/4 && @warn "Packing fraction $packing_fraction exceeds the maximum for spherical particles (π/4 ≈ 0.785). Consider reducing the number of particles or increasing the box size to avoid unphysical overlaps."
end

# Info printing on shell and file
datestamp = Dates.format(now(), "YYYYmmdd-HHMMSS")

infos = @sprintf "%s\nBox shape: %s\nNumber of particles = %i\nNumber density = %s μm⁻²\nR=%.1f μm \nT = %.1f (K)\nv=%s (μm/s) \nω=%s (rad/s)\nCharacteristic lengths: (a=%.1f b=%.1f) μm\npf=%s\nIntegration step: dt=%.0e s \nSimulation downsampling: %i\nNumber of steps: Nt=%.1e\nTotal simulated time T_tot = %.2e s\n\nInteraction function: %s with parameters: %s\nRange: %.1f μm\nOffcenter: %s" datestamp box_shape Np density R T v ω a b packing_fraction δt measevery  Nt T_tot int_func int_params intrange offcenter

println(infos)

for i in 1:ICS
    pathf= joinpath(simulation_folder, "run$i")
        if !isdir(pathf)
            mkdir(pathf)
        else
            @warn "Directory $pathf already exists. Contents will be overwritten."
            rm(pathf, recursive=true, force = true)
            mkdir(pathf)
        end
        filename= "simulation"*"_run$i"
        filepath = joinpath(pathf, filename)

        datafname = filepath*".txt"

        # Simulation and file storage
        open(datafname, "w") do infile
            writedlm(infile, ["N" "Time" "xpos" "ypos" "orientation" "omega"], ",")
        end

        start = now()
        @info "$(start) Started simulation #$i"


    ABPE, matrices = initABPE( Np, L, R, T, v, ω, int_func, offcenter, intrange, int_params...,)
    if int_func == lennard_jones
        offcenter_nosuperpose!(ABPE, δt, offcenter, 2R+1e-3, int_func, int_params...)
    end
    for nt in 0:Nt
        if nt % measevery == 0
            pnumber = collect(1:Np)
            time = fill(nt, Np)
            #creating Data
            data = DataFrame(
                N= pnumber,
                Time= time,
                xpos= ABPE.x,
                ypos= ABPE.y,
                orientation=ABPE.θ,
                omega=ABPE.ω,
            )  
            CSV.write(datafname, data, append = true)
        end
        ABPE =update_heun(ABPE,matrices,δt, offcenter, intrange, int_func , int_params...)
        if nt % (Nt÷100) == 0
            elapsed = Dates.canonicalize(Dates.round((now()-start), Dates.Second))
            print("$((100*nt÷Nt))%... Step $nt, total elapsed time $(elapsed)\r")
        end
    end
    @info "$(now()) Simulation and file writing finished"
end
