using Distributed

total_cores = 16 
threads_per_sim = 2 
num_workers = total_cores ÷ threads_per_sim
addprocs(num_workers)

const PKGDIR   = @__DIR__
const JULIABIN = joinpath(Sys.BINDIR, Base.julia_exename())
const SCRIPT   = joinpath(PKGDIR, "src", "run_simulation.jl")

# Runs one initial condition of one simulation folder in its own process
@everywhere function run_one(folder, ic, threads, juliabin, pkgdir, script)
    # The log stays outside run$ic, which run_simulation.jl empties when it starts
    logfile = joinpath(folder, "output_run$ic.log")
    cmd = Cmd(`$juliabin --project=$pkgdir --threads $threads $script $folder --ic $ic`, dir = pkgdir)
    open(logfile, "w") do log
        run(pipeline(cmd, stdout = log, stderr = log))
    end
    return nothing
end

# Redo only the run02 folders left unfinished (per check_progress.sh).
unfinished = Dict(
    "run02_voronoi" => [
        "scaled_attractive_pf=1e-01_offcenter=0e+00",
        "scaled_attractive_pf=2e-01_offcenter=0e+00",
        "scaled_attractive_pf=2e-01_offcenter=-5e-01",
        "scaled_attractive_pf=2e-01_offcenter=5e-01",
        "scaled_repulsive_pf=2e-01_offcenter=-5e-01",
    ],
    "run02_range" => [
        "scaled_attractive_pf=2e-01_offcenter=0e+00",
        "scaled_repulsive_pf=2e-01_offcenter=-5e-01",
    ],
)

# Initial conditions to run in some of the folders above, e.g. "run02_range/scaled_repulsive_pf=2e-01_offcenter=-5e-01" => [3, 5].
# The folders not listed here run all the ICS of their settings.jl.
only_ics = Dict{String, Vector{Int}}(
)

# Number of initial conditions in a settings.jl
read_ics(folder) = parse(Int, match(r"^\s*ICS\s*=\s*(\d+)"m, read(joinpath(folder, "settings.jl"), String))[1])

# One job per initial condition of each folder
jobs = Tuple{String, Int}[]
for (rundir, names) in unfinished, name in names
    folder = joinpath(PKGDIR, "..", "simulations", rundir, name)
    @assert isfile(joinpath(folder, "settings.jl")) "missing settings.jl in $folder"
    for ic in get(only_ics, joinpath(rundir, name), 1:read_ics(folder))
        push!(jobs, (folder, ic))
    end
end

# A failed job is reported at the end without stopping the others
results = pmap(jobs; on_error = e -> sprint(showerror, e)) do (folder, ic)
    run_one(folder, ic, threads_per_sim, JULIABIN, PKGDIR, SCRIPT)
end
failed = [job for (job, result) in zip(jobs, results) if result isa String]
isempty(failed) || @warn "$(length(failed)) of $(length(jobs)) jobs failed, see their output_run<ic>.log" failed

