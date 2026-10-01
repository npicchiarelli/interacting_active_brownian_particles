using Distributed

total_cores = 16 
threads_per_sim = 2 
num_workers = total_cores ÷ threads_per_sim
addprocs(num_workers)

const PKGDIR   = @__DIR__
const JULIABIN = joinpath(Sys.BINDIR, Base.julia_exename())
const SCRIPT   = joinpath(PKGDIR, "src", "run_simulation.jl")

@everywhere function run_one(folder, threads, juliabin, pkgdir, script)
    logfile = joinpath(folder, "output.log")
    cmd = Cmd(`$juliabin --project=$pkgdir --threads $threads $script $folder`, dir = pkgdir)
    run(pipeline(cmd,
        stdout = open(logfile, "w"),
        stderr = open(logfile, "w")))
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

folders = [joinpath(PKGDIR, "..", "simulations", run, name)
           for (run, names) in unfinished for name in names]
for d in folders
    @assert isfile(joinpath(d, "settings.jl")) "missing settings.jl in $d"
end

# println(folders, length(folders))
pmap((folder) -> run_one(folder, threads_per_sim, JULIABIN, PKGDIR, SCRIPT), folders)

