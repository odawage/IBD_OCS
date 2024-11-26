
using Distributed, ClusterManagers
using Base.Threads
using Revise
using Pkg
using xyBnG

# make sure to use the fuuuull fucking path
# ENV["PYTHON"] = "/mnt/users/odwa/.conda/envs/simulation/bin/python"
# Pkg.build("PyCall")
# PyCall.libpython
@info "  we're in Julia"
 
#addprocs_slurm(parse(Int, ENV["SLURM_NTASKS"]))

@info "Number of threads: $(Threads.nthreads())"


using PyCall,CSV,DataFrames, Mmap,Serialization,VCFTools,xyBnG 
include("/mnt/users/odwa/paper-2/xyBnG_phased/IBD_seg.jl")    

#lock = ReentrantLock()
#include("/mnt/users/odwa/paper-2/xyBnG_phased/phased_matrix.jl")

#@time IBD_seg(data="/mnt/SCRATCH/odwa/xyBnG_sim/rst")

@time IBD_seg(data=".")


# np = parse(Int, ENV["SLURM_NTASKS"])
# addprocs(SlurmManager(np); exeflags="--project")

# #results = Vector{Int}(undef, 10)

# @distributed for i in 1:5
#     filename = "$(i).txt"  # Create filename based on the iteration index
#     open(filename, "w") do file
#         write(file, "test")
#     end
# end 

#df = DataFrame(Index=1:10, Value=results)
