           
using Distributed, ClusterManagers
using Base.Threads
using Revise
using Pkg
using xyBnG
using PyCall,CSV,DataFrames, Mmap,Serialization,VCFTools 

# make sure to use the fuuuull fucking path
# ENV["PYTHON"] = "/mnt/users/odwa/.conda/envs/simulation/bin/python"
# Pkg.build("PyCall")
# PyCall.libpython
@info "  we're in Julia"
 
#addprocs_slurm(parse(Int, ENV["SLURM_NTASKS"]))

@info "Number of threads: $(Threads.nthreads())"


include("/mnt/users/odwa/IBD_OCS/xyBnG_phased/ALT_IBD_seg.jl")    

NE = ARGS[1]
# REMEMBER IS IT A ANON OR COW SIM!!!

@time IBD_seg(data=".",
    baseDir = "NE_$NE",
    testDir = "cattle_NE_$NE",
    NE = NE ,
    animal = "Anon",
    nchp = 100_000)

