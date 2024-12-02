           
using Distributed, ClusterManagers
using Base.Threads
using Revise
using Pkg
using xyBnG


@info "Number of threads: $(Threads.nthreads())"


using PyCall,CSV,DataFrames, Mmap,Serialization,VCFTools,xyBnG 
include("/mnt/users/odwa/IBD_OCS/xyBnG_phased/IBD_seg.jl")    

using DataFrames
using Distributed
using LinearAlgebra
using Markdown
using Mmap
using Random
using Base.Threads
using Serialization
using Statistics
using xyBnG.XY
using xyBnG.Conn
using xyBnG.Sum
import xyBnG.xps: agocs,ggocs, iiocs, igocs, iiocs,initPop, randbrd 
import xyBnG.xyTypes: Cattle, Species, Trait, Plan, name
import xyBnG.Breeding: phenotype!, Predict!, Select, reproduce!, TM1997, TM2024
import xyBnG.Founder: ts_base, uniq, sample_xy, macs_base
import xyBnG.RS: nrm, irm, grm, xirm

include("/mnt/users/odwa/IBD_OCS/xyBnG_phased/phased_matrix_v2.jl")

test = "."
bar = "3-pgocs_10"
foo =  "$tag-rand"

ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
cp("$test/$foo.xy", xy, force=true)

lmp = deserialize("3-founder.lmp")
F0 = nothing
ε = 1e-6

species = Cattle(5_000)
    trait = Trait("growth", 0.25, 10_000)
    nchp = 50_000
    nref = 10_000
    nrng = 5
    nsel = 5
    plan = Plan(100, 100, 400) #Plan(25, 50, 200) Plan(12, 24, 96)
    fixed = ["grt"]
    dF = 0.005 #0.022 0.011
    nrpt = 5
    keep = true


for ign in 1:ngn
    rader = nrow(ped)
    @info " -  number of animals for generation $ign is : $rader "
    #print("Generation  $ign")
    ids = view(ped, ped.grt .== ped.grt[end], :id)
    phenotype!(ids, ped, trait)
    G = grm(xy, lmp.chip, lmp.frq) + ε * I
    giv = inv(G)
    Predict!(ids, ped, fixed, giv, trait)
   
    # subset the xy file 
    r,c = XY.dim(xy)  
    XY.sub("$test/$bar.xy", 1:r, c-plan.noff*2+1:c, "$xy-sub")
    

    @info "  - Generating segments"
    
    g22 =  prm(test, bar, "$xy-sub", lmp, l)

    if ign == 1 
        F0 =  mean(diag(g22) ) -1
        @info "  - Phased estimated F0 for length $l is $F0 "
    end 

    ng = Select(ids, plan, ped, g22, trait, dF, ign; F0 = F0)

    reproduce!(ng, ped, xy, lmp, trait)
    GC.gc()