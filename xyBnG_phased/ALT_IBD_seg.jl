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


# REMEMBER IS IT A ANON OR COW SIM!!!

include("/mnt/users/odwa/IBD_OCS/xyBnG_phased/phased_matrix_v2.jl")
include("/mnt/users/odwa/IBD_OCS/xyBnG_phased/base_filer.jl")
include("/mnt/users/odwa/IBD_OCS/xyBnG_phased/browing_prm.jl")
include("/mnt/users/odwa/IBD_OCS/xyBnG_phased/zooroh_prm.jl")

function IBD_seg(;
    data = "rst",
    baseDir = "tskit_50",
    testDir = "cattle_50",
    species = Cattle(6_000),
    animal = "cow",
    NE = 10000,
    trait = Trait("growth", 0.25, 10_000),
    nchp = 50_000,
    nref = 10_000,
    nrng = 10,
    nsel = 10,
    plan = Plan(100,100,400), #Plan(25, 50, 200) Plan(12, 24, 96) plan(100,100,400)
    fixed = ["grt"],
    dF = 0.005, #0.022 0.011  0.005
    nrpt = 50,
    keep = true
    )

    @info "  - Running $nchp "
    # Scenario recording
    base, test = "$data/$baseDir", "$data/$testDir"
    #base, test = "$baseDir", "$testDir"
    isdir("$test") || mkpath("$test")
    scenario = (Data=data,
                BaseDir = baseDir,
                TestDir = testDir,
                Species = species,
                Animal = animal,
                Trait = trait,
                Nchp = nchp, 
                Nref = nref,
                Nrng = nrng,
                Nsel = nsel,
                Plan = plan,
                Fixed = fixed,
                ΔF = dF, 
                Nrpt = nrpt,)

    xyBnG.Sum.savepar(scenario, "$test/scenario.par")
    isfile("$test/summary.ser") && rm("$test/summary.ser", force=true)

    # Prepare a base population
    println(pwd())

    
    maf = 0.0050
    NE = parse(Int64,NE)

    ALT_chkbase(base, NE, animal ,species ) # Prepare/verify a base population
    fxy, fmp = "$base/$(species.name).xy", "$base/$(species.name).lmp"

    hdr = xyBnG.xyTypes.header(fxy)
    hdr.type = 13
    xyBnG.XY.header!(fxy, hdr)

#    fix_header(fxy)
     # "agocs","ggocs","pgocs_0.5","pgocs_01","pgocs_05"
    schemes = ["iiocs", "ggocs","pgocs_01","pgocs_05","pgocs_08","pgocs_10"] # 
    #schemes = ["zoogocs", "iiocs", "pgocs_05"] # "rgocs"


    # Simulations
    #1:nrpt
    for irpt in 1:nrpt
        tag = lpad(irpt, ndigits(nrpt), '0')
        println()
        @info "==========> Repetition: $tag / $nrpt <=========="
        @info "  - Prepare a founder population"


        lmp, F0 = initPop(fxy, fmp, test, plan, maf, nchp, nref, nrng, trait, tag)
        @info "  - F0 for repetition $tag is $F0 "
        
        for scheme in schemes
            foo, bar = "$tag-rand", "$tag-$scheme"

            if occursin("pgocs",bar)

                l_cm = parse(Float64, split(bar, "_")[2])

                pgocs(test, foo, bar, lmp, nrng, nsel, trait, fixed, plan, dF, l_cm)

            elseif occursin("rgocs",bar)
                rgocs(test, foo, bar, lmp, nrng, nsel, trait, fixed, plan, dF, animal)
            elseif occursin("zoogocs",bar)
                zoogocs(test, foo, bar, lmp, nrng, nsel, trait, fixed, plan, dF)
            else 

                scheme_2 =  @eval $(Symbol(scheme)) 
                # re-add F0
                ped = deserialize("$test/$foo.ped")
                ids = view(ped, ped.grt .== ped.grt[end], :id)
                xy = "$test/$foo.xy"
                if first(scheme) == 'a' 
                    G = nrm(ped)
                    F0 = xyBnG.xps.meanoffd(G[ids,ids])
                elseif first(scheme) == 'g'
                    G = grm(xy, lmp.chip, lmp.frq) 
                    F0 = xyBnG.xps.meanoffd(G[ids,ids])
                    @info "  - GRM F0 for repetition $tag is $F0 "
                else
                    mid = size(ped, 1)
                    G = irm(xy, lmp.chip, mid+1-length(ids):mid) 
                    F0 = xyBnG.xps.meanoffd(G)
                     @info "  - IRM F0 for repetition $tag is $F0 "
                end 
                ε = 1e-6
                scheme_2(test, foo, bar, lmp, nsel, trait, fixed, plan, dF, F0; ε = ε)
            end
            summary = Sum.xysum("$test/$bar.ped", "$test/$bar.xy", lmp, trait)
            Sum.savesum("$test/summary.ser", summary)
        end

        if !keep
            for f in readdir("$test")
                occursin(Regex("^$(tag)"), f) && rm("$test/$f", force=true)
            end
        end

        cmd ="for file in *$testDir/$tag-*; do pigz \$file; done" 
        run(`bash -c $cmd`)

    end
    open("$test/scenario.par", "a") do io
        println(io, "Ended: ", time())
        
    end
end
