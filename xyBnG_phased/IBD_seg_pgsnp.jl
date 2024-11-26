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

include("/mnt/users/odwa/paper-2/xyBnG_phased/phased_matrix_v2.jl")


function IBD_seg(;
    data = "rst",
    baseDir = "tskit_50",
    testDir = "cattle_test",
    species = Cattle(5_000),
    trait = Trait("growth", 0.25, 10_000),
    nchp = 50_000,
    nref = 10_000,
    nrng = 5,
    nsel = 5,
    plan = Plan(100, 100, 400), #Plan(25, 50, 200) Plan(12, 24, 96)
    fixed = ["grt"],
    dF = 0.005, #0.022 0.011
    nrpt = 5,
    keep = true
    )

    @info "  - OK leeeeets GO"
    # Scenario recording
    base, test = "$data/$baseDir", "$data/$testDir"
    #base, test = "$baseDir", "$testDir"
    isdir("$test") || mkpath("$test")
    scenario = (Data=data,
                BaseDir = baseDir,
                TestDir = testDir,
                Species = species,
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

    # Prepare a base population
    # if !isfile("$base/desc.txt")
    #     ts_base(species, base)
    #     desc = readlines("$base/desc.txt")
    # else
    #     desc = readlines("$base/desc.txt")
    #     parse(Int, desc[2]) < species.nid ||
    #         desc[1] ≠ species.name && ts_base(cattle, baseDir)
    # end

    # sname = desc[1]
    # fxy, fmp, maf = "$base/$sname.xy", "$base/$sname.lmp", 0.0
    
    maf = 0.0
    
    xyBnG.xps.chkbase(base, species) # Prepare/verify a base population
    fxy, fmp = "$base/$(species.name).xy", "$base/$(species.name).lmp"

    fix_header(fxy)
     # "agocs","ggocs","pgocs_0.5","pgocs_01","pgocs_05"
    schemes = ["pgocs_10"] # 

    # toxy = Conn.TS.toxy
    # isfile("$base/$sname.xy") || toxy(base)

    # Simulations
    #1:nrpt
    for irpt in 1:nrpt
        tag = lpad(irpt, ndigits(nrpt), '0')
        println()
        @info "==========> Repetition: $tag / $nrpt <=========="
        @info "  - Prepare a founder population"

        #sample_xy(fxy, fmp, test, plan.noff, maf, nchp, nref, trait)
        
        #mv("$test/founder.xy", "$test/snp.xy", force=true)
        
        #uniq("$test/snp.xy", "$test/founder.xy")
        #lmp = deserialize("$test/founder.lmp")
        
        # The starting point: random selection
        # randbrd(test, "founder", "$tag-rand", lmp, nrng, trait, plan; ibd=true)
        lmp, F0 = initPop(fxy, fmp, test, plan, maf, nchp, nref, nrng, trait, tag)
        @info "  - F0 for repetition $tag is $F0 "
        
        for scheme in schemes
            foo, bar = "$tag-rand", "$tag-$scheme"

            ifoccursin("pgocs",bar)

                l_cm = parse(Float64, split(bar, "_")[2])

                pgocs(test, foo, bar, lmp, nrng, nsel, trait, fixed, plan, dF, l_cm)
            else 
                scheme_2 =  @eval $(Symbol(scheme)) 
                # re-add F0
                ped = deserialize("$test/$foo.ped")
                ids = view(ped, ped.grt .== ped.grt[end], :id)
                xy = "$test/$foo.xy"
                if first(scheme) == 'a' 
                    G = nrm(ped)
                    F0 = mean(diag(G[ids,ids]))
                elseif first(scheme) == 'g'
                    G = grm(xy, lmp.chip, lmp.frq) 
                    F0 = mean(diag(G[ids,ids]))
                else
                    mid = size(ped, 1)
                    G = irm(xy, lmp.chip, mid+1-length(ids):mid) 
                    F0 = mean(diag(G))-1

                end 
                
                scheme_2(test, foo, bar, lmp, nsel, trait, fixed, plan, dF, F0)
            end
            summary = Sum.xysum("$test/$bar.ped", "$test/$bar.xy", lmp, trait)
            Sum.savesum("$test/summary.ser", summary)
        end
        if !keep
            for f in readdir("$test")
                occursin(Regex("^$(tag)"), f) && rm("$test/$f", force=true)
            end
        end
    end
    open("$test/scenario.par", "a") do io
        println(io, "Ended: ", time())
    end
end
