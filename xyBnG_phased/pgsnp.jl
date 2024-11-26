using DataFrames
using LinearAlgebra
using Serialization
using xyBnG
import xyBnG.Sum: xysum, savesum, cormat, savepar
import xyBnG.xyTypes: Plan
import xyBnG.xps: initPop, chkbase, ablup, gblup, iblup, aaocs, ggocs, iiocs, agocs, igocs, riocs

"""
    pgsnp(; data = "pgsnp", rst = "rst", nrpt = 1)
# Using `pgsnp` results as base population.

## Manule

Suppose you are in a working directory ``dir``. This file ``pgsnp.jl`` is in
``script/``. You want to save the results in ``rst``

- Date: 2024-08-27, works ok with xyBnG v1.2.4
"""
function pgsnp(chr;
    nchp = 600,
    nref = 1_000,
    rst = "1-chr-pgsnp",
    trait = Trait("growth", 0.25, 1_000),
    nrpt = 1,
    nsel = 30,
    nrng = 5,
    ne = 200,
    maf = 0.0,
    dF = 0.011,
    hist = 5000,
    mr = 4.0,
)
    sdir = @__DIR__
    if !isfile("$sdir/pgsnp")
        @info "Compile pgsnp.cpp"
        run(`g++ -Wall -O3 --std=c++11 -pthread -o $sdir/pgsnp $sdir/pgsnp.cpp`)
    end
    species = Cattle(ne)
    plan, plnb, fixed = Plan(25, 50, 200), Plan(50, 50, 200), ["grt"]
    CULLS = (gblup, ablup, iblup)
    OCSS = (aaocs, iiocs, ggocs, agocs, igocs)
    rst = abspath(rst)

    isdir(rst) && rm(rst, force = true, recursive = true)
    mkpath(rst)
    scenario = (
        Result = rst,
        Species = species,
        History = hist,
        MutationRate = mr, # per 1e8 bp per meiosis
        Trait = trait,
        Chr_len = chr,
        Nchp = nchp,
        Nref = nref,
        Nrng = nrng,
        Nsel = nsel,
        Plan = plan,
        Fixed = fixed,
        MAF = maf,
        ΔF = dF,
        Nrpt = nrpt,
        Schemes = union(CULLS, OCSS),
    )
    savepar(scenario, "$rst/scenario.par")

    npd = ndigits(nrpt)
    open("$rst/desc.txt", "w") do io
        println(io, "BosTau")
        println(io, species.nid)
    end
    for irpt = 1:nrpt
        tag = lpad(irpt, npd, '0')
        @info "==========> Repeat: $tag / $nrpt <=========="
        @info "  - Prepare a founder population"
        cmd = `$sdir/pgsnp $(species.nid) $hist $chr $mr`
        run(pipeline(cmd, stdout = "$rst/chr.1", stderr = devnull))
        xyBnG.Conn.PG.toxy("$rst")
        fxy = "$rst/$(species.name).xy"
        fmp = "$rst/$(species.name).lmp"
        lmp, F0 = initPop(fxy, fmp, rst, plan, maf, nchp, nref, nrng, trait, tag)
        lmp.chip = lmp.chip .&& .!lmp[!, trait.name] .&& .!lmp[!, "dark"]
        lmp.dark = lmp.dark .&& .!lmp[!, trait.name]

        for scheme in CULLS
            foo, bar = "$tag-rand", tag * '-' * string(scheme)
            scheme(rst, foo, bar, lmp, nsel, trait, fixed, plan)
            summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
            savesum("$rst/summary.ser", summary)
        end
        for scheme in OCSS
            foo, bar = "$tag-rand", tag * '-' * string(scheme)
            scheme(rst, foo, bar, lmp, nsel, trait, fixed, plnb, dF, F0)
            summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
            savesum("$rst/summary.ser", summary)
        end
    end
    open("$rst/scenario.par", "a") do io
        println(io, "Ended: ", time())
    end
end