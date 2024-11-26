using DataFrames
using LinearAlgebra
using Serialization
using xyBnG
import xyBnG.Sum: xysum, savesum, cormat, savepar
import xyBnG.xyTypes: Plan
import xyBnG.xps: initPop, chkbase, ggocs, iiocs, aaocs, riocs


sdir = @__DIR__
if !isfile("$sdir/pgsnp")
    @info "Compile pgsnp.cpp"
    run(`g++ -Wall -O3 --std=c++11 -pthread -o $sdir/pgsnp $sdir/pgsnp.cpp`)
end
chr = [1,2,3,4,5]
nchr = length([1,2,3,4,5])

nrng, nsel, maf, dF, hist, mr = 5, 30, 0.0, 0.011, 5000, 4.0
plan = Plan(25, 50, 200)
species = Cattle(plan.noff)

Threads.@threads for i = 1:nchr
    cmd = `$sdir/pgsnp $(species.nid) $hist $(chr[i]) $mr`
    run(pipeline(cmd, stdout = "$rst/chr.$i", stderr = devnull))
end
xyBnG.Conn.PG.toxy("$rst")