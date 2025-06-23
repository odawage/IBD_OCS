using DataFrames
using Mmap
using Statistics
import xyBnG.Util: commas
using xyBnG.XY

"""
    function meanoffd(a)
Possible the fastest way to calculate the mean off-diagonals of a square matrix.
"""

# function meanoffd(a)
#     m, n = size(a)
#     m == n || error("Not square")
#     (sum(a) - sum(diag(a))) / (n * (n - 1))
# end



function cow_genome(chr::Any)
    # [Chromosome length in bp](https://www.ncbi.nlm.nih.gov/projects/r_gencoll/ftp_service/nph-gc-ftp-service.cgi/?HistoryId=MCID_642d5d40ceff2e2c64293c60&QueryKey=1&ReleaseType=RefSeq&FileType=GENOME_FASTA&Flat=true)
    chrs = [
    1.58534110,
    1.36231102,
    1.21005158,
    1.20000601,
    1.20089316,
    1.17806340,
    1.10682743,
    1.13319770,
    1.05454467,
    1.03308737,
    1.06982474,
    0.87216183,
    0.83472345,
    0.82403003,
    0.85007780,
    0.81013979,
    0.73167244,
    0.65820629,
    0.63449741,
    0.71974595,
    0.69862954,
    0.60773035,
    0.52498615,
    0.62317253,
    0.42350435,
    0.51992305,
    0.45612108,
    0.45940150,
    0.51098607] 
    
    return chrs[chr]
end 



"""
    ts_base_alt(pop::Cattle, dir::AbstractString)

Simulate a cattle population of name `pop.name`, and `pop.nid` ID in `dir`.
with ne 
Note:
- This is a simulation with the coancestor/backward simulator msprime.
- Needs to have `tskit`, `msprime`, `scipy` and `stdpopsim` installed.
"""
function ts_base_alt(pop::Species, dir::AbstractString, NE::Any)
    '~' ∈ dir && error("Character '~' is forbidden in dir")
    isdir(dir) || mkpath(dir)
    @info "Simulating a cattle population of name $(pop.name), and $(pop.nid) ID in $dir
      Note: This is a simulation with the coancestor/backward simulator msprime."
    Threads.@threads for chr = 1:29
        print(" $chr")
        cmd = `python /mnt/users/odwa/IBD_OCS/xyBnG_phased/msprime/constant_ne.py "$dir" $chr $(pop.nid) $NE`
        run(pipeline(cmd, stderr = devnull))
       
    end
    open("$dir/desc.txt", "w") do io
        println(io, pop.name)
        println(io, pop.nid)
        println(io, NE)
    end
end


"""
    ts_base_NoSpecies(pop::Cattle, dir::AbstractString)

Simulate a population of name `pop.name`, and `pop.nid` ID in `dir`.
with ne 
Note:
- This is a simulation with the coancestor/backward simulator msprime.
- Needs to have `tskit`, `msprime`, `scipy` and `stdpopsim` installed.
"""
function ts_base_NoSpecies(pop::Species, dir::AbstractString, NE::Any)
    '~' ∈ dir && error("Character '~' is forbidden in dir")
    isdir(dir) || mkpath(dir)
    @info "Simulating a cattle population of name $(pop.name), and $(pop.nid) ID in $dir
      Note: This is a simulation with the coancestor/backward simulator msprime."
    Threads.@threads for chr = 1:30
        print(" $chr")
        cmd = `python /mnt/users/odwa/IBD_OCS/xyBnG_phased/msprime/Constant_ne_anon_species.py "$dir" $chr $(pop.nid) $NE`
        run(pipeline(cmd, stderr = devnull))
       
    end
    open("$dir/desc.txt", "w") do io
        println(io, pop.name)
        println(io, pop.nid)
        println(io, NE)
    end
end



"""
    chkbase(dir::AbstractString, pop::Species; base = ts_base)
Check if the base population `pop` is already created. If not, create it.
"""
function ALT_chkbase(dir::AbstractString,NE::Int64, animal::AbstractString, pop::Species; base = ts_base)
    if isfile("$dir/desc.txt")
        desc = readlines("$dir/desc.txt")
        desc[1] == pop.name && parse(Int, desc[2]) ≥ pop.nid && return
    end
    if base == ts_base
        if animal == "cow"
        ts_base_alt(pop, dir,NE)
        cmd = "bcftools concat $dir/{1..29}.vcf.gz -Oz -o $dir/Merged.vcf.gz"
        run(`bash -c $cmd`)
    
        Conn.vcf.toxy("$dir/Merged.vcf.gz", "$dir/BosTau")
        #Conn.TS.toxy(dir)
        else 
        ts_base_NoSpecies(pop, dir,NE)
        cmd = "bcftools concat $dir/{1..30}.vcf.gz -Oz -o $dir/Merged.vcf.gz"
        run(`bash -c $cmd`)
    
        Conn.vcf.toxy("$dir/Merged.vcf.gz", "$dir/BosTau")    
        end 

    elseif base == macs_base
        macs_base(pop, dir)
        Conn.MaCS.toxy(dir)
    else
        error("Unknown base population creator: $base")
    end
end


# for i in {1..29}; do
#     bcftools index ${i}.vcf.gz
# done

# bcftools concat \
#     {1..29}.vcf.gz \
#     -o merged_genome.vcf.gz \
#     -Oz     # Output in compressed format

# # Index the final merged VCF
# bcftools index merged_genome.vcf.gz

# Replace PRJXXXXXX with your project accession


