using Base.Threads
using ClusterManagers
using CSV
using DataFrames
using Distributed
using LinearAlgebra
using Markdown
using Mmap
using PyCall
using Random
using Revise
using Base.Threads
using Serialization
using Statistics
using VCFTools
using xyBnG
using xyBnG.XY
using xyBnG.Conn
using xyBnG.Sum
import xyBnG.xps: agocs,ggocs, iiocs, igocs, iiocs,initPop, randbrd 
import xyBnG.xyTypes: Cattle, Species, Trait, Plan, name
import xyBnG.Breeding: phenotype!, Predict!, Select, reproduce!, TM1997, TM2024
import xyBnG.Founder: ts_base, uniq, sample_xy, macs_base
import xyBnG.RS: nrm, irm, grm, xirm

include("/mnt/users/odwa/paper-2/xyBnG_phased/phased_matrix.jl")


"""
The goal for this script is to write a function for the whole segment pipeline. 
And a function or loop that calculates the rel-mat for all the replicas of the simulation. 

In the all replicas loop I'am going to use the final .xy file to end up with one giant matrix per replica.
    If needed I can later make within generation subsets. 
I might also have to write a function that converts a matrix to a long form dataframe to be able to plot it. 

"""


function seg3mat(path::AbstractString, spath::AbstractString)
    genome_l = 0 
    for chr in 1:29 
        first = read(pipeline(ignorestatus(`grep -v '^#' $spath/$chr.vcf`), `head -n 1`, `awk '{print $2}'`), String) |> strip
        first = parse(Int, first)
        
        last = read(pipeline(ignorestatus(`grep -v '^#' $spath/$chr.vcf`), `tail -n 1`, `awk '{print $2}'`), String) |> strip
        last = parse(Int, last) 
        genome_l = genome_l + (last - first) 
    end
    

    samples = VCFTools.sampleID("$spath/1.vcf") 
    mat = zeros(length(samples), length(samples))
    
    seg = DataFrame()
    for chr in 1:29
        isfile("$spath/$chr.seg") || error("Segment file $chr not found.")
        t_seg = deserialize("$spath/$chr.seg")
        seg = vcat(seg, t_seg)
    end 

    seg.len = seg.stop_bp - seg.start_bp
    grouped_df = groupby(seg, [:id1, :id2])
    result_df = combine(grouped_df, :len => (x -> sum(x) ) => :sum)
    result_df.rel = [row.id1 == row.id2 ? 1 + (row.sum / genome_l) : (row.sum / (genome_l * 4))*2 for row in eachrow(result_df)]
    
    for row in eachrow(result_df)
        x_idx = row.id1 + 1
        y_idx = row.id2 + 1
        mat[x_idx, y_idx] = mat[y_idx, x_idx] = row.rel
    end
    mat

end


"""
inpath:     where the .xy and the .lmp files are 
outpath:    where all the inbetween and final files will reside
xy:         full name of the xy file 
lmp:        Full name of the lmp file. 
"""

function prm(inpath::AbstractString, outpath::AbstractString, xy::AbstractString, lmp::AbstractString, vcf::AbstractString, length::Int64)
    Conn.xy.tovcf("$xy", "$inpath/$lmp", "$outpath/$vcf")

    file_lock = ReentrantLock()
        Threads.@threads for chr in 1:29 
            lock(file_lock) do
                generate_seg(chr, "$vcf", "$outpath", length)
            end
                Splice("$outpath/segtemp/$chr.csv","$outpath/segtemp/$chr.seg")  
        end
        G = seg3mat(outpath, "$outpath/segtemp")
    G
end

# G = irm(xy, chp, id) + ε * I
# G = grm(xy, lmp.chip, lmp.frq) + ε * I

function GenRm(inpath::AbstractString, outpath::AbstractString, rep::Any, ped::AbstractString, xy::AbstractString, lmp::AbstractString, vcf::AbstractString, lengths::Any, E::Any )
    
    lmp_s = lmp
    
    ped, xy, lmp = deserialize("$inpath/$ped"), "$inpath/$xy", deserialize("$inpath/$lmp")

    ε = E
    # subset the xy file 
    r,c = XY.dim(xy)  
    XY.sub(xy, 1:r, c-400*10*2+1:c, "$inpath/xy-sub")
    xy_sub = "$inpath/xy-sub"

    #ped = Base.filter(row -> row[:grt] in Vector(11:20), ped)s
    ped = ped[Int(c/2)-(400*10)+1:Int(c/2),:]
    id = 1:size(ped, 1)

    G = irm(xy_sub, lmp.chip, id) + ε * I
    serialize("$outpath/$rep-IRM.mat", G)

    G =  grm(xy_sub, lmp.chip, lmp.frq) + ε * I
    serialize("$outpath/$rep-GRM.mat", G)


    for l in lengths
        G = prm(inpath, outpath, "$inpath/xy-sub", lmp_s, vcf, l) + ε * I
        serialize("$outpath/$rep-$l-PRM.mat", G)
        
    end 
end 


"""
Aaaaand then we run all the functions 
"""

wd = "/mnt/SCRATCH/odwa/xyBnG_sim/rst"
outpath = "$wd/cattle_50/sum_cat"
isdir("$outpath") || mkdir("$outpath")

nrpt = 10


for irpt in 1:nrpt 
    i = lpad(irpt, ndigits(nrpt), '0')
    GenRm("$wd/cattle_50", "$outpath", i, "$i-iiocs.ped", "$i-iiocs.xy","$i-founder.lmp","$i-iiocs.vcf",[01,05,10],1e-6 )
    GC.gc()

end


function matrix_to_dataframe(mat::Matrix)
    rows, cols = size(mat)
    df = DataFrame(id1 = Int[], id2 = Int[], value = Float64[])
    
    for r in 1:rows
        for c in 1:cols
            push!(df, (r, c, mat[r, c]))
        end
    end
    
    return df
end




df_total = DataFrame(id1 = Int[], id2 = Int[], value = Float64[], age1 = Int[], rep = Int[])
schemes = ["IRM","GRM","1-PRM","5-PRM","10-PRM"] #,"5-PRM","10-PRM"
nrpt = 10
n= 400
for scheme in schemes

    df_rep = DataFrame()
    for irpt in 1:nrpt
    irpt = lpad(irpt,2,"0")
    mat = deserialize("$outpath/$irpt-$scheme.mat")
        # if !occursin("PRM",scheme)
        #     r,c = size(mat)  
        #     mat = mat[1:r, c-n*10*2+1:c]
        # end 
    df = matrix_to_dataframe(mat)
    rename!(df, :value => Symbol("$(scheme)"))
    df.age1 = ceil.(Int, df.id1 ./ n)
    df.age2 = ceil.(Int, df.id2 ./ n)
    df[!,:rep] = fill(irpt,nrow(df))
    df_rep = vcat(df_rep, df)
    end 
    # If df_total is empty, initialize it with the first df_rep
    if nrow(df_total) == 0
        global df_total = df_rep
    else
        # Join df_total with df_rep
        global df_total = innerjoin(df_total, df_rep[!,Not([:age1,:age2])], on = [:id1, :id2, :rep], makeunique=true)
    end
end 


serialize("$outpath/df_total", df_total)
