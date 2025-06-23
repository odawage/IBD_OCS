"""
    phasedibd(env::AbstractString, vcf::AbstractString, lmp::AbstractString)

Calculate IBD relationships with genotypes in `vcf` and (plink formatted)
linkage map `lmp`, using Python package `phasedibd`. The trick here is that the
`phasedibd` package is old. It is using Cython language level 2 to compile. It
also needs older python environment. 

In my test, I used conda configured Python 3.7.16.

```bash
conda create -n py37 python=3.7
git clone https://github.com/23andMe/phasedibd
cd phasedibd
make
python setup.py install
```

You need to tell this Julia `phasedibd` function the environment name, in which
that Python `phasedibd` package was installed. For example, `py37`, as shown
above. Then the calculation is done by run ``Cmd`` `conda run -n py37 python
phasedibd.py`.`

"""

using CodecZlib
using CSV
using DataFrames
using Mmap
using Serialization
using VCFTools
using PyCall

"""
Calls the phasedibd function to find the IBD segments 
    chr = 
    vcf = 

"""
function generate_seg(chr::Int64, vcf::AbstractString, path::AbstractString, l_cm::Float64)
    isdir("$path/segtemp") || mkdir("$path/segtemp")
    VCFTools.filter_chr("$path/$vcf", chr, des="$path/segtemp/$chr.vcf")

    cmd = `python /mnt/users/odwa/IBD_OCS/xyBnG_phased/phasedibd_julia.py "$path" $chr $l_cm`
    run(pipeline(cmd, stderr = devnull))
end 

#python /mnt/users/odwa/IBD_OCS/xyBnG_phased/phasedibd_julia.py "." 1 10

#  function generate_seg(chr::Int64, vcf::AbstractString, path::AbstractString, l_cm::Float64)
#     isdir("$path/segtemp") || mkdir("$path/segtemp")
#     VCFTools.filter_chr("$path/$vcf", chr, des="$path/segtemp/$chr.vcf")

#     py"""
#     import numpy as np
#     import pandas as pd
#     import os
#     import sys
#     import unittest
#     import phasedibd as ibd
    
#     l_cm = $l_cm
#     path = $path
#     chr = str($chr)

#     TEST_DATA_PATH = os.path.abspath(path + "/segtemp") 
#     #print(TEST_DATA_PATH)  # This will print the constructed path for verification
    
#     haplotypes = ibd.VcfHaplotypeAlignment(TEST_DATA_PATH +"/" + chr + ".vcf")
    
    
#     # This the template  basically decides how many error, phasing errors and missing you're ok with 
#     tpbwt = ibd.TPBWTAnalysis(template=[[1, 0, 1, 0],
#      [0, 1, 0, 1],
#      [1, 1, 0, 0],
#      [0, 0, 1, 1],
#      [1, 0, 0, 1],
#      [0, 1, 1, 0]])
    
#     # L_m is the minimim number of SNPs in a segment 
#     # L_F is the min length in cM  
#     # Use_phase_correction should technically fix phasing errors using a heuristic scan. 
#     # This means that technically if you have phased data you should not get phasing errors. 
#     # I've only That functinon on unphased data and  not found a difference in overlapped segments . 
    
#     ibd_results = tpbwt.compute_ibd(haplotypes, L_m=20, L_f=l_cm, segments_out_path=TEST_DATA_PATH, use_phase_correction=False, verbose=False)
    
#     #, use_phase_correction=TRUE
#     """

# end 

"""
    splice(csv::AbstractString, saved::AbstractString; maxid = 0)
Read GZipped CSV file `csv`, splice IBD haplotypes where they overlaps. The
spliced data is zipped and saved in `saved` in the same directory for future
reuse. `maxid` is the ID square dims that have been calculated.
"""

"""
V2 of the splice 
"""

function Splice(csv::AbstractString, saved::AbstractString, maxid = 0)
    
    segments = CSV.read(csv, DataFrame)
    rename!(segments, :end => :stop, :end_cm => :stop_cm, :end_bp => :stop_bp, :chromosome => :chr)
    # Sort the DataFrame by relevant columns
    sort!(segments, [:id1, :id2, :id1_hap, :id2_hap, :start_bp])
    

    merged_segments = DataFrame(
            id1 = Int[],
            id2 = Int[],
            id1_hap = Int[],
            id2_hap = Int[],
            start = Int[],
            stop = Int[],
            start_cm = Float64[],
            stop_cm = Float64[],
            start_bp = Int[],
            stop_bp = Int[],
            chr = Int[],
        )

    

    # Initialize the first row as the current segment to compare with
    current_segment = segments[1,:]

    # Iterate over the DataFrame starting from the second row
    for i in 2:nrow(segments)
        if  segments[i, :chr] == current_segment[ :chr] &&
            segments[i, :id1] == current_segment[ :id1] &&
            segments[i, :id2] == current_segment[ :id2] &&
            segments[i, :id1_hap] == current_segment[ :id1_hap] &&
            segments[i, :id2_hap] == current_segment[ :id2_hap] &&
            segments[i, :start_bp] <= current_segment[ :stop_bp]   # Handle overlapping or adjacent ranges
           
            # Merge segments: keep the smallest start_bp, start, start_cm, and the largest stop_bp, end_pos, stop_cm
            current_segment[ :stop] = max(current_segment[ :stop], segments[i, :stop])
            current_segment[ :stop_cm] = max(current_segment[ :stop_cm], segments[i, :stop_cm])
            current_segment[ :stop_bp] = max(current_segment[ :stop_bp], segments[i, :stop_bp])
        else
            # If no merge is possible, push the current_segment to the merged_segments and update current_segment
            push!(merged_segments, current_segment)
            current_segment = segments[i, :]
        end
        #push!(merged_segments, current_segment)
    end

    # Push the last segment after loop ends
    push!(merged_segments, current_segment)

    serialize(saved, merged_segments)
end



"""
    seg2mat(fseg::AbstractString, fmat::AbstractString)
Read the DataFrame of segments `fseg`, and write the result IBD relationship
matrix to `fmat`.
"""


function seg2mat(path::AbstractString, spath::AbstractString, max_chr::Any)
    genome_l = 0 
    for chr in 1:max_chr
        first = read(pipeline(ignorestatus(`grep -v '^#' $spath/$chr.vcf`), `head -n 1`, `awk '{print $2}'`), String) |> strip
        first = parse(Int, first)
        
        last = read(pipeline(ignorestatus(`grep -v '^#' $spath/$chr.vcf`), `tail -n 1`, `awk '{print $2}'`), String) |> strip
        last = parse(Int, last) 
        genome_l = genome_l + (last - first) 
    end
    genome_l

    samples = VCFTools.sampleID("$spath/1.vcf") 
    mat = zeros(length(samples), length(samples))
    
    seg = DataFrame()
    for chr in 1:max_chr
        isfile("$spath/$chr.seg") || error("Segment file $chr not found.")
        t_seg = deserialize("$spath/$chr.seg")
        seg = vcat(seg, t_seg)
    end 

    seg.len = seg.stop_bp - seg.start_bp
    grouped_df = groupby(seg, [:id1, :id2])
    result_df = combine(grouped_df, :len => (x -> sum(x) ) => :sum)
    result_df.rel = [row.id1 == row.id2 ? 1 + (row.sum / genome_l) : (row.sum / (genome_l * 4))*2 for row in eachrow(result_df)]
    #result_df.rel = [row.id1 == row.id2 ?  (row.sum / genome_l) : (row.sum / (genome_l * 4)) for row in eachrow(result_df)]

    for row in eachrow(result_df)
        x_idx = row.id1 + 1
        y_idx = row.id2 + 1
        mat[x_idx, y_idx] = mat[y_idx, x_idx] = row.rel
    end
    # for some of the combinations there is no observations for the diagonal. so we need to add the 1 when missing. 
    for i in 1:size(mat, 1)
        if mat[i, i] == 0.0
            mat[i, i] = 1.0
        end
    end
    mat

end

# prm(inpath, "tester", "$inpath/xy-sub", lmp, l)
  
function prm(test::AbstractString, bar::AbstractString, xy::AbstractString, lmp::DataFrame, l_cm::Float64)
    max_chr = maximum(lmp.chr)
    serialize("$test/$bar.lmp", lmp)
    Conn.xy.tovcf("$xy", "$test/$bar.lmp", "$test/$bar.vcf")

    # file_lock = ReentrantLock()
    GC.gc()

        Threads.@threads for chr in 1:max_chr 
            local_l_cm = Float64(l_cm) 

            #lock(file_lock) do
                generate_seg(chr, "$bar.vcf", "$test", local_l_cm)
            #end
                Splice("$test/segtemp/$chr.csv","$test/segtemp/$chr.seg")  
        end
    G = seg2mat(test, "$test/segtemp", max_chr)
    G
end



"""
pgocs(test, foo, bar, lmp, ngn, trait, fixed, plan)
OCS with `phased IBD` relationship matrix for constraint on `foo`.xy and
`foo`.ped in directory `test` for `ngn` generations. SNP linkage information are
in DataFrame `lmp`. The results are saved in `bar`.xy, `bar`.ped in directory
`test`. The selection is on a single trait `trait` with fixed effects `fixed`,
which is a column name vector in pedigree DataFrame. The selection is based on
the selection plan `plan`.

See also [`gblup`](@ref), [`ablup`](@ref).
"""

function pgocs(test, foo, bar, lmp, nrng, ngn, trait, fixed, plan, dF, l)
    @info "  - Directional selection PGOCS with $l cM segments, for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force=true)
    isdir("$test/segtemp") || mkdir("$test/segtemp")
    
    F0 = nothing

    for ign in 1:ngn
        rader = nrow(ped)
        @info " -  number of animals for generation $ign is : $rader "
        #print("Generation  $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        ε = 1e-6
        GC.gc()
        G = grm(xy, lmp.chip, lmp.frq)+ ε * I
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        G = nothing
        giv = nothing
        # subset the xy file 
        r,c = XY.dim(xy)  
        XY.sub("$test/$bar.xy", 1:r, c-plan.noff*2+1:c, "$xy-sub")
        

        @info "  - Generating segments"
        g22 =  prm(test, bar, "$xy-sub", lmp, l)

        if ign == 1 
            F0 =  xyBnG.xps.meanoffd(g22)
            @info "  - Phased estimated F0 for length $l is $F0 "
        end 

        ng = Select(ids, plan, ped, g22, trait, dF, ign; F0 = F0)

        reproduce!(ng, ped, xy, lmp, trait)
        GC.gc()
    end
    println()
    serialize("$test/$bar.ped", ped)
end

function pgocs_diag(test, foo, bar, lmp, nrng, ngn, trait, fixed, plan, dF, l)
    @info "  - Directional selection PGOCS with $l cM segments, for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force=true)
    isdir("$test/segtemp") || mkdir("$test/segtemp")
    
    F0 = nothing

    for ign in 1:ngn
        rader = nrow(ped)
        @info " -  number of animals for generation $ign is : $rader "
        #print("Generation  $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        ε = 1e-6
        G = grm(xy, lmp.chip, lmp.frq)+ ε * I
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
       
        # subset the xy file 
        r,c = XY.dim(xy)  
        XY.sub("$test/$bar.xy", 1:r, c-plan.noff*2+1:c, "$xy-sub")
        

        @info "  - Generating segments"
        
        g22 =  prm(test, bar, "$xy-sub", lmp, l)

        if ign == 1 
            F0 =  xyBnG.xps.meanoffd(g22)
            @info "  - Phased estimated F0 for length $l is $F0 "
        end 

        ng = Select(ids, plan, ped, g22, trait, dF, ign; F0 = F0)

        reproduce!(ng, ped, xy, lmp, trait)
        GC.gc()
    end
    println()
    serialize("$test/$bar.ped", ped)
end


function fix_header(file)
    hdr = xyBnG.XY.header(file)
    nlc, nhp = xyBnG.XY.dim(file)
    fz = filesize(file)
    if hdr.type ≠ 13
        sz = sizeof(xyBnG.xyTypes._type(hdr.type))
        if fz == sz * nlc * nhp + 24
            @info "This is a proper file, no fix needed"
            return
        end
        if abs(fz - nlc * nhp / 8 - 24) < 10
            @info "This is maybe of BitArray, header type will be changed ..."
            hdr.r, hdr.type = 1, 13
            xyBnG.XY.header!(file, hdr)
            @info "Done"
        end
    else
        if abs(fz - nlc * nhp / 8 - 24) < 10
            @info "This is a proper file, no fix needed"
        else
            @info "I don't know how to fix this file"
        end
    end
end
