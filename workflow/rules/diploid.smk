from common.config import (
    CoreLevel,
    parse_full_refkey_class,
    strip_full_refkey,
    flip_full_refkey,
)

# TODO don't hardcode minimap2 params (which might be changed if we move to a
# different asm)

dip = config.to_bed_dirs(CoreLevel.DIPLOID)


def all_diploid(ref_final_key, build_key):
    bd = config.to_build_data_full(ref_final_key, build_key)

    def go(name):
        return [
            Path(p)
            for p in expand(
                getattr(rules, name).output,
                allow_missing=True,
                merge_len=bd.build.include.hets,
            )
        ]

    return config.all_diploid(
        ref_final_key,
        build_key,
        all_xy(ref_final_key, build_key),
        go("merge_het_regions"),
        go("invert_het_regions"),
        go("merge_het_SNVorSV_regions"),
        go("invert_het_SNVorSV_regions"),
        Path(rules.diploid_readme.output[0]),
    )


def minimap_inputs(wildcards):
    rk = wildcards["ref_final_key"]
    other_rk = flip_full_refkey(rk)
    out = if_dip1_else(
        True, True, "filter_sort_split_ref", "filter_sort_ref", wildcards
    )
    fa = out["fa"]
    idx = out["index"]

    def expand_rk(path, rk):
        return expand(path, allow_missing=True, ref_final_key=rk)

    return {
        "this_hap": expand_rk(fa, rk),
        "_this_idx": expand_rk(idx, rk),
        "other_hap": expand_rk(fa, other_rk),
        "_other_idx": expand_rk(idx, other_rk),
    }


# Dipcall normally outputs a bed file that roughly corresponds to regions with
# no breaks. Since we are interested in large structural variation in
# addition to small variants, obtain this bed file using the same process and
# then flip it.
rule cross_align_breaks:
    input:
        unpack(minimap_inputs),
    output:
        dip.inter.postsort.data / "breaks_cross_align.paf.gz",
    threads: lambda w: config.thread_per_chromosome(w.ref_final_key, w.build_key, 8, True, True)
    log:
        dip.inter.postsort.log / "cross_align_breaks.log",
    benchmark:
        dip.inter.postsort.bench / "cross_align_breaks.txt"
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_final_key, w.build_key, lambda m: m.crossAlignBreaks
        ),
    conda:
        "../envs/quasi-dipcall.yml"
    shell:
        """
        minimap2 -c --paf-no-hit -t{threads} --cs -z200000,10000,200 -xasm5 \
          {input.this_hap} \
          {input.other_hap} \
          2> {log} | \
        bgzip -c > {output}
        """


rule breaks_cross_alignment_to_bed:
    input:
        paf=rules.cross_align_breaks.output,
        paftools_bin=rules.download_paftools.output,
    output:
        dip.inter.postsort.data / "breaks_cross_align.bed.gz",
    log:
        dip.inter.postsort.log / "breaks_cross_alignment_to_bed.log",
    conda:
        "../envs/quasi-dipcall.yml"
    shell:
        """
        gunzip -c {input.paf} | \
        sort -k6,6 -k8,8n | \
        k8 {input.paftools_bin} call - 2> {log} | \
        grep ^R | \
        cut -f2,3,4 | \
        gzip -c > {output}
        """


# These need to be sorted since the paftools call script needs a sorted input
# but for convenience we don't sort numerically (doing so would require
# threading complex python logic through the otherwise elegant pipeline above).
# We could also presort (rather than sort twice) but this would require lots of
# space to store the paf rather than a relatively small bed file (ie no sequence
# strings)
#
# TODO misleading name since this also complements after sorting
rule sort_breaks_bed:
    input:
        bed=rules.breaks_cross_alignment_to_bed.output[0],
        genome=lambda w: if_dip1_else(
            True, True, "filter_sort_split_ref", "filter_sort_ref", w
        )["genome"],
    output:
        dip.inter.postsort.data / "sorted_breaks_cross_align.bed.gz",
    log:
        dip.inter.postsort.log / "sort_breaks_bed.txt",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/diploid/sort_breaks.py"


# Get variants (large and small) between the two haplotypes using the same
# process dipcall uses to produce its *.pair.vcf.gz file.
rule cross_align_variants:
    input:
        unpack(minimap_inputs),
    output:
        dip.inter.postsort.data / "variant_cross_align.sam.gz",
    threads: lambda w: config.thread_per_chromosome(w.ref_final_key, w.build_key, 8, True, True)
    log:
        dip.inter.postsort.log / "cross_align_variants.log",
    # this is what dipcall does to produce the pair file in one step
    conda:
        "../envs/quasi-dipcall.yml"
    benchmark:
        dip.inter.postsort.bench / "cross_align_variants.txt"
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_final_key, w.build_key, lambda m: m.crossAlignVariants
        ),
    shell:
        """
        minimap2 -a -t{threads} --cs -z200000,10000,200 -xasm5 \
          {input.this_hap} \
          {input.other_hap} \
          2> {log} | \
        bgzip -c > {output}
        """


rule filter_sort_variant_cross_alignment:
    input:
        aux_bin=rules.download_dipcall_aux.output,
        sam=rules.cross_align_variants.output,
    output:
        dip.inter.postsort.data / "sorted_variant_cross_alignments.bam",
    benchmark:
        dip.inter.postsort.bench / "filter_sort_variant_cross_alignment.txt"
    log:
        dip.inter.postsort.log / "filter_sort_variant_cross_alignment.txt",
    conda:
        "../envs/quasi-dipcall.yml"
    threads: 4
    resources:
        mem_mb=lambda w, threads: config.buildkey_to_malloc(
            w.ref_final_key,
            w.build_key,
            lambda m: m.filterSortVariantCrossAlignment,
        ),
    params:
        mem_per_thread=lambda w, threads, resources: int(resources.mem_mb / threads),
    shell:
        """
        k8 {input.aux_bin} samflt {input.sam} | \
        samtools sort -m{params.mem_per_thread}M --threads {threads} \
        2> {log} > {output}
        """


# TODO not sure if I also need to filter out alignments where one chr is
# clearly not mapped to the other, something like this: awk '{if ((gensub("_PATERNAL", "", 1, $1) == gensub("_MATERNAL", "", 1, $3)) || ("@" == substr($0,0,1))) { print $0 }}'


rule variant_cross_alignment_to_bed:
    input:
        bam=rules.filter_sort_variant_cross_alignment.output,
        hap=lambda w: if_dip1_else(
            True, True, "filter_sort_split_ref", "filter_sort_ref", w
        )["fa"],
    output:
        dip.inter.postsort.data / "variant_cross_align.bed.gz",
    conda:
        "../envs/quasi-dipcall.yml"
    # Convert vcf tob bed using POS and length of REF for the region bounds, and
    # add a fourth column for the absolute value of the indel length (for
    # filtering later)
    shell:
        """
        htsbox pileup -q5 -evcf {input.hap} {input.bam} | \
        grep -v '^#' | \
        awk 'OFS="\t" {{
          ref = length($4);
          alt = length($5);
          print $1, $2-1, $2+length($4)-1, (ref>alt) ? ref-alt : alt-ref
        }}' | \
        gzip -c > {output}
        """


rule filter_SNVorSV:
    input:
        rules.variant_cross_alignment_to_bed.output,
    output:
        dip.inter.postsort.data / "SNVorSV_variants.bed.gz",
    params:
        ## TODO don't hard code this
        sv_cutoff=50,
    conda:
        "../envs/quasi-dipcall.yml"
    shell:
        """
        gunzip -c {input} | \
        awk '$4>={params.sv_cutoff} || ($4==0 && ($3-$2==1))' | \
        cut -f1,2,3 | \
        gzip -c > {output}
        """


# gzip here so we can concatenate later without crossing the beams
rule merge_all_hets:
    input:
        rules.variant_cross_alignment_to_bed.output,
        rules.sort_breaks_bed.output,
    output:
        dip.inter.postsort.data / "het_regions.bed.gz",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        multiIntersectBed -i {input} | \
        mergeBed -i stdin | \
        gzip -c > {output}
        """


use rule merge_all_hets as merge_SNVorSV_hets with:
    input:
        rules.filter_SNVorSV.output,
        rules.sort_breaks_bed.output,
    output:
        dip.inter.postsort.data / "SNV_or_SV_het_regions.bed.gz",


rule combine_dip1_hets:
    input:
        unpack(lambda w: combine_dip_inputs("merge_all_hets", w)),
    output:
        dip.inter.postsort.data / "combined_het_regions.bed.gz",
    shell:
        """
        cat {input.pat} {input.mat} > {output}
        """


use rule combine_dip1_hets as combine_SNVorSV_dip1_hets with:
    input:
        unpack(lambda w: combine_dip_inputs("merge_SNVorSV_hets", w)),
    output:
        dip.inter.postsort.data / "combined_SNV_SV_het_regions.bed.gz",


# This het/hom analysis only makes sense for diploid chromosomes, which X and Y
# are not with the exception of PAR. Therefore, make dummy bed for use in taking
# out the PAR. I don't necessarily want this to fail in the case of not
# providing PAR coordinates, so if X/Y are not provided, simply don't include
# that region. This obviously might result in a blank file, but this shouldn't
# matter when input to the second argument of substractBed (ie "subtract
# nothing").
rule concat_xy_nonpar:
    input:
        lambda w: all_xy(w.ref_final_key, w.build_key).sex.nonpar_paths,
        # lambda w: all_xy_nonPAR(w.ref_final_key, w.build_key),
    output:
        dip.inter.postsort.data / "xy_nonpar.bed",
    # /dev/null will keep cat from hanging in case there are no input files. The
    # -f flag is necessary to prevent gunzip from whining when it doesn't see a
    # header. Thus if we have no inputs, we get a blank file, which is what we
    # want.
    shell:
        "cat {input} /dev/null | gunzip -c -f > {output}"


rule merge_het_regions:
    input:
        bed=lambda w: if_dip1_else(
            False, True, "combine_dip1_hets", "merge_all_hets", w
        ),
        valid=rules.get_valid_regions.output.auto,
        genome=rules.filter_sort_ref.output["genome"],
        xy_nonpar=rules.concat_xy_nonpar.output,
    output:
        dip.final("het_regions_{merge_len}"),
    conda:
        "../envs/bedtools.yml"
    params:
        size=lambda w: int(w["merge_len"]),
    wildcard_constraints:
        merge_len=f"\d+",
    # NOTE: The initial gunzip is very necessary despite the fact that mergeBed
    # can read gzip'ed files. For the dip1 case, the incoming bed will be two
    # concat'ed gzip'ed beds, and for some reason bedtools will only read the
    # first bed (probably because there's a header in the middle of the file?).
    # Unzipping first to a text stream will solve this issue.
    shell:
        """
        gunzip -c {input.bed} | \
        slopBed -i stdin -g {input.genome} -b {params.size} | \
        mergeBed -i stdin | \
        intersectBed -a stdin -b {input.valid} -sorted -g {input.genome} | \
        subtractBed -a stdin -b {input.xy_nonpar} -sorted | \
        bgzip -c > {output}
        """


use rule merge_het_regions as merge_het_SNVorSV_regions with:
    input:
        bed=lambda w: if_dip1_else(
            False, True, "combine_SNVorSV_dip1_hets", "merge_SNVorSV_hets", w
        ),
        valid=rules.get_valid_regions.output.auto,
        genome=rules.filter_sort_ref.output["genome"],
        xy_nonpar=rules.concat_xy_nonpar.output,
    output:
        dip.final("het_SNVorSV_regions_{merge_len}"),
    wildcard_constraints:
        merge_len=f"\d+",


rule invert_het_regions:
    input:
        bed=rules.merge_het_regions.output,
        valid=rules.get_valid_regions.output.auto,
        genome=rules.filter_sort_ref.output["genome"],
        xy_nonpar=rules.concat_xy_nonpar.output,
    output:
        dip.final("hom_regions_{merge_len}"),
    conda:
        "../envs/bedtools.yml"
    wildcard_constraints:
        merge_len=f"\d+",
    shell:
        """
        complementBed -i {input.bed} -g {input.genome} | \
        intersectBed -a stdin -b {input.valid} -sorted -g {input.genome} | \
        subtractBed -a stdin -b {input.xy_nonpar} -sorted | \
        bgzip -c > {output}
        """


use rule invert_het_regions as invert_het_SNVorSV_regions with:
    input:
        bed=rules.merge_het_SNVorSV_regions.output,
        valid=rules.get_valid_regions.output.auto,
        genome=rules.filter_sort_ref.output["genome"],
        xy_nonpar=rules.concat_xy_nonpar.output,
    output:
        dip.final("hom_SNVorSV_regions_{merge_len}"),
    wildcard_constraints:
        merge_len=f"\d+",


rule diploid_readme:
    input:
        common="workflow/templates/common.j2",
        description="workflow/templates/diploid_description.j2",
        methods="workflow/templates/diploid_methods.j2",
        bedtools_env="workflow/envs/bedtools.yml",
        diploid_env="workflow/envs/quasi-dipcall.yml",
    params:
        paths=lambda w: all_diploid(w["ref_final_key"], w["build_key"]),
    output:
        dip.readme,
    conda:
        "../envs/templates.yml"
    localrule: True
    script:
        "../scripts/python/templates/format_readme/format_diploid.py"
