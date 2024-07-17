from common.config import (
    si_to_gaps,
    bd_to_bench_bed,
    bd_to_bench_vcf,
    bd_to_query_vcf,
    refkey_config_to_prefix,
    ChrIndex,
)

ref = config.ref_dirs

# reference


rule download_ref:
    output:
        ref.src.reference.data / "ref.fna.gz",
    params:
        src=lambda w: config.refsrckey_to_ref_src(w.ref_src_key),
    log:
        ref.src.reference.log / "download_ref.log",
    conda:
        "../envs/bedtools.yml"
    localrule: True
    script:
        "../scripts/python/bedtools/misc/get_ref.py"


# note this is only needed for mappability where we need all chromosome names;
# for anything else we can get an index for the filtered and sorted FASTA
rule index_full_ref:
    input:
        lambda w: expand_final_to_src(rules.download_ref.output, w),
    output:
        ref.inter.prebuild.data / "ref.fna.fai",
    conda:
        "../envs/utils.yml"
    log:
        ref.inter.prebuild.log / "index_ref.log",
    shell:
        """
        samtools faidx {input} -o - 2> {log} > {output}
        """


# filtered/sorted reference + index + genome
#
# For each build, select the chromosomes we actually want and sort them in in
# numerical order. If we do this here, nearly all downstream steps won't need to
# be sorted, since they will query the reference in the order presented within.
# Additionally, make and index and genome file corresponding to each of these
# sorted/filtered references.
#
# Note that in the case of dip1 the reference can be split into two separate
# haplotypes or not. When split, there are also instances where we want to
# enforce diploid-only references (so either dip1 or dip2) and disallow hap
# references. Encode these three possibilities here


def filter_sort_ref_outputs(split, nohap):
    prefix = refkey_config_to_prefix(split, nohap)
    fa = f"{prefix}_filtered_ref.fa"
    parent = ref.inter.build.data
    return {
        "fa": parent / fa,
        "index": parent / (fa + ".fai"),
        "genome": parent / f"{prefix}_genome.txt",
    }


rule filter_sort_ref:
    input:
        lambda w: expand_final_to_src(rules.download_ref.output, w),
    output:
        **filter_sort_ref_outputs(False, False),
    conda:
        "../envs/utils.yml"
    log:
        ref.inter.build.log / "filter_sort_ref.log",
    script:
        "../scripts/python/bedtools/ref/filter_sort_ref.py"


use rule filter_sort_ref as filter_sort_split_ref with:
    input:
        lambda w: expand(
            rules.download_ref.output,
            allow_missing=True,
            ref_src_key=config.refkey_strip_if_dip1(w["ref_final_key"], False),
        ),
    output:
        **filter_sort_ref_outputs(True, False),
    log:
        ref.inter.build.log / "filter_sort_split_ref.log",


use rule filter_sort_ref as filter_sort_split_ref_nohap with:
    input:
        lambda w: expand(
            rules.download_ref.output,
            allow_missing=True,
            ref_src_key=config.refkey_strip_if_dip1(w["ref_final_key"], True),
        ),
    output:
        **filter_sort_ref_outputs(True, True),
    log:
        ref.inter.build.log / "filter_sort_split_ref_nohap.log",


rule download_gaps:
    input:
        rules.build_kent.output.bb2bed,
    output:
        ref.src.reference.data / "gap.bed.gz",
    params:
        src=lambda w: to_bed_src(si_to_gaps, w),
    localrule: True
    log:
        ref.src.reference.log / "download_gaps.log",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/misc/get_bedlike.py"


def valid_input(wildcards):
    rk = wildcards.ref_final_key
    bk = wildcards.build_key

    bd = config.to_build_data_full(rk, bk)
    sex_chrs = config.buildkey_to_wanted_xy(rk, bk)

    gaps_target = (
        expand(
            rules.download_gaps.output,
            ref_src_key=config.refkey_to_bed_refsrckeys_smk(si_to_gaps, rk),
        )
        if bd.have_gaps
        else None
    )

    par_target = (
        expand(
            rules.write_PAR_intermediate.output[0],
            allow_missing=True,
            sex_chr="Y",
        )[0]
        if bd.have_y_PAR and ChrIndex.CHRY in sex_chrs
        else None
    )

    return {
        k: v
        for k, v in [
            ("gaps", gaps_target),
            ("parY", par_target),
        ]
        if v is not None
    }


rule get_valid_regions:
    input:
        unpack(valid_input),
        genome=rules.filter_sort_ref.output["genome"],
    output:
        auto=ref.inter.build.data / "valid_regions_auto.bed.gz",
        parY=ref.inter.build.data / "valid_regions_parY.bed.gz",
    log:
        ref.inter.build.log / "valid_regions.txt",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/ref/get_valid_regions.py"


# benchmark


use rule download_gaps as download_bench_vcf with:
    output:
        ref.src.benchmark.data / "bench.vcf.gz",
    params:
        src=lambda w: config.buildkey_to_vcf_src(
            bd_to_bench_vcf,
            w.ref_src_key,
            w.build_key,
        ),
    localrule: True
    log:
        ref.src.benchmark.log / "download_bench_vcf.log",


use rule download_gaps as download_bench_bed with:
    output:
        ref.src.benchmark.data / "bench.bed.gz",
    params:
        src=lambda w: config.buildkey_to_bed_src(
            bd_to_bench_bed,
            w.ref_src_key,
            w.build_key,
        ),
    localrule: True
    log:
        ref.src.benchmark.log / "download_bench_bed.log",


use rule download_gaps as download_query_vcf with:
    output:
        ref.src.benchmark.data / "query.vcf.gz",
    params:
        src=lambda w: config.buildkey_to_vcf_src(
            bd_to_query_vcf,
            w.ref_src_key,
            w.build_key,
        ),
    localrule: True
    log:
        ref.src.benchmark.log / "download_query_vcf.log",


checkpoint normalize_bench_bed:
    input:
        lambda w: bed_src_bd_inputs(rules.download_bench_bed.output, bd_to_bench_bed, w),
    output:
        ref.inter.filtersort.data / "bench_filtered.json",
    conda:
        "../envs/bedtools.yml"
    params:
        output_pattern=lambda w: expand(
            ref.inter.filtersort.subbed / "bench_filtered.bed.gz",
            build_key=w.build_key,
        )[0],
    script:
        "../scripts/python/bedtools/ref/normalize_bench_bed.py"


# misc


# lots of things depend on PAR which is why this isn't part of the XY module
rule write_PAR_intermediate:
    output:
        ref.inter.build.data / "chr{sex_chr}_PAR.bed.gz",
    conda:
        "../envs/bedtools.yml"
    localrule: True
    script:
        "../scripts/python/bedtools/xy/write_par.py"


# skeleton rules (meant to be overridden elsewhere)


def invert_region_inputs(bed_path):
    return {
        "bed": bed_path,
        "genome": rules.filter_sort_ref.output["genome"],
        "valid": rules.get_valid_regions.output.auto,
    }


rule _invert_autosomal_regions:
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        complementBed -i {input.bed} -g {input.genome} | \
        intersectBed -a stdin -b {input.valid} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """
