from os.path import splitext, basename
from pathlib import Path
from common.config import CoreLevel, MutualPathPair

mlty = config.to_bed_dirs(CoreLevel.MAPPABILITY)

gem_wc_constraints = {
    "l": "\d+",
    "m": "\d+",
    "e": "\d+",
}


def single_nonunique_inputs(ref_final_key, build_key):
    c = checkpoints.merge_nonunique.get(
        ref_final_key=ref_final_key, build_key=build_key
    )
    with c.output.single_lowmap.open() as f:
        return json.load(f)


def all_mappability(ref_final_key, build_key):
    # guard to prevent checkpoint from firing when we call this function
    bd = config.to_build_data_full(ref_final_key, build_key)
    if not bd.want_mappability:
        return None

    ss = single_nonunique_inputs(ref_final_key, build_key)
    return config.all_lowmap(
        ref_final_key,
        build_key,
        MutualPathPair(
            Path(rules.merge_nonunique.output.all_lowmap),
            Path(rules.invert_merged_nonunique.output[0]),
        ),
        [Path(p) for p in ss],
        Path(rules.mappability_readme.output[0]),
    )


################################################################################
# index/align


def plaid_mode(wildcards):
    t = config.thread_per_chromosome(
        wildcards.ref_final_key,
        wildcards.build_key,
        4,
        True,
        False,
    )
    return t * 1.5


def filter_mappability_ref_inputs(wildcards):
    rk = config.refkey_strip_if_dip1(wildcards["ref_final_key"], False)
    return {
        "fa": expand(
            rules.download_ref.output,
            allow_missing=True,
            ref_src_key=rk,
        )[0],
        "idx": expand(
            rules.index_full_ref.output,
            allow_missing=True,
            ref_final_key=rk,
        )[0],
    }


rule filter_mappability_ref:
    input:
        unpack(filter_mappability_ref_inputs),
    output:
        mlty.inter.postsort.data / "ref.fa",
    log:
        mlty.inter.postsort.log / "filter_mappability_ref.txt",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/mappability/filter_ref.py"


rule gem_index:
    input:
        fa=rules.filter_mappability_ref.output,
        bin=rules.unpack_gem.output.indexer,
    output:
        mlty.inter.postsort.data / "index.gem",
    params:
        base=lambda wildcards, output: splitext(output[0])[0],
    threads: plaid_mode
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_final_key, w.build_key, lambda m: m.gemIndex
        ),
    log:
        mlty.inter.postsort.log / "index.log",
    benchmark:
        mlty.inter.postsort.bench / "index.txt"
    shell:
        """
        PATH={config.tools_bin_dir}:$PATH
        {input.bin} \
        --complement emulate \
        -T {threads} \
        -i {input.fa} \
        -o {params.base} > {log} 2>&1
        """


rule gem_mappability:
    input:
        idx=rules.gem_index.output[0],
        bin=rules.unpack_gem.output.mappability,
    output:
        mlty.inter.postsort.data / "unique_l{l}_m{m}_e{e}.mappability",
    params:
        base=lambda wildcards, output: splitext(output[0])[0],
    threads: plaid_mode
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_final_key, w.build_key, lambda m: m.gemMappability
        ),
    log:
        mlty.inter.postsort.data / "mappability_l{l}_m{m}_e{e}.log",
    benchmark:
        mlty.inter.postsort.bench / "mappability_l{l}_m{m}_e{e}.txt"
    wildcard_constraints:
        **gem_wc_constraints,
    shell:
        """
        {input.bin} \
        -m {wildcards.m} \
        -e {wildcards.e} \
        -l {wildcards.l} \
        -T {threads} \
        -I {input.idx} \
        -o {params.base} > {log} 2>&1
        """


rule gem_to_wig:
    input:
        idx=rules.gem_index.output,
        map=rules.gem_mappability.output[0],
        bin=rules.unpack_gem.output.gem2wig,
    output:
        mlty.inter.postsort.data / "unique_l{l}_m{m}_e{e}.wig",
    params:
        base=lambda wildcards, output: splitext(output[0])[0],
    log:
        mlty.inter.postsort.log / "gem2wig_l{l}_m{m}_e{e}.log",
    benchmark:
        mlty.inter.postsort.bench / "gem2wig_l{l}_m{m}_e{e}.txt"
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_final_key, w.build_key, lambda m: m.gemToWig
        ),
    wildcard_constraints:
        **gem_wc_constraints,
    shell:
        """
        {input.bin} \
        -I {input.idx} \
        -i {input.map} \
        -o {params.base} 2> {log}
        """


# NOTE: -d option to not sort and save some speed/memory (we shall sort later,
# because this command does not do it the way I want)
rule wig_to_bed:
    input:
        rules.gem_to_wig.output,
    output:
        mlty.inter.postsort.data / "unique_l{l}_m{m}_e{e}.bed.gz",
    conda:
        "../envs/map.yml"
    wildcard_constraints:
        **gem_wc_constraints,
    shell:
        """
        sed 's/ AC//' {input} | \
        wig2bed -d | \
        awk '$5>0.9' | \
        cut -f1-3 | \
        gzip -c > \
        {output}
        """


rule combine_dip1_nonunique_beds:
    input:
        unpack(lambda w: combine_dip_inputs("wig_to_bed", w)),
    output:
        mlty.inter.postsort.data / "combined_unique_l{l}_m{m}_e{e}.bed.gz",
    shell:
        """
        cat {input.pat} {input.mat} > {output}
        """


################################################################################
# create stratifications


def merge_nonunique_inputs(wildcards):
    rk = wildcards["ref_final_key"]
    bk = wildcards["build_key"]
    l, m, e = config.to_build_data_full(rk, bk).mappability_params
    out = if_dip1_else(
        False, False, "combine_dip1_nonunique_beds", "wig_to_bed", wildcards
    )
    return expand(out, zip, allow_missing=True, l=l, m=m, e=e)


checkpoint merge_nonunique:
    input:
        bed=merge_nonunique_inputs,
        gapless=rules.get_gapless.output.auto,
        genome=rules.filter_sort_ref.output["genome"],
    output:
        single_lowmap=mlty.inter.postsort.data / "single_lowmap.json",
        all_lowmap=mlty.final("lowmappabilityall"),
    log:
        mlty.inter.postsort.log / "nonunique_output.txt",
    params:
        path_pattern=lambda w: expand(
            mlty.final("{{}}"),
            ref_final_key=w.ref_final_key,
            build_key=w.build_key,
        )[0],
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_final_key, w.build_key, lambda m: m.mergeNonunique
        ),
    benchmark:
        mlty.inter.postsort.bench / "merge_nonunique.txt"
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/mappability/merge_nonunique.py"


use rule _invert_autosomal_regions as invert_merged_nonunique with:
    input:
        **invert_region_inputs(rules.merge_nonunique.output.all_lowmap),
    output:
        mlty.final("notinlowmappabilityall"),


# TODO what is the >0.9 thing in awk above?
rule mappability_readme:
    input:
        common="workflow/templates/common.j2",
        description="workflow/templates/mappability_description.j2",
        methods="workflow/templates/mappability_methods.j2",
        map_env="workflow/envs/map.yml",
        bedtools_env="workflow/envs/bedtools.yml",
    params:
        paths=lambda w: all_mappability(w["ref_final_key"], w["build_key"]),
    output:
        mlty.readme,
    conda:
        "../envs/templates.yml"
    localrule: True
    script:
        "../scripts/python/templates/format_readme/format_mappability.py"
