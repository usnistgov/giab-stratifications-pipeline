from common.config import (
    CoreLevel,
    si_to_cds,
    si_to_mhc,
    si_to_kir,
    si_to_vdj,
    bd_to_other_bed,
    strip_full_refkey,
    MutualPathPair,
    BUILTIN_OTHER,
)

# NOTE this contains rules that generate stratifications in Functional,
# OtherDifficult, and others without explicit rules (Ancestry et al)

func = config.to_bed_dirs(CoreLevel.FUNCTIONAL)
odiff = config.to_bed_dirs(CoreLevel.OTHER_DIFFICULT)

# TODO don't hardcode this...
all_region_types = ["cds", "vdj", "kir", "mhc"]


def all_otherdifficult_sources(ref_key, build_key):
    return config.all_otherdifficult_sources(
        ref_key,
        build_key,
        Path(rules.download_gaps.output[0]),
        Path(rules.download_cds.output[0]),
        Path(rules.download_kir.output[0]),
        Path(rules.download_mhc.output[0]),
        Path(rules.download_vdj.output[0]),
        Path(rules.download_other.output[0]),
    )


def all_otherdifficult_sources_wc(wildcards):
    return all_otherdifficult_sources(wildcards["ref_key"], wildcards["build_key"])


def all_functional(ref_final_key, build_key):
    return config.all_functional(
        ref_final_key,
        build_key,
        all_otherdifficult_sources(strip_full_refkey(ref_final_key), build_key),
        MutualPathPair(
            Path(rules.merge_cds.output[0]),
            Path(rules.invert_cds.output[0]),
        ),
        Path(rules.functional_readme.output[0]),
    )


def all_otherdifficult(ref_final_key, build_key):
    return config.all_otherdifficult(
        ref_final_key,
        build_key,
        all_otherdifficult_sources(strip_full_refkey(ref_final_key), build_key),
        Path(rules.get_gaps.output[0]),
        Path(rules.remove_kir_gaps.output[0]),
        Path(rules.remove_mhc_gaps.output[0]),
        Path(rules.remove_vdj_gaps.output[0]),
        Path(rules.remove_gaps_other.output[0]),
        Path(rules.otherdifficult_readme.output[0]),
    )


def all_misc(other_level_key, ref_final_key, build_key):
    return config.all_misc(
        ref_final_key,
        build_key,
        other_level_key,
        all_otherdifficult_sources(strip_full_refkey(ref_final_key), build_key),
        Path(rules.remove_gaps_other.output[0]),
        Path(rules.misc_readme.output[0]),
    )


################################################################################
# gaps


rule get_gaps:
    input:
        gapless=rules.get_gapless.output.auto,
        genome=rules.filter_sort_ref.output["genome"],
    output:
        odiff.final("gaps_slop15kb"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        complementBed -i {input.gapless} -g {input.genome} | \
        slopBed -i stdin -b 15000 -g {input.genome} | \
        mergeBed -i stdin | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


################################################################################
# CDS/VDJ/MHC/KIR


use rule download_gaps as download_cds with:
    output:
        func.src.data / "gff.txt.gz",
    params:
        src=lambda w: to_bed_src(si_to_cds, w),
    localrule: True
    log:
        func.src.log / "gff.log",


use rule download_gaps as download_vdj with:
    output:
        func.src.data / "vdj.bed.gz",
    params:
        src=lambda w: to_bed_src(si_to_vdj, w),
    localrule: True
    log:
        func.src.log / "vdj.log",


use rule download_gaps as download_mhc with:
    output:
        func.src.data / "mhc.bed.gz",
    params:
        src=lambda w: to_bed_src(si_to_mhc, w),
    localrule: True
    log:
        func.src.log / "mhc.log",


use rule download_gaps as download_kir with:
    output:
        func.src.data / "kir.bed.gz",
    params:
        src=lambda w: to_bed_src(si_to_kir, w),
    localrule: True
    log:
        func.src.log / "kir.log",


def functional_inputs(wildcards):
    src = all_otherdifficult_sources_wc(wildcards)
    return {
        "cds": src.refseq_sources,
        "mhc": src.mhc_sources,
        "kir": src.kir_sources,
        "vdj": src.vdj_sources,
    }


checkpoint normalize_cds:
    input:
        unpack(functional_inputs),
    output:
        **{k: func.inter.filtersort.data / f"{k}.json" for k in all_region_types},
    params:
        cds=lambda w: to_output_pattern(func, "cds", w),
        vdj=lambda w: to_output_pattern(func, "vdj", w),
        mhc=lambda w: to_output_pattern(func, "mhc", w),
        kir=lambda w: to_output_pattern(func, "kir", w),
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_key, w.build_key, lambda m: m.normalizeCds
        ),
    benchmark:
        func.inter.filtersort.bench / "normalize_cds.txt"
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/functional/normalize_cds.py"


rule merge_cds:
    input:
        bed=lambda w: read_named_checkpoint("normalize_cds", "cds", w),
        genome=rules.filter_sort_ref.output["genome"],
        gapless=rules.get_gapless.output.auto,
    output:
        func.final("refseq_cds"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        mergeBed -i {input.bed} | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


use rule _invert_autosomal_regions as invert_cds with:
    input:
        **invert_region_inputs(rules.merge_cds.output),
    output:
        func.final("notinrefseq_cds"),


rule remove_vdj_gaps:
    input:
        bed=lambda w: read_named_checkpoint("normalize_cds", "vdj", w),
        genome=rules.filter_sort_ref.output["genome"],
        gapless=rules.get_gapless.output.auto,
    output:
        odiff.final("VDJ"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        mergeBed -i {input.bed} | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


use rule remove_vdj_gaps as remove_mhc_gaps with:
    input:
        bed=lambda w: read_named_checkpoint("normalize_cds", "mhc", w),
        genome=rules.filter_sort_ref.output["genome"],
        gapless=rules.get_gapless.output.auto,
    output:
        odiff.final("MHC"),


use rule remove_vdj_gaps as remove_kir_gaps with:
    input:
        bed=lambda w: read_named_checkpoint("normalize_cds", "kir", w),
        genome=rules.filter_sort_ref.output["genome"],
        gapless=rules.get_gapless.output.auto,
    output:
        odiff.final("KIR"),


################################################################################
# other

# make sure the wildcards here can match everything except the "built-in" other-
# difficult beds (listed here)
#
# NOTE: this is necessary because while I have a validator in the config to
# prevent built-in bed files to be added, snakemake doesn't know about these
# and will give an ambiguous file exception without a constraint.

other_constraints = {
    "other_level_key": f"({'|'.join(config.all_other_levels)})",
    "other_strat_key": f"(?!({'|'.join(BUILTIN_OTHER)}))[A-Za-z0-9-._]+",
}


use rule download_gaps as download_other with:
    output:
        config.ref_src_dir
        / "{build_key}"
        / "other"
        / "{other_level_key}_{other_strat_key}.bed.gz",
    log:
        config.ref_src_dir
        / "{build_key}"
        / "other"
        / "{other_level_key}_{other_strat_key}.log",
    params:
        src=lambda w: config.buildkey_to_bed_src(
            lambda bd: bd_to_other_bed(w.other_level_key, w.other_strat_key, bd),
            w.ref_src_key,
            w.build_key,
        ),
    wildcard_constraints:
        **other_constraints,
    localrule: True


checkpoint normalize_other:
    input:
        lambda w: all_otherdifficult_sources_wc(w).other_source(
            w["other_level_key"], w["other_strat_key"]
        ),
    output:
        config.intermediate_build_hapless_dir
        / "other"
        / "{other_level_key}"
        / "{other_strat_key}.json",
    params:
        # TODO this is sloppy
        output_pattern=lambda w: expand(
            config.intermediate_build_dir
            / "other"
            / w.other_level_key
            / f"{w.other_strat_key}.bed.gz",
            ref_final_key="%s",
            build_key=w.build_key,
        )[0],
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_key, w.build_key, lambda m: m.normalizeOther
        ),
    conda:
        "../envs/bedtools.yml"
    wildcard_constraints:
        **other_constraints,
    script:
        "../scripts/python/bedtools/functional/normalize_other.py"


rule remove_gaps_other:
    input:
        bed=lambda w: read_checkpoint(
            "normalize_other",
            w,
            ["other_level_key", "other_strat_key"],
        ),
        genome=rules.filter_sort_ref.output["genome"],
        gapless=rules.get_gapless.output.auto,
    output:
        config.build_final_strat_path("{other_level_key}", "{other_strat_key}"),
    conda:
        "../envs/bedtools.yml"
    wildcard_constraints:
        **other_constraints,
    shell:
        """
        intersectBed -a {input.bed} -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


################################################################################
# readmes


rule functional_readme:
    input:
        common="workflow/templates/common.j2",
        description="workflow/templates/functional_description.j2",
        methods="workflow/templates/functional_methods.j2",
        _sources=lambda w: all_functional(
            w["ref_final_key"], w["build_key"]
        ).all_sources,
        bedtools_env="workflow/envs/bedtools.yml",
    params:
        paths=lambda w: all_functional(w["ref_final_key"], w["build_key"]),
    output:
        func.readme,
    conda:
        "../envs/templates.yml"
    localrule: True
    script:
        "../scripts/python/templates/format_readme/format_functional.py"


rule otherdifficult_readme:
    input:
        common="workflow/templates/common.j2",
        description="workflow/templates/otherdifficult_description.j2",
        methods="workflow/templates/otherdifficult_methods.j2",
        bedtools_env="workflow/envs/bedtools.yml",
        _sources=lambda w: all_otherdifficult(
            w["ref_final_key"], w["build_key"]
        ).all_sources,
    params:
        paths=lambda w: all_otherdifficult(w["ref_final_key"], w["build_key"]),
    output:
        odiff.readme,
    conda:
        "../envs/templates.yml"
    localrule: True
    script:
        "../scripts/python/templates/format_readme/format_otherdifficult.py"


rule misc_readme:
    input:
        common="workflow/templates/common.j2",
        description="workflow/templates/misc_description.j2",
        methods="workflow/templates/misc_methods.j2",
        bedtools_env="workflow/envs/bedtools.yml",
        _sources=lambda w: all_misc(
            w["other_level_key"],
            w["ref_final_key"],
            w["build_key"],
        ).all_sources,
    params:
        paths=lambda w: all_misc(
            w["other_level_key"],
            w["ref_final_key"],
            w["build_key"],
        ),
    output:
        config.build_final_readme_path("{other_level_key}"),
    conda:
        "../envs/templates.yml"
    localrule: True
    script:
        "../scripts/python/templates/format_readme/format_misc.py"
