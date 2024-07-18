from functools import partial
from common.config import CoreLevel
from bedtools.gc.helpers import seqtk_args, range_bounds
from more_itertools import unzip
import json

gc = config.to_bed_dirs(CoreLevel.GC)


def gc_inputs(ref_final_key, build_key):
    c = checkpoints.intersect_gc_ranges.get(
        ref_final_key=ref_final_key,
        build_key=build_key,
    )
    with c.output.all_gc.open() as f:
        return json.load(f)


def all_gc(ref_final_key, build_key):
    # guard to prevent checkpoint from firing if we don't need GC
    bd = config.to_build_data_full(ref_final_key, build_key)
    if not bd.want_gc:
        return None

    res = gc_inputs(ref_final_key, build_key)
    return config.all_gc(
        ref_final_key,
        build_key,
        [Path(p) for p in res["gc_ranges"]],
        [Path(p) for p in [*res["other_extremes"], res["widest_extreme"]]],
        Path(rules.gc_readme.output[0]),
    )


rule find_gc_content:
    input:
        ref=rules.filter_sort_ref.output["fa"],
        genome=rules.filter_sort_ref.output["genome"],
        valid=rules.get_valid_regions.output.auto,
    output:
        gc.inter.postsort.data / "gc{frac}.bed.gz",
    params:
        args=lambda w: seqtk_args(config, w["ref_final_key"], w["build_key"], w["frac"]),
    conda:
        "../envs/seqtk.yml"
    wildcard_constraints:
        frac="\d+",
    shell:
        """
        seqtk gc {params.args} -l 100 {input.ref} | \
        cut -f1-3 | \
        slopBed -i stdin -g {input.genome} -b 50 | \
        mergeBed -i stdin | \
        intersectBed -a stdin -b {input.valid} -sorted -g {input.genome} | \
        bgzip -c \
        > {output}
        """


use rule find_gc_content as find_gc_content_final with:
    output:
        gc.final("gc{frac}_slop50"),


def range_inputs(wildcards):
    def _expand(p, frac):
        return expand(p, allow_missing=True, frac=frac)

    def expand_inter(frac):
        return _expand(rules.find_gc_content.output, frac)

    def expand_final(frac):
        return _expand(rules.find_gc_content_final.output, frac)

    lowest, lower, higher, highest = range_bounds(
        config,
        wildcards["ref_final_key"],
        wildcards["build_key"],
    )

    return {
        "low": expand_final(lowest) + expand_inter(lower),
        "high": expand_inter(higher) + expand_final(highest),
    }


checkpoint intersect_gc_ranges:
    input:
        unpack(range_inputs),
        genome=rules.filter_sort_ref.output["genome"],
        valid=rules.get_valid_regions.output.auto,
    # 'widest' is a symlink to the widest extreme gc range bed which should
    # always be present if this rule is run. This minor complexifier allows
    # other rules to refer to this particular file without invoking the
    # checkpoint and possibly causing a global nuclear war.
    output:
        all_gc=gc.inter.postsort.data / "intersect_output.json",
        widest_extreme=gc.inter.postsort.data / "widest.bed.gz",
    log:
        gc.inter.postsort.log / "intersect_ranges.txt",
    conda:
        "../envs/bedtools.yml"
    # hack together a format pattern that will be used for output
    params:
        path_pattern=lambda w: expand(
            gc.final("{{}}"),
            ref_final_key=w.ref_final_key,
            build_key=w.build_key,
        )[0],
    script:
        "../scripts/python/bedtools/gc/intersect_ranges.py"


rule gc_readme:
    input:
        common="workflow/templates/common.j2",
        description="workflow/templates/gc_description.j2",
        methods="workflow/templates/gc_methods.j2",
        bedtools_env="workflow/envs/bedtools.yml",
        seqtk_env="workflow/envs/seqtk.yml",
    params:
        paths=lambda w: all_gc(w["ref_final_key"], w["build_key"]),
    output:
        gc.readme,
    conda:
        "../envs/templates.yml"
    localrule: True
    script:
        "../scripts/python/templates/format_readme/format_gc.py"
