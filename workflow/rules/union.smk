from common.config import CoreLevel, MutualPathPair, strip_full_refkey

uni = config.to_bed_dirs(CoreLevel.UNION)


def all_union(ref_final_key, build_key):
    return config.all_union(
        ref_final_key,
        build_key,
        all_segdups(ref_final_key, build_key),
        Path(rules.merge_nonunique.output.all_lowmap),
        Path(rules.intersect_gc_ranges.output.widest_extreme),
        all_low_complexity(ref_final_key, build_key),
        all_xy(ref_final_key, build_key),
        MutualPathPair(
            Path(rules.intersect_segdup_and_map.output[0]),
            Path(rules.invert_segdup_and_map.output[0]),
        ),
        MutualPathPair(
            Path(rules.intersect_alldifficult.output[0]),
            Path(rules.invert_alldifficult.output[0]),
        ),
        Path(rules.union_readme.output[0]),
    )


rule intersect_segdup_and_map:
    input:
        lambda w: all_union(w["ref_final_key"], w["build_key"]).segdup_lowmap_inputs,
    output:
        uni.final("alllowmapandsegdupregions"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        multiIntersectBed -i {input} | \
        mergeBed -i stdin | \
        bgzip -c \
        > {output}
        """


use rule _invert_autosomal_regions as invert_segdup_and_map with:
    input:
        **invert_region_inputs(rules.intersect_segdup_and_map.output),
    output:
        uni.final("notinalllowmapandsegdupregions"),


# TODO include other difficult in here as well
use rule intersect_segdup_and_map as intersect_alldifficult with:
    input:
        lambda w: all_union(w["ref_final_key"], w["build_key"]).all_difficult_inputs,
    output:
        uni.final("alldifficultregions"),


use rule invert_segdup_and_map as invert_alldifficult with:
    input:
        **invert_region_inputs(rules.intersect_alldifficult.output),
    output:
        uni.final("notinalldifficultregions"),


rule union_readme:
    input:
        common="workflow/templates/common.j2",
        description="workflow/templates/union_description.j2",
        methods="workflow/templates/union_methods.j2",
        bedtools_env="workflow/envs/bedtools.yml",
        _inputs=lambda w: all_union(w["ref_final_key"], w["build_key"]).all_inputs,
    params:
        paths=lambda w: all_union(w["ref_final_key"], w["build_key"]),
    output:
        uni.readme,
    conda:
        "../envs/templates.yml"
    localrule: True
    script:
        "../scripts/python/templates/format_readme/format_union.py"
