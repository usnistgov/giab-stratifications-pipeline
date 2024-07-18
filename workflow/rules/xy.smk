from common.config import CoreLevel, MutualPathPair

xy = config.to_bed_dirs(CoreLevel.XY)


def all_xy(ref_final_key, build_key):
    return config.all_xy(
        ref_final_key,
        build_key,
        Path(rules.download_genome_features_bed.output[0]),
        Path(rules.filter_XTR_features.output[0]),
        Path(rules.filter_ampliconic_features.output[0]),
        MutualPathPair(
            Path(rules.write_PAR_final.output[0]),
            Path(rules.invert_PAR.output[0]),
        ),
        Path(rules.filter_autosomes.output[0]),
        Path(rules.xy_readme.output[0]),
    )


use rule download_gaps as download_genome_features_bed with:
    output:
        xy.src.data / "genome_features_{sex_chr}.bed.gz",
    log:
        xy.src.log / "genome_features_{sex_chr}.log",
    params:
        src=lambda w: (
            config.refsrckey_to_y_features_src
            if w.sex_chr == "Y"
            else config.refsrckey_to_x_features_src
        )(w.ref_src_key),
    localrule: True


rule write_PAR_final:
    input:
        bed=rules.write_PAR_intermediate.output,
        valid=rules.get_valid_regions.output.parY,
        genome=rules.filter_sort_ref.output["genome"],
    output:
        xy.final("chr{sex_chr}_PAR"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        intersectBed -a {input.bed} -b {input.valid} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


feature_inputs = {
    "bed": lambda w: expand_final_to_src(rules.download_genome_features_bed.output, w)[
        0
    ],
    "valid": rules.get_valid_regions.output.parY,
    "genome": rules.filter_sort_ref.output["genome"],
}


rule filter_XTR_features:
    input:
        **feature_inputs,
    output:
        xy.final("chr{sex_chr}_XTR"),
    log:
        xy.inter.postsort.log / "{sex_chr}_filter_XTR_features.txt",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/xy/filter_sort_xtr.py"


rule filter_ampliconic_features:
    input:
        **feature_inputs,
    output:
        xy.final("chr{sex_chr}_ampliconic"),
    log:
        xy.inter.postsort.log / "{sex_chr}_filter_ampliconic_features.txt",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/xy/filter_sort_ampliconic.py"


rule invert_PAR:
    input:
        bed=rules.write_PAR_final.output,
        genome=rules.filter_sort_ref.output["genome"],
        valid=rules.get_valid_regions.output.parY,
    output:
        xy.final("chr{sex_chr}_nonPAR"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        complementBed -i {input.bed} -g {input.genome} | \
        grep {wildcards.sex_chr} | \
        intersectBed -a stdin -b {input.valid} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


rule filter_autosomes:
    input:
        rules.get_valid_regions.output.auto,
    output:
        xy.final("AllAutosomes"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        gunzip -c {input} | \
        grep -v \"X\|Y\" | \
        bgzip -c > {output}
        """


rule xy_readme:
    input:
        common="workflow/templates/common.j2",
        description="workflow/templates/xy_description.j2",
        methods="workflow/templates/xy_methods.j2",
        bedtools_env="workflow/envs/bedtools.yml",
        _sources=lambda w: all_xy(w["ref_final_key"], w["build_key"]).all_sources,
    params:
        paths=lambda w: all_xy(w["ref_final_key"], w["build_key"]),
    output:
        xy.readme,
    conda:
        "../envs/templates.yml"
    localrule: True
    script:
        "../scripts/python/templates/format_readme/format_xy.py"
