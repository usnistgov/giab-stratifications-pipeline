from os.path import dirname, basename
from more_itertools import unique_everseen, unzip
from os import scandir
from common.config import CoreLevel
from common.functional import DesignError

post_inter_dir = config.intermediate_build_dir / "postprocess"
post_log_dir = config.log_build_dir / "postprocess"
post_bench_dir = config.bench_build_dir / "postprocess"
validation_dir = config.final_root_dir / "validation"


# split these functions apart so I can use non-checkpoint paths without angering
# the snakemake gods
ALL_NON_CP_TARGETS = [
    all_diploid,
    all_functional,
    all_otherdifficult,
    all_low_complexity,
    all_segdups,
    all_telomeres,
    all_union,
    all_xy,
] + [
    lambda rk, bk, lk=lk: all_misc(lk, rk, bk) for lk in config.other_level_keys
]  # semi-official hack to get lambda to make closure around lk (if lk=lk isn't
# present then each lambda will only use the last value of lk in the loop,
# see https://docs.python.org/3.4/faq/programming.html#id10)

ALL_TARGETS = [
    all_gc,
    all_mappability,
    *ALL_NON_CP_TARGETS,
]


def all_strat_targets(wildcards):
    rk = wildcards["ref_final_key"]
    bk = wildcards["build_key"]
    return [t for f in ALL_TARGETS if (r := f(rk, bk)) is not None for t in r.all_final]


def all_readme_targets(wildcards):
    rk = wildcards["ref_final_key"]
    bk = wildcards["build_key"]
    return [r.readme for f in ALL_TARGETS if (r := f(rk, bk)) is not None]


rule list_all_strats:
    input:
        all_strat_targets,
    output:
        post_inter_dir / "all_strats.txt",
    localrule: True
    script:
        "../scripts/python/bedtools/postprocess/list_strats.py"


rule list_all_bb_strats:
    input:
        bed2bb=rules.build_kent.output.bed2bb,
        targets=all_strat_targets,
        genome=rules.filter_sort_ref.output.genome,
    output:
        post_inter_dir / "all_bb_strats.txt",
    log:
        post_log_dir / "all_bb_strats.log",
    script:
        "../scripts/python/bedtools/postprocess/list_convert_bb.py"


rule generate_md5sums:
    input:
        rules.list_all_strats.output,
    output:
        config.final_build_dir / "{ref_final_key}-genome-stratifications-md5s.txt",
    conda:
        "../envs/bedtools.yml"
    localrule: True
    script:
        "../scripts/python/bedtools/postprocess/list_md5.py"


use rule generate_md5sums as generate_bb_md5sums with:
    input:
        rules.list_all_bb_strats.output,
    output:
        config.final_build_dir / "{ref_final_key}-genome-stratifications-bb-md5s.txt",
    localrule: True


rule generate_tsv_list:
    input:
        rules.list_all_strats.output,
    output:
        config.final_build_dir / "{ref_final_key}-all-stratifications.tsv",
    params:
        suffix=".bed.gz",
    localrule: True
    script:
        "../scripts/python/bedtools/postprocess/generate_tsv.py"


use rule generate_tsv_list as generate_bb_tsv_list with:
    input:
        rules.list_all_bb_strats.output,
    output:
        config.final_build_dir / "{ref_final_key}-all-bb-stratifications.tsv",
    params:
        suffix=".bb",
    localrule: True


rule unit_test_strats:
    input:
        strats=rules.list_all_strats.output[0],
        valid_auto=rules.get_valid_regions.output.auto,
        valid_parY=rules.get_valid_regions.output.parY,
        genome=rules.filter_sort_ref.output["genome"],
        strat_list=rules.generate_tsv_list.output[0],
        checksums=rules.generate_md5sums.output[0],
    output:
        touch(post_inter_dir / "unit_tests.done"),
    log:
        error=post_log_dir / "unit_tests_errors.log",
        failed=post_log_dir / "failed_tests.log",
    benchmark:
        post_bench_dir / "unit_test_strats.txt"
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/postprocess/run_unit_tests.py"


rule get_coverage_table:
    input:
        _test=rules.unit_test_strats.output,
        bedlist=rules.list_all_strats.output[0],
        valid=rules.get_valid_regions.output.auto,
        genome=rules.filter_sort_ref.output["genome"],
    # TODO don't hardcode in the future
    params:
        window_size=int(1e6),
    output:
        full=touch(post_inter_dir / "coverage_full.tsv.gz"),
        window=touch(post_inter_dir / "coverage_window.tsv.gz"),
    benchmark:
        post_bench_dir / "get_coverage_table.txt"
    script:
        "../scripts/python/bedtools/postprocess/get_coverage_table.py"


rule download_comparison_strat_tarball:
    output:
        directory(config.resources_dir / "comparisons" / "{compare_key}"),
    conda:
        "../envs/utils.yml"
    # ASSUME this is not None (we control for this when building the "all" rule)
    params:
        url=lambda w: config.comparison_strats[w.compare_key],
        tar=lambda w: f"/tmp/comparison_{w.compare_key}.tar.gz",
    localrule: True
    shell:
        """
        curl -fsSq -o {params.tar} {params.url} && \
        mkdir {output} && \
        tar xzf {params.tar} --directory {output} --strip-components=1
        """


rule compare_strats:
    input:
        # ensure tests are run before this rule
        _test=rules.unit_test_strats.output,
        # NOTE: in the case id dip2 configurations, each comparison key will
        # correspond to one half of the diploid dataset (each tarball is one
        # haplotype)
        old=lambda w: expand(
            rules.download_comparison_strat_tarball.output,
            compare_key=config.compare_key(w.ref_final_key, w.build_key),
        )[0],
        # use this to target a specific rule to satisfy the snakemake scheduler,
        # the thing I actually need here is the parent directory
        new_list=rules.generate_tsv_list.output[0],
    output:
        validation_dir / "{ref_final_key}@{build_key}" / "diagnostics.tsv",
    log:
        post_log_dir / "comparison.log",
    # ASSUME this env will have pandas and bedtools (only real deps for this)
    conda:
        "../envs/bedtools.yml"
    threads: 8
    script:
        "../scripts/python/bedtools/postprocess/diff_previous_strats.py"


rule all_comparisons:
    input:
        [
            expand(rules.compare_strats.output, ref_final_key=rk, build_key=bk)[0]
            for rk, bk in zip(*config.all_full_build_keys)
            if config.compare_key(rk, bk) is not None
        ],
    localrule: True


rule make_coverage_plots:
    input:
        # this first input isn't actually used, but ensures the unit tests pass
        # before running the rmd script
        **{
            "_test": expand(
                rules.unit_test_strats.output,
                zip,
                ref_final_key=(t := config.all_full_build_keys)[0],
                build_key=t[1],
            ),
            "coverage": expand(
                rules.get_coverage_table.output.full,
                zip,
                ref_final_key=t[0],
                build_key=t[1],
            ),
        },
    output:
        validation_dir / "coverage_plots.html",
    conda:
        "../envs/rmarkdown.yml"
    params:
        core_levels=[c.value for c in CoreLevel],
        other_levels=[l.key for l in config.other_levels],
    script:
        "../scripts/rmarkdown/rmarkdown/coverage_plots.Rmd"


rule make_window_coverage_plots:
    input:
        _tests=rules.unit_test_strats.output,
        coverage=rules.get_coverage_table.output.window,
    output:
        validation_dir / "window_coverages" / "{ref_final_key}@{build_key}.html",
    conda:
        "../envs/rmarkdown.yml"
    params:
        core_levels=[c.value for c in CoreLevel],
        other_levels=[l.key for l in config.other_levels],
    script:
        "../scripts/rmarkdown/rmarkdown/window_coverage_plots.Rmd"


rule all_window_coverage_plots:
    input:
        [
            expand(
                rules.make_window_coverage_plots.output,
                ref_final_key=rk,
                build_key=bk,
            )[0]
            for rk, bk in zip(*config.all_full_build_keys)
        ],
    localrule: True


rule run_happy:
    input:
        _test=rules.unit_test_strats.output,
        ref=rules.filter_sort_ref.output["fa"],
        bench_vcf=lambda w: expand_final_to_src(rules.download_bench_vcf.output, w)[0],
        bench_bed=lambda w: read_checkpoint("normalize_bench_bed", w),
        query_vcf=lambda w: expand_final_to_src(rules.download_query_vcf.output, w)[0],
        strats=rules.generate_tsv_list.output[0],
    output:
        post_inter_dir / "happy" / "happy.extended.csv",
    params:
        prefix=lambda _, output: str(output[0]).replace(".extended.csv", ""),
    conda:
        "../envs/happy.yml"
    log:
        post_log_dir / "happy" / "happy.log",
    threads: 8
    benchmark:
        post_bench_dir / "run_happy.txt"
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_final_key, w.build_key, lambda m: m.runHappy
        ),
    script:
        "../scripts/python/bedtools/postprocess/run_happy.py"


rule summarize_happy:
    input:
        [
            expand(
                rules.run_happy.output,
                ref_final_key=rk,
                build_key=bk,
            )
            for rk, bk in zip(*config.all_build_keys)
            if config.to_build_data(rk, bk).have_benchmark
        ],
    output:
        validation_dir / "benchmark_summary.html",
    params:
        subsets=config.benchmark_subsets,
    conda:
        "../envs/rmarkdown.yml"
    script:
        "../scripts/rmarkdown/rmarkdown/benchmark.Rmd"


rule copy_strat_background:
    input:
        "workflow/files/BACKGROUND.md",
    output:
        config.final_build_dir / "{ref_final_key}-BACKGROUND.md",
    localrule: True
    shell:
        """
        cp {input} {output}
        """


rule copy_validation_README:
    input:
        "workflow/files/README_validation.md",
    output:
        validation_dir / "README.md",
    localrule: True
    shell:
        """
        cp {input} {output}
        """


rule build_strat_README:
    input:
        readme="workflow/templates/main.j2",
        ref=lambda w: expand_final_to_src(rules.download_ref.output, w)[0],
    output:
        config.final_build_dir / "{ref_final_key}-README.md",
    params:
        segdups=lambda w: all_segdups(w["ref_final_key"], w["build_key"]),
        functional=lambda w: all_functional(w["ref_final_key"], w["build_key"]),
        otherdiff=lambda w: all_otherdifficult(w["ref_final_key"], w["build_key"]),
        union=lambda w: all_union(w["ref_final_key"], w["build_key"]),
        xy=lambda w: all_xy(w["ref_final_key"], w["build_key"]),
    conda:
        "../envs/templates.yml"
    localrule: True
    script:
        "../scripts/python/templates/format_readme/format_main.py"


# include bed/bb in both of these since tar will complain of the bed->bb step is
# changing the directory contents while it is running
def bb_inputs(wildcards):
    rk = wildcards["ref_final_key"]
    bk = wildcards["build_key"]
    t = config.to_build_data_full(rk, bk).want_bb
    return {
        "_bb_checksums": (rules.generate_bb_md5sums.output if t else []),
        "_all_bb": rules.generate_bb_tsv_list.output if t else [],
    }


rule generate_tarballs:
    input:
        unpack(bb_inputs),
        all_bed=rules.generate_tsv_list.output,
        background=rules.copy_strat_background.output,
        readme=rules.build_strat_README.output,
        _bed_checksums=rules.generate_md5sums.output,
        _strat_readmes=all_readme_targets,
    output:
        config.final_root_dir
        / "genome-stratifications-{ref_final_key}@{build_key}.tar.gz",
    params:
        parent=lambda _, input: Path(input.all_bed[0]).parent.parent,
        target=lambda _, input: Path(input.all_bed[0]).parent.name,
    localrule: True
    shell:
        """
        tar czf {output} -C {params.parent} \
          --exclude=*.bb \
          --exclude=*-all-bb-stratifications.tsv \
          --exclude=*-genome-stratifications-bb-md5s.txt \
          {params.target} 
        """


rule generate_bb_tarballs:
    input:
        # since double-star apparently doesn't apply to functions
        unpack(bb_inputs),
        **rules.generate_tarballs.input,
    output:
        config.final_root_dir
        / "genome-stratifications-bb-{ref_final_key}@{build_key}.tar.gz",
    params:
        parent=lambda _, input: Path(input.all_bed[0]).parent.parent,
        target=lambda _, input: Path(input.all_bed[0]).parent.name,
    localrule: True
    shell:
        """
        tar czf {output} -C {params.parent} \
          --exclude=*.bed.gz \
          --exclude=*-all-stratifications.tsv \
          --exclude=*-genome-stratifications-md5s.txt \
          {params.target} 
        """


rule format_global_README:
    input:
        readme="workflow/templates/README.j2",
    output:
        config.final_root_dir / "README.md",
    conda:
        "../envs/templates.yml"
    localrule: True
    script:
        "../scripts/python/templates/format_readme/format_global.py"


rule checksum_everything:
    input:
        rules.format_global_README.output,
        rules.copy_validation_README.output,
        rules.make_coverage_plots.output,
        rules.summarize_happy.output,
        rules.all_comparisons.input,
        rules.all_window_coverage_plots.input,
        [
            expand(rules.generate_tarballs.output, ref_final_key=rk, build_key=bk)
            for rk, bk in zip(*config.all_full_build_keys)
        ],
        [
            expand(rules.generate_bb_tarballs.output, ref_final_key=rk, build_key=bk)
            for rk, bk in zip(*config.all_full_build_keys)
            if config.to_build_data_full(rk, bk).want_bb
        ],
    output:
        config.final_root_dir / "genome-stratifications-md5s.txt",
    params:
        root=config.final_root_dir,
    localrule: True
    shell:
        """
        find {params.root} -type f -exec md5sum {{}} + | \
        sed 's|{params.root}|\.|' | \
        grep -v "\./genome-stratifications-md5s\.txt" | \
        sort -k2,2 > {output} && \
        cd {params.root} && \
        md5sum -c --quiet genome-stratifications-md5s.txt
        """
