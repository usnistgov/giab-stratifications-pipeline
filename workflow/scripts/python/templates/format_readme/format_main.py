from typing import Any, TypeVar, Type
import common.config as cfg
from common.functional import DesignError, fmap_maybe, fmap_maybe_def
import template_utils as tu

X = TypeVar("X")

GAPS_DESC = "gaps in the reference"
KIR_DESC = (
    "the killer-cell immunoglobulin-like receptor (KIR) region (highly polymorphic)"
)
MHC_DESC = "the major histocompatibility complex (MHC) region (highly polymorphic)"
VDJ_DESC = "the VDJ loci (highly polymorphic with somatic rearrangements)"


def if_source(s: str, x: cfg.SourceOutputPaths | None) -> str | None:
    return s if fmap_maybe(lambda y: y.source, x) is not None else None


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    def get_paths(name: str, t: Type[X]) -> X | None:
        x = smk.params[name]
        if x is not None and not isinstance(x, t):
            raise DesignError()
        return x

    def have_final_outputs(name: str, t: Type[cfg._HasFinalBeds]) -> bool:
        x = get_paths(name, t)
        return fmap_maybe_def(False, lambda x: len(x.all_final) > 0, x)

    ws: dict[str, str] = smk.wildcards

    readme = cfg.smk_to_input_name(smk, "readme")
    refpath = cfg.smk_to_input_name(smk, "ref")

    rfk = cfg.wc_to_reffinalkey(ws)
    bk = cfg.wc_to_buildkey(ws)

    bd = sconf.to_build_data_full(rfk, bk)

    otherdiff = fmap_maybe(
        lambda o: (
            [
                cfg.readme_fill(x)
                for x in [
                    if_source(GAPS_DESC, o.gaps),
                    if_source(KIR_DESC, o.kir),
                    if_source(MHC_DESC, o.mhc),
                    if_source(VDJ_DESC, o.vdj),
                ]
                if x is not None
            ],
            len(o.other) > 0,
        ),
        get_paths("otherdiff", cfg.OtherDifficultPaths),
    )

    # TODO this seems like something I should check everywhere
    # ultra paranoid, make sure the keys we want are actually populated
    other_keys = [lk for lk, lv in bd.build.other_strats.items() if len(lv) > 0]

    ref_src = sconf.with_ref_data_full(
        rfk,
        lambda rd: rd.ref.src.elem,
        lambda rd: rd.ref.src.elem,
        lambda hap, rd: rd.ref.src.choose(hap),
    )

    ref_src_txt = cfg.readme_fill(
        cfg.format_src(ref_src.documentation, refpath, "The reference FASTA")
    )

    compare_txt = fmap_maybe(
        lambda c: cfg.readme_fill(
            "This version of the stratifications was compared against "
            f"a previous version obtained from {c}"
        ),
        (
            sconf.comparison_strats[ck]
            if (ck := sconf.compare_key(rfk, bk)) is not None
            and ck in sconf.comparison_strats
            else None
        ),
    )

    txt = tu.load_template_path(readme).render(
        refname=rfk,
        haplotypes=[h.name for h in sconf.refkey_haplotypes(rfk)],
        have_diploid=bd.want_diploid,
        have_low_complexity=bd.want_low_complexity,
        have_segdups=have_final_outputs("segdups", cfg.SegdupPaths),
        have_sex=have_final_outputs("xy", cfg.SexPaths),
        have_functional=have_final_outputs("functional", cfg.FunctionalPaths),
        have_gc=bd.want_gc,
        have_mappability=bd.want_mappability,
        otherdiff=otherdiff,
        have_telomeres=bd.want_telomeres,
        have_union=have_final_outputs("union", cfg.UnionPaths),
        other_levels=[o for o in sconf.other_levels if o.key in other_keys],
        ref_src=ref_src_txt,
        pipeline_repo=sconf.docs.pipeline_repo,
        config_repo=sconf.docs.config_repo,
        compare_txt=compare_txt,
    )

    out = cfg.smk_to_output(smk)

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore
