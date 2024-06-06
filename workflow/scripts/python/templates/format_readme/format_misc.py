import jinja2 as j2
from typing import Any
import common.config as cfg
from common.functional import DesignError, fmap_maybe_def, fmap_maybe
import template_utils as tu


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    rfk = cfg.wc_to_reffinalkey(ws)
    bk = cfg.wc_to_buildkey(ws)

    paths = smk.params["paths"]

    if not isinstance(paths, cfg.MiscPaths):
        raise DesignError()

    lk = paths.desc.key

    desc = tu.fmt_other_descriptions(sconf, rfk, bk, lk, paths.paths)
    src = tu.fmt_other_srcs(sconf, rfk, bk, lk, paths.paths)

    bedtools_env_path = cfg.smk_to_input_name(smk, "bedtools_env")
    bedtools_deps = tu.env_dependencies(bedtools_env_path, {"bedtools", "samtools"})

    def render_description(t: j2.Template) -> str:
        return t.render(desc=desc)

    def render_methods(t: j2.Template) -> str:
        return t.render(
            src=src,
            deps=bedtools_deps,
        )

    txt = tu.render_readme(
        smk,
        render_description,
        render_methods,
        paths.desc.desc,
        paths.desc.key,
        sconf.refkey_haplotypes(rfk),
    )

    out = cfg.smk_to_output(smk)

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore
