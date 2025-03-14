from typing import Any
import common.config as cfg
import template_utils as tu


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    readme = cfg.smk_to_input_name(smk, "readme")

    txt = tu.load_template_path(readme).render(
        version=sconf.version,
    )

    out = cfg.smk_to_output(smk)

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore
