from typing import Any
import common.config as cfg
from helpers import filter_xy_feature


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    filter_xy_feature(smk, sconf, lambda fs: fs.xtr)


main(snakemake, snakemake.config)  # type: ignore
