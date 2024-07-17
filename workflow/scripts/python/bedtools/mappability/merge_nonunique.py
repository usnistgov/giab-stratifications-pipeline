import re
import pandas as pd
from typing import Any, IO
from pathlib import Path
import common.config as cfg
import subprocess as sp
from common.bed import (
    read_bed_default,
    filter_sort_bed,
    bed_to_stream,
    mergeBed,
    intersectBed,
    complementBed,
    multiIntersectBed,
)
from common.io import bgzip_file, check_processes
import json


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    inputs = cfg.smk_to_inputs_name(smk, "bed")
    genome = cfg.smk_to_input_name(smk, "genome")
    valid = cfg.smk_to_input_name(smk, "valid")
    log = cfg.smk_to_log(smk)
    out_single = cfg.smk_to_output_name(smk, "single_lowmap")
    out_all = cfg.smk_to_output_name(smk, "all_lowmap")
    path_pattern = cfg.smk_to_param_str(smk, "path_pattern")

    im, fm = sconf.buildkey_to_ref_mappers(
        cfg.wc_to_reffinalkey(ws),
        cfg.wc_to_buildkey(ws),
    )

    def final_path(name: str) -> Path:
        p = Path(path_pattern.format(name))
        p.parent.mkdir(exist_ok=True, parents=True)
        return p

    def to_single_output(basename: str) -> Path:
        m = re.match(".*(l\\d+_m\\d+_e\\d+).*", basename)
        assert m is not None, f"malformed mappability file name: {basename}"
        return final_path(f"nonunique_{m[1]}")

    def read_sort_bed(p: Path) -> pd.DataFrame:
        # Sort here because we can't assume wig2bed sorts its output. Also,
        # filtering is necessary because the output should have unplaced contigs
        # in it that we don't want.
        return filter_sort_bed(im, fm, read_bed_default(p))

    def merge_bed(bed: IO[bytes], out: Path) -> tuple[sp.Popen[bytes], sp.Popen[bytes]]:
        p1, o1 = mergeBed(bed, ["-d", "100"])
        p2, o2 = intersectBed(o1, valid, genome)
        bgzip_file(o2, out)
        return p1, p2

    def merge_single(i: Path, o: Path) -> None:
        bed = read_sort_bed(i)
        with bed_to_stream(bed) as s:
            p1, o1 = complementBed(s, genome)
            p2, p3 = merge_bed(o1, o)
            check_processes([p1, p2, p3], log)

    # all_lowmap = final_path("lowmappabilityall")
    single_lowmap = []

    # If there is only one input, merge this to make one "all_nonunique" bed.
    # Otherwise, merge each individual input and then combine these with
    # multi-intersect to make the "all_nonunique" bed.
    if len(inputs) == 1:
        merge_single(Path(inputs[0]), out_all)
    else:
        # first read/sort/write all single bed files
        single_lowmap = [to_single_output(Path(i).name) for i in inputs]
        for i, o in zip(inputs, single_lowmap):
            merge_single(i, o)
        # once all the single files are on disk, stream them together; this
        # allows us to avoid keeping multiple dataframes in memory at once
        p1, mi_out = multiIntersectBed(single_lowmap)
        p2, p3 = merge_bed(mi_out, out_all)
        check_processes([p1, p2, p3], log)

    with open(out_single, "w") as f:
        obj = [str(p) for p in single_lowmap]
        json.dump(obj, f)


main(snakemake, snakemake.config)  # type: ignore
