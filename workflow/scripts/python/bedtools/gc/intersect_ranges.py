from pathlib import Path
from typing import Any, NamedTuple, Callable
import common.config as cfg
from common.functional import DesignError
from common.io import bgzip_file, gunzip, check_processes
from common.bed import (
    complementBed,
    subtractBed,
    intersectBed,
    multiIntersectBed,
    mergeBed,
)
import json


class GCInput(NamedTuple):
    bed: Path
    fraction: int
    is_range_bound: bool


def write_simple_range_beds(
    final_path: Callable[[str], Path],
    gs: list[GCInput],
    genome: Path,
    is_low: bool,
    log: Path,
) -> list[str]:
    def fmt_out(bigger_frac: int, smaller_frac: int) -> Path:
        lower_frac, upper_frac = (
            (smaller_frac, bigger_frac) if is_low else (bigger_frac, smaller_frac)
        )
        return final_path(f"gc{lower_frac}to{upper_frac}_slop50")

    pairs = zip(gs[1:], gs[:-1]) if is_low else zip(gs[:-1], gs[1:])
    torun = [
        (fmt_out(bigger.fraction, smaller.fraction), bigger, smaller)
        for bigger, smaller in pairs
    ]
    for out, bigger, smaller in torun:
        p1, o1 = gunzip(bigger.bed)
        p2, o2 = subtractBed(o1, smaller.bed, genome)
        bgzip_file(o2, out)
        check_processes([p1, p2], log)
    return [str(t[0]) for t in torun]


def write_middle_range_bed(
    final_path: Callable[[str], Path],
    lower: GCInput,
    upper: GCInput,
    genome: Path,
    valid_regions: Path,
    log: Path,
) -> str:
    out = final_path(f"gc{lower.fraction}to{upper.fraction}_slop50")
    with open(upper.bed, "rb") as i:
        p1, o0 = complementBed(i, genome)
        p2, o1 = subtractBed(o0, lower.bed, genome)
        p3, o2 = intersectBed(o1, valid_regions, genome)
        bgzip_file(o2, out)
        check_processes([p1, p2, p3], log)
    return str(out)


# ASSUME low/high are non-empty and the same length subset to range bounds
def write_intersected_range_beds(
    final_path: Callable[[str], Path],
    low: list[GCInput],
    high: list[GCInput],
    log: Path,
) -> tuple[Path, list[str]]:
    pairs = zip(
        [x for x in low if x.is_range_bound],
        [x for x in reversed(high) if x.is_range_bound],
    )
    torun = [
        (final_path(f"gclt{b1.fraction}orgt{b2.fraction}_slop50"), b1, b2)
        for b1, b2 in pairs
    ]
    for bed_out, b1, b2 in torun:
        p1, o1 = multiIntersectBed([b1.bed, b2.bed])
        p2, o2 = mergeBed(o1, [])
        bgzip_file(o2, bed_out)
        check_processes([p1, p2], log)
    out = [t[0] for t in torun]
    return (out[0], [str(p) for p in out[1:]])


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    path_pattern = cfg.smk_to_param_str(smk, "path_pattern")

    low_paths = cfg.smk_to_inputs_name(smk, "low")
    high_paths = cfg.smk_to_inputs_name(smk, "high")
    genome_path = cfg.smk_to_input_name(smk, "genome")
    valid_path = cfg.smk_to_input_name(smk, "valid")
    log_path = cfg.smk_to_log(smk)
    out_widest = cfg.smk_to_output_name(smk, "widest_extreme")
    out_all = cfg.smk_to_output_name(smk, "all_gc")

    rfk = cfg.wc_to_reffinalkey(ws)
    bk = cfg.wc_to_buildkey(ws)
    bd = sconf.to_build_data(cfg.strip_full_refkey(rfk), bk)
    gps = bd.build.include.gc
    if gps is None:
        raise DesignError
    # ASSUME both of these input lists are sorted by GC fraction
    low = [GCInput(p, f, r) for p, (f, r) in zip(low_paths, gps.low)]
    high = [GCInput(p, f, r) for p, (f, r) in zip(high_paths, gps.high)]

    def final_path(name: str) -> Path:
        p = Path(path_pattern.format(name))
        p.parent.mkdir(exist_ok=True, parents=True)
        return p

    low_strats = write_simple_range_beds(
        final_path,
        low,
        genome_path,
        True,
        log_path,
    )
    high_strats = write_simple_range_beds(
        final_path,
        high,
        genome_path,
        False,
        log_path,
    )
    range_strat = write_middle_range_bed(
        final_path,
        low[-1],
        high[0],
        genome_path,
        valid_path,
        log_path,
    )
    widest_extreme, other_extremes = write_intersected_range_beds(
        final_path,
        low,
        high,
        log_path,
    )

    # ASSUME there is one extreme denoted in the config (otherwise we won't have
    # something to put in "widest extreme"
    with open(out_all, "w") as f:
        # put the first low and last high input here since these are already
        # in the final directory
        ranges: list[str] = [
            str(low[0].bed),
            *low_strats,
            range_strat,
            *high_strats,
            str(high[-1].bed),
        ]
        obj = {
            "gc_ranges": ranges,
            "widest_extreme": str(widest_extreme),
            "other_extremes": other_extremes,
        }
        json.dump(obj, f, indent=4)

    out_widest.symlink_to(widest_extreme.resolve())


main(snakemake, snakemake.config)  # type: ignore
