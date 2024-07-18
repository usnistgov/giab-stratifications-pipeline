import gzip
from math import floor
import re
from pathlib import Path
from typing import Any
import common.config as cfg
from common.functional import not_none_unsafe, DesignError
from common.bed import ChrName, InternalChrIndex

ChrMapper = dict[ChrName, tuple[cfg.ShortChrName, cfg.Haplotype, InternalChrIndex]]


def get_chr_paired_mapper(
    sconf: cfg.GiabStrats,
    rk: cfg.RefKeyFullS,
    bk: cfg.BuildKey,
) -> ChrMapper:
    cis = sconf.to_build_data(cfg.strip_full_refkey(rk), bk).build_chrs
    xs = sconf.with_ref_data_full(
        rk,
        lambda rd: rd.ref.chr_pattern.to_chr_data(cis, cfg.Haplotype.PAT),
        lambda rd: rd.ref.chr_pattern.to_chr_data(cis),
        lambda hap, rd: rd.ref.chr_pattern.choose(hap).to_chr_data(cis, hap),
    )
    return {x.name: (x.shortname, x.haplotype, x.idx) for x in xs}


def get_distribution(
    ipath: Path,
    genome: dict[ChrName, int],
    k: int,
) -> dict[ChrName, list[tuple[float, float, float]]]:
    acc: dict[ChrName, dict[int, int]] = {}

    with gzip.open(ipath, "rb") as f:
        for x in f:
            s = x.split(b"\t")
            chrom = ChrName(s[0].decode())
            start = int(s[1])
            end = int(s[2])
            cur = start
            if chrom not in acc:
                acc[chrom] = {}
            while cur < end:
                idx = floor(cur / k) * k
                if idx not in acc[chrom]:
                    acc[chrom][idx] = 0
                next_idx = idx + k
                acc[chrom][idx] += (next_idx if end >= next_idx else end) - cur
                cur = next_idx

    return {
        chrom: [(i / (t := genome[chrom]), c / k, k / t) for i, c in ks.items()]
        for chrom, ks in acc.items()
    }


def read_genome(path: Path) -> dict[ChrName, int]:
    acc = {}
    with open(path, "r") as f:
        for x in f:
            s = x.split("\t")
            acc[ChrName(s[0])] = int(s[1])
    return acc


def sum_bed_file(path: Path) -> dict[ChrName, int]:
    # ASSUME bed file has been tested for validity
    acc = {}
    with gzip.open(path, "r") as f:
        for i in f:
            cells = i.split(b"\t")
            chrom = ChrName(cells[0].decode())
            if chrom not in acc:
                acc[chrom] = 0
            acc[chrom] += int(cells[2]) - int(cells[1])

    return acc


def bedpath_to_strat_names(bp: Path) -> tuple[str, str]:
    strat_name = not_none_unsafe(
        re.match("^[^_]+_(.+).bed.gz", bp.name),
        lambda x: x[1],
    )

    return (bp.parent.name, strat_name)


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    valid_path = cfg.smk_to_input_name(smk, "valid")
    genome_path = cfg.smk_to_input_name(smk, "genome")
    bedlist = cfg.smk_to_input_name(smk, "bedlist")
    full_out = cfg.smk_to_output_name(smk, "full")
    window_out = cfg.smk_to_output_name(smk, "window")

    window_size = smk.params["window_size"]

    if not isinstance(window_size, int):
        raise DesignError("window size is not an int")

    rk = cfg.wc_to_reffinalkey(ws)
    bk = cfg.wc_to_buildkey(ws)

    chr_hap_mapper = get_chr_paired_mapper(sconf, rk, bk)
    valid_sum = sum_bed_file(valid_path)
    genome = read_genome(genome_path)

    with open(bedlist, "r") as i, gzip.open(full_out, "w") as of, gzip.open(
        window_out, "w"
    ) as ow:
        for bedpath in i:
            bp = Path(bedpath.strip())
            level_name, strat_name = bedpath_to_strat_names(bp)
            # TODO I could make these names a bit shorter by only appending
            # the haplotype to dip2 references since that's the only place
            # we need to do this
            for chrom, bedtotal in sum_bed_file(bp).items():
                shortChromName, hap, idx = chr_hap_mapper[chrom]
                fraction = bedtotal / valid_sum[chrom]
                fullkey = f"{rk}@{bk}-{hap.value + 1}"
                newline = "\t".join(
                    [
                        fullkey,
                        level_name,
                        strat_name,
                        shortChromName,
                        str(idx - hap.value * 24),
                        str(fraction),
                    ]
                )
                of.write((newline + "\n").encode())

            for chrom, rest in get_distribution(bp, genome, window_size).items():
                shortChromName, hap, idx = chr_hap_mapper[chrom]
                for pos, coverage, width in rest:
                    newline = "\t".join(
                        [
                            str(hap.value + 1),
                            level_name,
                            strat_name,
                            shortChromName,
                            str(idx - hap.value * 24),
                            str(pos),
                            str(coverage),
                            str(width),
                        ]
                    )
                    ow.write((newline + "\n").encode())


main(snakemake, snakemake.config)  # type: ignore
