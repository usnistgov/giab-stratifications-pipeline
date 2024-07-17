import gzip
from typing import Any, NamedTuple
from pathlib import Path
from os.path import dirname, basename
from os import scandir
import common.config as cfg
from common.io import is_bgzip, gunzip
from common.bed import InternalChrIndex, subtractBed, ChrName
import subprocess as sp

RevMapper = dict[ChrName, InternalChrIndex]


class PathsEnv(NamedTuple):
    auto: Path
    parY: Path
    genome: Path
    error_log: Path


def test_bgzip(strat_file: Path) -> list[str]:
    """Each bed file should be bgzip'ed (not just gzip'ed)."""
    return [] if is_bgzip(strat_file) else ["is not bgzip file"]


def test_bed_format(strat_file: Path, reverse_map: RevMapper) -> str | None:
    with gzip.open(strat_file) as f:
        prevChrom = ""
        prevChromIdx = -1
        prevStart = -1
        prevEnd = -1
        for i in f:
            # should have three columns separated by tabs
            cells = i.split(b"\t")
            if len(cells) != 3:
                return "bed file has wrong number of columns"

            # types of each column should be str/int/int
            chrom = cells[0].decode()
            try:
                start = int(cells[1])
                end = int(cells[2])
            except ValueError as e:
                return str(e)

            # chrom should be valid
            try:
                chromIdx = reverse_map[ChrName(chrom)]
            except KeyError:
                return f"invalid chr: {chrom}"

            # chrom column should be sorted in ascending order
            if not prevChromIdx <= chromIdx:
                return "chrom column not sorted"

            # end should be greater than start
            if not start < end:
                return f"invalid region: {chrom} {start} {end}"

            # regions should be separated by at least one bp
            if not ((prevEnd < start) or (prevChromIdx < chromIdx)):
                return (
                    "non-disjoint regions: "
                    f"{prevChrom} {prevStart} {prevEnd}"
                    f"and {chrom} {start} {end}"
                )

            prevChrom = chrom
            prevChromIdx = chromIdx
            prevStart = start
            prevEnd = end

    # we made it, this file is legit :)
    return None


def test_bed_valid_regions(strat_file: Path, env: PathsEnv) -> bool:
    # ignore the gaps stratification since this is out of the valid regions by
    # definition
    if "gaps_slop15kb" in basename(strat_file):
        return True
    else:
        isXY = basename(dirname(strat_file)) == "XY"
        valid_path = env.parY if isXY else env.auto
        # gunzip needed here because subtract bed will output gibberish if we
        # give it an empty gzip file
        p1, o = gunzip(strat_file)
        p2, _ = subtractBed(o, valid_path, env.genome)
        out, err = p2.communicate()
        p1.wait()
        with open(env.error_log, "a") as f:
            haserror = False
            if p1.returncode != 0:
                f.write("gunzip error\n")
                haserror = True
            if p2.returncode != 0:
                f.write(err.decode() + "\n")
                haserror = True
            if haserror:
                exit(1)
        return len(out) == 0


def test_bed(strat_file: Path, reverse_map: RevMapper, env: PathsEnv) -> list[str]:
    if (format_error := test_bed_format(strat_file, reverse_map)) is not None:
        return [format_error]
    else:
        if test_bed_valid_regions(strat_file, env):
            return []
        else:
            return ["has gaps"]


def test_checksums(checksums: Path) -> list[str]:
    """The md5sum utility should pass when invoked."""
    res = sp.run(
        ["md5sum", "-c", "--strict", "--quiet", checksums.name],
        cwd=checksums.parent,
        capture_output=True,
        text=True,
    )
    errors = [] if (s := res.stdout) == "" else s.strip().split("\n")
    return [f"checksum error: {e}" for e in errors]


def test_tsv_list(tsv_list: Path) -> list[str]:
    """The final stratifications list should match final beds exactly.

    We are running our hotel on very tight margins; no extra or missing beds
    allowed.
    """

    def strat_set(root: Path, sub: Path) -> set[str]:
        if root.is_dir():
            deeper = (Path(p.path) for p in scandir(root))
            return {str(q) for d in deeper for q in strat_set(d, sub / d.name)}
        elif root.is_file() and root.name.endswith(".bed.gz"):
            return {str(sub)}
        else:
            return set()

    current = strat_set(tsv_list.parent, Path("./"))
    with open(tsv_list, "r") as f:
        listed = {line.strip().split("\t")[1] for line in f}
        missing = [f"not in final directory: {p}" for p in listed - current]
        extra = [f"not in final list: {p}" for p in current - listed]
        return missing + extra


def run_all_tests(
    strat_file: Path,
    reverse_map: RevMapper,
    env: PathsEnv,
) -> list[tuple[Path, str]]:
    return [
        (strat_file, msg)
        for msg in test_bgzip(strat_file) + test_bed(strat_file, reverse_map, env)
    ]


def strat_files(path: Path) -> list[Path]:
    with open(path, "r") as f:
        return [Path(s.strip()) for s in f]


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    # check global stuff first (since this is faster)

    strat_list_path = cfg.smk_to_input_name(smk, "strat_list")
    checksums_path = cfg.smk_to_input_name(smk, "checksums")
    failed_tests_path = cfg.smk_to_log_name(smk, "failed")

    global_failures = test_checksums(checksums_path) + test_tsv_list(strat_list_path)

    with open(failed_tests_path, "w") as f:
        for j in global_failures:
            f.write(j + "\n")

    if len(global_failures) > 0:
        exit(1)

    # check individual stratification files

    env = PathsEnv(
        auto=cfg.smk_to_input_name(smk, "valid_auto"),
        parY=cfg.smk_to_input_name(smk, "valid_parY"),
        genome=cfg.smk_to_input_name(smk, "genome"),
        error_log=cfg.smk_to_log_name(smk, "error"),
    )

    strats_path = cfg.smk_to_input_name(smk, "strats")

    fm = sconf.buildkey_to_ref_mappers(
        cfg.wc_to_reffinalkey(ws),
        cfg.wc_to_buildkey(ws),
    )[1]
    reverse_map = {v: k for k, v in fm.items()}

    strat_failures = [
        res
        for p in strat_files(strats_path)
        for res in run_all_tests(p, reverse_map, env)
    ]

    with open(failed_tests_path, "a") as f:
        for i in strat_failures:
            f.write("%s: %s" % i + "\n")

    if len(strat_failures) > 0:
        exit(1)


main(snakemake, snakemake.config)  # type: ignore
