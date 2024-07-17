from typing import Any, Callable
import common.config as cfg
from common.bed import filter_sort_bed, bed_to_stream, intersectBed
from common.io import bgzip_file, check_processes
from common.functional import DesignError


def filter_xy_feature(
    smk: Any,
    sconf: cfg.GiabStrats,
    get_level: Callable[[cfg.XYFeatures], str | None],
) -> None:
    ws: dict[str, str] = smk.wildcards
    bed_input = cfg.smk_to_input_name(smk, "bed")
    valid_path = cfg.smk_to_input_name(smk, "valid")
    genome_path = cfg.smk_to_input_name(smk, "genome")
    bed_output = cfg.smk_to_output(smk)
    log = cfg.smk_to_log(smk)

    i = cfg.ChrIndex.from_name_unsafe(ws["sex_chr"])

    rfk = cfg.wc_to_reffinalkey(ws)
    rk = cfg.strip_full_refkey(rfk)

    pat = sconf.refkey_to_xy_ref_chr_pattern(rfk, i)

    bd = sconf.to_build_data(rk, cfg.wc_to_buildkey(ws))
    si = bd.refdata.strat_inputs
    if (
        si.xy is None
        or (fs := si.xy.features) is None
        or (level := get_level(fs)) is None
    ):
        raise DesignError()
    bf = si.xy_feature_bed_unsafe(i)
    conv = cfg.HapToHapChrConversion(bf.bed.chr_pattern, pat, bd.build_chrs)

    df = bf.read(bed_input)
    df_sorted = filter_sort_bed(conv.init_mapper, conv.final_mapper, df)
    level_mask = df_sorted[bf.level_col].str.contains(level)
    df_filtered = df_sorted[level_mask].drop(columns=[bf.level_col])
    # TODO put this in its own rule to simplify script?
    with bed_to_stream(df_filtered) as s:
        p, o = intersectBed(s, valid_path, genome_path)
        bgzip_file(o, bed_output)
        check_processes([p], log)
