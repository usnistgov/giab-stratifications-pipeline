# Generate valid regions bed files
#
# The long version - These define the regions that are allowed to be in the
# final stratification bed files (the unit tests enforce this). Historically,
# these are generated from the gaps bed file which accompanies hg19/38 and
# denotes places in the genome that are filled with Ns, where gaps are not
# included in the final stratifications. With gapless genomes such as CHM13,
# these gaps files don't exist (obviously) so its valid regions includes
# everything (almost, more on that next).
#
# To further complicate matters, any haploid genome which includes a Y
# chromosome normally masks the PAR since this region (by definition) behaves
# like an autosomal (ie diploid) chromosome, and it is redundant and confusing
# to include a diploid region in an otherwise haploid reference. Therefore, it
# makes sense to remove this region from the final stratifications, with the
# exception (just to be more complicated) of the Y PAR stratification itself.
# Full diploid genomes such as the HG002 Q100 assembly don't have this problem
# since all autosomes have two haplotype and the X and Y can then be naturally
# grouped with their respective haplotype.
#
# To deal with this masking issue, here we make two versions of the valid
# regions denoted "auto" and "parY". For haploid references, "auto" does not
# include the Y PAR which reflects the fact that this is masked. "parY" includes
# the Y PAR and is used in the XY stratifications to include the PAR region
# itself. For diploid genomes, "auto" and "parY" are identical and do not
# exclude the Y PAR. In all cases, the valid regions do not include gaps in the
# reference.


from typing import Any
from pathlib import Path
import common.config as cfg
from common.bed import (
    write_bed,
    bed_to_stream,
    complementBed,
    mergeBed,
    subtractBed,
)
from common.io import check_processes, tee, bgzip_file
from common.functional import DesignError
import pandas as pd


# convert genome to bed file (where each region is just the length of one
# chromosome)
def read_genome_bed(p: Path) -> pd.DataFrame:
    df = pd.read_table(
        p,
        header=None,
        names=["chrom", "end"],
        dtype={"chrom": str, "end": int},
    )
    df["start"] = 0
    return df[["chrom", "start", "end"]].copy()


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    genome_path = cfg.smk_to_input_name(smk, "genome")
    auto_out = cfg.smk_to_output_name(smk, "auto")
    parY_out = cfg.smk_to_output_name(smk, "parY")
    log = cfg.smk_to_log(smk)
    ws: dict[str, Any] = smk.wildcards
    rk = cfg.wc_to_reffinalkey(ws)

    is_haploid = len(sconf.refkey_haplotypes(rk)) == 0

    try:
        gap_inputs = cfg.smk_to_inputs_name(smk, "gaps")
    except DesignError:
        gap_inputs = None

    try:
        parY_path = cfg.smk_to_input_name(smk, "parY")
    except DesignError:
        parY_path = None

    # If we have gap input, invert this to get the gapless regions. These will
    # processed further depending on if we have a Y PAR and/or a haploid
    # reference.
    if gap_inputs is not None:
        gaps_df = sconf.with_build_data_and_bed_full(
            rk,
            cfg.wc_to_buildkey(ws),
            cfg.bd_to_gaps,
            lambda bd, b: cfg.with_inputs_hap(
                b,
                gap_inputs,
                lambda i, bf: cfg.read_filter_sort_hap_bed(bd, bf, i),
                lambda bc: cfg.build_hap_coords_df(bd, bc),
            ),
            lambda bd, b: cfg.with_inputs_dip1(
                b,
                gap_inputs,
                lambda i, bf: cfg.read_filter_sort_dip1to1_bed(bd, bf, i),
                lambda bc: cfg.build_dip1to1_coords_df(bd, bc),
            ),
            lambda hap, bd, b: cfg.with_inputs_dip1(
                b,
                gap_inputs,
                lambda i, bf: cfg.read_filter_sort_dip1to2_bed(bd, bf, i),
                lambda bc: cfg.build_dip1to2_coords_df(bd, bc),
            ).choose(hap),
            lambda bd, b: cfg.with_inputs_dip2(
                b,
                gap_inputs,
                lambda i, bf: cfg.read_filter_sort_dip2to1_bed(bd, bf, i),
                lambda bc: cfg.build_dip2to1_coords_df(bd, bc),
            ),
            lambda hap, bd, b: cfg.with_inputs_dip2(
                b,
                gap_inputs,
                lambda i, bf: cfg.read_filter_sort_dip2to2_bed(
                    bd, bf, i.choose(hap), hap
                ),
                lambda bc: cfg.build_dip2to2_coords_df(hap, bd, bc),
            ),
        )

        with bed_to_stream(gaps_df) as s:
            p1, o1 = mergeBed(s, ["-d", "100"])
            p2, o2 = complementBed(o1, genome_path)

            # If the Y PAR is given and this is a haploid reference, subtract
            # parY from the gapless regions. Otherwise just link them since we
            # either don't know what to subtract or have a diploid reference
            # which doesn't require this auto/parY distinction.
            if parY_path is not None and is_haploid:
                p3, o3, o4 = tee(o2)
                p4, o5 = subtractBed(o3, parY_path, genome_path)

                bgzip_file(o4, parY_out)
                bgzip_file(o5, auto_out)

                check_processes([p1, p2, p3, p4], log)
            else:
                bgzip_file(o2, auto_out)
                parY_out.symlink_to(auto_out.resolve())

                check_processes([p1, p2], log)

    # If we don't have a gaps bed file, use the genome bed file as the gapless
    # regions (which means everything is gapless). Then subtract off the Y
    # region if we have a haploid reference and we know the Y PAR coordinates.
    else:
        genome_bed = read_genome_bed(genome_path)
        write_bed(parY_out, genome_bed)

        if parY_path is not None and is_haploid:
            with open(parY_out, "rb") as f:
                p, o = subtractBed(f, parY_path, genome_path)
                bgzip_file(o, auto_out)
                check_processes([p], log)
        else:
            auto_out.symlink_to(parY_out.resolve())


main(snakemake, snakemake.config)  # type: ignore
