import pandas as pd
import common.config as cfg
from typing import Callable


def sort_bed_numerically(df: pd.DataFrame) -> pd.DataFrame:
    cols = df.columns.tolist()
    return df.sort_values(
        by=[cols[0], cols[1], cols[2]],
        axis=0,
        ignore_index=True,
    )


def filter_sort_bed(
    sconf: cfg.GiabStrats,
    f: Callable[[cfg.Stratification], cfg.BedFile],
    rk: cfg.RefKey,
    bk: cfg.BuildKey,
    df: pd.DataFrame,
) -> pd.DataFrame:
    from_map = sconf.buildkey_to_init_chr_mapping(f, rk, bk)
    to_map = sconf.buildkey_to_final_chr_mapping(rk, bk)
    df[0] = df[0].map(from_map)
    df = sort_bed_numerically(df.dropna(subset=[0]))
    df[0] = df[0].map(to_map)
    return df
