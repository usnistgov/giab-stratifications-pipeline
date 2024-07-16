"""Master configuration definition for the entire pipeline :)

Overview:

The entire pipeline configuration is defined here as a pydantic model, which
validates and types all data prior to downstream processing in snakemake rules
or in scripts called by snakemake. Since pydantic interfaces with mypy, and
since mypy can be used to lint python code for correctness, the logic of this
pipeline is heavily biased toward python code rather than snakemake.

In general, snakemake is used to handle build dependencies, and almost all other
logic is handled via python scripts which are type-checked "at the boundaries.

Data hierarchy and flow:

Each target reference is index by a key called a ref key ('ref_key'). Each
reference can have one of more builds indexed by a build key ('build_key') which
describes what should be included (chromosome numbers and various stratification
levels) within a given stratification 'build' for the target reference under
which it is located. The number of builds is found by the Cartesian product of
all ref keys and build keys.

Additionally, ref/build keys are grouped under one of three categories
corresponding to the way its haplotypes are configured: haploid, diploid1, and
diploid2 (see next section). In the case of diploid2, each haplotype is spread
across two files. For this reason (and others), we further distinguish the ref
key into a "full ref key" which may have a haplotype associated with it. A full
ref key associated with the reference is also called a "final ref key" and one
associated with an input is called a "source ref key."

Any given reference (indexed by the ref key) can have input files associated
with it; these will be used to make some of the stratifications. In general
these are bed files or files that can be coerced into being bed files by
selecting certain columns. Note that in the case of diploid assemblies, these
inputs may have to be merged or split depending on if the target reference and
the input file have each haplotype in one file or two files (see next section).

Haploid vs Diploid:

This pipeline is meant to generate stratifications for human genomes which may
either be haploid or diploid. In the case of diploid references, there are two
cases to consider: dip1 or dip2 (see below for terminology). In the case of dip1
we only output one final directory since the two haplotypes are in one file (and
presumably will be consumed as such). In the case of dip2, we make two output
directories (one per haplotype) since the haplotypes are spread across two
reference files. Thus from the reference perspective (ie the fasta for which the
stratifications must be built) we have 3 cases to consider: dip1, dip2, or hap.

To further complicate matters, the inputs for the dip1 and dip2 cases may
themselves be some combination of dip1 or dip2 (ie a reference might be dip1 but
a bed file for the reference might be dip2). Thus for input files, we have 5
cases to consider:
* hap -> hap: input and reference are both haploid
* dip1 -> dip1: input and reference are both dip1
* dip1 -> dip2: input is dip1, reference is dip2
* dip2 -> dip2: input is dip2, reference is dip1
* dip2 -> dip2: input and reference are both dip2

In the case of dip1 -> dip2, the input must be split. In the case of
dip2 -> dip1, the inputs (two) must be combined. For the other three, each input
can be used more or less as-is since the cardinality of the inputs and reference
matches. To make the dip2 -> dip1 case even more complex, the chromosome names
might be identical across the two input files, but will need to be distinguished
when combined by adding a suffix to them.

To make this process as simple as possible, input bed files (and related) will
be downloaded and processed immediately into bed output that correspond directly
to their target reference. This post-download processing step is called
"normalization" and generally consists of the following:
* filters desired chromosomes
* sorts all chromosomes by numeric index (with X and Y being 23 and 24, and the
  paternal haplotype sorted prior to the maternal haplotype when applicable)
* splits or combines the input as described above when applicable
* renames all chromosomes to match the target reference
* rearranges all columns into proper bed format

Each normalization is also generally a checkpoint, since it is not known until
runtime how many files it needs to consume or produce. (Note that this is a
design choice to reduce repeated code. It is technically feasible to produce
snakemake rules and scripts that individual handle each of the 5 cases in a way
that doesn't require checkpoints, but these 5 rules would need to be repeated
for each input file which would be more error prone than the solution chosen
here)

Terminology:
* haploid (or "hap"): half a diploid dataset, ie one haplotype
* diploid1 (or "dip1"): a diploid dataset that is all in one file
* diploid2 (or "dip2"): a diploid dataset that is in two files

Conventions:
* 'DesignError' exceptions are those that should not happen; if they do the code
  is incorrect
* bed files are gzip'ed, fasta files are bgzip'ed, and all files are compressed
  with one of these after downloading
* paternal sorts before maternal
* chromosomes are numbered and sorted 1-24 where X and Y are 23/24 respectively
"""

from __future__ import annotations
import sys
import json
import pandas as pd
import re
from textwrap import fill
from pathlib import Path
from pydantic import BaseModel as BaseModel_
from pydantic.generics import GenericModel as GenericModel_
from pydantic.generics import GenericModelT
from pydantic import validator, HttpUrl, FilePath, NonNegativeInt, Field
from dataclasses import dataclass
from enum import unique, Enum
from typing import (
    IO,
    Union,
    NewType,
    Any,
    Callable,
    TypeVar,
    Type,
    NamedTuple,
    cast,
    Annotated,
    Generic,
    TypeGuard,
    Protocol,
)
from typing_extensions import Self, assert_never
from functools import reduce
from itertools import chain
from more_itertools import duplicates_everseen, flatten
from common.functional import (
    maybe2,
    from_maybe,
    fmap_maybe,
    fmap_maybe_def,
    maybe_to_list,
    with_first,
    DesignError,
    not_none_unsafe,
    none_unsafe,
    unzip2,
    unzip3,
    noop,
    raise_inline,
)
from common.io import is_gzip, is_bgzip, get_md5
import common.bed as bed


################################################################################
# Type aliases

Percent = Annotated[int, Field(ge=0, le=100)]

RefKey = NewType("RefKey", str)
# full refkey represented as a string (see below for class)
RefKeyFullS = NewType("RefKeyFullS", str)
BuildKey = NewType("BuildKey", str)
CompareKey = NewType("CompareKey", str)
OtherLevelKey = NewType("OtherLevelKey", str)
OtherStratKey = NewType("OtherStratKey", str)
HaplotypeName = NewType("HaplotypeName", str)

# helper type representing snakemake wildcards
SmkWildcards = dict[str, Any]

GCBound = tuple[Percent, bool]

# a bare chromosome name (like "1" or "X")
ShortChrName = NewType("ShortChrName", str)

# the set chromosomes desired per build
BuildChrs = NewType("BuildChrs", "set[ChrIndex]")

# the set of chromosomes specific to a haplotype
HapChrs = NewType("HapChrs", "set[ChrIndex]")

# an ordered list of chromosomes specific to a haplotype (sorted numerically)
OrderedHapChrs = NewType("OrderedHapChrs", "list[ChrIndex]")

# an ordered list of chromosome names specific to a haplotype (sorted numerically)
OrderedHapChrNames = NewType("OrderedHapChrNames", list[bed.ChrName])


################################################################################
# Type variables

W = TypeVar("W")
X = TypeVar("X")
Y = TypeVar("Y")
Z = TypeVar("Z")


################################################################################
# Helper functions


def sub_wildcard(s: str, wc: str, rep: str) -> str:
    return s.replace("{" + wc + "}", rep)


def sub_wildcards(s: str, wc: dict[str, str]) -> str:
    return reduce(lambda acc, nxt: sub_wildcard(acc, *nxt), wc.items(), s)


def sub_wildcard_path(s: Path, wc: str, rep: str) -> Path:
    return Path(sub_wildcard(str(s), wc, rep))


def sub_wildcards_path(s: Path, wc: dict[str, str]) -> Path:
    return Path(sub_wildcards(str(s), wc))


def flip_hap(h: Haplotype) -> Haplotype:
    return h.choose(Haplotype.MAT, Haplotype.PAT)


def parse_full_refkey_class(s: RefKeyFullS) -> RefKeyFull:
    m = re.match("(.+)\\.([mp]at)", s)
    # ASSUME this will never fail due to the pat/mat permitted match pattern
    rk, hap = (s, None) if m is None else (m[1], Haplotype.from_name(m[2]))
    return RefKeyFull(RefKey(rk), hap)


def parse_full_refkey(s: RefKeyFullS) -> tuple[RefKey, Haplotype | None]:
    return parse_full_refkey_class(s).as_tuple


def strip_full_refkey(s: RefKeyFullS) -> RefKey:
    return parse_full_refkey(s)[0]


def flip_full_refkey_class(r: RefKeyFull) -> RefKeyFull:
    return RefKeyFull(r.key, fmap_maybe(flip_hap, r.hap))


def flip_full_refkey(s: RefKeyFullS) -> RefKeyFullS:
    return flip_full_refkey_class(parse_full_refkey_class(s)).name


def choose_xy_unsafe(c: ChrIndex, x_res: X, y_res: X) -> X:
    if c is ChrIndex.CHRX:
        return x_res
    elif c is ChrIndex.CHRY:
        return y_res
    else:
        raise DesignError(f"I am not an X or Y, I am a {c}")


def sort_chr_indices(cs: HapChrs) -> OrderedHapChrs:
    return OrderedHapChrs(
        [x for x, _ in sorted([(c, c.value) for c in cs], key=lambda x: x[1])]
    )


def refkey_config_to_prefix(split: bool, nohap: bool) -> str:
    sp = "split" if split else "nosplit"
    hp = "nohap" if nohap else "withhap"
    return f"{sp}_{hp}"


def prefix_to_refkey_config(s: str) -> tuple[bool, bool]:
    m = re.match("^([^_]+)_([^_]+)", s)

    if m is None:
        raise DesignError(f"could not parse refkeys config: {s}")

    match m[1]:
        case "split":
            split = True
        case "nosplit":
            split = False
        case _ as e:
            raise DesignError(f"unknown split {e}")

    match m[2]:
        case "withhap":
            nohap = False
        case "nohap":
            nohap = True
        case _ as e:
            raise DesignError(f"unknown split {e}")

    return split, nohap


def make_double(f: Callable[[Haplotype], X]) -> Double[X]:
    return Double(pat=f(Haplotype.PAT), mat=f(Haplotype.MAT))


def to_refkeys(x: SingleOrDouble[X], rk: RefKey) -> RefKeyFull1or2:
    if isinstance(x, Single):
        return x.key(rk)
    elif isinstance(x, Double):
        return x.keys(rk)
    else:
        assert_never(x)


def to_str_refkeys(x: SingleOrDouble[X], rk: RefKey) -> RefKeyFullS1or2:
    return to_refkeys(x, rk).map(lambda k: k.name)


def match1_unsafe(xs: list[X], f: Callable[[X], Y], msg: None | str = None) -> Y:
    """Call function with the one value from a singleton list.

    Error if list is not a singleton.
    """
    match xs:
        case [x]:
            return f(x)
        case _:
            raise DesignError(
                msg if msg is not None else f"One input expected, got {len(xs)}"
            )


def match2_unsafe(
    xs: list[X], f: Callable[[Double[X]], Y], msg: None | str = None
) -> Y:
    """Call function with the twos value from a 2-ary list.

    Error if list does not have two members.
    """
    match xs:
        case [x1, x2]:
            return f(Double(x1, x2))
        case _:
            raise DesignError(
                msg if msg is not None else f"Two inputs expected, got {len(xs)}"
            )


def match12_unsafe(
    xs: list[X],
    f1: Callable[[X], Y],
    f2: Callable[[Double[X]], Y],
    msg: None | str = None,
) -> Y:
    """Combination of `match1_unsafe` and `match2_unsafe` with two functions
    for each case. Error if input list is does not have one or two elements.
    """
    match xs:
        case [x1]:
            return f1(x1)
        case [x1, x2]:
            return f2(Double(x1, x2))
        case _:
            raise DesignError(
                msg if msg is not None else f"Two inputs expected, got {len(xs)}"
            )


# type helpers


# TODO mypy for some reason doesn't understand how to narrow a
# Something[Union[X, Y]] to a Something[X] using 'isinstance'
def is_dip1_bed(
    x: BedFile[Dip1Src[S] | Dip2Src[S]] | DipBedCoords,
) -> TypeGuard[BedFile[Dip1Src[S]] | Dip1BedCoords]:
    return isinstance(x, Dip1BedCoords) or (
        isinstance(x, BedFile) and isinstance(x.bed, Dip1Src)
    )


def is_dip2_bed(
    x: BedFile[Dip1Src[S] | Dip2Src[S]] | DipBedCoords,
) -> TypeGuard[BedFile[Dip2Src[S]] | Dip2BedCoords]:
    return isinstance(x, Dip2BedCoords) or (
        isinstance(x, BedFile) and isinstance(x.bed, Dip2Src)
    )


# union detanglers


def with_dip_bedfile(
    bf: BedFile[Dip1Src[S] | Dip2Src[S]] | DipBedCoords,
    dip1: Callable[[BedFile[Dip1Src[S]] | Dip1BedCoords], Y],
    dip2: Callable[[BedFile[Dip2Src[S]] | Dip2BedCoords], Y],
) -> Y:
    if is_dip1_bed(bf):
        return dip1(bf)
    elif is_dip2_bed(bf):
        return dip2(bf)
    else:
        # TODO this is a mypy bug, I should be able to use assert_never here
        raise DesignError("not a dip1 or dip2")
        # assert_never(bf)


def match_single_unsafe(f: Callable[[X], Y], x: SingleOrDouble[X]) -> Y:
    if isinstance(x, Single):
        return f(x.elem)
    else:
        raise DesignError()


def match_double_unsafe(f: Callable[[Double[X]], Y], x: SingleOrDouble[X]) -> Y:
    if isinstance(x, Double):
        return f(x)
    else:
        raise DesignError()


def with_single_or_double(
    single_f: Callable[[Single[X]], Y],
    double_f: Callable[[Double[X]], Y],
    x: SingleOrDouble[X],
) -> Y:
    if isinstance(x, Single):
        return single_f(x)
    if isinstance(x, Double):
        return double_f(x)
    else:
        assert_never(x)


def from_single_or_double(x: Single[X] | Double[X], hap: Haplotype | None) -> X:
    return with_single_or_double(
        lambda x: none_unsafe(hap, x.elem),
        lambda x: not_none_unsafe(hap, lambda h: x.choose(h)),
        x,
    )


def with_ref_data(
    rd: AnyRefData,
    hap_f: Callable[[HapRefData], X],
    dip1_f: Callable[[Dip1RefData], X],
    dip2_f: Callable[[Dip2RefData], X],
) -> X:
    if isinstance(rd.ref, HapSrc):
        return hap_f(rd)
    elif isinstance(rd.ref, Dip1Src):
        return dip1_f(rd)
    elif isinstance(rd.ref, Dip2Src):
        return dip2_f(rd)
    else:
        assert_never(rd)


def with_build_data(
    bd: AnyBuildData,
    hap_f: Callable[[HapBuildData], X],
    dip1_f: Callable[[Dip1BuildData], X],
    dip2_f: Callable[[Dip2BuildData], X],
) -> X:
    if isinstance(bd.refdata.ref, HapSrc):
        return hap_f(bd)
    elif isinstance(bd.refdata.ref, Dip1Src):
        return dip1_f(bd)
    elif isinstance(bd.refdata.ref, Dip2Src):
        return dip2_f(bd)
    else:
        assert_never(bd)


def with_inputs_null_hap(
    b: BedFile[HapSrc[S]] | HapBedCoords | None,
    inputs: list[X],
    src_f: Callable[[X, BedFile[HapSrc[S]]], Y],
    coords_f: Callable[[HapBedCoords], Y],
    none_f: Callable[[], Y],
) -> Y:
    if isinstance(b, BedFile):
        return match1_unsafe(inputs, lambda i: src_f(i, b))
    elif isinstance(b, HapBedCoords):
        if len(inputs) > 0:
            raise DesignError()
        return coords_f(b)
    elif b is None:
        return none_f()
    else:
        assert_never(b)


def with_inputs_hap(
    b: BedFile[HapSrc[S]] | HapBedCoords | None,
    inputs: list[X],
    src_f: Callable[[X, BedFile[HapSrc[S]]], Y],
    coords_f: Callable[[HapBedCoords], Y],
) -> Y:
    return with_inputs_null_hap(b, inputs, src_f, coords_f, raise_inline)


def with_inputs_null_dip1(
    b: Dip1BedFileOrCoords,
    inputs: list[X],
    src_f: Callable[[X, Dip1BedFile], Y],
    coords_f: Callable[[Dip1BedCoords], Y],
    none_f: Callable[[], Y],
) -> Y:
    if isinstance(b, BedFile):
        return match1_unsafe(inputs, lambda i: src_f(i, b))
    elif isinstance(b, Dip1BedCoords) or isinstance(b, Dip2BedCoords):
        if len(inputs) > 0:
            raise DesignError()
        return coords_f(b)
    elif b is not None:
        return coords_f()
    else:
        assert_never(b)


def with_inputs_dip1(
    b: Dip1BedFileOrCoords,
    inputs: list[X],
    src_f: Callable[[X, Dip1BedFile], Y],
    coords_f: Callable[[Dip1BedCoords], Y],
) -> Y:
    return with_inputs_null_dip1(b, inputs, src_f, coords_f, raise_inline)


def with_inputs_null_dip2(
    b: Dip2BedFileOrCoords,
    inputs: list[X],
    src_f: Callable[[Double[X], Dip2BedFile], Y],
    coords_f: Callable[[Dip2BedCoords], Y],
    none_f: Callable[[], Y],
) -> Y:
    if isinstance(b, BedFile):
        return match2_unsafe(inputs, lambda i: src_f(i, b))
    elif isinstance(b, Dip1BedCoords) or isinstance(b, Dip2BedCoords):
        if len(inputs) > 0:
            raise DesignError()
        return coords_f(b)
    elif b is not None:
        return coords_f()
    else:
        assert_never(b)


def with_inputs_dip2(
    b: Dip2BedFileOrCoords,
    inputs: list[X],
    src_f: Callable[[Double[X], Dip2BedFile], Y],
    coords_f: Callable[[Dip2BedCoords], Y],
) -> Y:
    return with_inputs_null_dip2(b, inputs, src_f, coords_f, raise_inline)


# noop conversion getters


def hap_noop_conversion(bd: HapBuildData) -> HapToHapChrConversion:
    return bd.refdata.ref.noop_conversion(bd.build_chrs)


def dip1_noop_conversion(bd: Dip1BuildData) -> DipToDipChrConversion:
    return bd.refdata.ref.noop_conversion(bd.build_chrs)


def dip1_split_noop_conversion(
    h: Haplotype, bd: Dip1BuildData
) -> HapToHapChrConversion:
    return bd.refdata.ref.noop_conversion(bd.build_chrs).split(h)


def dip2_noop_conversion(h: Haplotype, bd: Dip2BuildData) -> HapToHapChrConversion:
    return bd.refdata.ref.noop_conversion(bd.build_chrs).choose(h)


# functions for dealing with 'dict[RefKey, X]' type things


def to_ref_data_unsafe(
    xs: dict[
        RefKey,
        Stratification[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
    ],
    rk: RefKey,
) -> RefData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT]:
    try:
        s = xs[rk]
        return RefData_(rk, s.ref, s.strat_inputs, s.builds)
    except KeyError:
        raise DesignError(f"Could not get ref data for key '{rk}'")


def all_ref_data(
    xs: dict[
        RefKey,
        Stratification[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
    ],
) -> list[RefData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT]]:
    return [to_ref_data_unsafe(xs, rk) for rk in xs]


def all_ref_refsrckeys(
    xs: dict[
        RefKey,
        Stratification[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
    ],
) -> list[RefKeyFullS]:
    return [s for k, v in xs.items() for s in to_str_refkeys(v.ref.src, k).as_list]


def all_build_data(
    xs: dict[
        RefKey,
        Stratification[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
    ],
) -> list[BuildData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT]]:
    return [r.to_build_data_unsafe(b) for r in all_ref_data(xs) for b in r.builds]


def all_bed_build_and_refsrckeys(
    xs: dict[
        RefKey,
        Stratification[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
    ],
    f: BuildDataToSrc,
) -> list[tuple[RefKeyFullS, BuildKey]]:
    return [
        (rk, b.buildkey)
        for b in all_build_data(xs)
        if (src := f(b)) is not None
        for rk in to_str_refkeys(src, b.refdata.refkey).as_list
    ]


def all_bed_refsrckeys(
    xs: dict[
        RefKey,
        Stratification[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
    ],
    f: BuildDataToSrc,
) -> list[RefKeyFullS]:
    return [rk for rk, _ in all_bed_build_and_refsrckeys(xs, f)]


def all_build_keys(
    xs: dict[
        RefKey,
        Stratification[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
    ],
) -> list[tuple[RefKey, BuildKey]]:
    return [(r.refdata.refkey, r.buildkey) for r in all_build_data(xs)]


def all_ref_build_keys(
    xs: dict[
        RefKey,
        Stratification[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
    ],
) -> list[tuple[RefKeyFullS, BuildKey]]:
    return [
        (rk, r.buildkey)
        for r in all_build_data(xs)
        for rk in to_str_refkeys(r.refdata.ref.src, r.refdata.refkey).as_list
    ]


# path formatters


def sub_output_path(pat: str, rk: RefKeyFull) -> Path:
    if "{" in pat or "}" in pat:
        raise DesignError(f"not all wildcards replaced in pattern {pat}")
    return Path(pat.replace("%s", rk.name))


def prepare_output_path(path: Path) -> Path:
    return Path(str(path).replace("{ref_key}", "%s"))


# itty bitty accessor functions


def bd_to_si(
    f: StratInputToBed,
    x: BuildData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
) -> BedFile[BedSrcT] | BedCoordsT | None:
    return f(x.refdata.strat_inputs)


def si_to_cds(x: StratInputs[BedSrcT, BedCoordsT]) -> CDS[BedSrcT] | BedCoordsT | None:
    return x.functional.cds


def si_to_mhc(
    x: StratInputs[BedSrcT, BedCoordsT]
) -> BedFile[BedSrcT] | BedCoordsT | None:
    return x.functional.mhc


def si_to_kir(
    x: StratInputs[BedSrcT, BedCoordsT]
) -> BedFile[BedSrcT] | BedCoordsT | None:
    return x.functional.kir


def si_to_vdj(
    x: StratInputs[BedSrcT, BedCoordsT]
) -> BedFile[BedSrcT] | BedCoordsT | None:
    return x.functional.vdj


def bd_to_cds(
    x: BuildData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
) -> CDS[BedSrcT] | BedCoordsT | None:
    return si_to_cds(x.refdata.strat_inputs)


def bd_to_mhc(
    x: BuildData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
) -> BedFile[BedSrcT] | BedCoordsT | None:
    return si_to_mhc(x.refdata.strat_inputs)


def bd_to_kir(
    x: BuildData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
) -> BedFile[BedSrcT] | BedCoordsT | None:
    return si_to_kir(x.refdata.strat_inputs)


def bd_to_vdj(
    x: BuildData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
) -> BedFile[BedSrcT] | BedCoordsT | None:
    return si_to_vdj(x.refdata.strat_inputs)


def si_to_simreps(
    x: StratInputs[BedSrcT, BedCoordsT]
) -> BedFile[BedSrcT] | BedCoordsT | None:
    return x.low_complexity.simreps


def bd_to_simreps(
    x: BuildData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
) -> BedFile[BedSrcT] | BedCoordsT | None:
    return si_to_simreps(x.refdata.strat_inputs)


def si_to_rmsk(
    x: StratInputs[BedSrcT, BedCoordsT]
) -> RMSKFile[BedSrcT] | BedCoordsT | None:
    return x.low_complexity.rmsk


# TODO take the boolean switch out of here
def bd_to_rmsk(
    x: BuildData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
) -> RMSKFile[BedSrcT] | BedCoordsT | None:
    return si_to_rmsk(x.refdata.strat_inputs)


def si_to_satellites(
    x: StratInputs[BedSrcT, BedCoordsT]
) -> SatFile[BedSrcT] | BedCoordsT | None:
    return x.low_complexity.satellites


def bd_to_satellites(
    x: BuildData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
) -> SatFile[BedSrcT] | BedCoordsT | None:
    return si_to_satellites(x.refdata.strat_inputs)


def si_to_superdups(
    x: StratInputs[BedSrcT, BedCoordsT]
) -> BedFile[BedSrcT] | BedCoordsT | None:
    return x.segdups.superdups


def bd_to_superdups(
    x: BuildData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
) -> BedFile[BedSrcT] | BedCoordsT | None:
    return si_to_superdups(x.refdata.strat_inputs)


def si_to_gaps(
    x: StratInputs[BedSrcT, BedCoordsT]
) -> BedFile[BedSrcT] | BedCoordsT | None:
    return x.gap


def bd_to_gaps(
    x: BuildData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
) -> BedFile[BedSrcT] | BedCoordsT | None:
    return si_to_gaps(x.refdata.strat_inputs)


def bd_to_other(
    lk: OtherLevelKey,
    sk: OtherStratKey,
    x: BuildData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
) -> OtherBedFile[BedSrcT, BedCoordsT]:
    return x.build.other_strats[lk][sk]


def bd_to_other_bed(
    lk: OtherLevelKey,
    sk: OtherStratKey,
    x: BuildData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
) -> BedFile[BedSrcT] | BedCoordsT | None:
    return x.build.other_strats[lk][sk].data


def bd_to_bench_bed(
    x: BuildData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
) -> BedFile[BedSrcT] | None:
    return fmap_maybe(lambda y: y.bench_bed, x.build.bench)


def bd_to_bench_vcf(
    x: BuildData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
) -> VCFFile[VcfSrcT] | None:
    return fmap_maybe(lambda y: y.bench_vcf, x.build.bench)


def bd_to_query_vcf(
    x: BuildData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT],
) -> VCFFile[VcfSrcT] | None:
    return fmap_maybe(lambda y: y.query_vcf, x.build.bench)


# snakemake helpers


def wc_lookup(ws: SmkWildcards, k: str) -> Any:
    try:
        return ws[k]
    except KeyError:
        raise DesignError(f"Could not find {k} in wildcards")


def wc_to_refkey(ws: SmkWildcards) -> RefKey:
    return RefKey(wc_lookup(ws, "ref_key"))


def wc_to_buildkey(ws: SmkWildcards) -> BuildKey:
    return BuildKey(wc_lookup(ws, "build_key"))


def wc_to_reffinalkey(ws: SmkWildcards) -> RefKeyFullS:
    return RefKeyFullS(wc_lookup(ws, "ref_final_key"))


def smk_to_param_str(smk: Any, name: str) -> str:
    ps = smk.params
    if hasattr(ps, name):
        x = ps[name]
        if isinstance(x, str):
            return x
        else:
            raise DesignError(f"{name} in params is not a string")
    else:
        raise DesignError(f"Params does not have {name}")


def smk_to_param_strs(smk: Any, name: str) -> list[str]:
    ps = smk.params
    if hasattr(ps, name):
        x = ps[name]
        if isinstance(x, list) and all([isinstance(y, str) for y in x]):
            return x
        else:
            raise DesignError(f"Params for {name} are not a list")
    else:
        raise DesignError(f"Params does not have {name}")


def smk_to_param_path(smk: Any, name: str) -> Path:
    return Path(smk_to_param_str(smk, name))


def smk_to_param_paths(smk: Any, name: str) -> list[Path]:
    return [Path(p) for p in smk_to_param_strs(smk, name)]


def smk_to_output(smk: Any, n: int = 0) -> Path:
    try:
        return Path(smk.output[n])
    except IndexError:
        raise DesignError(f"No output file for index {n}")


def smk_to_output_name(smk: Any, name: str) -> Path:
    # TODO not DRY
    ps = smk.output
    if hasattr(ps, name):
        p = ps[name]
        if isinstance(p, str):
            return Path(p)
        else:
            raise DesignError(f"Output file for {name} is not a string")
    else:
        raise DesignError(f"Output files do not have name {name}")


def smk_to_log(smk: Any, n: int = 0) -> Path:
    try:
        return Path(smk.log[n])
    except IndexError:
        raise DesignError(f"No log file for index {n}")


def smk_to_log_name(smk: Any, name: str) -> Path:
    ps = smk.log
    if hasattr(ps, name):
        p = ps[name]
        if isinstance(p, str):
            return Path(p)
        else:
            raise DesignError(f"Log file for {name} is not a string")
    else:
        raise DesignError(f"Log files do not have name {name}")


def smk_to_input(smk: Any, n: int = 0) -> Path:
    i = smk.input[n]
    if isinstance(i, str):
        return Path(i)
    else:
        raise DesignError(f"Input files for {i} are a list")


def smk_to_inputs(smk: Any, n: int = 0, allow_empty: bool = False) -> list[Path]:
    i = smk.input[n]
    if isinstance(i, str):
        raise DesignError(f"Input files for {i} are not a list")
    else:
        ps = [Path(p) for p in i]
        if len(i) == 0 and not allow_empty:
            raise DesignError(f"Input files for {i} is an empty list")
        return ps


def smk_to_inputs_all(smk: Any, allow_empty: bool = False) -> list[Path]:
    i = smk.input
    if isinstance(i, str):
        raise DesignError(f"Input files for {i} are not a list")
    else:
        ps = [Path(p) for p in i]
        if len(i) == 0 and not allow_empty:
            raise DesignError(f"Input files for {i} is an empty list")
        return ps


def smk_to_input_name(smk: Any, name: str) -> Path:
    i = smk.input
    if hasattr(i, name):
        x = i[name]
        if isinstance(x, str):
            return Path(x)
        else:
            raise DesignError(f"Input files for {name} are a list")
    else:
        raise DesignError(f"Input files do not have name {name}")


def smk_to_inputs_name(smk: Any, name: str, allow_empty: bool = False) -> list[Path]:
    i = smk.input
    if hasattr(i, name):
        x = i[name]
        if isinstance(x, str):
            raise DesignError(f"Input files for {name} are not a list")
        else:
            ps = [Path(p) for p in x]
            if len(ps) == 0 and not allow_empty:
                raise DesignError(f"Input files for {name} is an empty list")
            return ps

    else:
        raise DesignError(f"Input files do not have name {name}")


# IO functions for processing bed files of various flavors


def read_filter_sort_hap_bed(
    bd: HapBuildData, bf: HapBedFile, ipath: Path
) -> pd.DataFrame:
    """Read a haploid bed file, sort it, and write it in bgzip format."""
    conv = bd.refdata.ref.chr_conversion(bf.bed.chr_pattern, bd.build_chrs)
    df = bf.read(ipath)
    return bed.filter_sort_bed(conv.init_mapper, conv.final_mapper, df)


def read_write_filter_sort_hap_bed(
    ipath: Path,
    opath: Path,
    bd: HapBuildData,
    bf: HapBedFile,
    g: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
) -> None:
    """Read a haploid bed file, sort it, and write it in bgzip format."""
    df = read_filter_sort_hap_bed(bd, bf, ipath)
    bed.write_bed(opath, g(df))


def read_filter_sort_dip1to1_bed(
    bd: Dip1BuildData,
    bf: Dip1BedFile,
    ipath: Path,
) -> pd.DataFrame:
    """Read a diploid bed file, sort it, and write it in bgzip format."""
    conv = bd.refdata.ref.dip_chr_conversion(bf.bed.chr_pattern, bd.build_chrs)
    df = bf.read(ipath)
    return bed.filter_sort_bed(conv.init_mapper, conv.final_mapper, df)


def read_write_filter_sort_dip1to1_bed(
    ipath: Path,
    opath: Path,
    bd: Dip1BuildData,
    bf: Dip1BedFile,
    g: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
) -> None:
    """Read a haploid bed file, sort it, and write it in bgzip format."""
    df = read_filter_sort_dip1to1_bed(bd, bf, ipath)
    bed.write_bed(opath, g(df))


def read_filter_sort_dip2to1_bed(
    bd: Dip1BuildData,
    bf: Dip2BedFile,
    ipath: Double[Path],
) -> pd.DataFrame:
    """Read two haploid bed files, combine and sort them as diploid, and write
    it in bgzip format.
    """

    conv = bd.refdata.ref.hap_chr_conversion(bf.bed.chr_pattern, bd.build_chrs)
    imap = conv.init_mapper
    fmap = conv.final_mapper

    return pd.concat(
        imap.both(
            lambda i, hap: bed.filter_sort_bed(i, fmap, bf.read(ipath.choose(hap)))
        ).as_list
    )


def read_write_filter_sort_dip2to1_bed(
    ipath: Double[Path],
    opath: Path,
    bd: Dip1BuildData,
    bf: Dip2BedFile,
    g: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
) -> None:
    """Read a haploid bed file, sort it, and write it in bgzip format."""
    df = read_filter_sort_dip2to1_bed(bd, bf, ipath)
    bed.write_bed(opath, g(df))


def read_filter_sort_dip1to2_bed(
    bd: Dip2BuildData,
    bf: Dip1BedFile,
    ipath: Path,
) -> Double[pd.DataFrame]:
    conv = bd.refdata.ref.dip_chr_conversion(bf.bed.chr_pattern, bd.build_chrs)
    imap, splitter = conv.init_mapper
    fmap = conv.final_mapper

    # TODO small micro-optimization for when I feel like it; splitting the bed
    # will also filter it, so we need only sort it after vs filtering and
    # sorting
    def go(df: pd.DataFrame, fmap: bed.FinalMapper) -> pd.DataFrame:
        return bed.filter_sort_bed(imap, fmap, df)

    df = bf.read(ipath)
    df0, df1 = bed.split_bed(splitter, df)
    return Double(go(df0, fmap.pat), go(df1, fmap.mat))


def read_write_filter_sort_dip1to2_bed(
    ipath: Path,
    opath: Double[Path],
    bd: Dip2BuildData,
    bf: Dip1BedFile,
    g0: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
    g1: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
) -> None:
    """Read a haploid bed file, sort it, and write it in bgzip format."""
    res = read_filter_sort_dip1to2_bed(bd, bf, ipath).sum(opath)
    res.both(lambda r, hap: bed.write_bed(r[1], hap.choose(g0, g1)(r[0])))


def read_filter_sort_dip2to2_bed(
    bd: Dip2BuildData,
    bf: Dip2BedFile,
    ipath: Path,
    hap: Haplotype,
) -> pd.DataFrame:
    conv = bd.refdata.ref.hap_chr_conversion(bf.bed.chr_pattern, bd.build_chrs)
    df = bf.read(ipath)
    conv_ = conv.choose(hap)
    return bed.filter_sort_bed(conv_.init_mapper, conv_.final_mapper, df)


def read_write_filter_sort_dip2to2_bed(
    ipath: Path,
    opath: Path,
    hap: Haplotype,
    bd: Dip2BuildData,
    bf: Dip2BedFile,
    g: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
) -> None:
    df = read_filter_sort_dip2to2_bed(bd, bf, ipath, hap)
    bed.write_bed(opath, g(df))


def build_hap_coords_df(bd: HapBuildData, bf: HapBedCoords) -> pd.DataFrame:
    to_map = bd.refdata.ref.chr_pattern.final_mapper(bd.build_chrs, Haplotype.PAT)
    return bed.indexed_bedlines_to_df(bf.lines(Haplotype.PAT).elem, to_map)


def build_dip1to1_coords_df(bd: Dip1BuildData, bf: Dip1BedCoords) -> pd.DataFrame:
    to_map = bd.refdata.ref.chr_pattern.final_mapper(bd.build_chrs)
    return bed.indexed_bedlines_to_df(bf.lines.elem, to_map)


def build_dip1to2_coords_df(
    bd: Dip2BuildData, bf: Dip1BedCoords
) -> Double[pd.DataFrame]:
    lines = bf.lines.elem
    pat = [x for x in lines if x.chr < 24]
    mat = [x for x in lines if x.chr >= 24]
    return bd.refdata.ref.chr_pattern.both(
        lambda p, h: bed.indexed_bedlines_to_df(
            h.choose(pat, mat),
            p.final_mapper(bd.build_chrs, h),
        )
    )


def build_dip2to1_coords_df(bd: Dip1BuildData, bf: Dip2BedCoords) -> pd.DataFrame:
    to_map = bd.refdata.ref.chr_pattern.final_mapper(bd.build_chrs)
    return bed.indexed_bedlines_to_df([*chain(*bf.lines.as_tuple)], to_map)


def build_dip2to2_coords_df(
    hap: Haplotype, bd: Dip2BuildData, bf: Dip2BedCoords
) -> pd.DataFrame:
    pattern = bd.refdata.ref.chr_pattern.choose(hap)
    to_map = pattern.final_mapper(bd.build_chrs, hap)
    lines = bf.lines.choose(hap)
    return bed.indexed_bedlines_to_df(lines, to_map)


def filter_sort_bed_main_inner(
    sconf: GiabStrats,
    rk: RefKey,
    bk: BuildKey,
    inputs: list[Path],
    output: Path,
    output_pattern: str,
    f: BuildDataToBed,
    g: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
) -> list[Path]:
    # if no inputs given, assume that the bed file is specified in yaml and
    # we don't need to download/read anything from disk
    if len(inputs) == 0:

        def hap(o: Path, bd: HapBuildData, bf: HapBedCoords) -> list[Path]:
            bed.write_bed(o, build_hap_coords_df(bd, bf))
            return [o]

        def dip1to1(o: Path, bd: Dip1BuildData, bf: Dip1BedCoords) -> list[Path]:
            bed.write_bed(o, build_dip1to1_coords_df(bd, bf))
            return [o]

        def dip1to2(
            o: Double[Path], bd: Dip2BuildData, bf: Dip1BedCoords
        ) -> list[Path]:
            build_dip1to2_coords_df(bd, bf).both(
                lambda df, hap: bed.write_bed(o.choose(hap), df)
            )
            return o.as_list

        def dip2to1(o: Path, bd: Dip1BuildData, bf: Dip2BedCoords) -> list[Path]:
            bed.write_bed(o, build_dip2to1_coords_df(bd, bf))
            return [o]

        def dip2to2(
            o: Path, hap: Haplotype, bd: Dip2BuildData, bf: Dip2BedCoords
        ) -> list[Path]:
            bed.write_bed(o, build_dip2to2_coords_df(hap, bd, bf))
            return [o]

        return sconf.with_build_data_and_bed_o(
            rk,
            bk,
            output,
            output_pattern,
            f,
            hap,
            dip1to1,
            dip1to2,
            dip2to1,
            dip2to2,
        )
    else:

        def _hap(i: Path, o: Path, bd: HapBuildData, bf: HapBedFile) -> list[Path]:
            read_write_filter_sort_hap_bed(i, o, bd, bf, g)
            return [o]

        def _dip1to1(
            i: Path, o: Path, bd: Dip1BuildData, bf: Dip1BedFile
        ) -> list[Path]:
            read_write_filter_sort_dip1to1_bed(i, o, bd, bf, g)
            return [o]

        def _dip1to2(
            i: Path, o: Double[Path], bd: Dip2BuildData, bf: Dip1BedFile
        ) -> list[Path]:
            read_write_filter_sort_dip1to2_bed(i, o, bd, bf, g)
            return o.as_list

        def _dip2to1(
            i: Double[Path], o: Path, bd: Dip1BuildData, bf: Dip2BedFile
        ) -> list[Path]:
            read_write_filter_sort_dip2to1_bed(i, o, bd, bf, g)
            return [o]

        def _dip2to2(
            i: Path, o: Path, hap: Haplotype, bd: Dip2BuildData, bf: Dip2BedFile
        ) -> list[Path]:
            read_write_filter_sort_dip2to2_bed(i, o, hap, bd, bf, g)
            return [o]

        return sconf.with_build_data_and_bed_io(
            rk,
            bk,
            inputs,
            output,
            output_pattern,
            f,
            _hap,
            _dip1to1,
            _dip1to2,
            _dip2to1,
            _dip2to2,
        )


def filter_sort_bed_main(
    f: BuildDataToBed,
    smk: Any,
    g: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
) -> None:
    """Read a bed and filter/sort it appropriately.

    This is meant to be called in snakemake scripts, as this operations is
    very common. 'smk' is the snakemake object and 'f' is a function to
    retrieve the bed configuration from the config instance (which will be
    obtained from the snakemake object).
    """
    sconf: GiabStrats = smk.config
    ws: SmkWildcards = smk.wildcards

    if not isinstance((ins := smk.input), list) and not all(
        isinstance(x, str) for x in ins
    ):
        raise DesignError(f"Inputs must be a list of strings, got {ins}")

    if not isinstance(output_pattern := smk.params["output_pattern"], str):
        raise DesignError(f"Output pattern must be a string, got {output_pattern}")

    filter_sort_bed_main_inner(
        sconf,
        wc_to_refkey(ws),
        wc_to_buildkey(ws),
        [Path(i) for i in ins],
        smk.output[0],
        output_pattern,
        f,
        g,
    )


# doc formatters


def format_bed_params(ps: BedFileParams) -> str:
    cols = ps.bed_cols
    chrom = cols.chr + 1
    start = cols.start + 1
    end = cols.end + 1

    columns_txt = " ".join(
        [
            f"To construct the bed file, columns {chrom}, {start}, and {end}",
            "were selected as the 'chrom', 'start', and 'end' columns",
            "respectively.",
        ]
    )

    skip_txt = (
        None if ps.skip_lines == 0 else f"The first {ps.skip_lines} were skipped."
    )

    offset = (
        "Each coordinate was reduced by 1 since the incoming file was 1-indexed."
        if ps.one_indexed
        else None
    )

    return " ".join([x for x in [columns_txt, skip_txt, offset] if x is not None])


def format_md5(p: Path) -> str:
    h = get_md5(p)
    return f"The md5 hash of the uncompressed file was {h}."


def format_this(this: str | None) -> str:
    return from_maybe("This file", this)


def format_url_src(src: UrlSrcDoc, p: Path, this: str | None) -> str:
    url_txt = f"{format_this(this)} was downloaded from {src.url}."
    md5_txt = format_md5(p)
    return " ".join([x for x in [url_txt, md5_txt, src.comment] if x is not None])


def format_local_src(src: LocalSrcDoc, p: Path, this: str | None) -> str:
    local_txt = (
        f"{format_this(this)} was saved locally on the machine running the pipeline."
    )
    md5_txt = format_md5(p)
    return " ".join([x for x in [local_txt, md5_txt, src.comment] if x is not None])


def format_src(src: SrcDoc, p: Path, this: str | None) -> str:
    if isinstance(src, UrlSrcDoc):
        return format_url_src(src, p, this)
    elif isinstance(src, LocalSrcDoc):
        return format_local_src(src, p, this)
    else:
        assert_never(src)


def readme_fill(s: str) -> str:
    return fill(s, width=80, break_long_words=False)


################################################################################
# Helper classes


class Null(Generic[X]):
    """Dummy to represent "none" of X (in contrast to Single and Double)."""

    @property
    def as_list(self) -> list[X]:
        return []

    def map(self, f: Callable[[X], Y]) -> Null[Y]:
        return Null()


@dataclass(frozen=True)
class Single(Generic[X]):
    """Wrapper class for one thing "X" which generally means a haploid or
    diploid2 source.

    This is mostly to pair with 'Double' (see below) so that we can check the
    union of "Single" and "Double" and operate on either one or two "X"s
    accordingly.

    Contains overloaded key functions so that functions that possess this
    object can get the correct refkey corresponding to "X"
    """

    elem: X

    @property
    def as_list(self) -> list[X]:
        return [self.elem]

    def map(self, f: Callable[[X], Y]) -> Single[Y]:
        return Single(elem=f(self.elem))

    def key(self, rk: RefKey) -> Single[RefKeyFull]:
        return Single(elem=RefKeyFull(rk, None))


@dataclass(frozen=True)
class Double(Generic[X]):
    """Wrapper class for two things "X" which are generally diploid2 sources.

    The two elements can correspond to each haplotype, and methods exist to
    operate on each accordingly.
    """

    pat: X
    mat: X

    @property
    def as_tuple(self) -> tuple[X, X]:
        return (self.pat, self.mat)

    @property
    def as_list(self) -> list[X]:
        return [*self.as_tuple]

    def keys(self, rk: RefKey) -> Double[RefKeyFull]:
        return make_double(lambda h: RefKeyFull(rk, h))

    def map(self, f: Callable[[X], Y]) -> Double[Y]:
        return Double(pat=f(self.pat), mat=f(self.mat))

    def sum(self, y: Double[Y]) -> Double[tuple[X, Y]]:
        return Double((self.pat, y.pat), (self.mat, y.mat))

    def sum2(self, y: Double[Y], z: Double[Z]) -> Double[tuple[X, Y, Z]]:
        return Double((self.pat, y.pat, z.pat), (self.mat, y.mat, z.mat))

    def choose(self, hap: Haplotype) -> X:
        return hap.choose(self.pat, self.mat)

    def both(self, f: Callable[[X, Haplotype], Y]) -> Double[Y]:
        return Double(f(self.pat, Haplotype.PAT), f(self.mat, Haplotype.MAT))

    def both_(self, f: Callable[[X], Y]) -> Double[Y]:
        return self.both(lambda x, _: f(x))


SingleOrDouble = Single[X] | Double[X]
NullOrSingleOrDouble = Null[X] | SingleOrDouble[X]


@dataclass(frozen=True)
class RefKeyFull:
    """Ref key which may or may not have a haplotype appended to it."""

    key: RefKey
    hap: Haplotype | None

    @property
    def strip(self) -> RefKey:
        return RefKey(self.key)

    @property
    def has_hap(self) -> bool:
        return self.hap is not None

    @property
    def as_tuple(self) -> tuple[RefKey, Haplotype | None]:
        return (self.key, self.hap)

    @property
    def name(self) -> RefKeyFullS:
        k, h = self.as_tuple
        return RefKeyFullS(f"{k}.{h.name}" if h is not None else k)


class Haplotype(Enum):
    "One of the human diploid haplotypes. 0 = Paternal, 1 = Maternal"
    PAT: int = 0
    MAT: int = 1

    @classmethod
    def from_name(cls, n: str) -> Self:
        "Build haplotype from a string. Must be exactly 'pat' or 'mat'."
        try:
            return next(i for i in cls if i.name == n)
        except StopIteration:
            raise ValueError(f"could not make haplotype from name '{n}'")

    @property
    def name(self) -> HaplotypeName:
        return HaplotypeName(self.choose("pat", "mat"))

    def choose(self, left: X, right: X) -> X:
        "Do either left (pat) or right (mat) depending on the haplotype."
        if self is Haplotype.PAT:
            return left
        elif self is Haplotype.MAT:
            return right
        else:
            assert_never(self)


@unique
class ChrIndex(Enum):
    """Represents a valid chromosome index.

    Chromosomes are numbered by integers 1-24 (23 and 24 being X and Y
    respectively). These integers reflect the sort order in output bed files.
    """

    # NOTE: these start at 1 not 0 to conincide with the names of (most)
    # the chromosomes
    CHR1: int = 1
    CHR2: int = 2
    CHR3: int = 3
    CHR4: int = 4
    CHR5: int = 5
    CHR6: int = 6
    CHR7: int = 7
    CHR8: int = 8
    CHR9: int = 9
    CHR10: int = 10
    CHR11: int = 11
    CHR12: int = 12
    CHR13: int = 13
    CHR14: int = 14
    CHR15: int = 15
    CHR16: int = 16
    CHR17: int = 17
    CHR18: int = 18
    CHR19: int = 19
    CHR20: int = 20
    CHR21: int = 21
    CHR22: int = 22
    CHRX: int = 23
    CHRY: int = 24

    @classmethod
    def from_name(cls, n: str) -> Self:
        "Build chr index from a string. Must be a valid digit or 'X' or 'Y'"
        try:
            return next(i for i in cls if i.chr_name == n)
        except StopIteration:
            raise ValueError(f"could make chr index from name '{n}'")

    @classmethod
    def from_name_unsafe(cls, n: str) -> Self:
        "Like 'from_name' but raises DesignError"
        try:
            return cls.from_name(n)
        except ValueError as e:
            raise DesignError(e)

    def __init__(self, i: int) -> None:
        "Build chr index from an integer (which must be in [1,24])"
        self.chr_name = ShortChrName("X" if i == 23 else ("Y" if i == 24 else str(i)))

    def to_internal_index(self, hap: Haplotype) -> bed.InternalChrIndex:
        "Convert this index into an integer corresponding to sort order"
        return bed.InternalChrIndex(hap.value * 24 + self.value - 1)

    # TODO this obviously only makes sense for males
    @property
    def xy_to_hap_unsafe(self) -> Haplotype:
        """
        Convert this index to a haplotype given it is either X or Y.

        Throw DesignError if not X or Y.
        """
        return choose_xy_unsafe(self, Haplotype.MAT, Haplotype.PAT)


@unique
class CoreLevel(Enum):
    """A stratification level (eg "GCcontent" or "mappability")

    These are the only "built-in" levels contained within the pipeline.
    Users may add other levels if they wish to include other source
    files, but these must be specified manually (see below).
    """

    FUNCTIONAL = "Functional"
    LOWCOMPLEXITY = "LowComplexity"
    GC = "GCcontent"
    MAPPABILITY = "Mappability"
    SEGDUPS = "SegmentalDuplications"
    UNION = "Union"
    TELOMERES = "Telomere"
    XY = "XY"
    # overlaps with "other" strat categories, needed because this is where
    # the gaps strat will go
    OTHER_DIFFICULT = "OtherDifficult"
    DIPLOID = "Diploid"


# chromosome name conversions


class _NonDivergentConversion:
    """A chromosome name conversion that doesn't involve a split or merge"""

    @property
    def init_mapper(self) -> bed.InitMapper:
        return NotImplemented

    @property
    def final_mapper(self) -> bed.FinalMapper:
        return NotImplemented


@dataclass(frozen=True)
class HapToHapChrConversion(_NonDivergentConversion):
    fromPattern: HapChrPattern
    toPattern: HapChrPattern
    indices: BuildChrs

    # NOTE dummy haplotype used here, the only reason we chose PAT is because
    # it is numerically zero and thus makes downstream calculations work.
    # This is obviously meaningless for haploid case
    @property
    def init_mapper(self) -> bed.InitMapper:
        return self.fromPattern.init_mapper(self.indices, Haplotype.PAT)

    @property
    def final_mapper(self) -> bed.FinalMapper:
        return self.toPattern.final_mapper(self.indices, Haplotype.PAT)


@dataclass(frozen=True)
class DipToDipChrConversion(_NonDivergentConversion):
    fromPattern: DipChrPattern
    toPattern: DipChrPattern
    indices: BuildChrs

    @property
    def init_mapper(self) -> bed.InitMapper:
        return self.fromPattern.init_mapper(self.indices)

    @property
    def final_mapper(self) -> bed.FinalMapper:
        return self.toPattern.final_mapper(self.indices)

    def split(self, h: Haplotype) -> HapToHapChrConversion:
        return HapToHapChrConversion(
            self.fromPattern.to_hap_pattern(h),
            self.toPattern.to_hap_pattern(h),
            self.indices,
        )


@dataclass(frozen=True)
class HapToDipChrConversion:
    fromPattern: Double[HapChrPattern]
    toPattern: DipChrPattern
    indices: BuildChrs

    @property
    def init_mapper(self) -> Double[bed.InitMapper]:
        return self.fromPattern.both(lambda p, h: p.init_mapper(self.indices, h))

    @property
    def final_mapper(self) -> bed.FinalMapper:
        return self.toPattern.final_mapper(self.indices)


@dataclass(frozen=True)
class DipToHapChrConversion:
    fromPattern: DipChrPattern
    toPattern: Double[HapChrPattern]
    indices: BuildChrs

    @property
    def init_mapper(self) -> tuple[bed.InitMapper, bed.SplitMapper]:
        im = self.fromPattern.init_mapper(self.indices)
        fm0 = self.toPattern.pat.final_mapper(self.indices, Haplotype.PAT)
        return (im, bed.make_split_mapper(im, fm0))

    @property
    def final_mapper(self) -> Double[bed.FinalMapper]:
        return self.toPattern.both(lambda p, h: p.final_mapper(self.indices, h))


class ChrData(NamedTuple):
    idx: bed.InternalChrIndex
    name: bed.ChrName
    shortname: ShortChrName
    haplotype: Haplotype


# tuples representing file paths for the pipeline


class DataLogDirs(NamedTuple):
    data: Path
    log: Path


class DataLogBenchDirs(NamedTuple):
    data: Path
    log: Path
    bench: Path


class FilterSortDirs(NamedTuple):
    data: Path
    bench: Path
    log: Path
    subbed: Path


class BedInterDirs(NamedTuple):
    filtersort: FilterSortDirs
    postsort: DataLogBenchDirs


class BedDirs(NamedTuple):
    src: DataLogDirs
    inter: BedInterDirs
    final: Callable[[str], Path]
    readme: Path


class RefInterDirs(NamedTuple):
    prebuild: DataLogBenchDirs
    filtersort: FilterSortDirs
    build: DataLogBenchDirs


class RefSrcDirs(NamedTuple):
    benchmark: DataLogDirs
    reference: DataLogDirs


class RefDirs(NamedTuple):
    src: RefSrcDirs
    inter: RefInterDirs


# Documentation classes


class LocalSrcDoc(NamedTuple):
    comment: str


class UrlSrcDoc(NamedTuple):
    comment: str | None
    url: str


SrcDoc = LocalSrcDoc | UrlSrcDoc


class BedDoc(NamedTuple):
    params: BedFileParams
    bed: Single[SrcDoc] | Double[SrcDoc]


# Rule aggregation classes


@dataclass(frozen=True)
class _HasFinalBeds:
    readme: Path

    @property
    def all_final(self) -> list[Path]:
        return NotImplemented

    @property
    def all_final_str(self) -> list[str]:
        return [str(x) for x in self.all_final]


class _HasSources:
    @property
    def all_sources(self) -> list[Path]:
        return NotImplemented


@dataclass(frozen=True)
class MutualPathPair:
    positive: Path
    negative: Path

    def both(self, f: Callable[[Path], Path]) -> MutualPathPair:
        return MutualPathPair(f(self.positive), f(self.negative))

    @property
    def paths(self) -> list[Path]:
        return [self.positive, self.negative]


@dataclass(frozen=True)
class UniformRepeatPaths:
    perfect: list[Path]  # ASSUME this will always be non-empty
    imperfect: list[Path]  # ASSUME this will always be non-empty
    homopolymers: MutualPathPair

    @property
    def _all_final(self) -> list[Path]:
        return [*self.perfect, *self.imperfect, *self.homopolymers.paths]


@dataclass(frozen=True)
class RepeatsPaths(_HasSources):
    trf_src: Path0or1or2
    rmsk_src: Path0or1or2

    filtered_trs: list[Path]  # ASSUME this is non-empty
    all_trs: MutualPathPair
    all_repeats: MutualPathPair

    @property
    def all_final(self) -> list[Path]:
        return [
            *self.filtered_trs,
            *self.all_trs.paths,
            *self.all_repeats.paths,
        ]


@dataclass(frozen=True)
class SatellitesPaths:
    sat_src: Path0or1or2

    sats: MutualPathPair

    used_censat: bool
    all_repeats: RepeatsPaths | None

    @property
    def _all_final(self) -> list[Path]:
        return [
            *self.sats.paths,
            *(r.all_final if (r := self.all_repeats) is not None else []),
        ]


@dataclass(frozen=True)
class LowComplexityPaths(_HasSources, _HasFinalBeds):
    uniform_repeats: UniformRepeatPaths
    satellites: SatellitesPaths | None

    @property
    def all_sources(self) -> list[Path]:
        return fmap_maybe_def(
            [],
            lambda s: fmap_maybe_def(
                s.sat_src.as_list,
                lambda r: r.trf_src.as_list
                + s.sat_src.as_list
                + (r.rmsk_src.as_list if s.used_censat else []),
                s.all_repeats,
            ),
            self.satellites,
        )

    @property
    def all_final(self) -> list[Path]:
        return [
            x
            for x in (
                self.uniform_repeats._all_final
                + (self.satellites._all_final if self.satellites is not None else [])
            )
        ]


@dataclass(frozen=True)
class XYFeaturePaths:
    src: Path

    bed: XYFile
    xtr_path: Path | None
    ampliconic_path: Path | None
    xtr: str | None
    ampliconic: str | None


@dataclass(frozen=True)
class PARPaths:
    path: MutualPathPair
    doc: str


@dataclass(frozen=True)
class SubSexPaths:
    par: PARPaths | None
    features: XYFeaturePaths | None

    @property
    def all_inputs(self) -> list[Path]:
        return fmap_maybe_def([], lambda z: [z.src], self.features)

    @property
    def par_path(self) -> Path | None:
        return fmap_maybe(lambda z: z.path.positive, self.par)

    @property
    def nonpar_path(self) -> Path | None:
        return fmap_maybe(lambda z: z.path.negative, self.par)

    @property
    def xtr_path(self) -> Path | None:
        return fmap_maybe(lambda z: z.xtr_path, self.features)

    @property
    def ampliconic_path(self) -> Path | None:
        return fmap_maybe(lambda z: z.ampliconic_path, self.features)

    @property
    def all_output_paths(self) -> list[Path]:
        return [
            x
            for x in [
                self.par_path,
                self.nonpar_path,
                self.xtr_path,
                self.ampliconic_path,
            ]
            if x is not None
        ]


@dataclass(frozen=True)
class _SexPaths:
    @property
    def all_paths(self) -> list[SubSexPaths]:
        return NotImplemented

    @property
    def par_paths(self) -> list[Path]:
        return [p for x in self.all_paths if (p := x.par_path) is not None]

    @property
    def nonpar_paths(self) -> list[Path]:
        return [p for x in self.all_paths if (p := x.nonpar_path) is not None]

    @property
    def xtr_paths(self) -> list[Path]:
        return [p for x in self.all_paths if (p := x.xtr_path) is not None]

    @property
    def ampliconic_paths(self) -> list[Path]:
        return [p for x in self.all_paths if (p := x.ampliconic_path) is not None]

    @property
    def all_features(self) -> list[Path]:
        return self.xtr_paths + self.ampliconic_paths

    @property
    def all_inputs(self) -> list[Path]:
        return [i for x in self.all_paths for i in x.all_inputs]

    @property
    def all_output_paths(self) -> list[Path]:
        return [o for x in self.all_paths for o in x.all_output_paths]


@dataclass(frozen=True)
class MaleHapSexPaths(_SexPaths):
    x: SubSexPaths | None
    y: SubSexPaths | None

    @property
    def all_paths(self) -> list[SubSexPaths]:
        return [x for x in [self.x, self.y] if x is not None]


# TODO when we need it...
# class FemaleHapSexPaths(NamedTuple):
#     x: SubSexPaths


@dataclass(frozen=True)
class Dip1SexPaths(_SexPaths):
    sex1: SubSexPaths | None  # X
    sex2: SubSexPaths | None  # Y if male

    @property
    def all_paths(self) -> list[SubSexPaths]:
        return [x for x in [self.sex1, self.sex2] if x is not None]


@dataclass(frozen=True)
class Dip2SexPaths(_SexPaths):
    paths: SubSexPaths | None
    hap: Haplotype

    @property
    def all_paths(self) -> list[SubSexPaths]:
        return fmap_maybe_def([], lambda z: [z], self.paths)


AnySexPaths = MaleHapSexPaths | Dip1SexPaths | Dip2SexPaths


@dataclass(frozen=True)
class SexPaths(_HasSources, _HasFinalBeds):
    sex: AnySexPaths
    auto: Path | None

    @property
    def all_sources(self) -> list[Path]:
        return self.sex.all_inputs

    @property
    def all_final(self) -> list[Path]:
        return self.sex.all_output_paths + maybe_to_list(self.auto)


@dataclass(frozen=True)
class LowComplexitySources:
    trf: Path0or1or2 | None
    sat: Path0or1or2 | None
    rmsk: Path0or1or2 | None

    @property
    def trf_sources(self) -> list[Path]:
        if self.trf is None:
            raise DesignError()
        else:
            return self.trf.as_list

    @property
    def sat_sources(self) -> list[Path]:
        if self.sat is None:
            raise DesignError()
        else:
            return self.sat.as_list

    @property
    def rmsk_sources(self) -> list[Path]:
        if self.rmsk is None:
            raise DesignError()
        else:
            return self.rmsk.as_list


Path1or2 = SingleOrDouble[Path]
Path0or1or2 = NullOrSingleOrDouble[Path]

RefKeyFull1or2 = SingleOrDouble[RefKeyFull]
RefKeyFullS1or2 = SingleOrDouble[RefKeyFullS]


@dataclass(frozen=True)
class OtherDifficultSources:
    gaps: Path0or1or2 | None
    refseq: Path0or1or2 | None
    vdj: Path0or1or2 | None
    kir: Path0or1or2 | None
    mhc: Path0or1or2 | None
    other: dict[OtherLevelKey, dict[OtherStratKey, Path0or1or2]]

    # NOTE return empty lists here to avoid failure when calling in rules
    @property
    def gaps_sources(self) -> list[Path]:
        return self.gaps.as_list if self.gaps is not None else []

    @property
    def refseq_sources(self) -> list[Path]:
        return self.refseq.as_list if self.refseq is not None else []

    @property
    def vdj_sources(self) -> list[Path]:
        return self.vdj.as_list if self.vdj is not None else []

    @property
    def kir_sources(self) -> list[Path]:
        return self.kir.as_list if self.kir is not None else []

    @property
    def mhc_sources(self) -> list[Path]:
        return self.mhc.as_list if self.mhc is not None else []

    @property
    def other_sources(self) -> list[Path]:
        return [i for x in self.other.values() for y in x.values() for i in y.as_list]

    def other_source(self, lk: OtherLevelKey, sk: OtherStratKey) -> list[Path]:
        try:
            return self.other[lk][sk].as_list
        except KeyError:
            return []


@dataclass(frozen=True)
class FunctionalPaths(_HasSources, _HasFinalBeds):
    cds_source: Path0or1or2
    cds_output: MutualPathPair | None

    @property
    def all_sources(self) -> list[Path]:
        return self.cds_source.as_list

    @property
    def all_final(self) -> list[Path]:
        return fmap_maybe_def([], lambda z: z.paths, self.cds_output)


@dataclass(frozen=True)
class SourceOutputPaths:
    source: Path0or1or2
    output: Path


@dataclass(frozen=True)
class ImmunoPaths(SourceOutputPaths):
    source: Path0or1or2
    output: Path
    source_is_refseq: bool


@dataclass(frozen=True)
class OtherDifficultPaths(_HasSources, _HasFinalBeds):
    gaps: SourceOutputPaths | None
    vdj: ImmunoPaths | None
    mhc: ImmunoPaths | None
    kir: ImmunoPaths | None
    other: dict[OtherStratKey, SourceOutputPaths]

    @property
    def all_sources(self) -> list[Path]:
        return [
            p
            for x in [self.gaps, self.vdj, self.mhc, self.kir, *self.other.values()]
            if x is not None
            for p in x.source.as_list
        ]

    @property
    def all_final(self) -> list[Path]:
        return [
            p.output
            for p in [
                self.gaps,
                self.vdj,
                self.mhc,
                self.kir,
                *self.other.values(),
            ]
            if p is not None
        ]


@dataclass(frozen=True)
class MiscPaths(_HasSources, _HasFinalBeds):
    paths: dict[OtherStratKey, SourceOutputPaths]
    desc: OtherLevelDescription

    @property
    def all_sources(self) -> list[Path]:
        return [p for x in self.paths.values() for p in x.source.as_list]

    @property
    def all_final(self) -> list[Path]:
        return [x.output for x in self.paths.values()]


@dataclass(frozen=True)
class LowmapPaths(_HasFinalBeds):
    union: MutualPathPair
    single: list[Path]
    params: list[LowMapParams]  # ASSUME non empty

    @property
    def all_final(self) -> list[Path]:
        return self.union.paths + self.single


@dataclass(frozen=True)
class GCPaths(_HasFinalBeds):
    lowGC: Path
    middleGC: list[Path]
    highGC: Path
    extremes: list[Path]
    params: GCParams

    @property
    def all_final(self) -> list[Path]:
        return [self.lowGC, *self.middleGC, self.highGC, *self.extremes]


@dataclass(frozen=True)
class TelomerePaths(_HasFinalBeds, _HasSources):
    telomeres: Path

    # dummy property
    @property
    def all_sources(self) -> list[Path]:
        return []

    @property
    def all_final(self) -> list[Path]:
        return [self.telomeres]


@dataclass(frozen=True)
class SegdupSources(_HasSources):
    superdup: Path0or1or2 | None

    @property
    def all_sources(self) -> list[Path]:
        return self.superdup.as_list if self.superdup is not None else []


@dataclass(frozen=True)
class SegdupPaths(_HasSources, _HasFinalBeds):
    superdups: Path0or1or2
    all_segdups: MutualPathPair
    long_segdups: MutualPathPair

    @property
    def all_sources(self) -> list[Path]:
        return self.superdups.as_list

    @property
    def all_final(self) -> list[Path]:
        return self.all_segdups.paths + self.long_segdups.paths


@dataclass(frozen=True)
class DiploidPaths(_HasFinalBeds, _HasSources):
    hets: list[Path]
    SNVorSV_hets: list[Path]
    homs: list[Path]
    SNVorSV_homs: list[Path]

    nonpar: list[Path]

    # dummy property
    @property
    def all_sources(self) -> list[Path]:
        return []

    @property
    def non_par(self) -> list[Path]:
        return self.nonpar

    @property
    def all_final(self) -> list[Path]:
        return self.hets + self.homs + self.SNVorSV_hets + self.SNVorSV_homs


@dataclass(frozen=True)
class SegdupLowmapPaths:
    segdup_input: SegdupPaths
    lowmap_input: Path
    output: MutualPathPair

    @property
    def _all_inputs(self) -> list[Path]:
        return [self.lowmap_input, self.segdup_input.all_segdups.positive]

    @property
    def _all_final(self) -> list[Path]:
        return self.output.paths


# ASSUME at least two of the sources are non-null
@dataclass(frozen=True)
class AllDifficultPaths:
    gc_input: Path | None
    repeat_input: Path | None
    xy_inputs: list[Path]
    output: MutualPathPair

    @property
    def _all_inputs(self) -> list[Path]:
        return [
            p for p in [self.gc_input, self.repeat_input] if p is not None
        ] + self.xy_inputs

    @property
    def _all_final(self) -> list[Path]:
        return self.output.paths


@dataclass(frozen=True)
class UnionPaths(_HasFinalBeds, _HasSources):
    segdup_lowmap: SegdupLowmapPaths
    all_difficult: AllDifficultPaths | None

    # dummy property
    @property
    def all_sources(self) -> list[Path]:
        return []

    @property
    def segdup_lowmap_inputs(self) -> list[Path]:
        return self.segdup_lowmap._all_inputs

    @property
    def all_difficult_inputs(self) -> list[Path]:
        empty: list[Path] = []
        return fmap_maybe_def(empty, lambda x: x._all_inputs, self.all_difficult) + [
            self.segdup_lowmap.output.positive
        ]

    @property
    def all_inputs(self) -> list[Path]:
        return self.segdup_lowmap_inputs + self.all_difficult_inputs

    @property
    def segdup_lowmap_final(self) -> list[Path]:
        return self.segdup_lowmap.output.paths

    @property
    def all_difficult_final(self) -> list[Path]:
        empty: list[Path] = []
        return fmap_maybe_def(empty, lambda x: x.output.paths, self.all_difficult)

    @property
    def all_final(self) -> list[Path]:
        return self.segdup_lowmap_final + self.all_difficult_final


################################################################################
# Constants

CHR_INDEX_PLACEHOLDER = "%i"
CHR_HAP_PLACEHOLDER = "%h"

MHC_CHR = ChrIndex(6)
KIR_CHR = ChrIndex(19)

VDJ_CHRS = {ChrIndex(i) for i in [2, 7, 14, 22]}

# strats in "OtherDifficult" that are built-in and should not be included
# manually using the "other_strats" directive in "build"
BUILTIN_OTHER = {"VDJ", "KIR", "MHC", "gaps_slop15kb"}

OTHERDIFF_KEY = OtherLevelKey("OtherDifficult")


################################################################################
# Snakemake configuration model


class BaseModel(BaseModel_):
    class Config:
        frozen = True
        extra = "forbid"


class GenericModel(GenericModel_):
    class Config:
        frozen = True
        extra = "forbid"

    # dirty hack to get pickling to work for generic model types; see
    # https://github.com/pydantic/pydantic/issues/1667
    #
    # it seems this was fixed, but there might be some issue with getting
    # snakemake to recognize to get its paths correct
    def __class_getitem__(
        cls: Type[GenericModelT], params: Union[Type[Any], tuple[Type[Any], ...]]
    ) -> Type[Any]:
        created_class = super().__class_getitem__(params)
        setattr(
            sys.modules[created_class.__module__], created_class.__name__, created_class
        )
        return created_class


class Documentable:
    """Means to provide documentation for a file source"""

    @property
    def documentation(self) -> SrcDoc:
        return NotImplemented


class Documentable1:
    """Means to provide documentation for one source file"""

    @property
    def documentation(self) -> Single[SrcDoc]:
        return NotImplemented


class Documentable2:
    """Means to provide documentation for two source files"""

    @property
    def documentation(self) -> Double[SrcDoc]:
        return NotImplemented


class BaseDocumentable(BaseModel, Documentable):
    """BaseModel class which has a source and is documented"""

    pass


class GenericDocumentable1(GenericModel, Documentable1):
    """BaseModel class which has one source file and is documented"""

    pass


class GenericDocumentable2(GenericModel, Documentable2):
    """BaseModel class which has two source files and is documented"""

    pass


S = TypeVar("S", bound=BaseDocumentable)

Q = TypeVar(
    "Q",
    bound=GenericDocumentable1 | GenericDocumentable2,
)


class Diploid(GenericModel, Generic[X]):
    """A diploid thing"""

    pat: X
    mat: X

    @property
    def double(self) -> Double[X]:
        return Double(pat=self.pat, mat=self.mat)


class ChrPattern:
    """A general chromosome pattern providing interface to convert indices to
    names."""

    def to_names(self, cs: BuildChrs) -> OrderedHapChrNames:
        return NotImplemented


class HapChrPattern(BaseModel, ChrPattern):
    """Chromosome pattern for a haploid file.

    'template' contains a placeholder for the chromosome index.
    """

    template: str = "chr%i"
    special: dict[ChrIndex, bed.ChrName] = {}
    exclusions: set[ChrIndex] = set()

    @validator("template")
    def is_valid_template(cls, v: str) -> str:
        assert v.count(CHR_INDEX_PLACEHOLDER) == 1, "chr template must have '%i' in it"
        return v

    def _is_excluded(self, i: ChrIndex) -> bool:
        return i in self.exclusions

    def filter_indices(self, cs: BuildChrs) -> HapChrs:
        return HapChrs({i for i in cs if not self._is_excluded(i)})

    def to_chr_name(self, i: ChrIndex) -> bed.ChrName | None:
        if self._is_excluded(i):
            return None
        elif i in self.special:
            return self.special[i]
        else:
            return bed.ChrName(
                self.template.replace(
                    CHR_INDEX_PLACEHOLDER,
                    str(i.chr_name),
                )
            )

    def to_chr_data(self, cs: BuildChrs, h: Haplotype) -> list[ChrData]:
        return [
            ChrData(c.to_internal_index(h), n, c.chr_name, h)
            for c in sort_chr_indices(self.filter_indices(cs))
            if (n := self.to_chr_name(c)) is not None
        ]

    def to_names(self, cs: BuildChrs) -> OrderedHapChrNames:
        # NOTE: the haplotype argument is doing nothing since it is only
        # used to make the index which I remove before returning here
        return OrderedHapChrNames(
            [bed.ChrName(x[1]) for x in self.to_chr_data(cs, Haplotype.PAT)]
        )

    def init_mapper(self, cs: BuildChrs, hap: Haplotype) -> bed.InitMapper:
        return {d.name: d.idx for d in self.to_chr_data(cs, hap)}

    def final_mapper(self, cs: BuildChrs, hap: Haplotype) -> bed.FinalMapper:
        return {d.idx: d.name for d in self.to_chr_data(cs, hap)}


class DipChrPattern(BaseModel, ChrPattern):
    """Chromosome pattern for a haploid file.

    'template' contains placeholders for both the bare chromosome name/index
    and the haplotype, which maps to a specific haplotype index via 'hapnames'.
    """

    template: str = "chr%i_%h"
    special: dict[ChrIndex, bed.ChrName] = {}
    hapnames: Diploid[HaplotypeName] = Diploid(
        pat=HaplotypeName("PATERNAL"),
        mat=HaplotypeName("MATERNAL"),
    )
    # By default, paternal doesn't have X and maternal doesn't have Y
    exclusions: Diploid[set[ChrIndex]] = Diploid(
        pat={ChrIndex.CHRX},
        mat={ChrIndex.CHRY},
    )

    @validator("template")
    def is_valid_template(cls, v: str) -> str:
        assert (
            v.count(CHR_INDEX_PLACEHOLDER) == 1 and v.count(CHR_HAP_PLACEHOLDER) == 1
        ), "chr template must have '%i' and '%h' in it"
        return v

    @validator("hapnames")
    def is_valid_hapname(cls, v: Diploid[HaplotypeName]) -> Diploid[HaplotypeName]:
        def is_valid(n: HaplotypeName, h: Haplotype) -> None:
            t = CHR_INDEX_PLACEHOLDER in n or CHR_HAP_PLACEHOLDER in n
            assert not t, f"name for {h.name} must not have '%i' and '%h' in it"

        v.double.both(is_valid)
        return v

    def _is_excluded(self, i: ChrIndex, h: Haplotype) -> bool:
        return i in self.exclusions.double.choose(h)

    def filter_indices(self, cis: BuildChrs, h: Haplotype) -> HapChrs:
        return HapChrs({i for i in cis if not self._is_excluded(i, h)})

    def to_chr_name(self, i: ChrIndex, h: Haplotype) -> bed.ChrName | None:
        if self._is_excluded(i, h):
            return None
        elif i in self.special:
            return self.special[i]
        else:
            name = self.hapnames.double.choose(h)
            return bed.ChrName(
                self.template.replace(
                    CHR_INDEX_PLACEHOLDER,
                    str(i.chr_name),
                ).replace(
                    CHR_HAP_PLACEHOLDER,
                    name,
                )
            )

    def to_chr_data(self, cs: BuildChrs) -> list[ChrData]:
        # order is really important here; we want to iterate through the first
        # haplotype before the second so that the chromosome order is like
        # chr1_mat, chr2_mat ... chr1_pat, chr2_pat rather than chr1_mat,
        # chr1_pat ... etc
        return [
            ChrData(c.to_internal_index(h), n, c.chr_name, h)
            for h in Haplotype
            for c in sort_chr_indices(self.filter_indices(cs, h))
            if (n := self.to_chr_name(c, h)) is not None
        ]

    def to_names(self, cs: BuildChrs) -> OrderedHapChrNames:
        return OrderedHapChrNames([bed.ChrName(x[1]) for x in self.to_chr_data(cs)])

    def init_mapper(self, cs: BuildChrs) -> bed.InitMapper:
        return {d.name: d.idx for d in self.to_chr_data(cs)}

    def final_mapper(self, cs: BuildChrs) -> bed.FinalMapper:
        return {d.idx: d.name for d in self.to_chr_data(cs)}

    def to_hap_pattern(self, hap: Haplotype) -> HapChrPattern:
        hs = self.hapnames.double.choose(hap)
        return HapChrPattern(
            template=self.template.replace(CHR_HAP_PLACEHOLDER, hs),
            special=self.special,
            exclusions=self.exclusions.double.choose(hap),
        )


class HapSrc(GenericDocumentable1, Generic[S]):
    """Specification for a haploid source file."""

    chr_pattern_: HapChrPattern = Field(HapChrPattern(), alias="chr_pattern")
    hap: S

    def chr_conversion(
        self, fromChr: HapChrPattern, cis: BuildChrs
    ) -> HapToHapChrConversion:
        """Create a chromosome names conversion corresponding to 'fromChr'.

        'cis' is the list of chromosome indices that the conversion will
        consider not excluded.
        """
        return HapToHapChrConversion(fromChr, self.chr_pattern, cis)

    def noop_conversion(self, cis: BuildChrs) -> HapToHapChrConversion:
        """Create a chromosome conversion for this source itself.

        Useful for anything generated via the reference itself, since those
        will have chromosome names corresponding directly to the reference.
        """
        return self.chr_conversion(self.chr_pattern, cis)

    def hap_chrs(self, cis: BuildChrs) -> HapChrs:
        return self.chr_pattern.filter_indices(cis)

    @property
    def documentation(self) -> Single[SrcDoc]:
        return Single(elem=self.hap.documentation)

    @property
    def src(self) -> Single[S]:
        return Single(elem=self.hap)

    @property
    def chr_pattern(self) -> HapChrPattern:
        return self.chr_pattern_


class Dip1Src(GenericDocumentable1, Generic[S]):
    """Specification for a combined (diploid1) source file.

    The 'src' is assumed to have all chromosomes for both haplotypes in one
    file, which implies they are labeled so as to distinguish the haps. The
    pattern will match both the chromosome number and the haplotype within the
    chromosome name.
    """

    chr_pattern_: DipChrPattern = Field(DipChrPattern(), alias="chr_pattern")
    dip: S

    @property
    def documentation(self) -> Single[SrcDoc]:
        return Single(elem=self.dip.documentation)

    @property
    def src(self) -> Single[S]:
        return Single(elem=self.dip)

    @property
    def chr_pattern(self) -> DipChrPattern:
        return self.chr_pattern_

    def hap_chr_conversion(
        self,
        fromChr: Double[HapChrPattern],
        cis: BuildChrs,
    ) -> HapToDipChrConversion:
        """Create a dip2->dip1 conversion corresponding to 'fromChr'.

        'cis' is the list of chromosome indices that the conversion will
        consider not excluded.
        """
        return HapToDipChrConversion(fromChr, self.chr_pattern, cis)

    def dip_chr_conversion(
        self,
        fromChr: DipChrPattern,
        cis: BuildChrs,
    ) -> DipToDipChrConversion:
        """Create a dip1->dip1 conversion corresponding to 'fromChr'.

        'cis' is the list of chromosome indices that the conversion will
        consider not excluded.
        """
        return DipToDipChrConversion(fromChr, self.chr_pattern, cis)

    def noop_conversion(self, cis: BuildChrs) -> DipToDipChrConversion:
        """Create a chromosome conversion for this source itself.

        Useful for anything generated via the reference itself, since those
        will have chromosome names corresponding directly to the reference.
        """
        return self.dip_chr_conversion(self.chr_pattern, cis)

    def hap_chrs(self, cis: BuildChrs, h: Haplotype) -> HapChrs:
        return self.chr_pattern.filter_indices(cis, h)

    def all_chrs(self, cis: BuildChrs) -> HapChrs:
        return HapChrs(
            self.hap_chrs(cis, Haplotype.PAT) | self.hap_chrs(cis, Haplotype.MAT)
        )


class Dip2Src(GenericDocumentable2, Generic[S]):
    """Specification for split (diploid2) source files.

    Each source may or may not have each haplotype labeled; the identity of each
    haplotype in either source file is determined based on the configuration key
    under which it appears (hap1 or hap2) and the chromosome names for each are
    matched according to its corresponding entry in `chr_pattern`.
    """

    # TODO this could be cleaner (don't make one hap nested and the other flat)
    chr_pattern_: Diploid[HapChrPattern] = Field(
        Diploid(
            pat=HapChrPattern(
                template="chr%i_PATERNAL",
                exclusions=[ChrIndex.CHRX],
            ),
            mat=HapChrPattern(
                template="chr%i_MATERNAL",
                exclusions=[ChrIndex.CHRY],
            ),
        ),
        alias="chr_pattern",
    )
    pat: S
    mat: S

    @property
    def documentation(self) -> Double[SrcDoc]:
        return self.src.map(lambda x: x.documentation)

    @property
    def src(self) -> Double[S]:
        return Double(pat=self.pat, mat=self.mat)

    @property
    def chr_pattern(self) -> Double[HapChrPattern]:
        return self.chr_pattern_.double

    def hap_chr_conversion(
        self,
        fromChr: Double[HapChrPattern],
        cis: BuildChrs,
    ) -> Double[HapToHapChrConversion]:
        """Create a dip2->dip2 conversion corresponding to 'fromChr'.

        'cis' is the list of chromosome indices that the conversion will
        consider not excluded.
        """
        toChr = self.chr_pattern
        return Double(
            HapToHapChrConversion(fromChr.pat, toChr.pat, cis),
            HapToHapChrConversion(fromChr.mat, toChr.mat, cis),
        )

    def dip_chr_conversion(
        self,
        fromChr: DipChrPattern,
        cis: BuildChrs,
    ) -> DipToHapChrConversion:
        """Create a dip1->dip2 conversion corresponding to 'fromChr'.

        'cis' is the list of chromosome indices that the conversion will
        consider not excluded.
        """
        return DipToHapChrConversion(fromChr, self.chr_pattern, cis)

    def noop_conversion(self, cis: BuildChrs) -> Double[HapToHapChrConversion]:
        """Create a chromosome conversion for this source itself.

        Useful for anything generated via the reference itself, since those
        will have chromosome names corresponding directly to the reference.
        """
        return self.hap_chr_conversion(self.chr_pattern, cis)

    def hap_chrs(self, cis: BuildChrs, h: Haplotype) -> HapChrs:
        return self.chr_pattern.choose(h).filter_indices(cis)


class BedCoord(BaseModel):
    """Represents one line in a bed file, encoded directly in the yaml config.

    'chr', 'start', and 'end' represent the mandatory bed coordindates. 'more'
    is a list of additional columns should they be necessary. Note that when
    referencing these columns later, 'chr', 'start' and 'end' will always be
    columns 0-2 and 'more' will start at column 3 and increase in the order
    presented in the list.
    """

    chr: ChrIndex
    start: NonNegativeInt
    end: NonNegativeInt

    @validator("end")
    def positive_region(cls, end: int, values: dict[str, Any]) -> int:
        try:
            start: int = values["start"]
            assert end > start, "End must be greater than start"
        except KeyError:
            pass
        return end


class HapBedCoord(BedCoord):
    def line(self, h: Haplotype) -> bed.IndexedBedLine:
        return bed.IndexedBedLine(
            self.chr.to_internal_index(h),
            self.start,
            self.end,
        )


class DipBedCoord(BedCoord):
    hap: Haplotype

    @property
    def line(self) -> bed.IndexedBedLine:
        return bed.IndexedBedLine(
            self.chr.to_internal_index(self.hap),
            self.start,
            self.end,
        )


BedCoordT = TypeVar("BedCoordT", HapBedCoord, DipBedCoord)


class BedLines(GenericModel, Generic[BedCoordT], Documentable):
    """Bed file encoded directly in the configuration file."""

    lines: list[BedCoordT]
    provenance: str = "This bed file was specified manually in the yaml config."


HapBedLines = BedLines[HapBedCoord]
DipBedLines = BedLines[DipBedCoord]


class HapBedCoords(BaseModel):
    """Coordinates for haploid bed file."""

    hap: HapBedLines

    def lines(self, h: Haplotype) -> Single[list[bed.IndexedBedLine]]:
        return Single(elem=[x.line(h) for x in self.hap.lines])

    @property
    def lines_nohap(self) -> Single[list[bed.IndexedBedLine]]:
        return self.lines(Haplotype.PAT)


class Dip1BedCoords(BaseModel):
    """Coordinates for combined diploid bed file."""

    dip: DipBedLines

    @property
    def lines(self) -> Single[list[bed.IndexedBedLine]]:
        return Single(elem=[x.line for x in self.dip.lines])


class Dip2BedCoords(Diploid[HapBedLines]):
    """Coordinates for paternal and maternal bed files (seperate)."""

    pat: HapBedLines
    mat: HapBedLines

    @property
    def lines(self) -> Double[list[bed.IndexedBedLine]]:
        return self.double.both(lambda p, hap: [c.line(hap) for c in p.lines])


# TODO make this mandatory
class HashedSrc_(BaseDocumentable):
    """A source that may be hashed to verify its integrity"""

    md5: str | None = None


class FileSrc_(HashedSrc_):
    """Source for a local file"""

    filepath: FilePath
    comment: str = "This file was local on the filesystem when the pipeline was run."

    @property
    def documentation(self) -> SrcDoc:
        return LocalSrcDoc(comment=self.comment)


class BedLocalSrc(FileSrc_):
    """Source for a local bed file in gzip format"""

    @validator("filepath")
    def is_gzip(cls, v: FilePath) -> FilePath:
        assert is_gzip(v), "not in gzip format"
        return v


class RefFileSrc(FileSrc_):
    """Source for a local fasta file in bgzip format"""

    @validator("filepath")
    def path_is_bgzip(cls, v: FilePath) -> FilePath:
        assert is_bgzip(v), "not in bgzip format"
        return v


class HttpSrc_(HashedSrc_):
    """Source for a remote file to be downloaded via HTTP"""

    url: HttpUrl
    comment: str | None = None

    @property
    def documentation(self) -> SrcDoc:
        return UrlSrcDoc(url=self.url, comment=self.comment)


class BedHttpSrc(HttpSrc_):
    """Source for a remote bed file to be downloaded via HTTP."""

    pass


class RefHttpSrc(HttpSrc_):
    """Source for a remote fasta file to be downloaded via HTTP."""

    pass


class Paths(BaseModel):
    """Local build paths for snakemake."""

    resources: Path = Path("resources")
    results: Path = Path("results")


# TODO add checksum for GEM
class Tools(BaseModel):
    """Urls for tools to download/build/use in the pipeline."""

    kent: HttpUrl = "https://github.com/ucscGenomeBrowser/kent/archive/refs/tags/v462_base.tar.gz"  # type: ignore
    repseq: HttpUrl = "https://github.com/usnistgov/giab-repseq/archive/refs/tags/v1.1.0.tar.gz"  # type: ignore
    paftools: HttpUrl = "https://raw.githubusercontent.com/lh3/minimap2/e28a55be86b298708a2a67c924d665a00b8d829c/misc/paftools.js"  # type: ignore
    dipcall_aux: HttpUrl = "https://raw.githubusercontent.com/lh3/dipcall/6bd5d7724699491f215aeb5fb628490ebf2cc3ae/dipcall-aux.js"  # type: ignore
    gemlib: HttpUrl = "https://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%203/GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2/download"  # type: ignore
    gemhash: str = "d0ae279b249d675d53aa23d2b7d60e16"


class BedColumns(BaseModel):
    """Denotes coordinate columns in a bed file (0-indexed)."""

    chr: NonNegativeInt = 0
    start: NonNegativeInt = 1
    end: NonNegativeInt = 2

    @validator("start")
    def start_different(
        cls,
        v: NonNegativeInt,
        values: dict[str, int],
    ) -> NonNegativeInt:
        try:
            assert values["chr"] != v, "Bed columns must be different"
        except KeyError:
            pass
        return v

    @validator("end")
    def end_different(
        cls,
        v: NonNegativeInt,
        values: dict[str, int],
    ) -> NonNegativeInt:
        try:
            assert (
                values["chr"] != v and values["start"] != v
            ), "Bed columns must be different"
        except KeyError:
            pass
        return v

    def assert_different(self, x: int) -> None:
        assert (
            self.chr != x and self.start != x and self.end != x
        ), "Column must be different index"

    @property
    def columns(self) -> bed.BedColumns:
        return (self.chr, self.start, self.end)


class BedFileParams(BaseModel):
    """Parameters decribing how to parse a bed-like file.

    Members:
    bed_cols - the columns for the bed coordinates
    skip_lines - how many input lines to skip
    sep - column separator regexp (for "beds" with spaces instead of tabs)
    one_indexed - starting index is considered 0 if False, and 1 if True
    """

    bed_cols: BedColumns = BedColumns()
    skip_lines: NonNegativeInt = 0
    sep: str = "\t"
    one_indexed: bool = False


class BedFile(GenericModel, Generic[Q]):
    """Import specs for a bed-like file

    Members:
    bed - source for (means to obtain) the bed file
    params - parameters to apply to bed file when downloading/formatting
    """

    bed: Q
    params: BedFileParams = BedFileParams()

    def read(self, path: Path) -> pd.DataFrame:
        """Read bed file with params from 'path'."""
        return self._read(path, [])

    def _read(
        self,
        path: Path,
        more: list[int] = [],
        comment: str | None = None,
    ) -> pd.DataFrame:
        "Read bed file with params from 'path', optionally with 'more' columns."
        p = self.params
        return bed.read_bed(
            path,
            p.bed_cols.columns,
            p.skip_lines,
            p.sep,
            p.one_indexed,
            more,
            comment,
        )


# TODO clean this up with real polymorphism when mypy catches up with Haskell
# 98, see https://github.com/python/typing/issues/548

RefSrc = RefFileSrc | RefHttpSrc

HapRefFile = HapSrc[RefSrc]
Dip1RefFile = Dip1Src[RefSrc]
Dip2RefFile = Dip2Src[RefSrc]

RefSrcT = TypeVar("RefSrcT", HapRefFile, Dip1RefFile, Dip2RefFile)

# bed-like files may be remote, local, or specified manually in the config
BedSrc = BedLocalSrc | BedHttpSrc

HapBedSrc = HapSrc[BedSrc]
DipBedSrc = Dip1Src[BedSrc] | Dip2Src[BedSrc]
Dip1BedSrc = Dip1Src[BedSrc]
Dip2BedSrc = Dip2Src[BedSrc]

BedSrcT = TypeVar("BedSrcT", HapBedSrc, Dip1BedSrc | Dip2BedSrc)

HapBedFile = BedFile[HapBedSrc]
Dip1BedFile = BedFile[Dip1BedSrc]
Dip2BedFile = BedFile[Dip2BedSrc]

HapBedFileOrCoords = HapBedFile | HapBedCoords
Dip1BedFileOrCoords = Dip1BedFile | Dip1BedCoords
Dip2BedFileOrCoords = Dip2BedFile | Dip2BedCoords

DipBedCoords = Dip1BedCoords | Dip2BedCoords

BedCoordsT = TypeVar("BedCoordsT", HapBedCoords, Dip1BedCoords | Dip2BedCoords)

# vcf files may only be remote or local, and unlike bed files, there is no
# option to use a dip1 bed file for a dip2 reference and vice versa
HapVcfSrc = HapSrc[BedSrc]
Dip1VcfSrc = Dip1Src[BedSrc]
Dip2VcfSrc = Dip2Src[BedSrc]

VcfSrcT = TypeVar("VcfSrcT", HapVcfSrc, Dip1VcfSrc, Dip2VcfSrc)
VcfSrc = HapVcfSrc | Dip1VcfSrc | Dip2VcfSrc


class BuildCompare1(BaseModel):
    """Configuration for comparing generated strats to previous versions."""

    other: CompareKey
    path_mapper: dict[Path, Path] = {}
    replacements: list[tuple[str, str]] = []
    ignore_other: list[str] = []
    ignore_generated: list[str] = []


BuildCompare2 = Diploid[BuildCompare1]

BuildCompareT = TypeVar("BuildCompareT", BuildCompare1, BuildCompare2)


class VCFFile(GenericModel, Generic[VcfSrcT]):
    """Inport specs for a vcf file."""

    vcf: VcfSrcT


class RMSKFile(BedFile[Q], Generic[Q]):
    """Input file for repeat masker stratification."""

    bed: Q  # type narrowing won't work without this redfinition
    class_col: NonNegativeInt

    @validator("class_col")
    def end_different(
        cls,
        v: NonNegativeInt,
        values: dict[str, BedColumns],
    ) -> NonNegativeInt:
        try:
            values["bed_cols"].assert_different(v)
        except KeyError:
            pass
        return v

    def read(self, path: Path) -> pd.DataFrame:
        """Read a bed file at 'path' on disk and return dataframe"""
        return super()._read(path, [self.class_col])


class SatFile(BedFile[Q], Generic[Q]):
    """Configuration for a satellites file."""

    bed: Q
    sat_col: NonNegativeInt

    def read(self, path: Path) -> pd.DataFrame:
        return super()._read(path, [self.sat_col])


class LowComplexity(GenericModel, Generic[Q, BedCoordsT]):
    """Configuration for low complexity stratification."""

    rmsk: RMSKFile[Q] | BedCoordsT | None = None
    simreps: BedFile[Q] | BedCoordsT | None = None
    satellites: SatFile[Q] | BedCoordsT | None = None


class XYFile(HapBedFile):
    """Bed file input for XY features."""

    level_col: NonNegativeInt

    @validator("level_col")
    def level_different(
        cls,
        v: NonNegativeInt,
        values: dict[str, BedColumns],
    ) -> NonNegativeInt:
        try:
            values["bed_cols"].assert_different(v)
        except KeyError:
            pass
        return v

    def read(self, path: Path) -> pd.DataFrame:
        return super()._read(path, [self.level_col])


# TODO what if the reference is XX?
# TODO what if I want to specify this file manually?
class XYFeatures(BaseModel):
    """Configuration for XY features stratifications."""

    x_bed: XYFile
    y_bed: XYFile
    ampliconic: str | None = "Ampliconic"
    xtr: str | None = "XTR"


class XYPar(BaseModel):
    """Regions for the PARs on the X/Y chromosomes."""

    start: tuple[NonNegativeInt, NonNegativeInt]
    end: tuple[NonNegativeInt, NonNegativeInt]
    comment: str = (
        "The PAR coordinates were specified manually in the pipeline configuration."
    )

    @validator("start", "end")
    def positive_region(cls, v: tuple[int, int]) -> tuple[int, int]:
        assert v[1] > v[0], "End must be greater than start"
        return v

    def fmt(self, i: ChrIndex, pattern: HapChrPattern) -> str:
        # TODO this smells like something I'll be doing alot
        c = pattern.to_chr_name(i)
        return "\n".join(
            [
                f"{c}\t{self.start[0]}\t{self.start[1]}",
                f"{c}\t{self.end[0]}\t{self.end[1]}",
            ]
        )


class XY(BaseModel):
    """Configuration for the XY stratification."""

    features: XYFeatures | None = None
    x_par: XYPar | None = None
    y_par: XYPar | None = None

    def fmt_x_par(self, pattern: HapChrPattern) -> str | None:
        return fmap_maybe(lambda x: x.fmt(ChrIndex.CHRX, pattern), self.x_par)

    def fmt_x_par_unsafe(self, pattern: HapChrPattern) -> str:
        s = self.fmt_x_par(pattern)
        if s is None:
            raise DesignError("X PAR does not exist")
        return s

    def fmt_y_par(self, pattern: HapChrPattern) -> str | None:
        return fmap_maybe(lambda x: x.fmt(ChrIndex.CHRY, pattern), self.y_par)

    def fmt_y_par_unsafe(self, pattern: HapChrPattern) -> str:
        s = self.fmt_y_par(pattern)
        if s is None:
            raise DesignError("Y PAR does not exist")
        return s


class Mappability(BaseModel):
    """Configuration for Mappability stratification.

    members:
    - unplaced_chr_patterns: a list of regexps that will be used to identify
      non-primary chromosomes in the reference to be included in mappability
      evaluation.
    """

    unplaced_chr_patterns: list[str]


class SegDups(GenericModel, Generic[Q, BedCoordsT]):
    """Configuration for Segdup stratifications."""

    superdups: BedFile[Q] | BedCoordsT | None = None


class LowMapParams(BaseModel):
    """Parameters for a single mappability bed file."""

    length: NonNegativeInt
    mismatches: NonNegativeInt
    indels: NonNegativeInt


class GCParams(BaseModel):
    """The params by which to generate GC stratifications.

    Members:
    low: the lower boundaries to use; for instance, a list like [X, Y, Z] will
         correspond to bed files with GC content <X, X-Y, and Y-Z
    high: reverse of 'low'

    The second part of the bound corresponds to whether the boundary should be
    used to create combined boundary (True means yes). The number of True's in
    each list must equal, and will be matched in inverse order in each list.
    Thus something like low = [(X1, True) (X2, True)] and
    high = [(Y1, True), (Y2, True)] will correspond to two bed files with GC
    content <X1 & >Y2 and <X2 & >Y1.
    """

    low: list[GCBound] = [
        (15, False),
        (20, False),
        (25, True),
        (30, True),
    ]
    high: list[GCBound] = [
        (55, True),
        (60, False),
        (65, True),
        (70, False),
        (75, False),
        (80, False),
        (85, False),
    ]

    @validator("low", "high")
    def non_empty_range(cls, rng: list[GCBound]) -> list[GCBound]:
        try:
            assert len(rng) > 0, "GC low and high must be non-empty lists"
        except KeyError:
            pass
        return rng

    @validator("high")
    def has_balanced_nonempty_ranges(
        cls,
        high: list[GCBound],
        values: dict[str, Any],
    ) -> list[GCBound]:
        try:
            low = cast(list[GCBound], values["low"])
            low_len = len([x for x in low if x[1]])
            high_len = len([x for x in high if x[1]])
            assert (
                (low_len == high_len) and low_len > 0 and high_len > 0
            ), "GC low/high must have at least one and the same number of range boundaries"
        except KeyError:
            pass
        return high

    @property
    def low_sorted(self) -> list[GCBound]:
        return sorted(self.low, key=lambda x: x[0])

    @property
    def high_sorted(self) -> list[GCBound]:
        return sorted(self.high, key=lambda x: x[0])

    @property
    def low_fractions(self) -> list[int]:
        return [x[0] for x in self.low_sorted]

    @property
    def high_fractions(self) -> list[int]:
        return [x[0] for x in self.high_sorted]

    # NOTE these assume that the low/high lists are non-empty
    @property
    def low_bounds(self) -> tuple[int, list[int]]:
        bounds = self.low_fractions
        return (bounds[0], bounds[1:])

    @property
    def high_bounds(self) -> tuple[int, list[int]]:
        bounds = self.high_fractions
        return (bounds[-1], bounds[:-1])


class Include(BaseModel):
    """Flags to control which stratification levels are included."""

    low_complexity: bool = True
    xy: bool = True
    segdups: bool = True
    union: bool = True
    telomeres: bool = True
    cds: bool = True
    vdj: bool = True
    mhc: bool = False  # default to false since this isn't implemented yet
    kir: bool = False  # ditto
    mappability: set[LowMapParams] = {
        LowMapParams(length=250, mismatches=0, indels=0),
        LowMapParams(length=100, mismatches=2, indels=1),
    }
    gc: GCParams | None = GCParams()
    # NOTE: This is crude but it should a) work, b) provide a decent user xp
    # and c) typecheck nicely without requiring me to use a zillionth typevar
    # in all my signatures.
    #
    # This parameter is only used for diploid configurations. While it will be
    # present for all configurations, it will only be read when needed.
    # Defaulting to a finite set means that the user never needs to specify it
    # if they want this for diploid (which they probably do) and they don't need
    # to care in haploid cases. The only issue would be if the user specified
    # this in the haploid case; it technically should be a validation error
    # since it makes no sense in the case of haploid, but here it is setup to
    # not hurt anything.
    hets: set[int] = {100, 250, 500, 1000, 5000, 10000, 25000, 50000, 100000}


class OtherBedFile(GenericModel, Generic[BedSrcT, BedCoordsT]):
    """A bed file that is imported with minimal processing and included as-is
    in a given stratification package. Useful for one-off bed files made in a
    somewhat hacky (but documented) manner that I don't feel like enshrining
    via code here.

    The only processing done on these files is gap removal if desired. They are
    still checked for correctness.

    Attributes:
    - remove_gaps: If True, subtract that gaps bed if present.
    - description: Single sentence documenting the contents of the file. Will
      show up in the README.
    - data: the bed file source or yaml coordinates
    """

    description: str
    remove_gaps: bool = False
    data: BedFile[BedSrcT] | BedCoordsT


class Bench(GenericModel, Generic[BedSrcT, VcfSrcT]):
    """Configuration for benchmark to use when validating stratifications.

    Note: the two vcf files need to have haploid/diploid layouts that correspond
    to the target reference (ie if the reference is dip1, these two must also be
    dip1). This is a limitation of happy/vcfeval, which won't know what to do
    if we given them two files with different haplotypes.

    This restriction doesn't apply to the bed file since we can split/combine
    these are necessary to match the reference.
    """

    # TODO I could probably split/compine the VCFs as well...but that sounds
    # like too much work
    bench_vcf: VCFFile[VcfSrcT]
    query_vcf: VCFFile[VcfSrcT]
    bench_bed: BedFile[BedSrcT]


class Malloc(BaseModel):
    """Manual memory allocations for rules in Mb

    Note, this can obviously be done with snakemake itself, but I want to set
    these per-build if needed.

    These are only needed for "high" memory steps, which really only means
    things that involve sorting, hap.py, minimap, or GEM.
    """

    # mappability steps
    gemIndex: int = 16000  # GEM
    gemMappability: int = 12000  # GEM
    gemToWig: int = 8000  # GEM
    mergeNonunique: int = 4000  # sort

    # normalization steps (all of which involve a sort)
    normalizeRmsk: int = 4000
    normalizeSimreps: int = 4000
    normalizeCensat: int = 4000
    normalizeSuperdups: int = 4000
    normalizeCds: int = 4000
    normalizeOther: int = 4000

    # diploid steps
    crossAlignBreaks: int = 16000  # minimap2
    crossAlignVariants: int = 16000  # minimap2
    filterSortVariantCrossAlignment: int = 16000  # samtools sort

    # happy
    runHappy: int = 48000


OtherDict = dict[OtherLevelKey, dict[OtherStratKey, OtherBedFile[BedSrcT, BedCoordsT]]]


class Build(GenericModel, Generic[BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT]):
    chr_filter: set[ChrIndex]
    comparison: BuildCompareT | None = None
    bench: Bench[BedSrcT, VcfSrcT] | None = None
    other_strats: dict[
        OtherLevelKey, dict[OtherStratKey, OtherBedFile[BedSrcT, BedCoordsT]]
    ] = {}
    # TODO if I really want I could validate this such that the user would be
    # politely alerted in case they specify any diploid params for a haploid
    # config.
    include: Include = Include()
    malloc: Malloc | None = None
    bigbed: bool = False

    @validator("other_strats")
    def valid_other(
        cls,
        os: dict[OtherLevelKey, dict[OtherStratKey, OtherBedFile[BedSrcT, BedCoordsT]]],
    ) -> dict[OtherLevelKey, dict[OtherStratKey, OtherBedFile[BedSrcT, BedCoordsT]]]:
        for k, v in os.items():
            if k == CoreLevel.OTHER_DIFFICULT:
                bs = v.keys() & BUILTIN_OTHER
                assert len(bs) == 0, f"unallowed built-in: {', '.join(bs)}"
        return os


class OtherLevelDescription(BaseModel):
    key: OtherLevelKey
    desc: str


class CDSParams(BaseModel):
    # Defaults for a for a "normal" gff file with Refseq and CDS for source and
    # type columns respectively
    source_match: tuple[str, int] | None = (".*RefSeq", 1)
    type_match: tuple[str, int] | None = ("CDS", 2)
    attr_col: int = 8


class CDS(BedFile[Q], Generic[Q]):
    """Configuration for CDS stratifications.

    Note that the "bed" parameter actually expects a gff-ish formatted
    file by default as its "bed" file.
    """

    bed: Q
    cds_params: CDSParams = CDSParams()
    params: BedFileParams = BedFileParams(
        bed_cols=BedColumns(chr=0, start=3, end=4),
        one_indexed=True,
    )

    def read(self, path: Path) -> pd.DataFrame:
        """Read a bed file at 'path' on disk and return dataframe"""
        mempty: list[int] = []
        fps = self.cds_params
        r = fmap_maybe_def(mempty, lambda x: [x[1]], fps.source_match)
        c = fmap_maybe_def(mempty, lambda x: [x[1]], fps.type_match)
        # comment needed here since GFF files have ### at the end (womp)
        return super()._read(path, [fps.attr_col] + r + c, "#")


# TODO move these to the Functional directory given that the non-cds stuff
# describes "genes" anyways?
class Functional(GenericModel, Generic[Q, BedCoordsT]):
    """Configuration for all Functional-ish bed files.

    If not given, and the build has vdj/mhc/kir regions, attempt to get these
    regions from the bed files in 'cds' when parsing the functional genes.

    The 'mhc/kir/vdj' slots below are to override the above method in case it
    doesn't work (these regions are strange after all).
    """

    cds: CDS[Q] | BedCoordsT | None = None
    mhc: BedFile[Q] | BedCoordsT | None = None
    kir: BedFile[Q] | BedCoordsT | None = None
    vdj: BedFile[Q] | BedCoordsT | None = None


class StratInputs(GenericModel, Generic[BedSrcT, BedCoordsT]):
    gap: BedFile[BedSrcT] | BedCoordsT | None
    low_complexity: LowComplexity[BedSrcT, BedCoordsT] = LowComplexity()
    xy: XY = XY()
    mappability: Mappability | None
    segdups: SegDups[BedSrcT, BedCoordsT] = SegDups()
    functional: Functional[BedSrcT, BedCoordsT] = Functional()

    @property
    def xy_features_unsafe(self) -> XYFeatures:
        f = self.xy.features
        if f is None:
            raise DesignError("XY features does not exist")
        return f

    def xy_feature_bed_unsafe(self, i: ChrIndex) -> XYFile:
        f = self.xy_features_unsafe
        return choose_xy_unsafe(i, f.x_bed, f.y_bed)


HapStratInputs = StratInputs[HapBedSrc, HapBedCoords]
DipStratInputs = StratInputs[DipBedSrc, DipBedCoords]

StratInputT = TypeVar("StratInputT", HapStratInputs, DipStratInputs)


@dataclass(frozen=True)
class RefData_(Generic[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT]):
    """A helper class corresponding a given reference and its builds.

    This is primarily meant to provide a glue layer b/t the configuration
    structure and the functions that consume data from it. Because the logic of
    looking up a refkey and determining if it is hap/dip1/dip2 is tedious and
    annoying, this type will represent the results of such a lookup and provide
    an interface for downstream processing. It also is typed generically such
    that mypy can make inferences regarding its membership in hap/dip1/dip2.

    """

    refkey: RefKey
    ref: RefSrcT
    strat_inputs: StratInputs[BedSrcT, BedCoordsT]
    builds: dict[BuildKey, Build[BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT]]

    @property
    def ref_refkeys(self) -> RefKeyFull1or2:
        "The list of full refkeys for the reference (either one or two)"
        return to_refkeys(self.ref.src, self.refkey)

    @property
    def ref_str_refkeys(self) -> RefKeyFullS1or2:
        "Like 'ref_refkeys' but returns strings."
        return to_str_refkeys(self.ref.src, self.refkey)

    @property
    def mappability_patterns(self) -> list[str]:
        """List of mappability patterns for use in filtering extra contigs.

        Return an empty list if mappability is not given."""
        return fmap_maybe_def(
            [],
            lambda m: m.unplaced_chr_patterns,
            self.strat_inputs.mappability,
        )

    def to_build_data_unsafe(
        self,
        bk: BuildKey,
    ) -> "BuildData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT]":
        "Lookup a given build with a build key (and throw DesignError on fail)"
        bd = self.to_build_data(bk)
        if bd is None:
            raise DesignError(f"Could not create build data from key '{bk}'")
        return bd

    def to_build_data(
        self,
        bk: BuildKey,
    ) -> "BuildData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT] | None":
        "Lookup a given build with a build key"
        try:
            return BuildData_(self, bk, self.builds[bk])
        except KeyError:
            return None

    def get_refkeys(self, f: RefDataToSrc) -> RefKeyFullS1or2 | None:
        """
        Get the list of refkeys (either one or two) given a function
        that retrieves an input file
        """
        return fmap_maybe(
            lambda s: to_str_refkeys(s, self.refkey),
            f(self),
        )

    def get_si_refkeys(self, f: StratInputToSrc) -> RefKeyFullS1or2 | None:
        """
        Like 'get_refkeys' but 'f' takes strat_inputs and not a ref object.
        """
        return fmap_maybe(
            lambda s: to_str_refkeys(s, self.refkey),
            f(self.strat_inputs),
        )

    @property
    def has_low_complexity_rmsk(self) -> bool:
        """Return True if this reference has repeat masker specified."""
        return self.strat_inputs.low_complexity.rmsk is not None

    @property
    def has_low_complexity_simreps(self) -> bool:
        """Return True if this reference has simple repeats specified."""
        return self.strat_inputs.low_complexity.simreps is not None

    @property
    def has_low_complexity_censat(self) -> bool:
        """Return True if this reference has satellites specified."""
        return self.strat_inputs.low_complexity.satellites is not None


HapRefData = RefData_[HapRefFile, HapBedSrc, HapVcfSrc, HapBedCoords, BuildCompare1]
Dip1RefData = RefData_[Dip1RefFile, DipBedSrc, Dip1VcfSrc, DipBedCoords, BuildCompare1]
Dip2RefData = RefData_[Dip2RefFile, DipBedSrc, Dip2VcfSrc, DipBedCoords, BuildCompare2]

AnyRefData = HapRefData | Dip1RefData | Dip2RefData


@dataclass(frozen=True)
class BuildData_(Generic[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT]):
    """A helper class corresponding a given build.

    This follows a similar motivation as 'RefData_' above.
    """

    refdata: RefData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT]
    buildkey: BuildKey
    build: Build[BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT]

    @property
    def build_chrs(self) -> BuildChrs:
        """Return a set of all desired chromosomes for this build.

        NOTE: this is an unfiltered set, meaning that if the pattern for the
        reference excludes a set of chromosomes on a haplotype (ie the X on
        the paternal) this set will NOT reflect that exclusion.
        """
        cs = self.build.chr_filter
        return BuildChrs(set([x for x in ChrIndex]) if len(cs) == 0 else cs)

    @property
    def chr_indices(self) -> set[ChrIndex]:
        cs = self.build.chr_filter
        return set([x for x in ChrIndex]) if len(cs) == 0 else cs

    @property
    def want_bb(self) -> bool:
        return self.build.bigbed

    @property
    def want_diploid(self) -> bool:
        return len(self.build.include.hets) > 0

    @property
    def want_low_complexity(self) -> bool:
        return self.build.include.low_complexity

    @property
    def want_gc(self) -> bool:
        return self.build.include.gc is not None

    @property
    def want_telomeres(self) -> bool:
        return self.build.include.telomeres

    @property
    def want_segdups(self) -> bool:
        return self.build.include.segdups

    @property
    def want_union(self) -> bool:
        return self.build.include.union

    @property
    def have_gaps(self) -> bool:
        return self.refdata.strat_inputs.gap is not None

    @property
    def have_benchmark(self) -> bool:
        return self.build.bench is not None

    @property
    def want_hets(self) -> bool:
        r = self.refdata.ref
        if isinstance(r, HapSrc):
            return False
        elif isinstance(r, Dip1Src) or isinstance(r, Dip2Src):
            return len(self.build.include.hets) > 0
        else:
            assert_never(r)

    @property
    def mappability_params(
        self,
    ) -> tuple[list[int], list[int], list[int]]:
        ms = self.build.include.mappability
        return unzip3([(m.length, m.mismatches, m.indels) for m in ms])

    @property
    def want_mappability(self) -> bool:
        return len(self.build.include.mappability) > 0

    # this is the only xy-ish flag that's here since this won't depend on the
    # haplotype (just remove X and Y from whatever filter we set)
    #
    # TODO technically this isn't true because the autosomes could in theory be
    # excluded for each haplotype separately, but this should almost never happen
    # in real life (see vdj below)
    @property
    def want_xy_auto(self) -> bool:
        return len(self.build_chrs - set([ChrIndex.CHRX, ChrIndex.CHRY])) > 0

    # For each of these we could check if the X or Y chromosome(s) is/are
    # present in the chr filter. However, this would require
    # hap/dip1/dip2-specific subclasses which I don't feel like making. Instead
    # it is much easier to query the refkey and buildkey for the desired xy
    # chromosomes, inherit rules based on this, and then filter those rules
    # based on the following functions.
    @property
    def have_y_PAR(self) -> bool:
        return self.refdata.strat_inputs.xy.y_par is not None

    @property
    def want_cds(self) -> bool:
        return self.build.include.cds

    @property
    def want_mhc(self) -> bool:
        return self.build.include.mhc and MHC_CHR in self.build_chrs

    @property
    def want_kir(self) -> bool:
        return self.build.include.kir and KIR_CHR in self.build_chrs

    @property
    def want_vdj(self) -> bool:
        return self.build.include.vdj and len(VDJ_CHRS & self.build_chrs) > 0


class Stratification(
    GenericModel, Generic[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT]
):
    """Configuration for stratifications for a given reference."""

    ref: RefSrcT
    strat_inputs: StratInputs[BedSrcT, BedCoordsT]
    builds: dict[BuildKey, Build[BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT]]


HapBuildData = BuildData_[
    HapRefFile,
    HapBedSrc,
    HapVcfSrc,
    HapBedCoords,
    BuildCompare1,
]
Dip1BuildData = BuildData_[
    Dip1RefFile,
    DipBedSrc,
    Dip1VcfSrc,
    DipBedCoords,
    BuildCompare1,
]
Dip2BuildData = BuildData_[
    Dip2RefFile,
    DipBedSrc,
    Dip2VcfSrc,
    DipBedCoords,
    BuildCompare2,
]

AnyBuildData = HapBuildData | Dip1BuildData | Dip2BuildData

HapStrat = Stratification[
    HapRefFile,
    HapBedSrc,
    HapVcfSrc,
    HapBedCoords,
    BuildCompare1,
]
Dip1Strat = Stratification[
    Dip1RefFile,
    DipBedSrc,
    Dip1VcfSrc,
    DipBedCoords,
    BuildCompare1,
]
Dip2Strat = Stratification[
    Dip2RefFile,
    DipBedSrc,
    Dip2VcfSrc,
    DipBedCoords,
    BuildCompare2,
]

AnyStrat = HapStrat | Dip1Strat | Dip2Strat


class Documentation(BaseModel):
    pipeline_repo: HttpUrl = "https://github.com/usnistgov/giab-stratifications-pipeline"  # type: ignore
    config_repo: HttpUrl = (
        "https://github.com/usnistgov/giab-stratifications"  # type: ignore
    )


GENOME_SPECIFIC_DESC = (
    "difficult regions due to potentially difficult variation in a NIST/GIAB "
    "sample, including 1) regions containing putative compound heterozygous "
    "variants 2) small regions containing multiple phased variants, 3) regions "
    "with potential structural or copy number variation"
)

FUNCTIONAL_TECH_DESC = (
    "functional, or potentially functional, regions that are also likely to be "
    "technically difficult to sequences"
)


class GiabStrats(BaseModel):
    """Top level stratification object."""

    other_levels: list[OtherLevelDescription] = [
        OtherLevelDescription(
            key=OtherLevelKey("Ancestry"),
            desc="regions with inferred patterns of local ancestry",
        ),
        OtherLevelDescription(
            key=OtherLevelKey("FunctionalTechnicallyDifficult"),
            desc=FUNCTIONAL_TECH_DESC,
        ),
        OtherLevelDescription(
            key=OtherLevelKey("GenomeSpecific"),
            desc=GENOME_SPECIFIC_DESC,
        ),
    ]
    paths: Paths = Paths()
    tools: Tools = Tools()
    comparison_strats: dict[CompareKey, HttpUrl] = {}
    haploid_stratifications: dict[RefKey, HapStrat] = {}
    diploid1_stratifications: dict[RefKey, Dip1Strat] = {}
    diploid2_stratifications: dict[RefKey, Dip2Strat] = {}
    benchmark_subsets: list[str] = [
        "AllAutosomes",
        "AllTandemRepeats",
        "AllHomopolymers_ge7bp_imperfectge11bp_slop5",
        "SimpleRepeat_diTR_10to49_slop5",
        "SimpleRepeat_homopolymer_7to11_slop5",
        "SimpleRepeat_homopolymer_ge21_slop5",
        "SimpleRepeat_imperfecthomopolge11_slop5",
        "SimpleRepeat_imperfecthomopolge21_slop5",
        "SimpleRepeat_homopolymer_7to11_AT_slop5",
        "SimpleRepeat_homopolymer_ge21_AT_slop5",
        "SimpleRepeat_imperfecthomopolge11_AT_slop5",
        "SimpleRepeat_imperfecthomopolge21_AT_slop5",
        "SimpleRepeat_homopolymer_7to11_GC_slop5",
        "SimpleRepeat_homopolymer_ge21_GC_slop5",
        "SimpleRepeat_imperfecthomopolge11_GC_slop5",
        "SimpleRepeat_imperfecthomopolge21_GC_slop5",
        "alldifficultregions",
        "alllowmapandsegdupregions",
        "chrX_PAR",
        "chrX_XTR",
        "chrY_XTR",
        "notinalldifficultregions",
        "notinAllHomopolymers_ge7bp_imperfectge11bp_slop5",
        "notinAllTandemRepeatsandHomopolymers_slop5",
        "segdups",
    ]
    malloc: Malloc = Malloc()
    docs: Documentation = Documentation()

    @validator(
        "haploid_stratifications",
        "diploid1_stratifications",
        "diploid2_stratifications",
        each_item=True,
    )
    def builds_have_valid_existing(
        cls,
        v: AnyStrat,
        values: dict[str, Any],
    ) -> AnyStrat:
        try:
            levels: list[OtherLevelDescription] = values["other_levels"]
            assert OTHERDIFF_KEY not in [
                x.key for x in levels
            ], f"{OTHERDIFF_KEY} cannot be in other_levels"

            _levels = [OTHERDIFF_KEY, *[x.key for x in levels]]

            bad = [
                f"level='{lk}'; build='{bk}'"
                for bk, b in v.builds.items()
                for lk in b.other_strats
                if lk not in _levels
            ]
            if len(bad) > 0:
                assert (
                    False
                ), f"builds referencing invalid strat categories: {', '.join(bad)}"
        except KeyError:
            pass
        return v

    @validator(
        "haploid_stratifications",
        "diploid1_stratifications",
        "diploid2_stratifications",
        each_item=True,
    )
    def builds_have_valid_old_version(
        cls,
        v: AnyStrat,
        values: dict[str, Any],
    ) -> AnyStrat:
        try:
            prev: dict[CompareKey, HttpUrl] = values["comparison_strats"]
            bad = [
                f"version='{ck}'; build='{bk}'"
                for bk, b in v.builds.items()
                if (c := b.comparison) is not None
                for ck in (
                    [c.other]
                    if isinstance(c, BuildCompare1)
                    else [x.other for x in c.double.as_list]
                )
                if ck not in prev
            ]
            assert (
                len(bad) == 0
            ), f"builds referencing invalid previous version keys: {', '.join(bad)}"
        except KeyError:
            pass
        return v

    @validator("diploid2_stratifications")
    def no_overlapping_refkeys(
        cls,
        v: dict[RefKey, Dip2Strat],
        values: dict[str, Any],
    ) -> dict[RefKey, Dip2Strat]:
        try:
            hap: list[RefKey] = list(values["haploid_stratifications"])
            dip1: list[RefKey] = list(values["diploid2_stratifications"])
            dip2 = list(v)
            ds = list(duplicates_everseen(hap + dip1 + dip2))
            assert len(ds) == 0, f"duplicate refkeys: {', '.join(ds)}"
        except KeyError:
            pass
        return v

    # hack to make rmd scripts work with this (note this will totally kill
    # the config as it passes into an rmd script)
    def items(self) -> Any:
        return {}.items()

    # attributes

    @property
    def other_level_keys(self) -> list[OtherLevelKey]:
        return [x.key for x in self.other_levels]

    @property
    def all_other_levels(self) -> list[OtherLevelKey]:
        return [*self.other_level_keys, OTHERDIFF_KEY]

    # file paths

    @property
    def resources_dir(self) -> Path:
        return self.paths.resources

    @property
    def _tools_base_dir(self) -> Path:
        return self.resources_dir / "tools"

    @property
    def tools_src_dir(self) -> Path:
        return self._tools_base_dir / "src"

    @property
    def tools_make_dir(self) -> Path:
        return self._tools_base_dir / "make"

    @property
    def tools_bin_dir(self) -> Path:
        return self._tools_base_dir / "bin"

    @property
    def ref_src_dir(self) -> Path:
        return self.resources_dir / "{ref_src_key}"

    @property
    def results_dir(self) -> Path:
        return self.paths.results

    @property
    def final_root_dir(self) -> Path:
        return self.results_dir / "final"

    @property
    def final_build_dir(self) -> Path:
        return self.final_root_dir / "{ref_final_key}@{build_key}"

    @property
    def intermediate_root_dir(self) -> Path:
        return self.results_dir / "intermediates"

    @property
    def intermediate_build_hapless_dir(self) -> Path:
        return self.intermediate_root_dir / "{ref_key}@{build_key}"

    @property
    def intermediate_build_dir(self) -> Path:
        return self.intermediate_root_dir / "{ref_final_key}@{build_key}"

    @property
    def bench_root_dir(self) -> Path:
        return self.results_dir / "bench"

    @property
    def log_tools_dir(self) -> Path:
        return self.resources_dir / "log" / "tools"

    @property
    def log_src_dir(self) -> Path:
        return self.resources_dir / "log" / "{ref_src_key}"

    @property
    def log_results_dir(self) -> Path:
        return self.results_dir / "log" / "{ref_final_key}"

    @property
    def log_build_dir(self) -> Path:
        return self.log_results_dir / "builds" / "{ref_final_key}@{build_key}"

    @property
    def log_build_hapless_dir(self) -> Path:
        return self.log_results_dir / "builds" / "{ref_key}@{build_key}"

    @property
    def bench_build_dir(self) -> Path:
        return self.bench_root_dir / "{ref_final_key}@{build_key}"

    @property
    def bench_build_hapless_dir(self) -> Path:
        return self.bench_root_dir / "{ref_key}@{build_key}"

    def build_final_strat_path(self, level: str, name: str) -> Path:
        return self.final_build_dir / level / f"{{ref_final_key}}_{name}.bed.gz"

    def build_final_readme_path(self, level: str) -> Path:
        return self.final_build_dir / level / f"{{ref_final_key}}_{level}_README.md"

    def build_strat_path(self, level: CoreLevel, name: str) -> Path:
        return self.build_final_strat_path(level.value, name)

    def build_readme_path(self, level: CoreLevel) -> Path:
        return self.build_final_readme_path(level.value)

    @property
    def ref_dirs(self) -> RefDirs:
        return RefDirs(
            src=RefSrcDirs(
                reference=DataLogDirs(
                    data=self.ref_src_dir / "reference",
                    log=self.log_src_dir / "reference",
                ),
                benchmark=DataLogDirs(
                    data=self.ref_src_dir / "benchmark" / "{build_key}",
                    log=self.log_src_dir / "benchmark" / "{build_key}",
                ),
            ),
            inter=RefInterDirs(
                prebuild=DataLogBenchDirs(
                    data=self.intermediate_root_dir / "{ref_final_key}",
                    log=self.log_results_dir / "{ref_final_key}",
                    bench=self.bench_root_dir / "{ref_final_key}",
                ),
                filtersort=FilterSortDirs(
                    data=self.intermediate_build_hapless_dir / "ref",
                    subbed=prepare_output_path(
                        self.intermediate_build_hapless_dir / "ref"
                    ),
                    log=self.log_build_hapless_dir / "ref",
                    bench=self.bench_build_hapless_dir / "ref",
                ),
                build=DataLogBenchDirs(
                    data=self.intermediate_root_dir / "{ref_final_key}@{build_key}",
                    log=self.log_results_dir / "{ref_final_key}@{build_key}",
                    bench=self.bench_root_dir / "{ref_final_key}@{build_key}",
                ),
            ),
        )

    def to_bed_dirs(self, level: CoreLevel) -> BedDirs:
        v = level.value
        return BedDirs(
            src=DataLogDirs(
                self.ref_src_dir / v,
                self.log_src_dir / v,
            ),
            inter=BedInterDirs(
                filtersort=FilterSortDirs(
                    data=self.intermediate_build_hapless_dir / v,
                    log=self.log_build_hapless_dir / v,
                    bench=self.bench_build_hapless_dir / v,
                    subbed=prepare_output_path(self.intermediate_build_hapless_dir / v),
                ),
                postsort=DataLogBenchDirs(
                    data=self.intermediate_build_dir / v,
                    log=self.log_build_dir / v,
                    bench=self.bench_build_dir / v,
                ),
            ),
            final=lambda name: self.build_strat_path(level, name),
            readme=self.build_readme_path(level),
        )

    # because smk doesn't check these for existence yet:
    # https://github.com/snakemake/snakemake/issues/1657
    def _workflow_path(self, components: list[str]) -> Path:
        p = Path("workflow", *components)
        # except that it doesn't work too well in subworkflows...
        # assert p.exists(), f"{p} does not exist"
        return p

    def env_path(self, envname: str) -> Path:
        return self._workflow_path(["envs", f"{envname}.yml"])

    def _scripts_dir(self, rest: list[str]) -> Path:
        return self._workflow_path(["scripts", *rest])

    def python_script(self, basename: str) -> Path:
        return self._scripts_dir(["python", basename])

    def rmd_script(self, basename: str) -> Path:
        return self._scripts_dir(["rmarkdown", basename])

    # general accessors

    def buildkey_to_chrs(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
        split: bool,
        nohap: bool,
    ) -> HapChrs:
        return self.with_build_data_full_rconf(
            rk,
            bk,
            split,
            nohap,
            lambda bd: bd.refdata.ref.hap_chrs(bd.build_chrs),
            lambda bd: bd.refdata.ref.all_chrs(bd.build_chrs),
            lambda hap, bd: bd.refdata.ref.hap_chrs(bd.build_chrs, hap),
            lambda hap, bd: bd.refdata.ref.hap_chrs(bd.build_chrs, hap),
        )

    def buildkey_to_wanted_xy(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
    ) -> HapChrs:
        cis = self.buildkey_to_chrs(rk, bk, False, False)
        # NOTE order matters here since this will be used to create sorted bed
        # files
        return HapChrs({i for i in [ChrIndex.CHRX, ChrIndex.CHRY] if i in cis})

    def buildkey_to_wanted_xy_names(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
    ) -> set[ShortChrName]:
        return {c.chr_name for c in self.buildkey_to_wanted_xy(rk, bk)}

    def buildkey_to_malloc(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
        f: Callable[[Malloc], int],
    ) -> int:
        bd = self.to_build_data(strip_full_refkey(rk), bk)
        return max(
            fmap_maybe_def(f(self.malloc), lambda m: f(m), bd.build.malloc), 1000
        )

    def buildkey_to_ref_mappers(
        self, rk: RefKeyFullS, bk: BuildKey
    ) -> tuple[bed.InitMapper, bed.FinalMapper]:
        """Lookup a given build and return the init/final mappers
        corresponding to the reference chromosome names.

        This is useful for cases where the reference itself is used to
        generate a bed-like file which then needs to be sorted.
        """
        m = self.with_build_data_full(
            rk,
            bk,
            hap_noop_conversion,
            dip1_noop_conversion,
            dip2_noop_conversion,
        )
        return (m.init_mapper, m.final_mapper)

    def buildkey_to_ref_mappers_split(
        self, rk: RefKeyFullS, bk: BuildKey
    ) -> tuple[bed.InitMapper, bed.FinalMapper]:
        m = self.with_build_data_split_full(
            rk,
            bk,
            hap_noop_conversion,
            dip1_split_noop_conversion,
            dip2_noop_conversion,
        )
        return (m.init_mapper, m.final_mapper)

    def buildkey_to_ref_mappers_split_nohap(
        self, rk: RefKeyFullS, bk: BuildKey
    ) -> tuple[bed.InitMapper, bed.FinalMapper]:
        m = self.with_build_data_split_full_nohap(
            rk,
            bk,
            dip1_split_noop_conversion,
            dip2_noop_conversion,
        )
        return (m.init_mapper, m.final_mapper)

    def refsrckey_to_ref_src(self, rsk: RefKeyFullS) -> RefSrc:
        """Lookup a given reference and return its source object (haplotype
        specific)."""
        rk, hap = parse_full_refkey(rsk)
        src = self.to_ref_data(rk).ref.src
        return from_single_or_double(src, hap)

    def refkey_to_bed_refsrckeys(
        self, f: StratInputToBed, rk: RefKey
    ) -> RefKeyFullS1or2 | None:
        """Lookup a given reference and return the full refkeys for the
        bed file obtained with the given function.

        This is useful for bed file normalization rules which are all in terms
        of the bare refkey (ie not haplotype specific) but need to somehow
        get a list of inputs which are downloaded. Since there might be one or
        two inputs which may or may not have a haplotype associated with them,
        this function provides the full refkeys for obtaining said inputs.
        """
        # TODO this seems like a useful glue function (the labmda that is)
        return self.to_ref_data(rk).get_refkeys(
            lambda rd: (
                None
                if (bf := f(rd.strat_inputs)) is None
                else (bf.bed.src if isinstance(bf, BedFile) else None)
            )
        )

    def refkey_to_bed_refsrckeys_smk(
        self, f: StratInputToBed, rk: RefKey
    ) -> list[RefKeyFullS]:
        return (
            x.as_list
            if (x := self.refkey_to_bed_refsrckeys(f, rk)) is not None
            else raise_inline()
        )

    def _refkey_to_src(self, f: RefDataToSrc, rk: RefKeyFullS) -> BedSrc:
        rk_, hap = parse_full_refkey(rk)
        src = with_ref_data(
            self.to_ref_data(rk_), lambda rd: f(rd), lambda rd: f(rd), lambda rd: f(rd)
        )
        # TODO this should check for a single/double anybedsrc
        return from_single_or_double(src, hap) if src is not None else raise_inline()

    def refsrckey_to_bed_src(self, f: StratInputToBed, rk: RefKeyFullS) -> BedSrc:
        """Lookup a haplotype-specific bed file source with the given function."""
        return self._refkey_to_src(
            lambda rd: (
                b.bed.src
                if (b := f(rd.strat_inputs)) is not None and isinstance(b, BedFile)
                else raise_inline()
            ),
            rk,
        )

    def _refsrckey_to_xy_feature_src(self, rsk: RefKeyFullS, i: ChrIndex) -> BedSrc:
        return (
            self.to_ref_data(strip_full_refkey(rsk))
            .strat_inputs.xy_feature_bed_unsafe(i)
            .bed.src.elem
        )

    def refsrckey_to_x_features_src(self, rsk: RefKeyFullS) -> BedSrc:
        """Return the X features source file for a given reference."""
        return self._refsrckey_to_xy_feature_src(rsk, ChrIndex.CHRX)

    def refsrckey_to_y_features_src(self, rsk: RefKeyFullS) -> BedSrc:
        """Return the Y features source file for a given reference."""
        return self._refsrckey_to_xy_feature_src(rsk, ChrIndex.CHRY)

    def refkey_to_xy_ref_chr_pattern(
        self, rk: RefKeyFullS, i: ChrIndex
    ) -> HapChrPattern:
        """Return the XY chr pattern for a given reference and haplotype."""
        return self.with_ref_data_full(
            rk,
            lambda rd: rd.ref.chr_pattern,
            lambda rd: rd.ref.chr_pattern.to_hap_pattern(i.xy_to_hap_unsafe),
            lambda hap, rd: rd.ref.chr_pattern.choose(hap),
        )

    def buildkey_to_bed_refsrckeys(
        self, f: BuildDataToBed, rk: RefKey, bk: BuildKey
    ) -> RefKeyFullS1or2 | None:
        """Like 'refkey_to_bed_refsrckeys' but build-specific.

        Used for looking up benchmark files for each build.
        """
        # TODO this "update" function is not DRY
        # return self.refkey_to_bed_refsrckeys(lambda rd: f(rd.to_build_data(bk)), rk)
        return self.to_ref_data(rk).get_refkeys(
            lambda rd: (
                None
                if (bf := f(rd.to_build_data(bk))) is None
                else (bf.bed.src if isinstance(bf, BedFile) else None)
            )
        )

    def buildkey_to_bed_refsrckeys_smk(
        self, f: BuildDataToBed, rk: RefKey, bk: BuildKey
    ) -> list[RefKeyFullS]:
        return (
            x.as_list
            if (x := self.buildkey_to_bed_refsrckeys(f, rk, bk)) is not None
            else raise_inline()
        )

    def buildkey_to_bed_src(
        self, f: BuildDataToBed, rk: RefKeyFullS, bk: BuildKey
    ) -> BedSrc:
        """Like 'refsrckey_to_bed_src' but build-specific.

        Used for looking up benchmark sources for each build.
        """
        return self._refkey_to_src(
            lambda rd: (
                b.bed.src
                if (b := f(rd.to_build_data(bk))) is not None and isinstance(b, BedFile)
                else raise_inline()
            ),
            rk,
        )

    def buildkey_to_vcf_src(
        self, f: BuildDataToVCF, rk: RefKeyFullS, bk: BuildKey
    ) -> BedSrc:
        """Like 'buildkey_to_bed_src' but for benchmark VCF sources."""
        # TODO not DRY
        rk_, hap = parse_full_refkey(rk)
        bd = self.to_build_data(rk_, bk)
        src = with_build_data(bd, lambda bd: f(bd), lambda bd: f(bd), lambda bd: f(bd))
        if src is None:
            raise DesignError()
        return from_single_or_double(src.vcf.src, hap)

    def refkey_to_normalization_path(self, rk: RefKeyFullS, s: IO[bytes]) -> Path:
        """Return a list of paths for a given normalization checkpoint.

        This is assumed to specifically be called within a checkpoint that is
        downstream of a normalization rule. In this case, the rule spits out a
        JSON file with a list of paths that it produces, which either has one
        or two paths correspond to hap/dip1 or dip2 cases. In the latter case
        return the one path corresponding to the haplotype appended to 'rk'.

        's' is an open stream representing the output path from the checkpoint.
        """
        res = json.load(s)
        try:
            paths = [Path(p) for p in res]
        except TypeError:
            raise DesignError(f"Checkpoint does not have paths list, got {res}")

        return self.with_ref_data_full(
            rk,
            lambda _: match1_unsafe(paths, noop),
            lambda _: match1_unsafe(paths, noop),
            lambda hap, _: match2_unsafe(paths, lambda p: p.choose(hap)),
        )

    def to_ref_data(self, rk: RefKey) -> AnyRefData:
        """Lookup refdata object for a given refkey."""
        if rk in self.haploid_stratifications:
            return to_ref_data_unsafe(self.haploid_stratifications, rk)
        elif rk in self.diploid1_stratifications:
            return to_ref_data_unsafe(self.diploid1_stratifications, rk)
        elif rk in self.diploid2_stratifications:
            return to_ref_data_unsafe(self.diploid2_stratifications, rk)
        else:
            raise DesignError(f"invalid ref key: '{rk}'")

    def to_build_data(self, rk: RefKey, bk: BuildKey) -> AnyBuildData:
        """Lookup builddata object for a given refkey and build key."""

        def hap(rd: HapRefData) -> AnyBuildData:
            return rd.to_build_data_unsafe(bk)

        def dip1(rd: Dip1RefData) -> AnyBuildData:
            return rd.to_build_data_unsafe(bk)

        def dip2(rd: Dip2RefData) -> AnyBuildData:
            return rd.to_build_data_unsafe(bk)

        return with_ref_data(self.to_ref_data(rk), hap, dip1, dip2)

    def with_ref_data(
        self,
        rk: RefKey,
        hap_f: Callable[[HapRefData], X],
        dip1_f: Callable[[Dip1RefData], X],
        dip2_f: Callable[[Dip2RefData], X],
    ) -> X:
        """Lookup refdata and apply function depending on if it is hap, dip1,
        or dip2.
        """
        return with_ref_data(self.to_ref_data(rk), hap_f, dip1_f, dip2_f)

    def with_ref_data_full(
        self,
        rk: RefKeyFullS,
        hap_f: Callable[[HapRefData], X],
        dip1_f: Callable[[Dip1RefData], X],
        dip2_f: Callable[[Haplotype, Dip2RefData], X],
    ) -> X:
        """Like 'with_ref_data' but takes a full refkey and supplies the
        haplotype in the dip2 case.
        """
        rk_, hap = parse_full_refkey(rk)
        return self.with_ref_data(
            rk_,
            lambda rd: none_unsafe(hap, hap_f(rd)),
            lambda rd: none_unsafe(hap, dip1_f(rd)),
            lambda rd: not_none_unsafe(hap, lambda hap: dip2_f(hap, rd)),
        )

    def to_ref_data_full(self, rk: RefKeyFullS) -> AnyRefData:
        """Like 'to_ref_data' but takes a full refkey and does error checking."""

        def hap(rd: HapRefData) -> AnyRefData:
            return rd

        def dip1(rd: Dip1RefData) -> AnyRefData:
            return rd

        def dip2(_: Haplotype, rd: Dip2RefData) -> AnyRefData:
            return rd

        return self.with_ref_data_full(rk, hap, dip1, dip2)

    def with_ref_data_full_nohap(
        self,
        rk: RefKeyFullS,
        dip1_f: Callable[[Dip1RefData], X],
        dip2_f: Callable[[Haplotype, Dip2RefData], X],
    ) -> X:
        """Like 'with_ref_data_split_full' but forbids the hap1 case (ie
        for diploid only rules/scripts)
        """
        return self.with_ref_data_full(
            rk,
            lambda _: raise_inline(f"hap1 refkey not allowed: {rk}"),
            dip1_f,
            dip2_f,
        )

    def with_ref_data_split_full(
        self,
        rk: RefKeyFullS,
        hap_f: Callable[[HapRefData], X],
        dip1_f: Callable[[Haplotype, Dip1RefData], X],
        dip2_f: Callable[[Haplotype, Dip2RefData], X],
    ) -> X:
        """Like 'with_ref_data_full' but takes a refkey with a haplotype in
        the dip1 case (ie when "split")
        """
        rk_, hap = parse_full_refkey(rk)
        return self.with_ref_data(
            rk_,
            lambda rd: none_unsafe(hap, hap_f(rd)),
            lambda rd: not_none_unsafe(hap, lambda hap: dip1_f(hap, rd)),
            lambda rd: not_none_unsafe(hap, lambda hap: dip2_f(hap, rd)),
        )

    def with_ref_data_split_full_nohap(
        self,
        rk: RefKeyFullS,
        dip1_f: Callable[[Haplotype, Dip1RefData], X],
        dip2_f: Callable[[Haplotype, Dip2RefData], X],
    ) -> X:
        """Like 'with_ref_data_split_full' but forbids the hap1 case (ie
        for diploid only rules/scripts)
        """
        return self.with_ref_data_split_full(
            rk,
            lambda _: raise_inline(f"hap1 refkey not allowed: {rk}"),
            dip1_f,
            dip2_f,
        )

    def with_ref_data_full_rconf(
        self,
        rk: RefKeyFullS,
        split: bool,
        nohap: bool,
        hap_f: Callable[[HapRefData], X],
        dip1_f: Callable[[Dip1RefData], X],
        split_dip1_f: Callable[[Haplotype, Dip1RefData], X],
        dip2_f: Callable[[Haplotype, Dip2RefData], X],
    ) -> X:
        """Apply functions to ref data depending on if they are dip1/2/hap and
        depending on the refkey configuration (standard/split/nohap)
        """
        match (split, nohap):
            case (False, False):
                return self.with_ref_data_full(
                    rk,
                    hap_f,
                    dip1_f,
                    dip2_f,
                )
            case (False, True):
                return self.with_ref_data_full_nohap(
                    rk,
                    dip1_f,
                    dip2_f,
                )
            case (True, False):
                return self.with_ref_data_split_full(
                    rk,
                    hap_f,
                    split_dip1_f,
                    dip2_f,
                )
            case (True, True):
                return self.with_ref_data_split_full_nohap(
                    rk,
                    split_dip1_f,
                    dip2_f,
                )
            case _:
                # TODO indeed this should never happen, and mypy currently is
                # not smart enough to validate that this is the case:
                # https://github.com/python/mypy/issues/16722
                raise DesignError("This should never happen")

    def with_build_data(
        self,
        rk: RefKey,
        bk: BuildKey,
        hap_f: Callable[[HapBuildData], X],
        dip1_f: Callable[[Dip1BuildData], X],
        dip2_f: Callable[[Dip2BuildData], X],
    ) -> X:
        """Like 'with_ref_data' but for build data"""
        return with_build_data(self.to_build_data(rk, bk), hap_f, dip1_f, dip2_f)

    def with_build_data_full(
        self,
        rfk: RefKeyFullS,
        bk: BuildKey,
        hap_f: Callable[[HapBuildData], X],
        dip1_f: Callable[[Dip1BuildData], X],
        dip2_f: Callable[[Haplotype, Dip2BuildData], X],
    ) -> X:
        """Like 'with_ref_data_full' but for build data"""
        rk_, hap = parse_full_refkey(rfk)
        return self.with_build_data(
            rk_,
            bk,
            lambda bd: none_unsafe(hap, hap_f(bd)),
            lambda bd: none_unsafe(hap, dip1_f(bd)),
            lambda bd: not_none_unsafe(hap, lambda hap: dip2_f(hap, bd)),
        )

    def to_build_data_full(self, rk: RefKeyFullS, bk: BuildKey) -> AnyBuildData:
        """Like 'to_build_data' but takes a full refkey and does error checking."""

        def hap(rd: HapBuildData) -> AnyBuildData:
            return rd

        def dip1(rd: Dip1BuildData) -> AnyBuildData:
            return rd

        def dip2(_: Haplotype, rd: Dip2BuildData) -> AnyBuildData:
            return rd

        return self.with_build_data_full(rk, bk, hap, dip1, dip2)

    def with_build_data_full_nohap(
        self,
        rfk: RefKeyFullS,
        bk: BuildKey,
        dip1_f: Callable[[Dip1BuildData], X],
        dip2_f: Callable[[Haplotype, Dip2BuildData], X],
    ) -> X:
        """Like 'with_ref_data_full' but doesn't allow hap case"""
        return self.with_build_data_full(
            rfk,
            bk,
            lambda _: raise_inline(f"hap1 refkey not allowed: {rfk}"),
            dip1_f,
            dip2_f,
        )

    def with_build_data_split_full(
        self,
        rfk: RefKeyFullS,
        bk: BuildKey,
        hap_f: Callable[[HapBuildData], X],
        dip1_f: Callable[[Haplotype, Dip1BuildData], X],
        dip2_f: Callable[[Haplotype, Dip2BuildData], X],
    ) -> X:
        """Like 'with_ref_data_split_full' but for build data"""
        rk_, hap = parse_full_refkey(rfk)
        return self.with_build_data(
            rk_,
            bk,
            lambda rd: none_unsafe(hap, hap_f(rd)),
            lambda rd: not_none_unsafe(hap, lambda hap: dip1_f(hap, rd)),
            lambda rd: not_none_unsafe(hap, lambda hap: dip2_f(hap, rd)),
        )

    def with_build_data_split_full_nohap(
        self,
        rfk: RefKeyFullS,
        bk: BuildKey,
        dip1_f: Callable[[Haplotype, Dip1BuildData], X],
        dip2_f: Callable[[Haplotype, Dip2BuildData], X],
    ) -> X:
        """Like 'with_ref_data_split_full_nohap' but for build data"""
        return self.with_build_data_split_full(
            rfk,
            bk,
            lambda _: raise_inline(f"hap1 refkey not allowed: {rfk}"),
            dip1_f,
            dip2_f,
        )

    def with_build_data_full_rconf(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
        split: bool,
        nohap: bool,
        hap_f: Callable[[HapBuildData], X],
        dip1_f: Callable[[Dip1BuildData], X],
        split_dip1_f: Callable[[Haplotype, Dip1BuildData], X],
        dip2_f: Callable[[Haplotype, Dip2BuildData], X],
    ) -> X:
        """Apply functions to build data depending on if they are dip1/2/hap and
        depending on the refkey configuration (standard/split/nohap)
        """
        match (split, nohap):
            case (False, False):
                return self.with_build_data_full(
                    rk,
                    bk,
                    hap_f,
                    dip1_f,
                    dip2_f,
                )
            case (False, True):
                return self.with_build_data_full_nohap(
                    rk,
                    bk,
                    dip1_f,
                    dip2_f,
                )
            case (True, False):
                return self.with_build_data_split_full(
                    rk,
                    bk,
                    hap_f,
                    split_dip1_f,
                    dip2_f,
                )
            case (True, True):
                return self.with_build_data_split_full_nohap(
                    rk,
                    bk,
                    split_dip1_f,
                    dip2_f,
                )
            case _:
                # TODO indeed this should never happen, and mypy currently is
                # not smart enough to validate that this is the case:
                # https://github.com/python/mypy/issues/16722
                raise DesignError("This should never happen")

    def with_ref_data_and_bed(
        self,
        rk: RefKey,
        get_bed_f: RefDataToBed,
        hap_f: Callable[[HapRefData, HapBedFileOrCoords], X],
        dip_1to1_f: Callable[[Dip1RefData, Dip1BedFileOrCoords], X],
        dip_1to2_f: Callable[[Dip2RefData, Dip1BedFileOrCoords], X],
        dip_2to1_f: Callable[[Dip1RefData, Dip2BedFileOrCoords], X],
        dip_2to2_f: Callable[[Dip2RefData, Dip2BedFileOrCoords], X],
    ) -> X:
        """Lookup refdata and bed file according to the given function.
        Then apply function depending on if the bed is haploid or diploid and
        if the reference is hap/dip1/dip2.
        """
        return self.with_ref_data(
            rk,
            lambda rd: (
                raise_inline() if (b := get_bed_f(rd)) is None else hap_f(rd, b)
            ),
            lambda rd: with_dip_bedfile(
                raise_inline() if (b := get_bed_f(rd)) is None else b,
                lambda bf: dip_1to1_f(rd, bf),
                lambda bf: dip_2to1_f(rd, bf),
            ),
            lambda rd: with_dip_bedfile(
                raise_inline() if (b := get_bed_f(rd)) is None else b,
                lambda bf: dip_1to2_f(rd, bf),
                lambda bf: dip_2to2_f(rd, bf),
            ),
        )

    def with_ref_data_and_bed_full(
        self,
        rk: RefKeyFullS,
        get_bed_f: RefDataToBed,
        hap_f: Callable[[HapRefData, HapBedFileOrCoords], Z],
        dip_1to1_f: Callable[[Dip1RefData, Dip1BedFileOrCoords], Z],
        dip_1to2_f: Callable[[Haplotype, Dip2RefData, Dip1BedFileOrCoords], Z],
        dip_2to1_f: Callable[[Dip1RefData, Dip2BedFileOrCoords], Z],
        dip_2to2_f: Callable[[Haplotype, Dip2RefData, Dip2BedFileOrCoords], Z],
    ) -> Z:
        """Like 'with_ref_data_and_bed' but also takes a full refkey and
        supplies the haplotype to the dip2-reference case."""
        rk_, hap = parse_full_refkey(rk)
        return self.with_ref_data_and_bed(
            rk_,
            get_bed_f,
            hap_f,
            dip_1to1_f,
            lambda rd, bf: not_none_unsafe(hap, lambda h: dip_1to2_f(h, rd, bf)),
            dip_2to1_f,
            lambda rd, bf: not_none_unsafe(hap, lambda h: dip_2to2_f(h, rd, bf)),
        )

    def with_build_data_and_bed(
        self,
        rk: RefKey,
        bk: BuildKey,
        get_bed_f: BuildDataToBed,
        hap_f: Callable[[HapBuildData, HapBedFileOrCoords], Z],
        dip_1to1_f: Callable[[Dip1BuildData, Dip1BedFileOrCoords], Z],
        dip_1to2_f: Callable[[Dip2BuildData, Dip1BedFileOrCoords], Z],
        dip_2to1_f: Callable[[Dip1BuildData, Dip2BedFileOrCoords], Z],
        dip_2to2_f: Callable[[Dip2BuildData, Dip2BedFileOrCoords], Z],
    ) -> Z:
        """Like 'with_ref_data_and_bed' but for build data."""
        return self.with_ref_data_and_bed(
            rk,
            lambda rd: get_bed_f(rd.to_build_data_unsafe(bk)),
            lambda rd, bf: hap_f(rd.to_build_data_unsafe(bk), bf),
            lambda rd, bf: dip_1to1_f(rd.to_build_data_unsafe(bk), bf),
            lambda rd, bf: dip_1to2_f(rd.to_build_data_unsafe(bk), bf),
            lambda rd, bf: dip_2to1_f(rd.to_build_data_unsafe(bk), bf),
            lambda rd, bf: dip_2to2_f(rd.to_build_data_unsafe(bk), bf),
        )

    def with_build_data_and_bed_full(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
        get_bed_f: BuildDataToBed,
        hap_f: Callable[[HapBuildData, HapBedFileOrCoords], Z],
        dip_1to1_f: Callable[[Dip1BuildData, Dip1BedFileOrCoords], Z],
        dip_1to2_f: Callable[[Haplotype, Dip2BuildData, Dip1BedFileOrCoords], Z],
        dip_2to1_f: Callable[[Dip1BuildData, Dip2BedFileOrCoords], Z],
        dip_2to2_f: Callable[[Haplotype, Dip2BuildData, Dip2BedFileOrCoords], Z],
    ) -> Z:
        """Like 'with_ref_data_and_bed_full' but for build data."""
        return self.with_ref_data_and_bed_full(
            rk,
            lambda rd: get_bed_f(rd.to_build_data_unsafe(bk)),
            lambda rd, bf: hap_f(rd.to_build_data_unsafe(bk), bf),
            lambda rd, bf: dip_1to1_f(rd.to_build_data_unsafe(bk), bf),
            lambda hap, rd, bf: dip_1to2_f(hap, rd.to_build_data_unsafe(bk), bf),
            lambda rd, bf: dip_2to1_f(rd.to_build_data_unsafe(bk), bf),
            lambda hap, rd, bf: dip_2to2_f(hap, rd.to_build_data_unsafe(bk), bf),
        )

    def with_build_data_and_bed_i(
        self,
        rk: RefKey,
        bk: BuildKey,
        inputs: list[X],
        get_bed_f: BuildDataToBed,
        hap_f: Callable[[X, HapBuildData, HapBedFile], Z],
        dip_1to1_f: Callable[[X, Dip1BuildData, Dip1BedFile], Z],
        dip_1to2_f: Callable[[X, Dip2BuildData, Dip1BedFile], Z],
        dip_2to1_f: Callable[[Double[X], Dip1BuildData, Dip2BedFile], Z],
        dip_2to2_f: Callable[[Double[X], Dip2BuildData, Dip2BedFile], Z],
    ) -> Z:
        """Like 'with_build_data_and_bed' but also take a list of input files
        and supply it to one of the supplied higher-order functions.

        Throw DesignError if the list of inputs does not correspond to the
        expected number of inputs (either 1 or 2).

        """
        return self.with_build_data_and_bed(
            rk,
            bk,
            get_bed_f,
            lambda bd, bf: (
                match1_unsafe(inputs, lambda i: hap_f(i, bd, bf))
                if isinstance(bf, BedFile)
                else raise_inline()
            ),
            lambda bd, bf: (
                match1_unsafe(inputs, lambda i: dip_1to1_f(i, bd, bf))
                if isinstance(bf, BedFile)
                else raise_inline()
            ),
            lambda bd, bf: (
                match1_unsafe(inputs, lambda i: dip_1to2_f(i, bd, bf))
                if isinstance(bf, BedFile)
                else raise_inline()
            ),
            lambda bd, bf: match2_unsafe(
                inputs,
                lambda i: (
                    dip_2to1_f(i, bd, bf) if isinstance(bf, BedFile) else raise_inline()
                ),
            ),
            lambda bd, bf: match2_unsafe(
                inputs,
                lambda i: (
                    dip_2to2_f(i, bd, bf) if isinstance(bf, BedFile) else raise_inline()
                ),
            ),
        )

    def with_build_data_and_bed_o2(
        self,
        rk: RefKey,
        bk: BuildKey,
        output_f: Callable[[RefKeyFull], Y],
        write_outputs: Callable[[list[Y]], None],
        get_bed_f: BuildDataToBed,
        hap_f: Callable[[Y, HapBuildData, HapBedFileOrCoords], Z],
        dip_1to1_f: Callable[[Y, Dip1BuildData, Dip1BedFileOrCoords], Z],
        dip_1to2_f: Callable[[Double[Y], Dip2BuildData, Dip1BedFileOrCoords], Z],
        dip_2to1_f: Callable[[Y, Dip1BuildData, Dip2BedFileOrCoords], Z],
        dip_2to2_f: Callable[[Double[Y], Dip2BuildData, Dip2BedFileOrCoords], Z],
    ) -> Z:
        def out1(src: Single[RefSrc]) -> Y:
            return with_first(output_f(src.key(rk).elem), lambda o: write_outputs([o]))

        def out2(src: Double[RefSrc]) -> Double[Y]:
            return with_first(
                src.keys(rk).both_(output_f),
                lambda o: write_outputs(o.as_list),
            )

        return self.with_build_data_and_bed(
            rk,
            bk,
            get_bed_f,
            lambda bd, bf: hap_f(out1(bd.refdata.ref.src), bd, bf),
            lambda bd, bf: dip_1to1_f(out1(bd.refdata.ref.src), bd, bf),
            lambda bd, bf: dip_1to2_f(out2(bd.refdata.ref.src), bd, bf),
            lambda bd, bf: dip_2to1_f(out1(bd.refdata.ref.src), bd, bf),
            lambda bd, bf: dip_2to2_f(out2(bd.refdata.ref.src), bd, bf),
        )

    def with_build_data_and_bed_io2(
        self,
        rk: RefKey,
        bk: BuildKey,
        inputs: list[X],
        output_f: Callable[[RefKeyFull], Y],
        write_outputs: Callable[[list[Y]], None],
        get_bed_f: BuildDataToBed,
        hap_f: Callable[[X, Y, HapBuildData, HapBedFile], Z],
        dip_1to1_f: Callable[[X, Y, Dip1BuildData, Dip1BedFile], Z],
        dip_1to2_f: Callable[[X, Double[Y], Dip2BuildData, Dip1BedFile], Z],
        dip_2to1_f: Callable[[Double[X], Y, Dip1BuildData, Dip2BedFile], Z],
        dip_2to2_f: Callable[[Double[X], Double[Y], Dip2BuildData, Dip2BedFile], Z],
    ) -> Z:
        """Like '_with_build_data_and_bed_i' but also take a function that
        generates an output file path and function that writes the outputs paths
        (NOT the data itself) to disk. The five higher order functions
        corresponding to each of the haplotype configurations must take the
        correct number of output paths which they are assumed to use for
        writing.
        """

        def out1(src: Single[RefSrc]) -> Y:
            return with_first(output_f(src.key(rk).elem), lambda o: write_outputs([o]))

        def out2(src: Double[RefSrc]) -> Double[Y]:
            return with_first(
                src.keys(rk).both_(output_f),
                lambda o: write_outputs(o.as_list),
            )

        return self.with_build_data_and_bed_i(
            rk,
            bk,
            inputs,
            get_bed_f,
            lambda i, bd, bf: hap_f(i, out1(bd.refdata.ref.src), bd, bf),
            lambda i, bd, bf: dip_1to1_f(i, out1(bd.refdata.ref.src), bd, bf),
            lambda i, bd, bf: dip_1to2_f(i, out2(bd.refdata.ref.src), bd, bf),
            lambda i, bd, bf: dip_2to1_f(i, out1(bd.refdata.ref.src), bd, bf),
            lambda i, bd, bf: dip_2to2_f(i, out2(bd.refdata.ref.src), bd, bf),
        )

    def with_build_data_and_bed_o(
        self,
        rk: RefKey,
        bk: BuildKey,
        output: Path,
        output_pattern: str,
        get_bed_f: BuildDataToBed,
        hap_f: Callable[[Path, HapBuildData, HapBedCoords], list[Path]],
        dip_1to1_f: Callable[[Path, Dip1BuildData, Dip1BedCoords], list[Path]],
        dip_1to2_f: Callable[[Double[Path], Dip2BuildData, Dip1BedCoords], list[Path]],
        dip_2to1_f: Callable[[Path, Dip1BuildData, Dip2BedCoords], list[Path]],
        dip_2to2_f: Callable[
            [Path, Haplotype, Dip2BuildData, Dip2BedCoords], list[Path]
        ],
    ) -> list[Path]:
        def write_output(ps: list[Path]) -> None:
            with open(output, "w") as f:
                json.dump([str(p) for p in ps], f)

        return self.with_build_data_and_bed_o2(
            rk,
            bk,
            lambda rk: sub_output_path(output_pattern, rk),
            write_output,
            get_bed_f,
            lambda o, bd, bf: (
                hap_f(o, bd, bf) if not isinstance(bf, BedFile) else raise_inline()
            ),
            lambda o, bd, bf: (
                dip_1to1_f(o, bd, bf) if not isinstance(bf, BedFile) else raise_inline()
            ),
            lambda o, bd, bf: (
                dip_1to2_f(o, bd, bf) if not isinstance(bf, BedFile) else raise_inline()
            ),
            lambda o, bd, bf: (
                dip_2to1_f(o, bd, bf) if not isinstance(bf, BedFile) else raise_inline()
            ),
            lambda o, bd, bf: [
                y
                for x in o.both(
                    lambda o, hap: (
                        dip_2to2_f(o, hap, bd, bf)
                        if not isinstance(bf, BedFile)
                        else raise_inline()
                    )
                ).as_list
                for y in x
            ],
        )

    def with_build_data_and_bed_io(
        self,
        rk: RefKey,
        bk: BuildKey,
        inputs: list[X],
        output: Path,
        output_pattern: str,
        get_bed_f: BuildDataToBed,
        hap_f: Callable[[X, Path, HapBuildData, HapBedFile], list[Path]],
        dip_1to1_f: Callable[[X, Path, Dip1BuildData, Dip1BedFile], list[Path]],
        dip_1to2_f: Callable[[X, Double[Path], Dip2BuildData, Dip1BedFile], list[Path]],
        dip_2to1_f: Callable[[Double[X], Path, Dip1BuildData, Dip2BedFile], list[Path]],
        dip_2to2_f: Callable[
            [X, Path, Haplotype, Dip2BuildData, Dip2BedFile], list[Path]
        ],
    ) -> list[Path]:
        """Like 'with_build_data_and_bed_io2' with the following differences.

        * This function takes an input pattern to be used to generate the list
        of output files. This pattern must have a '%s' for where the full ref
        key will be subbed, and must not contain any unexpanded snakemake
        wildcards.

        * This functions takes a single path to which the output list will be
        written instead of a function. The output list will be written as a json
        dump.

        * here, the function corresponding to dip2->dip2 here only takes a
        single haplotype, which is useful since the haplotypes can be processed
        independently and thus the supplied function can be made less
        redundant.

        """

        def write_output(ps: list[Path]) -> None:
            with open(output, "w") as f:
                json.dump([str(p) for p in ps], f)

        return self.with_build_data_and_bed_io2(
            rk,
            bk,
            inputs,
            lambda rk: sub_output_path(output_pattern, rk),
            write_output,
            get_bed_f,
            hap_f,
            dip_1to1_f,
            dip_1to2_f,
            dip_2to1_f,
            lambda i, o, bd, bf: [
                y
                for x in i.both(
                    lambda i, hap: dip_2to2_f(i, o.choose(hap), hap, bd, bf)
                ).as_list
                for y in x
            ],
        )

    # TODO this really should be called "source_doc" or something
    def with_build_data_and_bed_doc(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
        inputs: Path0or1or2,
        get_bed_f: BuildDataToBed,
        this: str | None,
        level: int,
    ) -> str:
        """Format readme documentation for a given source.

        Input file(s) will be used to generate md5 hashes.

        'extra' will be appended as an additional paragraph to the end.

        Return a string of rmarkdown paragraphs formatted with 80-char wrap.
        """

        def hap(bf: HapBedFileOrCoords) -> list[str]:
            if isinstance(bf, BedFile) and not isinstance(inputs, Null):
                src_txt = match_single_unsafe(
                    lambda p: format_src(bf.bed.documentation.elem, p, this), inputs
                )
                params_txt = format_bed_params(bf.params)
                return [src_txt, params_txt]
            elif isinstance(bf, HapBedCoords) and isinstance(inputs, Null):
                return [bf.hap.provenance]
            else:
                raise DesignError()

        def dip1to1(bf: Dip1BedFileOrCoords) -> list[str]:
            if isinstance(bf, BedFile) and not isinstance(inputs, Null):
                dip_txt = "This source contained both haplotypes for this reference."
                src_txt = match_single_unsafe(
                    lambda p: format_src(bf.bed.documentation.elem, p, this), inputs
                )
                params_txt = format_bed_params(bf.params)
                return [dip_txt, src_txt, params_txt]
            elif isinstance(bf, Dip1BedCoords) and isinstance(inputs, Null):
                return [bf.dip.provenance]
            else:
                raise DesignError()

        def dip1to2(bf: Dip1BedFileOrCoords, hap: Haplotype) -> list[str]:
            dip_txt = " ".join(
                [
                    "This source contained both haplotypes for this reference,",
                    f"and only the {hap.name} haplotype was used to generate",
                    "this bed file.",
                ]
            )
            if isinstance(bf, BedFile) and not isinstance(inputs, Null):
                src_txt = match_single_unsafe(
                    lambda p: format_src(bf.bed.documentation.elem, p, this), inputs
                )
                params_txt = format_bed_params(bf.params)
                return [dip_txt, src_txt, params_txt]
            elif isinstance(bf, Dip1BedCoords) and isinstance(inputs, Null):
                return [dip_txt, bf.dip.provenance]
            else:
                raise DesignError()

        def dip2to1(bf: Dip2BedFileOrCoords) -> list[str]:
            indent = level * "#"
            dip_txt = " ".join(
                [
                    "The two haplotypes for this reference were obtained from",
                    "different source files",
                ]
            )

            def hap_header(h: Haplotype) -> str:
                t = h.choose("P", "M") + "ATERNAL"  # code golf ;)
                return f"{indent} {t} haplotype"

            if isinstance(bf, BedFile) and not isinstance(inputs, Null):

                def fmt_src(h: Haplotype) -> list[str]:
                    src_txt = match_double_unsafe(
                        lambda p: format_src(
                            bf.bed.documentation.choose(h), p.choose(h), this
                        ),
                        inputs,
                    )
                    return [hap_header(h), src_txt]

                src_paras = [*flatten(make_double(fmt_src).as_list)]
                params_txt = format_bed_params(bf.params)
                return [
                    dip_txt,
                    *src_paras,
                    f"{indent} Both haplotypes",
                    params_txt,
                ]
            elif isinstance(bf, Dip2BedCoords) and isinstance(inputs, Null):
                return [
                    dip_txt,
                    *flatten(
                        bf.double.both(
                            lambda x, hap: [hap_header(hap), x.provenance]
                        ).as_list
                    ),
                ]
            else:
                raise DesignError()

        def dip2to2(bf: Dip2BedFileOrCoords, hap: Haplotype) -> list[str]:
            dip_txt = " ".join(
                [
                    f"This source contained only the {hap.name} haplotype for",
                    "this reference, which was used to generate this bed file.",
                ]
            )
            if isinstance(bf, BedFile) and not isinstance(inputs, Null):
                src_txt = match_double_unsafe(
                    lambda p: format_src(
                        bf.bed.documentation.choose(hap), p.choose(hap), this
                    ),
                    inputs,
                )
                params_txt = format_bed_params(bf.params)
                return [dip_txt, src_txt, params_txt]
            elif isinstance(bf, Dip2BedCoords) and isinstance(inputs, Null):
                return [dip_txt, bf.double.choose(hap).provenance]
            else:
                raise DesignError()

        paragraphs = self.with_build_data_and_bed_full(
            rk,
            bk,
            get_bed_f,
            lambda _, bf: hap(bf),
            lambda _, bf: dip1to1(bf),
            lambda hap, _, bf: dip1to2(bf, hap),
            lambda _, bf: dip2to1(bf),
            lambda hap, _, bf: dip2to2(bf, hap),
        )

        return "\n\n".join([readme_fill(p) for p in paragraphs])

    # final refkey/buildkey lists (for the "all" target and related)

    @property
    def all_build_keys(self) -> tuple[list[RefKey], list[BuildKey]]:
        return unzip2(
            all_build_keys(self.haploid_stratifications)
            + all_build_keys(self.diploid1_stratifications)
            + all_build_keys(self.diploid2_stratifications)
        )

    @property
    def all_full_build_keys(self) -> tuple[list[RefKeyFullS], list[BuildKey]]:
        return unzip2(
            all_ref_build_keys(self.haploid_stratifications)
            + all_ref_build_keys(self.diploid1_stratifications)
            + all_ref_build_keys(self.diploid2_stratifications)
        )

    @property
    def all_full_ref_and_build_keys(self) -> list[tuple[RefKeyFullS, BuildKey]]:
        return (
            all_ref_build_keys(self.haploid_stratifications)
            + all_ref_build_keys(self.diploid1_stratifications)
            + all_ref_build_keys(self.diploid2_stratifications)
        )

    # source refkey/buildkey lists (for the "all resources" rule)

    @property
    def all_ref_refsrckeys(self) -> list[RefKeyFullS]:
        return (
            all_ref_refsrckeys(self.haploid_stratifications)
            + all_ref_refsrckeys(self.diploid1_stratifications)
            + all_ref_refsrckeys(self.diploid2_stratifications)
        )

    def _all_bed_build_and_refsrckeys(
        self, f: BuildDataToSrc
    ) -> list[tuple[RefKeyFullS, BuildKey]]:
        return (
            all_bed_build_and_refsrckeys(self.haploid_stratifications, f)
            + all_bed_build_and_refsrckeys(self.diploid1_stratifications, f)
            + all_bed_build_and_refsrckeys(self.diploid2_stratifications, f)
        )

    def _all_bed_refsrckeys(self, f: BuildDataToSrc) -> list[RefKeyFullS]:
        return [rk for rk, _ in self._all_bed_build_and_refsrckeys(f)]

    # @property
    # def all_refkey_gap(self) -> list[RefKeyFullS]:
    #     return self._all_bed_refsrckeys(
    #         lambda bd: fmap_maybe(
    #             lambda x: x.bed.src if isinstance(x, BedFile) else None, bd_to_gaps(bd)
    #         )
    #     )

    # @property
    # def all_refkey_rmsk(self) -> list[RefKeyFullS]:
    #     return self._all_bed_refsrckeys(
    #         lambda bd: fmap_maybe(
    #             lambda x: x.bed.src if isinstance(x, BedFile) else None, bd_to_rmsk(bd)
    #         )
    #     )

    # @property
    # def all_refkey_simreps(self) -> list[RefKeyFullS]:
    #     return self._all_bed_refsrckeys(
    #         lambda bd: fmap_maybe(
    #             lambda x: x.bed.src if isinstance(x, BedFile) else None,
    #             bd_to_simreps(bd),
    #         )
    #     )

    # @property
    # def all_refkey_censat(self) -> list[RefKeyFullS]:
    #     return self._all_bed_refsrckeys(
    #         lambda bd: fmap_maybe(
    #             lambda x: x.bed.src if isinstance(x, BedFile) else None,
    #             bd_to_satellites(bd),
    #         )
    #     )

    # @property
    # def all_refkey_segdups(self) -> list[RefKeyFullS]:
    #     return self._all_bed_refsrckeys(
    #         lambda bd: fmap_maybe(
    #             lambda x: x.bed.src if isinstance(x, BedFile) else None,
    #             bd_to_superdups(bd),
    #         )
    #     )

    @property
    def all_buildkey_bench(self) -> list[tuple[RefKeyFullS, BuildKey]]:
        return self._all_bed_build_and_refsrckeys(
            lambda bd: fmap_maybe(
                lambda x: x.bed.src if isinstance(x, BedFile) else None,
                bd_to_bench_bed(bd),
            )
        )

    # @property
    # def all_refkey_cds(self) -> list[RefKeyFullS]:
    #     return self._all_bed_refsrckeys(
    #         lambda bd: fmap_maybe(
    #             lambda x: x.bed.src if isinstance(x, BedFile) else None, bd_to_cds(bd)
    #         )
    #     )

    # @property
    # def all_refkey_mhc(self) -> list[RefKeyFullS]:
    #     return self._all_bed_refsrckeys(
    #         lambda bd: fmap_maybe(
    #             lambda x: x.bed.src if isinstance(x, BedFile) else None, bd_to_mhc(bd)
    #         )
    #     )

    # @property
    # def all_refkey_kir(self) -> list[RefKeyFullS]:
    #     return self._all_bed_refsrckeys(
    #         lambda bd: fmap_maybe(
    #             lambda x: x.bed.src if isinstance(x, BedFile) else None, bd_to_kir(bd)
    #         )
    #     )

    # @property
    # def all_refkey_vdj(self) -> list[RefKeyFullS]:
    #     return self._all_bed_refsrckeys(
    #         lambda bd: fmap_maybe(
    #             lambda x: x.bed.src if isinstance(x, BedFile) else None, bd_to_vdj(bd)
    #         )
    #     )

    # source and output functions

    def _test_if_final_path(
        self,
        p: Any,
        c: str,
        rk: RefKeyFullS,
        bk: BuildKey,
        sub: bool = False,
    ) -> None:
        if isinstance(p, Path):
            comp = self.final_build_dir / c
            _comp = (
                sub_wildcards_path(comp, {"ref_final_key": rk, "build_key": bk})
                if sub
                else comp
            )
            if _comp != p.parent:
                raise DesignError(f"{p} is not a within {_comp}")
        else:
            raise DesignError(f"{p} is not a path")

    def _test_if_final_paths(
        self,
        p: Any,
        c: str,
        rk: RefKeyFullS,
        bk: BuildKey,
        sub: bool = False,
    ) -> None:
        if isinstance(p, list):
            for x in p:
                self._test_if_final_path(x, c, rk, bk, sub)
        else:
            raise DesignError(f"Not a list {p}")

    def _test_if_final_mutual_path(
        self,
        p: Any,
        c: str,
        rk: RefKeyFullS,
        bk: BuildKey,
        sub: bool = False,
    ) -> None:
        if isinstance(p, MutualPathPair):
            self._test_if_final_path(p.positive, c, rk, bk, sub)
            self._test_if_final_path(p.negative, c, rk, bk, sub)
        else:
            raise DesignError(f"Not a mutual pair path: {p}")

    def _test_if_source_path(self, p: Any, rk: RefKey) -> None:
        if isinstance(p, Path):
            if not p.is_relative_to(self.resources_dir):
                raise DesignError(f"{p} is not a within {self.resources_dir}")
        else:
            raise DesignError(f"{p} is not a path")

    def _sub_rsk(
        self,
        p: Path,
        f: StratInputToBed,
        rk: RefKey,
    ) -> Path0or1or2 | None:
        self._test_if_source_path(p, rk)
        src = self.with_ref_data(
            rk,
            lambda rd: f(rd.strat_inputs),
            lambda rd: f(rd.strat_inputs),
            lambda rd: f(rd.strat_inputs),
        )
        if isinstance(src, BedFile):
            return (
                rsks.map(lambda s: sub_wildcards_path(p, {"ref_src_key": s}))
                if (rsks := self.refkey_to_bed_refsrckeys(f, rk)) is not None
                else None
            )
        elif src is None:
            return None
        else:
            return Null()

    def all_low_complexity_sources(
        self,
        rk: RefKey,
        rmsk: Path,
        censat: Path,
        trf: Path,
    ) -> LowComplexitySources:
        return LowComplexitySources(
            rmsk=self._sub_rsk(rmsk, si_to_rmsk, rk),
            sat=self._sub_rsk(censat, si_to_satellites, rk),
            trf=self._sub_rsk(trf, si_to_simreps, rk),
        )

    def all_low_complexity(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
        # sources
        src: LowComplexitySources,
        # uniform
        perfect: list[Path],
        imperfect: list[Path],
        homopolymers: MutualPathPair,
        # satellites
        sats: MutualPathPair,
        # all repeats
        filtered_trs: list[Path],
        all_trs: MutualPathPair,
        all_repeats: MutualPathPair,
        readme: Path,
    ) -> LowComplexityPaths | None:
        test_params = (CoreLevel.LOWCOMPLEXITY.value, rk, bk)
        self._test_if_final_paths(perfect, *test_params)
        self._test_if_final_paths(imperfect, *test_params)
        self._test_if_final_mutual_path(homopolymers, *test_params)
        self._test_if_final_mutual_path(sats, *test_params)
        self._test_if_final_paths(filtered_trs, *test_params)
        self._test_if_final_mutual_path(all_trs, *test_params)
        self._test_if_final_mutual_path(all_repeats, *test_params)
        self._test_if_final_path(readme, *test_params)

        bd = self.to_build_data_full(rk, bk)

        if not bd.want_low_complexity:
            return None
        # homopolymers and uniform repeats are included no matter what
        uniform = UniformRepeatPaths(
            perfect=perfect,
            imperfect=imperfect,
            homopolymers=homopolymers,
        )

        # include tandem repeats and merged output if we have rmsk/censat and simreps
        tm: tuple[Path0or1or2, Path0or1or2] | None = maybe2((src.trf, src.rmsk))
        if tm is None:
            repeats = None
        else:
            repeats = RepeatsPaths(
                trf_src=tm[0],
                rmsk_src=tm[1],
                filtered_trs=filtered_trs,
                all_trs=all_trs,
                all_repeats=all_repeats,
            )

        # include satellites only if we have rmsk or censat
        s: Path0or1or2 | None = from_maybe(src.rmsk, src.sat)
        if s is None:
            satpaths = None
        else:
            satpaths = SatellitesPaths(
                sat_src=s,
                sats=sats,
                used_censat=src.sat is not None,
                all_repeats=repeats,
            )

        return LowComplexityPaths(
            readme=readme,
            uniform_repeats=uniform,
            satellites=satpaths,
        )

    def all_otherdifficult_sources(
        self,
        rk: RefKey,
        bk: BuildKey,
        gaps: Path,
        refseq: Path,
        kir: Path,
        mhc: Path,
        vdj: Path,
        other: Path,
    ) -> OtherDifficultSources:
        self._test_if_source_path(other, rk)

        bd = self.to_build_data(rk, bk)

        return OtherDifficultSources(
            gaps=self._sub_rsk(gaps, si_to_gaps, rk),
            refseq=self._sub_rsk(refseq, si_to_cds, rk),
            mhc=self._sub_rsk(mhc, si_to_mhc, rk),
            kir=self._sub_rsk(kir, si_to_kir, rk),
            vdj=self._sub_rsk(vdj, si_to_vdj, rk),
            other={
                lk: {
                    sk: (
                        to_str_refkeys(v.data.bed.src, rk).map(
                            lambda r: sub_wildcards_path(
                                other,
                                {
                                    "ref_src_key": r,
                                    "other_level_key": lk,
                                    "other_strat_key": sk,
                                },
                            )
                        )
                        if isinstance(v.data, BedFile)
                        else Null()
                    )
                    for sk, v in o.items()
                }
                for lk, o in bd.build.other_strats.items()
            },
        )

    def all_functional(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
        src: OtherDifficultSources,
        cds: MutualPathPair,
        readme: Path,
    ) -> FunctionalPaths | None:
        test_params = (CoreLevel.FUNCTIONAL.value, rk, bk)
        self._test_if_final_mutual_path(cds, *test_params)
        self._test_if_final_path(readme, *test_params)

        bd = self.to_build_data_full(rk, bk)

        if not bd.want_cds:
            return None

        if src.refseq is None:
            return None

        return FunctionalPaths(
            readme=readme,
            cds_source=src.refseq,
            cds_output=cds,
        )

    def all_otherdifficult(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
        src: OtherDifficultSources,
        gaps: Path,
        kir: Path,
        mhc: Path,
        vdj: Path,
        other: Path,
        readme: Path,
    ) -> OtherDifficultPaths:
        test_params = (CoreLevel.OTHER_DIFFICULT.value, rk, bk)
        self._test_if_final_path(gaps, *test_params)
        self._test_if_final_path(kir, *test_params)
        self._test_if_final_path(mhc, *test_params)
        self._test_if_final_path(vdj, *test_params)
        self._test_if_final_path(
            sub_wildcard_path(other, "other_level_key", OTHERDIFF_KEY),
            *test_params,
        )
        self._test_if_final_path(readme, *test_params)

        bd = self.to_build_data_full(rk, bk)

        def immuno_paths(
            test: bool,
            immuno_src: Path0or1or2 | None,
            output: Path,
        ) -> ImmunoPaths | None:
            if not test:
                return None
            elif immuno_src is not None:
                isrc = immuno_src
                source_is_refseq = False
            elif src.refseq is not None:
                isrc = src.refseq
                source_is_refseq = True
            else:
                return None

            return ImmunoPaths(
                source=isrc,
                output=output,
                source_is_refseq=source_is_refseq,
            )

        other_diff = src.other[OTHERDIFF_KEY] if OTHERDIFF_KEY in src.other else {}
        other_paths = {
            k: SourceOutputPaths(
                source=v,
                output=sub_wildcards_path(
                    other,
                    {"other_level_key": OTHERDIFF_KEY, "other_strat_key": k},
                ),
            )
            for k, v in other_diff.items()
        }

        return OtherDifficultPaths(
            readme=readme,
            gaps=(
                SourceOutputPaths(source=src.gaps, output=gaps)
                if src.gaps is not None
                else None
            ),
            kir=immuno_paths(bd.want_kir, src.kir, kir),
            mhc=immuno_paths(bd.want_mhc, src.mhc, mhc),
            vdj=immuno_paths(bd.want_vdj, src.vdj, vdj),
            other=other_paths,
        )

    def all_misc(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
        lk: OtherLevelKey,
        src: OtherDifficultSources,
        other: Path,
        readme: Path,
    ) -> MiscPaths | None:
        def test(p: Path) -> None:
            _p = sub_wildcard_path(p, "other_level_key", lk)
            self._test_if_final_path(_p, lk, rk, bk)

        test(other)
        test(readme)

        _other = src.other[lk] if lk in src.other else {}

        if len(_other) == 0:
            return None

        desc = next((x for x in self.other_levels if x.key == lk), None)

        if desc is None:
            raise DesignError()

        return MiscPaths(
            readme=sub_wildcard_path(readme, "other_level_key", lk),
            desc=desc,
            paths={
                sk: SourceOutputPaths(
                    source=s,
                    output=sub_wildcards_path(
                        other,
                        {"other_level_key": lk, "other_strat_key": sk},
                    ),
                )
                for sk, s in _other.items()
            },
        )

    def all_union(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
        # sources
        segdups_src: SegdupPaths | None,
        lowmap: Path,
        gc_src: Path,
        lc: LowComplexityPaths | None,
        xy: SexPaths | None,
        # outputs
        segdup_lowmap_output: MutualPathPair,
        all_difficult: MutualPathPair,
        readme: Path,
    ) -> UnionPaths | None:
        test_params = (CoreLevel.UNION.value, rk, bk)
        self._test_if_final_mutual_path(segdup_lowmap_output, *test_params)
        self._test_if_final_mutual_path(all_difficult, *test_params)
        self._test_if_final_path(readme, *test_params)

        bd = self.to_build_data_full(rk, bk)

        if segdups_src is not None and not isinstance(segdups_src, SegdupPaths):
            raise DesignError()

        if lc is not None and not isinstance(lc, LowComplexityPaths):
            raise DesignError()

        if xy is not None and not isinstance(xy, SexPaths):
            raise DesignError()

        if not bd.want_union:
            return None

        # can't do anything if we don't have segdups or lowmap
        if segdups_src is None or not bd.want_mappability:
            return None

        sl = SegdupLowmapPaths(
            segdup_input=segdups_src,
            lowmap_input=lowmap,
            output=segdup_lowmap_output,
        )

        # make the alldifficult file if we have 2 of xy, repeats, or gc
        _gc_src = gc_src if bd.want_gc else None
        _repeat_src: Path | None = fmap_maybe(
            lambda x: fmap_maybe(
                lambda y: fmap_maybe(lambda z: z.all_repeats.positive, y.all_repeats),
                x.satellites,
            ),
            lc,
        )
        empty: list[Path] = []
        _xy_src = fmap_maybe_def(empty, lambda x: x.sex.all_features, xy)

        all_diff = (
            AllDifficultPaths(
                gc_input=_gc_src,
                repeat_input=_repeat_src,
                xy_inputs=_xy_src,
                output=all_difficult,
            )
            if ((_gc_src is not None) + (_repeat_src is not None) + (len(_xy_src) > 0))
            > 1
            else None
        )

        return UnionPaths(
            readme=readme,
            segdup_lowmap=sl,
            all_difficult=all_diff,
        )

    def all_gc(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
        ranges: list[Path],
        extremes: list[Path],
        readme: Path,
    ) -> GCPaths:
        # sub wildcards here during test since these paths will come from a
        # checkpoint which will automatically fill in the refkey/buildkey
        test_params = (CoreLevel.GC.value, rk, bk)
        self._test_if_final_paths(ranges, *test_params, sub=True)
        self._test_if_final_paths(extremes, *test_params, sub=True)
        self._test_if_final_path(readme, *test_params)

        bd = self.to_build_data_full(rk, bk)

        if not bd.build.include.gc:
            raise DesignError()

        return GCPaths(
            readme=readme,
            lowGC=ranges[0],
            middleGC=ranges[1:-1],
            highGC=ranges[-1],
            extremes=extremes,
            params=bd.build.include.gc,
        )

    def all_lowmap(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
        union: MutualPathPair,
        single: list[Path],
        readme: Path,
    ) -> LowmapPaths:
        test_params = (CoreLevel.MAPPABILITY.value, rk, bk)
        self._test_if_final_mutual_path(union, *test_params)
        # sub wildcards here during test since these paths will come from a
        # checkpoint which will automatically fill in the refkey/buildkey
        self._test_if_final_paths(single, *test_params, sub=True)
        self._test_if_final_path(readme, *test_params)

        bd = self.to_build_data_full(rk, bk)

        if not bd.want_mappability:
            raise DesignError()

        return LowmapPaths(
            readme=readme,
            union=union,
            single=single,
            params=[*bd.build.include.mappability],
        )

    def all_segdups_sources(self, rk: RefKeyFullS, superdups: Path) -> SegdupSources:
        sd = self._sub_rsk(superdups, si_to_superdups, strip_full_refkey(rk))
        return SegdupSources(superdup=sd)

    def all_segdups(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
        src: SegdupSources,
        segdups: MutualPathPair,
        long_segdups: MutualPathPair,
        readme: Path,
    ) -> SegdupPaths | None:
        test_params = (CoreLevel.SEGDUPS.value, rk, bk)
        self._test_if_final_mutual_path(segdups, *test_params)
        self._test_if_final_mutual_path(long_segdups, *test_params)
        self._test_if_final_path(readme, *test_params)

        bd = self.to_build_data_full(rk, bk)

        if bd.want_segdups and src is not None and src.superdup is not None:
            return SegdupPaths(
                readme=readme,
                superdups=src.superdup,
                all_segdups=segdups,
                long_segdups=long_segdups,
            )
        else:
            return None

    def all_diploid(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
        sex: SexPaths | None,
        hets: list[Path],
        homs: list[Path],
        SNVorSV_hets: list[Path],
        SNVorSV_homs: list[Path],
        readme: Path,
    ) -> DiploidPaths | None:
        test_params = (CoreLevel.DIPLOID.value, rk, bk)
        self._test_if_final_paths(hets, *test_params)
        self._test_if_final_paths(SNVorSV_hets, *test_params)
        self._test_if_final_paths(homs, *test_params)
        self._test_if_final_paths(SNVorSV_homs, *test_params)
        self._test_if_final_path(readme, *test_params)

        bd = self.to_build_data_full(rk, bk)

        if not bd.want_hets:
            return None

        # ASSUME merge_len has already been subbed here
        return DiploidPaths(
            readme=readme,
            hets=hets,
            homs=homs,
            SNVorSV_hets=SNVorSV_hets,
            SNVorSV_homs=SNVorSV_homs,
            nonpar=fmap_maybe_def([], lambda x: x.sex.nonpar_paths, sex),
        )

    def all_telomeres(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
        telomeres: Path,
        readme: Path,
    ) -> TelomerePaths | None:
        test_params = (CoreLevel.TELOMERES.value, rk, bk)
        self._test_if_final_path(telomeres, *test_params)
        self._test_if_final_path(readme, *test_params)

        bd = self.to_build_data_full(rk, bk)

        if not bd.want_telomeres:
            return None

        return TelomerePaths(readme=readme, telomeres=telomeres)

    def all_xy(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
        # sources
        features_src: Path,
        # outputs
        xtr: Path,
        ampliconic: Path,
        par: MutualPathPair,
        auto: Path,
        readme: Path,
    ) -> SexPaths:
        test_params = (CoreLevel.XY.value, rk, bk)
        self._test_if_final_path(xtr, *test_params)
        self._test_if_final_path(ampliconic, *test_params)
        self._test_if_final_mutual_path(par, *test_params)
        self._test_if_final_path(auto, *test_params)
        self._test_if_final_path(readme, *test_params)

        _features_src = sub_wildcard_path(
            features_src,
            "ref_src_key",
            strip_full_refkey(rk),
        )

        self._test_if_source_path(_features_src, strip_full_refkey(rk))

        def sub_sex(p: Path, sub_x: bool) -> Path:
            c = (ChrIndex.CHRX if sub_x else ChrIndex.CHRY).chr_name
            return sub_wildcard_path(p, "sex_chr", c)

        def to_features(use_x: bool, features: XYFeatures) -> XYFeaturePaths:
            return XYFeaturePaths(
                src=sub_sex(_features_src, use_x),
                bed=features.x_bed if use_x else features.y_bed,
                xtr_path=sub_sex(xtr, use_x) if features.xtr is not None else None,
                ampliconic_path=(
                    sub_sex(ampliconic, use_x)
                    if features.ampliconic is not None
                    else None
                ),
                xtr=features.xtr,
                ampliconic=features.ampliconic,
            )

        def to_par(use_x: bool, xy: XY) -> PARPaths | None:
            return fmap_maybe(
                lambda z: PARPaths(
                    path=par.both(lambda x: sub_sex(x, use_x)),
                    doc=z.comment,
                ),
                xy.x_par if use_x else xy.y_par,
            )

        def to_paths(use_x: bool, xy: XY) -> SubSexPaths:
            return SubSexPaths(
                features=fmap_maybe(lambda z: to_features(use_x, z), xy.features),
                par=to_par(use_x, xy),
            )

        def hap(bd: HapBuildData) -> AnySexPaths:
            cis = bd.refdata.ref.hap_chrs(bd.build_chrs)
            xy = bd.refdata.strat_inputs.xy
            return MaleHapSexPaths(
                x=to_paths(True, xy) if ChrIndex.CHRX in cis else None,
                y=to_paths(False, xy) if ChrIndex.CHRY in cis else None,
            )

        def dip1(bd: Dip1BuildData) -> AnySexPaths:
            cis = bd.refdata.ref.all_chrs(bd.build_chrs)
            xy = bd.refdata.strat_inputs.xy
            return Dip1SexPaths(
                sex1=to_paths(True, xy) if ChrIndex.CHRX in cis else None,
                sex2=to_paths(False, xy) if ChrIndex.CHRY in cis else None,
            )

        def dip2(hap: Haplotype, bd: Dip2BuildData) -> AnySexPaths:
            cis = bd.refdata.ref.hap_chrs(bd.build_chrs, hap)
            use_x = not (hap is Haplotype.PAT)
            xy = bd.refdata.strat_inputs.xy
            return Dip2SexPaths(
                paths=(
                    SubSexPaths(
                        features=fmap_maybe(
                            lambda z: to_features(use_x, z),
                            xy.features,
                        ),
                        par=to_par(use_x, xy),
                    )
                    if (use_x and ChrIndex.CHRX in cis) or ChrIndex.CHRY in cis
                    else None
                ),
                hap=hap,
            )

        bd = self.to_build_data_full(rk, bk)

        sex = self.with_build_data_full(
            rk,
            bk,
            hap,
            dip1,
            dip2,
        )
        return SexPaths(
            readme=readme,
            sex=sex,
            auto=auto if bd.want_xy_auto else None,
        )

    # other nice functions

    def get_comparison(self, rk: RefKeyFullS, bk: BuildKey) -> BuildCompare1 | None:
        return self.with_build_data_full(
            rk,
            bk,
            lambda bd: bd.build.comparison,
            lambda bd: bd.build.comparison,
            lambda hap, bd: fmap_maybe(
                lambda c: c.double.choose(hap),
                bd.build.comparison,
            ),
        )

    def compare_key(self, rk: RefKeyFullS, bk: BuildKey) -> CompareKey | None:
        return fmap_maybe(lambda x: x.other, self.get_comparison(rk, bk))

    def refkey_haplotypes(self, rk: RefKeyFullS) -> list[Haplotype]:
        """Test if refkey is dip1 or dip2.

        Return True if dip1, false if dip2.

        If split is True, require dip1 refkey to have a haplotype and error
        otherwise. The reverse is True if split is False.

        If nohap is True, throw error if refkey is hap. If False permit the hap
        case and return False.

        """
        return self.with_ref_data_full(
            rk,
            lambda _: [],
            lambda _: [h for h in Haplotype],
            lambda hap, _: [hap],
        )

    def refkey_is_dip1(self, rk: RefKeyFullS, split: bool, nohap: bool) -> bool:
        """Test if refkey is dip1 or dip2.

        Return True if dip1, false if dip2.

        If split is True, require dip1 refkey to have a haplotype and error
        otherwise. The reverse is True if split is False.

        If nohap is True, throw error if refkey is hap. If False permit the hap
        case and return False.

        """
        return self.with_ref_data_full_rconf(
            rk,
            split,
            nohap,
            lambda _: False,
            lambda _: True,
            lambda _, __: True,
            lambda _, __: False,
        )

    def refkey_strip_if_dip1(self, rk: RefKeyFullS, nohap: bool) -> RefKeyFullS:
        """Remove haplotype from refkey if dip1.

        Note this assumes a split refkey configuration.
        """
        rk_ = RefKeyFull(strip_full_refkey(rk), None).name
        if nohap:
            return self.with_ref_data_split_full_nohap(
                rk,
                lambda _, __: rk_,
                lambda _, __: rk,
            )
        else:
            return self.with_ref_data_split_full(
                rk,
                lambda _: rk,
                lambda _, __: rk_,
                lambda _, __: rk,
            )

    def refkey_append_if_dip1(self, rk: RefKeyFullS) -> list[RefKeyFull]:
        rk_ = strip_full_refkey(rk)
        return self.with_ref_data_full_nohap(
            rk,
            lambda _: [RefKeyFull(rk_, h) for h in Haplotype],
            lambda _, __: raise_inline("dip2 case not allowed"),
        )

    def thread_per_chromosome(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
        n: int,
        split: bool,
        nohap: bool,
    ) -> int:
        return max(n, len(self.buildkey_to_chrs(rk, bk, split, nohap)))


################################################################################
# protocols

# hacky rankN type mimicry


class RefDataToBed(Protocol):
    A = TypeVar("A", HapBedSrc, DipBedSrc)
    B = TypeVar("B", HapBedCoords, DipBedCoords)

    def __call__(
        self, __x: RefData_[RefSrcT, A, VcfSrcT, B, BuildCompareT]
    ) -> BedFile[A] | B | None:
        pass


class RefDataToSrc(Protocol):
    A = TypeVar("A", Single[BedSrc], Double[BedSrc])

    def __call__(
        self, __x: RefData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT]
    ) -> A | None:
        pass


class StratInputToBed(Protocol):
    A = TypeVar("A", HapBedSrc, DipBedSrc)
    B = TypeVar("B", HapBedCoords, DipBedCoords)

    def __call__(self, __x: StratInputs[A, B]) -> BedFile[A] | B | None:
        pass


# TODO not sure how to get manual text through this
class StratInputToSrc(Protocol):
    A = TypeVar("A", Single[BedSrc], Double[BedSrc])

    def __call__(self, __x: StratInputs[BedSrcT, BedCoordsT]) -> A | None:
        pass


class BuildDataToBed(Protocol):
    A = TypeVar("A", HapBedSrc, DipBedSrc)
    B = TypeVar("B", HapBedCoords, DipBedCoords)

    def __call__(
        self, __x: BuildData_[RefSrcT, A, VcfSrcT, B, BuildCompareT]
    ) -> BedFile[A] | B | None:
        pass


class BuildDataToVCF(Protocol):
    A = TypeVar("A", HapVcfSrc, Dip1VcfSrc, Dip2VcfSrc)

    def __call__(
        self, __x: BuildData_[RefSrcT, BedSrcT, A, BedCoordsT, BuildCompareT]
    ) -> VCFFile[A] | None:
        pass


class BuildDataToSrc(Protocol):
    A = TypeVar("A", Single[BedSrc], Double[BedSrc])

    def __call__(
        self, __x: BuildData_[RefSrcT, BedSrcT, VcfSrcT, BedCoordsT, BuildCompareT]
    ) -> A | None:
        pass
