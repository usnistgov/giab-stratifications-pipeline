"""
Conventions:
* Functions ending in "_unsafe" should never throw errors; if they do then the
  code is incorrect. This is in contrast with other errors which may happen due
  to network issues, invalid inputs, etc.
"""
import pandas as pd
import re
from pathlib import Path
from pydantic import BaseModel as BaseModel_
from pydantic import validator, HttpUrl, FilePath, NonNegativeInt, Field
from dataclasses import dataclass, astuple
from enum import Enum, unique
from typing import (
    NewType,
    Any,
    Callable,
    TypeVar,
    Type,
    NamedTuple,
    cast,
    Annotated,
    Generic,
    Union,
    TypeGuard,
    Protocol,
)
from typing_extensions import Self, assert_never
from more_itertools import unique_everseen
from common.io import (
    is_gzip,
    is_bgzip,
    DesignError,
    match1_unsafe,
    match2_unsafe,
    not_none_unsafe,
    none_unsafe,
)
import common.bed as bed

W = TypeVar("W")
X = TypeVar("X")
Y = TypeVar("Y")
Z = TypeVar("Z")


Percent = Annotated[int, Field(ge=0, le=100)]


class HapRefKey_(str):
    pass


class Dip1RefKey_(str):
    pass


class Dip2RefKey_(str):
    pass


RefKeyT_ = TypeVar("RefKeyT_", HapRefKey_, Dip1RefKey_, Dip2RefKey_)


class RefKey_(Generic[RefKeyT_]):
    def __init__(self, k: RefKeyT_) -> None:
        self.key: RefKeyT_ = k

    @property
    def hap1(self) -> "RefKeyFull[RefKeyT_]":
        return RefKeyFull(self.key, Haplotype.HAP1)

    @property
    def hap2(self) -> "RefKeyFull[RefKeyT_]":
        return RefKeyFull(self.key, Haplotype.HAP2)

    @property
    def haps(self) -> "list[RefKeyFull[RefKeyT_]]":
        return [self.hap1, self.hap2]


HapRefKey = RefKey_[HapRefKey_]
Dip1RefKey = RefKey_[Dip1RefKey_]
Dip2RefKey = RefKey_[Dip2RefKey_]
AnyRefKey = HapRefKey | Dip1RefKey | Dip2RefKey


RefKeyT = TypeVar("RefKeyT", HapRefKey, Dip1RefKey, Dip2RefKey)


@dataclass(frozen=True)
class RefKeyFull(Generic[RefKeyT_]):
    key: RefKeyT_
    hap: "Haplotype" | None

    @property
    def strip(self) -> RefKey_[RefKeyT_]:
        return RefKey_(self.key)

    @property
    def as_tuple(self) -> tuple[str, "Haplotype" | None]:
        return (str(self.key), self.hap)

    @property
    def name(self) -> str:
        k, h = self.as_tuple
        return f"{k}_{h}" if h is not None else k


# class wrapper so I can pattern match on them (which newtype won't allow)
class HapBuildKey(str):
    pass


class Dip1BuildKey(str):
    pass


class Dip2BuildKey(str):
    pass


class BuildPair_(NamedTuple, Generic[X, Y]):
    ref: X
    build: Y


BuildKey = HapBuildKey | Dip1BuildKey | Dip2BuildKey
RefKey = HapRefKey | Dip1RefKey | Dip2RefKey
HaploidBuildPair = BuildPair_[HapRefKey, HapBuildKey]
Diploid1BuildPair = BuildPair_[Dip1RefKey, Dip1BuildKey]
Diploid2BuildPair = BuildPair_[Dip2RefKey, Dip2BuildKey]
BuildPair = HaploidBuildPair | Diploid1BuildPair | Diploid2BuildPair
BuildPairT = TypeVar(
    "BuildPairT",
    HaploidBuildPair,
    Diploid1BuildPair,
    Diploid2BuildPair,
)


# def either_ref_key(
#     left: Callable[[HapRefKey], X], right: Callable[[Dip1RefKey], X], k: RefKey
# ) -> X:
#     if isinstance(k, HapRefKey):
#         return left(k)
#     elif isinstance(k, Dip1RefKey):
#         return right(k)
#     elif isinstance(k, Dip2RefKey):
#         # TODO
#         return right(k)
#     else:
#         assert_never(k)


# def either_build_pair(
#     left: Callable[[HaploidBuildPair], X],
#     right: Callable[[Diploid1BuildPair], X],
#     p: BuildPair,
# ) -> X:
#     if isinstance(p.ref, HapRefKey):
#         return left(p)
#     elif isinstance(p.ref, Dip1RefKey):
#         return right(p)
#     else:
#         assert_never(p.ref)


BuildKeyT = TypeVar("BuildKeyT", HapBuildKey, Dip1BuildKey, Dip2BuildKey)

CompareKey = NewType("CompareKey", str)
OtherLevelKey = NewType("OtherLevelKey", str)
OtherStratKey = NewType("OtherStratKey", str)
HaplotypeName = NewType("HaplotypeName", str)


CHR_INDEX_PLACEHOLDER = "%i"
CHR_HAP_PLACEHOLDER = "%h"


def fmap_maybe(
    f: Callable[[X], Y],
    x: X | None,
) -> None | Y:
    return None if x is None else f(x)


def fmap_maybe_def(
    default: Y,
    f: Callable[[X], Y],
    x: X | None,
) -> Y:
    return default if x is None else f(x)


def from_maybe(default: X, x: X | None) -> X:
    return default if x is None else x


class BaseModel(BaseModel_):
    class Config:
        frozen = True
        extra = "forbid"


class Haplotype(Enum):
    HAP1: int = 0
    HAP2: int = 1

    @classmethod
    def from_name(cls, n: str) -> Self:
        try:
            return next(i for i in cls if i.name == n)
        except StopIteration:
            raise ValueError(f"could make haplotype from name '{n}'")

    @classmethod
    def with_haps(cls, f: "Callable[[Haplotype], X]") -> tuple[X, X]:
        x0 = f(cls.HAP1)
        x1 = f(cls.HAP2)
        return (x0, x1)

    @property
    def name(self) -> HaplotypeName:
        return HaplotypeName(f"hap{self.value + 1}")

    def from_either(self, left: X, right: X) -> X:
        if self is Haplotype.HAP1:
            return left
        elif self is Haplotype.HAP2:
            return right
        else:
            assert_never(self)


class Diploid_(BaseModel, Generic[X]):
    hap1: X
    hap2: X

    def from_either(self, hap: Haplotype) -> X:
        return hap.from_either(self.hap1, self.hap2)


def to_haplotype(s: str) -> Haplotype | None:
    try:
        return Haplotype.from_name(s)
    except ValueError:
        return None


def parse_refkey(s: str) -> tuple[str, Haplotype | None]:
    m = re.match("(.+)_(hap[12])", s)
    # ASSUME this will never fail due to the hap1/2 permitted match pattern
    return (s, None) if m is None else (m[1], Haplotype.from_name(m[2]))


def strip_refkey(s: str) -> str:
    return parse_refkey(s)[0]


def insert_suffix(p: Path, i: str) -> Path:
    ss = p.suffixes
    ss_ = [i, *ss]
    return p.parent / p.name.replace("".join(ss), "".join(ss_))


# dummy Identity type to make higher-order types more consistent
HaploidOnly = Union[X]
HaploidOrDiploid = Union[X, Diploid_[X]]


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
        try:
            return next(i for i in cls if i.chr_name == n)
        except StopIteration:
            raise ValueError(f"could make chr index from name '{n}'")

    @classmethod
    def from_name_unsafe(cls, n: str) -> Self:
        try:
            return cls.from_name(n)
        except ValueError as e:
            raise DesignError(e)

    def __init__(self, i: int) -> None:
        self.chr_name: str = "X" if i == 23 else ("Y" if i == 24 else str(i))

    def chr_name_full(self, p: "HapChrPattern") -> str | None:
        return p.to_chr_name(self)

    def chr_name_full_dip(self, p: "DipChrPattern", hap: Haplotype) -> str | None:
        return p.to_chr_name(self, hap)

    def to_internal_index(self, hap: Haplotype) -> bed.InternalChrIndex:
        return bed.InternalChrIndex(hap.value * 24 + self.value - 1)

    def choose_xy_unsafe(self, x_res: X, y_res: X) -> X:
        if self.value is self.CHRX:
            return x_res
        elif self.value is self.CHRY:
            return y_res
        else:
            raise DesignError(f"I am not an X or Y, I am a {self}")

    @property
    def xy_to_hap_unsafe(self) -> Haplotype:
        return self.choose_xy_unsafe(Haplotype.HAP2, Haplotype.HAP1)


class ChrPattern:
    def to_names(self, cs: set[ChrIndex]) -> list[str]:
        return NotImplemented


class HapChrPattern(BaseModel, ChrPattern):
    template: str = "chr%i"
    special: dict[ChrIndex, str] = {}
    exclusions: list[ChrIndex] = []

    @validator("template")
    def is_valid_template(cls, v: str) -> str:
        assert v.count(CHR_INDEX_PLACEHOLDER) == 1, "chr template must have '%i' in it"
        return v

    def to_chr_name(self, i: ChrIndex) -> str | None:
        if i in self.exclusions:
            return None
        elif i in self.special:
            return self.special[i]
        else:
            return self.template.replace(CHR_INDEX_PLACEHOLDER, i.chr_name)

    def to_pairs(
        self,
        cs: set[ChrIndex],
        h: Haplotype,
    ) -> list[tuple[bed.InternalChrIndex, str]]:
        return [
            (c.to_internal_index(h), n)
            for c in cs
            if (n := self.to_chr_name(c)) is not None
        ]

    def to_names(self, cs: set[ChrIndex]) -> list[str]:
        # NOTE: the haplotype argument is doing nothing since it is only
        # used to make the index which I remove before returning here
        return [x[1] for x in self.to_pairs(cs, Haplotype.HAP1)]

    def init_mapper(self, cs: set[ChrIndex], hap: Haplotype) -> bed.InitMapper:
        return {n: i for i, n in self.to_pairs(cs, hap)}

    def final_mapper(self, cs: set[ChrIndex], hap: Haplotype) -> bed.FinalMapper:
        return {i: n for i, n in self.to_pairs(cs, hap)}


class DipChrPattern(BaseModel, ChrPattern):
    template: str = "chr%i_%h"
    special: dict[ChrIndex, str] = {}
    hapnames: Diploid_[HaplotypeName] = Diploid_(
        hap1=HaplotypeName("PATERNAL"),
        hap2=HaplotypeName("MATERNAL"),
    )
    exclusions: Diploid_[list[ChrIndex]] = Diploid_(
        hap1=[ChrIndex.CHRX],
        hap2=[ChrIndex.CHRY],
    )

    @validator("template")
    def is_valid_template(cls, v: str) -> str:
        assert (
            v.count(CHR_INDEX_PLACEHOLDER) == 1 and v.count(CHR_HAP_PLACEHOLDER) == 1
        ), "chr template must have '%i' and '%h' in it"
        return v

    def to_chr_name(self, i: ChrIndex, h: Haplotype) -> str | None:
        exc = self.exclusions.from_either(h)
        name = self.hapnames.from_either(h)
        if i in exc:
            return None
        elif i in self.special:
            return self.special[i]
        else:
            return self.template.replace(CHR_INDEX_PLACEHOLDER, i.chr_name).replace(
                CHR_HAP_PLACEHOLDER, name
            )

    def to_pairs(self, cs: set[ChrIndex]) -> list[tuple[bed.InternalChrIndex, str]]:
        return [
            (c.to_internal_index(h), n)
            for c in cs
            for h in Haplotype
            if (n := self.to_chr_name(c, h)) is not None
        ]

    def to_names(self, cs: set[ChrIndex]) -> list[str]:
        return [x[1] for x in self.to_pairs(cs)]

    def init_mapper(self, cs: set[ChrIndex]) -> bed.InitMapper:
        return {n: i for i, n in self.to_pairs(cs)}

    def final_mapper(self, cs: set[ChrIndex]) -> bed.FinalMapper:
        return {i: n for i, n in self.to_pairs(cs)}

    def to_hap_pattern(self, hap: Haplotype) -> HapChrPattern:
        hs = self.hapnames.from_either(hap)
        xs = self.exclusions.from_either(hap)
        return HapChrPattern(
            template=self.template.replace(CHR_HAP_PLACEHOLDER, hs),
            special=self.special,
            exclusions=xs,
        )


AnyChrPatternT = TypeVar("AnyChrPatternT", HapChrPattern, DipChrPattern)


class HapChrSource(BaseModel, Generic[X]):
    """Specification for a haploid source file."""

    src: X
    chr_pattern: HapChrPattern = HapChrPattern()


class DipChrSource1(BaseModel, Generic[X]):
    """Specification for a combined diploid source file.

    The 'src' is assumed to have all chromosomes for both haplotypes in one
    file, which implies they are labeled so as to distinguish the haps. The
    pattern will match both the chromosome number and the haplotype within the
    chromosome name.
    """

    src: X
    chr_pattern: DipChrPattern = DipChrPattern()


class DipChrSource2(BaseModel, Generic[X]):
    """Specification for split diploid source files.

    Each source may or may not have each haplotype labeled; the identity of each
    haplotype in either source file is determined based on the configuration key
    under which it appears (hap1 or hap2) and the chromosome names for each are
    matched according to its corresponding entry in `chr_pattern`.
    """

    src: Diploid_[X]
    chr_pattern: Diploid_[HapChrPattern] = Diploid_(
        hap1=HapChrPattern(),
        hap2=HapChrPattern(),
    )


DipChrSource = Union[DipChrSource1[X], DipChrSource2[X]]
AnyChrSource = Union[HapChrSource[X], DipChrSource[X]]


class HashedSrc_(BaseModel):
    md5: str | None = None


class FileSrc_(HashedSrc_):
    """Base class for local src files."""

    filepath: FilePath


class BedFileSrc(FileSrc_):
    """Filepath for bedfile."""

    @validator("filepath")
    def is_gzip(cls, v: FilePath) -> FilePath:
        assert is_gzip(v), "not in gzip format"
        return v


class RefFileSrc(FileSrc_):
    """Filepath for reference."""

    @validator("filepath")
    def path_is_bgzip(cls, v: FilePath) -> FilePath:
        assert is_bgzip(v), "not in bgzip format"
        return v


class HttpSrc_(HashedSrc_):
    """Base class for downloaded src files."""

    url: HttpUrl


class BedHttpSrc(HttpSrc_):
    """Url for bed file"""

    pass


class RefHttpSrc(HttpSrc_):
    """Url for reference"""

    pass


RefSrc = RefFileSrc | RefHttpSrc

# TODO this is for more than just "bed files" (right now it basically means "a
# file that is either not zipped or gzipped but not bgzipped")
BedSrc = BedFileSrc | BedHttpSrc


class BedSrcAndKey(NamedTuple):
    src: BedSrc
    key: list[str]  # TODO newtype here?


# TODO clean this up with real polymorphism when mypy catches up with Haskell
# 98, see https://github.com/python/typing/issues/548
AnyBedT = TypeVar(
    "AnyBedT",
    HapChrSource[BedSrc],
    DipChrSource[BedSrc],
)
AnyBedT_ = TypeVar(
    "AnyBedT_",
    HapChrSource[BedSrc],
    DipChrSource1[BedSrc],
    DipChrSource2[BedSrc],
)


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
    OtherDifficult = "OtherDifficult"


@dataclass
class HapToHapChrConversion:
    fromPattern: HapChrPattern
    toPattern: HapChrPattern
    indices: set[ChrIndex]

    @property
    def init_mapper(self) -> bed.InitMapper:
        return self.fromPattern.init_mapper(self.indices, Haplotype.HAP1)

    @property
    def final_mapper(self) -> bed.FinalMapper:
        return self.toPattern.final_mapper(self.indices, Haplotype.HAP1)


@dataclass
class DipToDipChrConversion:
    fromPattern: DipChrPattern
    toPattern: DipChrPattern
    indices: set[ChrIndex]

    @property
    def init_mapper(self) -> bed.InitMapper:
        return self.fromPattern.init_mapper(self.indices)

    @property
    def final_mapper(self) -> bed.FinalMapper:
        return self.toPattern.final_mapper(self.indices)


@dataclass
class HapToDipChrConversion:
    fromPattern: Diploid_[HapChrPattern]
    toPattern: DipChrPattern
    indices: set[ChrIndex]

    @property
    def init_mapper(self) -> tuple[bed.InitMapper, bed.InitMapper]:
        p = self.fromPattern
        i = self.indices
        return (
            p.hap1.init_mapper(i, Haplotype.HAP1),
            p.hap2.init_mapper(i, Haplotype.HAP2),
        )

    @property
    def final_mapper(self) -> bed.FinalMapper:
        return self.toPattern.final_mapper(self.indices)


@dataclass
class DipToHapChrConversion:
    fromPattern: DipChrPattern
    toPattern: Diploid_[HapChrPattern]
    indices: set[ChrIndex]

    @property
    def init_mapper(self) -> tuple[bed.InitMapper, bed.SplitMapper]:
        im = self.fromPattern.init_mapper(self.indices)
        fm0 = self.toPattern.hap1.final_mapper(self.indices, Haplotype.HAP1)
        return (im, bed.make_split_mapper(im, fm0))

    @property
    def final_mapper(self) -> tuple[bed.FinalMapper, bed.FinalMapper]:
        p = self.toPattern
        i = self.indices
        return (
            p.hap1.final_mapper(i, Haplotype.HAP1),
            p.hap2.final_mapper(i, Haplotype.HAP2),
        )


# For instances where we simply need to sort all the chromosomes and they are
# already named appropriately
# def fullset_hap_conv(p: HapChrPattern) -> HapToHapChrConversion:
#     return HapToHapChrConversion(p, p, set([i for i in ChrIndex]))


# def fullset_dip1_conv(p: DipChrPattern) -> DipToDipChrConversion:
#     return DipToDipChrConversion(
#         p, p, set([(h, i) for i in ChrIndex for h in Haplotype])
#     )


class Paths(BaseModel):
    """Local build paths for snakemake."""

    resources: Path = Path("resources")
    results: Path = Path("results")


class Tools(BaseModel):
    """Urls for tools to download/build/use in the pipeline."""

    repseq: HttpUrl = "https://github.com/ndwarshuis/repseq/archive/refs/tags/v1.1.0.tar.gz"  # type: ignore
    gemlib: HttpUrl = "https://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%203/GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2/download"  # type: ignore


# TODO non-negative ints which cannot equal each other
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
    def columns(self) -> dict[int, Type[int | str]]:
        return {self.chr: str, self.start: int, self.end: int}


class BedFileParams(BaseModel):
    """Parameters decribing how to parse a bed-like file.

    Members:
    chr_pattern - the pattern on the chromosomes; must include the special
      directive '%i' which will denote a standardized name (eg 1, 2 ...X, Y);
      the pattern is assumed to match the whole chromosome name.
    bed_cols - the columns for the bed coordinates
    skip_lines - how many input lines to skip
    sep - column separator regexp (for "beds" with spaces instead of tabs)
    """

    # chr_pattern: HapChrPattern = HapChrPattern()
    bed_cols: BedColumns = BedColumns()
    skip_lines: NonNegativeInt = 0
    sep: str = "\t"


class BedFile(BaseModel, Generic[X]):
    """Inport specs for a bed-like file."""

    data: X
    params: BedFileParams = BedFileParams()

    # @property
    # def src_list(self) -> list[BedSrc]:
    #     return diploid_to_list(self.data)

    def read(self, path: Path) -> pd.DataFrame:
        return self._read(path, [])

    def _read(self, path: Path, more: list[int] = []) -> pd.DataFrame:
        p = self.params
        return bed.read_bed(path, p.bed_cols.columns, p.skip_lines, p.sep, more)

    # def read_filter_sort_bed(
    #     self,
    #     ipath: Path,
    #     opath: Path,
    #     conv:
    #     more: list[int] = [],
    # ) -> None:
    #     """Read a haploid bed file, sort it, and write it in bgzip format."""
    #     # conv = sconf.haploid_stratifications.to_build_data_unsafe(
    #     #     rk, bk
    #     # ).chr_conversion(pat)
    #     im = None
    #     fm = None
    #     df = self._read(ipath, more)
    #     df_ = bed.filter_sort_bed(im, fm, df)
    #     bed.write_bed(opath, df_)


HapBedSrc = HapChrSource[BedSrc]
DipBedSrc = DipChrSource[BedSrc]
Dip1BedSrc = DipChrSource1[BedSrc]
Dip2BedSrc = DipChrSource2[BedSrc]

HapBedFile = BedFile[HapBedSrc]
DipBedFile = BedFile[DipBedSrc]
Dip1BedFile = BedFile[Dip1BedSrc]
Dip2BedFile = BedFile[Dip2BedSrc]
AnyBedFileT = BedFile[AnyBedT]


# TODO mypy for some reason doesn't understand how to narrow a
# Something[Union[X, Y]] to a Something[X] using 'isinstance'
def is_dip1_bed(x: BedFile[DipChrSource[X]]) -> TypeGuard[BedFile[DipChrSource1[X]]]:
    return isinstance(x.data, DipChrSource1)


def is_dip2_bed(x: BedFile[DipChrSource[X]]) -> TypeGuard[BedFile[DipChrSource2[X]]]:
    return isinstance(x.data, DipChrSource2)


# class VCFFile(BaseModel):
#     """Inport specs for a vcf file."""

#     src: BedSrc
#     chr_pattern: HapChrPattern = HapChrPattern()


class VCFFile(BedFile[X]):
    """Inport specs for a vcf file."""

    pass

    # src: BedSrc
    # chr_pattern: HapChrPattern = HapChrPattern()


class RMSKFile(BedFile[X]):
    """Input file for repeat masker stratification."""

    # TODO type narrowing won't work without this redfinition
    data: X
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
        return super()._read(path, [self.class_col])


class SatFile(BedFile[X]):
    """Configuration for a satellites file."""

    sat_col: NonNegativeInt

    def read(self, path: Path) -> pd.DataFrame:
        return super()._read(path, [self.sat_col])


# def either_diploid(left: Callable) -> list[X]:
#     s = x.src
#     if isinstance(s, Diploid_):
#         return [s.hap1, s.hap2]
#     else:
#         return [s]


def diploid_to_list(x: AnyChrSource[X]) -> list[X]:
    s = x.src
    if isinstance(s, Diploid_):
        return [s.hap1, s.hap2]
    else:
        return [s]


def map_diploid(f: Callable[[X], Y], x: AnyChrSource[X]) -> list[Y]:
    return [f(y) for y in diploid_to_list(x)]


class LowComplexity(BaseModel, Generic[X]):
    """Configuration for low complexity stratification."""

    rmsk: RMSKFile[X] | None
    simreps: BedFile[X] | None
    satellites: SatFile[X] | None


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
class XYFeatures(BaseModel):
    """Configuration for XY features stratifications."""

    x_bed: XYFile
    y_bed: XYFile
    ampliconic: bool
    xtr: bool


class XYPar(BaseModel):
    """Regions for the PARs on the X/Y chromosomes."""

    start: tuple[NonNegativeInt, NonNegativeInt]
    end: tuple[NonNegativeInt, NonNegativeInt]

    @validator("start", "end")
    def positive_region(cls, v: tuple[int, int]) -> tuple[int, int]:
        assert v[1] > v[0], "End must be greater than start"
        return v

    def fmt(self, i: ChrIndex, pattern: HapChrPattern) -> str:
        # TODO this smells like something I'll be doing alot
        c = i.chr_name_full(pattern)
        return "\n".join(
            [
                f"{c}\t{self.start[0]}\t{self.start[1]}",
                f"{c}\t{self.end[0]}\t{self.end[1]}",
            ]
        )


class XY(BaseModel):
    """Configuration for the XY stratification."""

    features: XYFeatures | None
    x_par: XYPar | None
    y_par: XYPar | None

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


class SegDups(BaseModel, Generic[X]):
    """Configuration for Segdup stratifications."""

    superdups: BedFile[X] | None


class LowMapParams(BaseModel):
    """Parameters for a single mappability bed file."""

    length: NonNegativeInt
    mismatches: NonNegativeInt
    indels: NonNegativeInt


GCBound = tuple[Percent, bool]


class GCParams(BaseModel):
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

    @validator("high")
    def has_balanced_ranges(
        cls,
        high: list[GCBound],
        values: dict[str, Any],
    ) -> list[tuple[Percent, bool]]:
        try:
            low = cast(list[GCBound], values["low"])
            assert len([x for x in low if x[1]]) == len(
                [x for x in high if x[1]]
            ), "GC low/high must have same number of range boundaries"
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

    @property
    def low_bounds(self) -> tuple[int, list[int]]:
        bounds = self.low_fractions
        return (bounds[0], bounds[1:])

    @property
    def high_bounds(self) -> tuple[int, list[int]]:
        bounds = self.high_fractions
        return (bounds[-1], bounds[:-1])


class HaploidInclude(BaseModel):
    """Flags to control which stratification levels are included."""

    low_complexity: bool = True
    xy: bool = True
    functional: bool = True
    segdups: bool = True
    union: bool = True
    telomeres: bool = True
    # TODO also add KIR and MHC since these should be derivable from refseq
    vdj: bool = True
    mappability: set[LowMapParams] = {
        LowMapParams(length=250, mismatches=0, indels=0),
        LowMapParams(length=100, mismatches=2, indels=1),
    }
    gc: GCParams | None = GCParams()


class DiploidInclude(HaploidInclude):
    """Flags to control which stratification levels are included."""

    hets: bool = True


IncludeT = TypeVar("IncludeT", HaploidInclude, DiploidInclude)


class OtherBedFile(BedFile[AnyBedT]):
    remove_gaps: bool = False


OtherStrats = dict[OtherLevelKey, dict[OtherStratKey, OtherBedFile[AnyBedT]]]


class Bench(BaseModel, Generic[AnyBedT, AnyBedT_]):
    """Configuration for benchmark to use when validating stratifications."""

    bench_vcf: VCFFile[AnyBedT_]
    query_vcf: VCFFile[AnyBedT_]
    bench_bed: BedFile[AnyBedT]


class BuildCompare(BaseModel):
    """Configuration for comparing generated strats to previous versions."""

    other: CompareKey
    path_mapper: dict[Path, Path] = {}
    replacements: list[tuple[str, str]] = []
    ignore_other: list[str] = []
    ignore_generated: list[str] = []


class Build_(BaseModel, Generic[AnyBedT, AnyBedT_, IncludeT]):
    chr_filter: set[ChrIndex]
    comparison: BuildCompare | None = None
    bench: Bench[AnyBedT, AnyBedT_] | None = None
    other_strats: OtherStrats[AnyBedT] = {}
    # TODO somehow get this to default to something, either with a validator
    # (maybe), hacky class method, or using a higher-order type var (which would
    # let me make the include type generic and definable from higher in the
    # stack)
    include: IncludeT


# # need class overrides here to get default values for include, which are
# # different b/t hap and dip cases
# class HaploidBuild(Build_[AnyBedT, AnyBedT_]):
#     """Spec for a stratification build."""

#     include: HaploidInclude = IncludeHaploid()


# # # AnyDipBedSrcT = TypeVar("AnyDipBedSrcT", DipChrSource1[BedSrc], DipChrSource2[BedSrc])


# class DiploidBuild(Build_[AnyBedT, AnyBedT_]):
#     """Spec for a stratification build."""

#     include: DiploidInclude = DiploidInclude()


# BenchT = TypeVar(
#     "BenchT",
#     Bench[HapChrSource[BedSrc]],
#     Bench[DipChrSource1[BedSrc]],
#     Bench[DipChrSource2[BedSrc]],
# )

# RefFile = AnyChrSource[RefSrc]


# class RefFile(BaseModel, Generic[X, Y]):
#     """Specification for a reference file."""

#     src: X
#     chr_pattern: Y


# class RefFileHaploid(BaseModel):
#     """Specification for a reference file."""

#     src: RefSrc
#     chr_pattern: HapChrPattern = HapChrPattern()


# RefFileDiploid = Diploid_[RefFileHaploid]


class Functional(BaseModel, Generic[X]):
    """Configuration for Functional stratifications."""

    ftbl_src: X
    gff_src: X


class StratInputs_(BaseModel, Generic[AnyBedT]):
    gap: BedFile[AnyBedT] | None
    low_complexity: LowComplexity[AnyBedT]
    xy: XY
    mappability: Mappability | None
    segdups: SegDups[AnyBedT]
    # TODO this type is wrong (specifically is has too much information)
    functional: Functional[AnyBedT] | None

    def _to_bed_src(
        self,
        f: Callable[[Self], BedFile[AnyBedT] | None],
    ) -> list[BedSrc]:
        x = f(self)
        if x is None:
            return []
        else:
            return diploid_to_list(x.data)

    def _to_hap_bed_src(
        self,
        h: Haplotype,
        f: Callable[[Self], BedFile[AnyBedT] | None],
    ) -> list[BedSrc]:
        x = f(self)
        if x is None:
            return []
        else:
            return diploid_to_list(x.data)

    @property
    def gap_src(self) -> list[BedSrc]:
        return self._to_bed_src(lambda x: x.gap)

    @property
    def rmsk_src(self) -> list[BedSrc]:
        return self._to_bed_src(lambda x: x.low_complexity.rmsk)

    @property
    def simreps_src(self) -> list[BedSrc]:
        return self._to_bed_src(lambda x: x.low_complexity.simreps)

    @property
    def satellites_src(self) -> list[BedSrc]:
        return self._to_bed_src(lambda x: x.low_complexity.satellites)

    @property
    def superdups_src(self) -> list[BedSrc]:
        return self._to_bed_src(lambda x: x.segdups.superdups)

    @property
    def ftbl_src(self) -> list[BedSrc]:
        return fmap_maybe_def(
            [], lambda x: diploid_to_list(x.ftbl_src), self.functional
        )

    @property
    def gff_src(self) -> list[BedSrc]:
        return fmap_maybe_def([], lambda x: diploid_to_list(x.gff_src), self.functional)

    @property
    def xy_features_unsafe(self) -> XYFeatures:
        f = self.xy.features
        if f is None:
            raise DesignError("XY features does not exist")
        return f

    def xy_feature_bed_unsafe(self, i: ChrIndex) -> XYFile:
        f = self.xy_features_unsafe
        return i.choose_xy_unsafe(f.x_bed, f.y_bed)


class StratInputToBed(Protocol):
    A = TypeVar("A", HapChrSource[BedSrc], DipChrSource[BedSrc])

    def __call__(self, __x: StratInputs_[A]) -> BedFile[A] | None:
        pass


HaploidStratInputs = StratInputs_[HapChrSource[BedSrc]]
DiploidStratInputs = StratInputs_[DipChrSource[BedSrc]]
AnyStratInputs = HaploidStratInputs | DiploidStratInputs

StratInputT = TypeVar("StratInputT", HaploidStratInputs, DiploidStratInputs)


# BuildT = TypeVar(
#     "BuildT",
#     HaploidBuild,
#     DiploidBuild[DipChrSource1[BedSrc]],
#     DiploidBuild[DipChrSource2[BedSrc]],
# )
RefSourceT = TypeVar(
    "RefSourceT",
    HapChrSource[RefSrc],
    DipChrSource1[RefSrc],
    DipChrSource2[RefSrc],
)


# HaploidRefData = RefData_[HapChrSource[RefSrc], HaploidStratInputs]
# Diploid1RefData = RefData_[DipChrSource1[RefSrc], DiploidStratInputs]
# Diploid2RefData = RefData_[DipChrSource2[RefSrc], DiploidStratInputs]
# AnyRefData = HaploidRefData | Diploid1RefData | Diploid2RefData


@dataclass
class BuildData_(Generic[RefKeyT_, BuildKeyT, RefSourceT, AnyBedT, AnyBedT_, IncludeT]):
    refkey: RefKey_[RefKeyT_]
    buildkey: BuildKeyT
    ref: RefSourceT
    strat_inputs: StratInputs_[AnyBedT]
    build: Build_[AnyBedT, AnyBedT_, IncludeT]

    # @property
    # def include(self) -> HaploidInclude:
    #     return self.build.include

    @property
    def mappability_patterns(self) -> list[str]:
        return fmap_maybe_def(
            [],
            lambda m: m.unplaced_chr_patterns,
            self.strat_inputs.mappability,
        )

    @property
    def chr_indices(self) -> set[ChrIndex]:
        cs = self.build.chr_filter
        return set([x for x in ChrIndex]) if len(cs) == 0 else cs

    # @property
    # def bench_vcf_src(self) -> BedSrc | None:
    #     return fmap_maybe(lambda x: x.bench_vcf.src, self.build.bench)

    # @property
    # def bench_bed_src(self) -> BedSrc | None:
    #     return fmap_maybe(lambda x: x.bench_bed.data.src, self.build.bench)

    # @property
    # def query_vcf_src(self) -> BedSrc | None:
    #     return fmap_maybe(lambda x: x.query_vcf.src, self.build.bench)

    @property
    def mappability_params(
        self,
    ) -> tuple[list[int], list[int], list[int]]:
        ms = self.build.include.mappability
        xs = [(m.length, m.mismatches, m.indels) for m in ms]
        # TODO mypy doesn't like unzip here for some reason :(
        return ([x[0] for x in xs], [x[1] for x in xs], [x[2] for x in xs])


@dataclass(frozen=True)
class RefData_(Generic[RefKeyT_, RefSourceT, AnyBedT, AnyBedT_, BuildKeyT, IncludeT]):
    refkey: RefKey_[RefKeyT_]
    ref: RefSourceT
    strat_inputs: StratInputs_[AnyBedT]
    builds: dict[BuildKeyT, Build_[AnyBedT, AnyBedT_, IncludeT]]

    def to_build_data_unsafe(
        self,
        bk: BuildKeyT,
    ) -> BuildData_[RefKeyT_, BuildKeyT, RefSourceT, AnyBedT, AnyBedT_, IncludeT]:
        bd = self.to_build_data(bk)
        if bd is None:
            raise DesignError(f"Could not create build data from key '{bk}'")
        return bd

    def to_build_data(
        self,
        bk: BuildKeyT,
    ) -> (
        BuildData_[RefKeyT_, BuildKeyT, RefSourceT, AnyBedT, AnyBedT_, IncludeT] | None
    ):
        try:
            return BuildData_(
                self.refkey,
                bk,
                self.ref,
                self.strat_inputs,
                self.builds[bk],
            )
        except KeyError:
            return None

    # def to_build_data_unsafe(
    #     self,
    #     rk: RefKeyT_,
    #     bk: BuildKeyT,
    # ) -> BuildData_[RefSourceT, AnyBedT, AnyBedT_, IncludeT]:
    #     return self.to_ref_data_unsafe(rk).to_build_data_unsafe(bk)

    # def to_build_data(
    #     self,
    #     rk: RefKeyT_,
    #     bk: BuildKeyT,
    # ) -> BuildData_[RefSourceT, AnyBedT, AnyBedT_, IncludeT] | None:
    #     try:
    #         return self.to_build_data_unsafe(rk, bk)
    #     except DesignError:
    #         return None

    @property
    def mappability_patterns(self) -> list[str]:
        return fmap_maybe_def(
            [],
            lambda m: m.unplaced_chr_patterns,
            self.strat_inputs.mappability,
        )


class BuildDataToBed(Protocol):
    A = TypeVar("A", HapBedSrc, DipBedSrc)

    def __call__(
        self, __x: BuildData_[RefKeyT_, BuildKeyT, RefSourceT, A, AnyBedT_, IncludeT]
    ) -> BedFile[A] | None:
        pass


class Stratification(
    BaseModel, Generic[RefSourceT, AnyBedT, AnyBedT_, BuildKeyT, IncludeT]
):
    """Configuration for stratifications for a given reference."""

    ref: RefSourceT
    strat_inputs: StratInputs_[AnyBedT]
    builds: dict[BuildKeyT, Build_[AnyBedT, AnyBedT_, IncludeT]]


# TODO wrap this in the instance classes for the chromosome converters
# def filter_sort_bed_conv(conv: HapToHapChrConversion, df: pd.DataFrame) -> pd.DataFrame:
#     """Filter and sort a bed file from a dataframe."""
#     from_map = conv.init_mapper
#     to_map = conv.final_mapper
#     return bed.filter_sort_bed(from_map, to_map, df)


@dataclass
class HaploidBuildData(
    BuildData_[
        HapRefKey_,
        HapBuildKey,
        HapChrSource[RefSrc],
        HapBedSrc,
        HapBedSrc,
        HaploidInclude,
    ]
):
    @property
    def final_mapper(self) -> bed.FinalMapper:
        return self.ref.chr_pattern.final_mapper(self.chr_indices, Haplotype.HAP1)

    @property
    def ref_chr_conversion(self) -> HapToHapChrConversion:
        p = self.ref.chr_pattern
        return HapToHapChrConversion(p, p, self.chr_indices)

    def chr_conversion(self, fromChr: HapChrPattern) -> HapToHapChrConversion:
        return HapToHapChrConversion(
            fromChr,
            self.ref.chr_pattern,
            self.chr_indices,
        )

    def read_filter_sort_hap_bed(self, bf: HapBedFile, ipath: Path) -> pd.DataFrame:
        """Read a haploid bed file, sort it, and write it in bgzip format."""
        conv = self.chr_conversion(bf.data.chr_pattern)
        df = bf.read(ipath)
        return bed.filter_sort_bed(conv.init_mapper, conv.final_mapper, df)

    def read_write_filter_sort_hap_bed(
        self,
        bf: HapBedFile,
        ipath: Path,
        opath: Path,
        g: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
    ) -> None:
        """Read a haploid bed file, sort it, and write it in bgzip format."""
        df = self.read_filter_sort_hap_bed(bf, ipath)
        bed.write_bed(opath, g(df))


@dataclass
class Diploid1BuildData(
    BuildData_[
        Dip1RefKey_,
        Dip1BuildKey,
        DipChrSource1[RefSrc],
        DipBedSrc,
        Dip1BedSrc,
        DiploidInclude,
    ]
):
    @property
    def final_mapper(self) -> bed.FinalMapper:
        return self.ref.chr_pattern.final_mapper(self.chr_indices)

    @property
    def ref_chr_conversion(self) -> DipToDipChrConversion:
        p = self.ref.chr_pattern
        return DipToDipChrConversion(p, p, self.chr_indices)

    def hap_chr_conversion(
        self,
        fromChr: Diploid_[HapChrPattern],
    ) -> HapToDipChrConversion:
        return HapToDipChrConversion(
            fromChr,
            self.ref.chr_pattern,
            self.chr_indices,
        )

    def dip_chr_conversion(
        self,
        fromChr: DipChrPattern,
    ) -> DipToDipChrConversion:
        return DipToDipChrConversion(
            fromChr,
            self.ref.chr_pattern,
            self.chr_indices,
        )

    def read_filter_sort_hap_bed(
        self,
        bf: Dip2BedFile,
        ipath: tuple[Path, Path],
    ) -> pd.DataFrame:
        """Read two haploid bed files, combine and sort them as diploid, and write
        it in bgzip format.
        """

        def go(b: Dip2BedFile, i: Path, imap: bed.InitMapper) -> pd.DataFrame:
            df = b.read(i)
            return bed.filter_sort_bed(imap, fmap, df)

        conv = self.hap_chr_conversion(bf.data.chr_pattern)
        imap1, imap2 = conv.init_mapper
        fmap = conv.final_mapper

        return pd.concat(
            [
                go(bf, *x)
                for x in [
                    (ipath[0], imap1),
                    (ipath[1], imap2),
                ]
            ]
        )

    def read_write_filter_sort_hap_bed(
        self,
        bf: Dip2BedFile,
        ipath: tuple[Path, Path],
        opath: Path,
        g: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
    ) -> None:
        """Read a haploid bed file, sort it, and write it in bgzip format."""
        df = self.read_filter_sort_hap_bed(bf, ipath)
        bed.write_bed(opath, g(df))

    def read_filter_sort_dip_bed(
        self,
        bf: Dip1BedFile,
        ipath: Path,
    ) -> pd.DataFrame:
        """Read a diploid bed file, sort it, and write it in bgzip format."""
        conv = self.dip_chr_conversion(bf.data.chr_pattern)
        df = bf.read(ipath)
        return bed.filter_sort_bed(conv.init_mapper, conv.final_mapper, df)

    def read_write_filter_sort_dip_bed(
        self,
        bf: Dip1BedFile,
        ipath: Path,
        opath: Path,
        g: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
    ) -> None:
        """Read a haploid bed file, sort it, and write it in bgzip format."""
        df = self.read_filter_sort_dip_bed(bf, ipath)
        bed.write_bed(opath, g(df))


@dataclass
class Diploid2BuildData(
    BuildData_[
        Dip1RefKey_,
        Dip2BuildKey,
        DipChrSource2[RefSrc],
        DipBedSrc,
        Dip2BedSrc,
        DiploidInclude,
    ]
):
    @property
    def final_mapper(self) -> tuple[bed.FinalMapper, bed.FinalMapper]:
        p = self.ref.chr_pattern
        i = self.chr_indices
        return (
            p.hap1.final_mapper(i, Haplotype.HAP1),
            p.hap2.final_mapper(i, Haplotype.HAP2),
        )

    @property
    def ref_chr_conversion(self) -> tuple[HapToHapChrConversion, HapToHapChrConversion]:
        def go(h: HapChrPattern) -> HapToHapChrConversion:
            return HapToHapChrConversion(h, h, self.chr_indices)

        p = self.ref.chr_pattern
        return (go(p.hap1), go(p.hap2))

    def hap_chr_conversion(
        self,
        fromChr: Diploid_[HapChrPattern],
    ) -> tuple[HapToHapChrConversion, HapToHapChrConversion]:
        toChr = self.ref.chr_pattern
        cis = self.chr_indices
        return (
            HapToHapChrConversion(fromChr.hap1, toChr.hap1, cis),
            HapToHapChrConversion(fromChr.hap2, toChr.hap2, cis),
        )

    def dip_chr_conversion(
        self,
        fromChr: DipChrPattern,
    ) -> DipToHapChrConversion:
        return DipToHapChrConversion(
            fromChr,
            self.ref.chr_pattern,
            self.chr_indices,
        )

    def read_filter_sort_hap_bed(
        self,
        bf: Dip2BedFile,
        ipath: Path,
        hap: Haplotype,
    ) -> pd.DataFrame:
        conv = self.hap_chr_conversion(bf.data.chr_pattern)
        df = bf.read(ipath)
        conv_ = hap.from_either(conv[0], conv[1])
        return bed.filter_sort_bed(conv_.init_mapper, conv_.final_mapper, df)

    def read_write_filter_sort_hap_bed(
        self,
        bf: Dip2BedFile,
        ipath: Path,
        opath: Path,
        hap: Haplotype,
        g: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
    ) -> None:
        df = self.read_filter_sort_hap_bed(bf, ipath, hap)
        bed.write_bed(opath, g(df))

    def read_filter_sort_dip_bed(
        self,
        bf: Dip1BedFile,
        ipath: Path,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        conv = self.dip_chr_conversion(bf.data.chr_pattern)
        imap, splitter = conv.init_mapper
        fmap0, fmap1 = conv.final_mapper

        def go(df: pd.DataFrame, fmap: bed.FinalMapper) -> pd.DataFrame:
            return bed.filter_sort_bed(imap, fmap, df)

        df = bf.read(ipath)
        df0, df1 = bed.split_bed(splitter, df)
        return (go(df0, fmap0), go(df1, fmap1))

    def read_write_filter_sort_dip_bed(
        self,
        bf: Dip1BedFile,
        ipath: Path,
        opath: tuple[Path, Path],
        g0: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
        g1: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
    ) -> None:
        """Read a haploid bed file, sort it, and write it in bgzip format."""
        df0, df1 = self.read_filter_sort_dip_bed(bf, ipath)
        bed.write_bed(opath[0], g0(df0))
        bed.write_bed(opath[1], g1(df1))


AnyBuildData = HaploidBuildData | Diploid1BuildData | Diploid2BuildData


# def apply_build_data(
#     bd: AnyBuildData,
#     f: Callable[[HaploidBuildData], X],
#     g: Callable[[Diploid1BuildData], X],
#     h: Callable[[Diploid2BuildData], X],
# ) -> X:
#     if isinstance(bd, HaploidBuildData):
#         return f(bd)
#     elif isinstance(bd, Diploid1BuildData):
#         return g(bd)
#     elif isinstance(bd, Diploid2BuildData):
#         return h(bd)
#     else:
#         assert_never(bd)


class StratDict_(
    dict[RefKeyT_, Stratification[RefSourceT, AnyBedT, AnyBedT_, BuildKeyT, IncludeT]],
    Generic[RefKeyT_, RefSourceT, AnyBedT, AnyBedT_, BuildKeyT, IncludeT],
):
    def to_ref_data_unsafe(
        self,
        rk: RefKeyT_,
    ) -> RefData_[RefKeyT_, RefSourceT, AnyBedT, AnyBedT_, BuildKeyT, IncludeT]:
        rd = self.to_ref_data(rk)
        if rd is None:
            raise DesignError(f"Could not get ref data for key '{rk}'")
        return rd

    def to_ref_data(
        self,
        rk: RefKeyT_,
    ) -> RefData_[RefKeyT_, RefSourceT, AnyBedT, AnyBedT_, BuildKeyT, IncludeT] | None:
        try:
            s = self[rk]
            return RefData_(RefKey_(rk), s.ref, s.strat_inputs, s.builds)
        except KeyError:
            return None


# class HaploidStratDict(
#     StratDict_[
#         HapRefKey,
#         HapChrSource[RefSrc],
#         HaploidStratInputs,
#         HapBuildKey,
#         HaploidBuild,
#     ]
# ):
#     def to_build_data_unsafe(
#         self,
#         rk: HapRefKey,
#         bk: HapBuildKey,
#     ) -> HaploidBuildData:
#         return HaploidBuildData(*astuple(super().to_build_data_unsafe(rk, bk)))

#     def to_build_data(
#         self,
#         rk: HapRefKey,
#         bk: HapBuildKey,
#     ) -> HaploidBuildData | None:
#         try:
#             return self.to_build_data_unsafe(rk, bk)
#         except KeyError:
#             return None


# class Diploid1StratDict(
#     StratDict_[
#         Dip1RefKey,
#         DipChrSource1[RefSrc],
#         DiploidStratInputs,
#         Dip1BuildKey,
#         DiploidBuild[DipChrSource1[BedSrc]],
#     ]
# ):
#     def to_build_data_unsafe(
#         self,
#         rk: Dip1RefKey,
#         bk: Dip1BuildKey,
#     ) -> Diploid1BuildData:
#         return Diploid1BuildData(*astuple(super().to_build_data_unsafe(rk, bk)))

#     def to_build_data(
#         self,
#         rk: Dip1RefKey,
#         bk: Dip1BuildKey,
#     ) -> Diploid1BuildData | None:
#         try:
#             return self.to_build_data_unsafe(rk, bk)
#         except KeyError:
#             return None


# class Diploid2StratDict(
#     StratDict_[
#         Dip2RefKey,
#         DipChrSource2[RefSrc],
#         DiploidStratInputs,
#         Dip2BuildKey,
#         DiploidBuild[DipChrSource2[BedSrc]],
#     ]
# ):
#     def to_build_data_unsafe(
#         self,
#         rk: Dip2RefKey,
#         bk: Dip2BuildKey,
#     ) -> Diploid2BuildData:
#         return Diploid2BuildData(*astuple(super().to_build_data_unsafe(rk, bk)))

#     def to_build_data(
#         self,
#         rk: Dip2RefKey,
#         bk: Dip2BuildKey,
#     ) -> Diploid2BuildData | None:
#         try:
#             return self.to_build_data_unsafe(rk, bk)
#         except KeyError:
#             return None

HaploidStratification = Stratification[
    HapChrSource[RefSrc],
    HapBedSrc,
    HapBedSrc,
    HapBuildKey,
    HaploidInclude,
]
Dip1Strat = Stratification[
    DipChrSource1[RefSrc],
    DipBedSrc,
    Dip1BedSrc,
    Dip1BuildKey,
    DiploidInclude,
]
Dip2Strat = Stratification[
    DipChrSource2[RefSrc],
    DipBedSrc,
    Dip2BedSrc,
    Dip2BuildKey,
    DiploidInclude,
]
AnyStratification = HaploidStratification | Dip1Strat | Dip2Strat

HaploidStratDict = StratDict_[
    HapRefKey_,
    HapChrSource[RefSrc],
    HapBedSrc,
    HapBedSrc,
    HapBuildKey,
    HaploidInclude,
]

Diploid1StratDict = StratDict_[
    Dip1RefKey_,
    DipChrSource1[RefSrc],
    DipBedSrc,
    Dip1BedSrc,
    Dip1BuildKey,
    DiploidInclude,
]

Diploid2StratDict = StratDict_[
    Dip2RefKey_,
    DipChrSource2[RefSrc],
    DipBedSrc,
    Dip2BedSrc,
    Dip2BuildKey,
    DiploidInclude,
]


# TODO add validator to ensure none of the keys in the strat/build dicts overlap
class GiabStrats(BaseModel):
    """Top level stratification object."""

    other_levels: list[OtherLevelKey] = [
        OtherLevelKey("Ancestry"),
        OtherLevelKey("FunctionalTechnicallyDifficult"),
        OtherLevelKey("GenomeSpecific"),
        OtherLevelKey("OtherDifficult"),
    ]
    paths: Paths = Paths()
    tools: Tools = Tools()
    haploid_stratifications: HaploidStratDict = HaploidStratDict()
    diploid1_stratifications: Diploid1StratDict = Diploid1StratDict()
    diploid2_stratifications: Diploid2StratDict = Diploid2StratDict()
    comparison_strats: dict[CompareKey, HttpUrl] = {}
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

    @validator("stratifications", each_item=True)
    def builds_have_valid_existing(
        cls,
        v: HaploidStratification,
        values: dict[str, Any],
    ) -> HaploidStratification:
        try:
            levels = cast(list[OtherLevelKey], values["other_levels"])
            bad = [
                f"level='{lk}'; build='{bk}'"
                for bk, b in v.builds.items()
                for lk in b.other_strats
                if lk not in levels
            ]
            if len(bad) > 0:
                assert (
                    False
                ), f"builds referencing invalid strat categories: {', '.join(bad)}"
        except KeyError:
            pass
        return v

    @validator("haploid_stratifications", each_item=True)
    def builds_have_valid_old_version(
        cls,
        v: HaploidStratification,
        values: dict[str, Any],
    ) -> HaploidStratification:
        try:
            prev = cast(dict[CompareKey, HttpUrl], values["previous"])
            bad = [
                f"version='{pk}'; build='{bk}'"
                for bk, b in v.builds.items()
                if b.comparison is not None
                for pk in b.comparison.other
                if pk not in prev
            ]
            if len(bad) > 0:
                assert (
                    False
                ), f"builds referencing invalid previous version keys: {', '.join(bad)}"
        except KeyError:
            pass
        return v

    # hack to make rmd scripts work with this (note this will totally kill
    # the config as it passes into an rmd script)
    def items(self) -> Any:
        return {}.items()

    # file paths

    @property
    def resources_dir(self) -> Path:
        return self.paths.resources

    @property
    def tools_src_dir(self) -> Path:
        return self.resources_dir / "tools"

    @property
    def ref_src_dir(self) -> Path:
        return self.paths.resources / "{ref_src_key}"

    @property
    def results_dir(self) -> Path:
        return self.paths.results

    @property
    def tools_make_dir(self) -> Path:
        return self.results_dir / "tools" / "make"

    @property
    def tools_bin_dir(self) -> Path:
        return self.results_dir / "tools" / "bin"

    @property
    def final_root_dir(self) -> Path:
        return self.results_dir / "final"

    @property
    def final_build_dir(self) -> Path:
        return self.final_root_dir / "{ref_key}@{build_key}"

    @property
    def intermediate_root_dir(self) -> Path:
        return self.results_dir / "intermediates"

    @property
    def intermediate_build_dir(self) -> Path:
        return self.intermediate_root_dir / "{ref_key}@{build_key}"

    @property
    def intermediate_ref_dir(self) -> Path:
        return self.intermediate_root_dir / "ref"

    @property
    def log_root_dir(self) -> Path:
        return self.results_dir / "log"

    @property
    def bench_root_dir(self) -> Path:
        return self.results_dir / "bench"

    @property
    def log_src_dir(self) -> Path:
        return self.log_root_dir / "resources" / "{ref_src_key}"

    @property
    def log_results_dir(self) -> Path:
        return self.log_root_dir / "results"

    @property
    def log_build_dir(self) -> Path:
        return self.log_results_dir / "builds" / "{ref_key}@{build_key}"

    @property
    def bench_build_dir(self) -> Path:
        return self.bench_root_dir / "{ref_key}@{build_key}"

    def build_final_strat_path(self, level: str, name: str) -> Path:
        return self.final_build_dir / level / f"{{ref_key}}_{name}.bed.gz"

    def build_strat_path(self, level: CoreLevel, name: str) -> Path:
        return self.build_final_strat_path(level.value, name)

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

    def refkey_to_ref_src(self, rk: str) -> RefSrc | None:
        return self.with_ref_data_unsafe(
            rk,
            lambda rd: rd.ref.src,
            lambda rd: rd.ref.src,
            lambda hap, rd: hap.from_either(rd.ref.src.hap1, rd.ref.src.hap2),
        )

    def refkey_to_bed_refsrckeys(self, f: StratInputToBed, rk: str) -> list[str]:
        return self.with_ref_data_and_bed(
            rk,
            f,
            lambda rk, _, __: [rk],
            lambda rk, _, __: [rk],
            lambda rk, _, __: [f"{rk}_{h.name}" for h in Haplotype],
            lambda rk, _, __: [rk],
            lambda rk, _, __: [f"{rk}_{h.name}" for h in Haplotype],
        )

    def refsrckey_to_bed_src(self, f: StratInputToBed, rsk: str) -> BedSrc:
        return self.with_ref_data_bed_unsafe(
            rsk,
            f,
            lambda _, __, bf: bf.data.src,
            lambda _, __, bf: bf.data.src,
            lambda _, __, bf: bf.data.src,
            lambda hap, _, __, bf: bf.data.src.from_either(hap),
            lambda hap, _, __, bf: bf.data.src.from_either(hap),
        )

    # def _refkey_to_bed(self, f: StratInputToBed, rk: str) -> BedSrcAndKey:
    #     return BedSrcAndKey(
    #         self._refkey_to_bed_src(f, rk),
    #         self._refkey_to_bed_refsrckey(f, rk),
    #     )

    # def refkey_to_gaps(self, rk: str) -> BedSrcAndKey:
    #     def go(si: StratInputs_[AnyBedT]) -> BedFile[AnyBedT] | None:
    #         return si.gap

    #     return self._refkey_to_bed(go, rk)

    # def refkey_to_trf(self, rk: str) -> BedSrcAndKey:
    #     def go(si: StratInputs_[AnyBedT]) -> BedFile[AnyBedT] | None:
    #         return si.low_complexity.simreps

    #     return self._refkey_to_bed(go, rk)

    # def refkey_to_rmsk(self, rk: str) -> BedSrcAndKey:
    #     def go(si: StratInputs_[AnyBedT]) -> BedFile[AnyBedT] | None:
    #         return si.low_complexity.rmsk

    #     return self._refkey_to_bed(go, rk)

    # def refkey_to_satellite(self, rk: str) -> BedSrcAndKey:
    #     def go(si: StratInputs_[AnyBedT]) -> BedFile[AnyBedT] | None:
    #         return si.low_complexity.satellites

    #     return self._refkey_to_bed(go, rk)

    # def refkey_to_strat(self, k: RefKey) -> AnyStratification:
    #     # def get_hap(x: HapRefKey) -> AnyStratification:
    #     #     return self.haploid_stratifications[x]

    #     # def get_dip(x: Dip1RefKey) -> AnyStratification:
    #     #     return self.diploid1_stratifications[x]

    #     # return either_ref_key(get_hap, get_dip, k)
    #     if isinstance(k, HapRefKey):
    #         return self.haploid_stratifications[k]
    #     elif isinstance(k, Dip1RefKey):
    #         return self.diploid1_stratifications[k]
    #     elif isinstance(k, Dip2RefKey):
    #         return self.diploid2_stratifications[k]
    #     else:
    #         assert_never(k)

    # def refkey_to_mappability_patterns(self, k: RefKey) -> list[str]:
    #     if (m := self.refkey_to_strat(k).strat_inputs.mappability) is None:
    #         return []
    #     else:
    #         return m.unplaced_chr_patterns

    # def buildkey_to_build(self, p: BuildPair) -> HaploidBuild | DiploidBuild:
    #     if isinstance(p.ref, HapRefKey):
    #         return self.haploid_stratifications[p.ref].builds[p.build]
    #     elif isinstance(p.ref, Dip1RefKey):
    #         return self.diploid1_stratifications[p.ref].builds[p.build]
    #     elif isinstance(p.ref, Dip2RefKey):
    #         return self.diploid2_stratifications[p.ref].builds[p.build]
    #     else:
    #         assert_never(p)

    # def buildkey_to_include(self, p: BuildPair) -> HaploidInclude:
    #     return self.buildkey_to_build(p).include

    # def buildkey_to_chr_indices(self, p: BuildPair) -> set[ChrIndex]:
    #     cs = self.buildkey_to_build(p).chr_filter
    #     return set([x for x in ChrIndex]) if len(cs) == 0 else cs

    # def otherkey_to_bed(
    #     self,
    #     p: BuildPair,
    #     lk: OtherLevelKey,
    #     sk: OtherStratKey,
    # ) -> OtherBedFile[HapChrSource[BedSrc]] | OtherBedFile[DipChrSource[BedSrc]]:
    #     return self.buildkey_to_build(p).other_strats[lk][sk]

    # # src getters (for use in downloading inputs)

    # def refkey_to_haploid_ref_src(self, k: HapRefKey) -> RefSrc:
    #     return self.haploid_stratifications[k].ref.src

    # def refkey_to_bench_vcf_src(self, p: BuildPair) -> BedSrc | None:
    #     return fmap_maybe(lambda x: x.bench_vcf.src, self.buildkey_to_build(p).bench)

    # def refkey_to_bench_bed_src(self, p: BuildPair) -> BedSrc | None:
    #     return fmap_maybe(
    #         lambda x: x.bench_bed.data.src, self.buildkey_to_build(p).bench
    #     )

    # def refkey_to_query_vcf_src(self, p: BuildPair) -> BedSrc | None:
    #     return fmap_maybe(lambda x: x.query_vcf.src, self.buildkey_to_build(p).bench)

    # # def _refkey_to_haploid_src(
    # #     self,
    # #     f: Callable[[HaploidStratInputs], BedFile | None],
    # #     k: HapRefKey,
    # # ) -> BedSrc | None:
    # #     i = self.haploid_stratifications[k].strat_inputs
    # #     return fmap_maybe(lambda x: x.src, f(i))

    # # def _refkey_to_diploid_src(
    # #     self,
    # #     f: Callable[[DiploidStratInputs], HaploidOrDiploid[BedFile] | None],
    # #     k: DiploidRefKey,
    # # ) -> list[BedSrc] | None:
    # #     i = f(self.diploid_stratifications[k].strat_inputs)
    # #     if isinstance(i, Diploid_):
    # #         return [i.hap1.src, i.hap2.src]
    # #     if isinstance(i, BedFile):
    # #         return [i.src]
    # #     else:
    # #         return None

    # # def _refkey_to_src(
    # #     self,
    # #     f: Callable[[HaploidStratInputs], BedFile | None],
    # #     g: Callable[[DiploidStratInputs], BedFile | Diploid_[BedFile] | None],
    # #     k: RefKey,
    # # ) -> list[BedSrc]:
    # #     if isinstance(k, HapRefKey):
    # #         return fmap_maybe_def(
    # #             lambda x: [x.src],
    # #             f(self.haploid_stratifications[k].strat_inputs),
    # #             [],
    # #         )
    # #     elif isinstance(k, DiploidRefKey):
    # #         b = g(self.diploid_stratifications[k].strat_inputs)
    # #         if b is None:
    # #             return []
    # #         elif isinstance(b, BedFile):
    # #             return [b.src]
    # #         elif isinstance(b, Diploid_):
    # #             return [b.hap1.src, b.hap2.src]
    # #         else:
    # #             assert_never(b)
    # #     else:
    # #         assert_never(k)

    # # def _refkey_to_src(
    # #     self,
    # #     f: Callable[[StratInputs_[AnyBedSrcT]], BedFile_[AnyBedSrcT] | None],
    # #     k: RefKey,
    # # ) -> list[BedSrc]:
    # #     return self.refkey_to_strat(k).strat_inputs.to_src(f)
    # #     # b = f(i)
    # #     # if b is None:
    # #     #     return []
    # #     # else:
    # #     #     return diploid_to_list(b.src)

    # def refkey_to_gap_src(self, k: RefKey) -> list[BedSrc]:
    #     return self.refkey_to_strat(k).strat_inputs.gap_src

    # # def refkey_to_x_features_src(self, k: RefKey) -> list[BedSrc]:
    # #     return fmap_maybe_def(
    # #         [],
    # #         lambda x: x.x_bed.src_list,
    # #         self.refkey_to_strat(k).strat_inputs.xy.features,
    # #     )

    # # def refkey_to_y_features_src(self, k: RefKey) -> list[BedSrc]:
    # #     return fmap_maybe_def(
    # #         [],
    # #         lambda x: x.y_bed.src_list,
    # #         self.refkey_to_strat(k).strat_inputs.xy.features,
    # #     )

    # def refkey_to_x_PAR(self, k: RefKey) -> XYPar | None:
    #     return self.refkey_to_strat(k).strat_inputs.xy.x_par

    # def refkey_to_y_PAR(self, k: RefKey) -> XYPar | None:
    #     return self.refkey_to_strat(k).strat_inputs.xy.y_par

    # def refkey_to_simreps_src(self, k: RefKey) -> list[BedSrc]:
    #     return self.refkey_to_strat(k).strat_inputs.simreps_src

    # def refkey_to_rmsk_src(self, k: RefKey) -> list[BedSrc]:
    #     return self.refkey_to_strat(k).strat_inputs.rmsk_src

    # def refkey_to_satellite_src(self, k: RefKey) -> list[BedSrc]:
    #     return self.refkey_to_strat(k).strat_inputs.satellites_src

    # def refkey_to_superdups_src(self, k: RefKey) -> list[BedSrc]:
    #     return self.refkey_to_strat(k).strat_inputs.superdups_src

    # def refkey_to_functional_ftbl_src(self, k: RefKey) -> list[BedSrc]:
    #     return self.refkey_to_strat(k).strat_inputs.ftbl_src

    # def refkey_to_functional_gff_src(self, k: RefKey) -> list[BedSrc]:
    #     return self.refkey_to_strat(k).strat_inputs.gff_src

    # def otherkey_to_src(
    #     self,
    #     p: BuildPair,
    #     lk: OtherLevelKey,
    #     sk: OtherStratKey,
    # ) -> list[BedSrc]:
    #     return self.otherkey_to_bed(p, lk, sk).src_list

    # chromosome standardization

    # def refkey_to_final_chr_pattern(self, k: RefKey) -> ChrPattern:
    #     return self.refkey_to_strat(k).ref.chr_pattern

    # def refkey_to_bench_chr_pattern(self, p: BuildPair) -> HapChrPattern | None:
    #     return fmap_maybe(
    #         lambda x: x.bench_vcf.chr_pattern,
    #         self.buildkey_to_build(p).bench,
    #     )

    # def refkey_to_query_chr_pattern(self, p: BuildPair) -> HapChrPattern | None:
    #     return fmap_maybe(
    #         lambda x: x.query_vcf.chr_pattern,
    #         self.buildkey_to_build(p).bench,
    #     )

    # # TODO this nomenclature sucks
    # def buildkey_to_hap_chr_conversion(
    #     self,
    #     p: HaploidBuildPair,
    #     fromChr: HapChrPattern,
    # ) -> HapToHapChrConversion:
    #     toChr = self.haploid_stratifications[p.ref].ref.chr_pattern
    #     cis = self.buildkey_to_chr_indices(p)
    #     return HapToHapChrConversion(fromChr, toChr, cis)

    # def buildkey_to_dip_chr_conversion1(
    #     self,
    #     p: Diploid1BuildPair,
    #     fromChr: Diploid_[HapChrPattern] | DipChrPattern,
    # ) -> DipToDipChrConversion:
    #     toChr = self.diploid1_stratifications[p.ref].ref.chr_pattern
    #     cis = self.buildkey_to_chr_indices(p)
    #     return DipToDipChrConversion(
    #         fromChr, toChr, set([(h, c) for c in cis for h in Haplotype])
    #     )

    # def buildkey_to_dip_chr_conversion2(
    #     self,
    #     p: Diploid2BuildPair,
    #     fromChr: Diploid_[HapChrPattern],
    #     h: Haplotype,
    # ) -> tuple[HapToHapChrConversion, HapToHapChrConversion]:
    #     toChr = self.diploid2_stratifications[p.ref].ref.chr_pattern
    #     cis = self.buildkey_to_chr_indices(p)
    #     return (
    #         HapToHapChrConversion(fromChr.hap1, toChr.hap1, cis),
    #         HapToHapChrConversion(fromChr.hap2, toChr.hap2, cis),
    #     )

    # def buildkey_to_dip_chr_conversion12(
    #     self,
    #     p: Diploid1BuildPair,
    #     fromChr: HapChrPattern,
    #     h: Haplotype,
    # ) -> HapToDipChrConversion:
    #     toChr = self.diploid1_stratifications[p.ref].ref.chr_pattern
    #     cis = self.buildkey_to_chr_indices(p)
    #     return HapToDipChrConversion(fromChr, toChr, cis, h)

    # def buildkey_to_chr_conversion(
    #     self,
    #     p: BuildPair,
    #     fromChr: HapChrPattern | ,
    # ) -> ChrConversion_:
    #     if isinstance(p.ref, HapBuildKey):
    #         return self.buildkey_to_hap_chr_conversion(p, fromChar)
    #     if isinstance(p.ref, DiploidBuildKey):
    #         return
    #     else:
    #         assert_never(p.ref)

    # def buildkey_to_final_chr_mapping(self, p: BuildPair) -> dict[int, str]:
    #     if isinstance(p.ref, HapRefKey):
    #         return self.buildkey_to_dip_chr_conversion
    #     elif isinstance(p.ref, DiploidRefKey):
    #         return
    #     else:
    #         assert_never(p.ref)
    #     # pat = self.refkey_to_final_chr_pattern(p.ref)
    #     # cs = self.buildkey_to_chr_indices(p)
    #     # return {c.value: n for c in cs if (n := pat.to_chr_names(c) is not None}

    # def buildkey_to_mappability(
    #     self, p: BuildPair
    # ) -> tuple[list[int], list[int], list[int]]:
    #     ms = self.buildkey_to_include(p).mappability
    #     l, m, e = unzip((m.length, m.mismatches, m.indels) for m in ms)
    #     return ([*l], [*m], [*e])

    # def buildkey_to_comparison(self, p: BuildPair) -> BuildCompare | None:
    #     return self.buildkey_to_build(p).comparison

    # def buildkey_to_comparekey(self, p: BuildPair) -> CompareKey | None:
    #     return fmap_maybe(lambda x: x.other, self.buildkey_to_comparison(p))

    # # include switches (for controlling which snakemake rules to activate)

    # def has_low_complexity_rmsk(self, rk: RefKey) -> bool:
    #     return self.refkey_to_strat(rk).strat_inputs.low_complexity.rmsk is not None

    # def has_low_complexity_simreps(self, rk: RefKey) -> bool:
    #     return self.refkey_to_strat(rk).strat_inputs.low_complexity.simreps is not None

    # def has_low_complexity_censat(self, rk: RefKey) -> bool:
    #     return (
    #         self.refkey_to_strat(rk).strat_inputs.low_complexity.satellites is not None
    #     )

    # def _want_chr_index(self, p: BuildPair, i: ChrIndex) -> bool:
    #     cis = self.buildkey_to_chr_indices(p)
    #     return i in cis

    # def want_xy_x(self, p: BuildPair) -> bool:
    #     return self._want_chr_index(p, ChrIndex.CHRX) and self.buildkey_to_include(p).xy

    # def want_xy_y(self, p: BuildPair) -> bool:
    #     return self._want_chr_index(p, ChrIndex.CHRY) and self.buildkey_to_include(p).xy

    # def wanted_xy_chr_names(self, p: BuildPair) -> list[str]:
    #     return [
    #         i.chr_name
    #         for i in [ChrIndex.CHRX, ChrIndex.CHRY]
    #         if self._want_chr_index(p, i)
    #     ]

    # # def want_xy_sex(self, rk: RefKey, bk: BuildKey) -> bool:
    # #     return self.want_xy_x(rk, bk) and self.want_xy_y(rk, bk)

    # def want_x_PAR(self, p: BuildPair) -> bool:
    #     return self.want_xy_x(p) and self.refkey_to_x_PAR(p.ref) is not None

    # def want_y_PAR(self, p: BuildPair) -> bool:
    #     return self.want_xy_y(p) and self.refkey_to_y_PAR(p.ref) is not None

    # def want_xy_auto(self, p: BuildPair) -> bool:
    #     cis = self.buildkey_to_chr_indices(p)
    #     return len(cis - set([ChrIndex.CHRX, ChrIndex.CHRY])) > 0

    # def want_xy_XTR(self, rk: RefKey) -> bool:
    #     f = self.refkey_to_strat(rk).strat_inputs.xy.features
    #     return f is not None and f.xtr

    # def want_xy_ampliconic(self, rk: RefKey) -> bool:
    #     f = self.refkey_to_strat(rk).strat_inputs.xy.features
    #     return f is not None and f.ampliconic

    # def want_low_complexity(self, p: BuildPair) -> bool:
    #     return self.buildkey_to_include(p).low_complexity

    # def want_gc(self, p: BuildPair) -> bool:
    #     return self.buildkey_to_include(p).gc is not None

    # def want_functional(self, p: BuildPair) -> bool:
    #     return (
    #         self.buildkey_to_include(p).functional
    #         and self.refkey_to_functional_ftbl_src(p.ref) is not None
    #         and self.refkey_to_functional_gff_src(p.ref) is not None
    #     )

    # def want_telomeres(self, p: BuildPair) -> bool:
    #     return self.buildkey_to_include(p).telomeres

    # def want_segdups(self, p: BuildPair) -> bool:
    #     return (
    #         self.refkey_to_superdups_src(p.ref) is not None
    #         and self.buildkey_to_include(p).segdups
    #     )

    # def _want_union(self, p: BuildPair) -> bool:
    #     return self.buildkey_to_include(p).union

    # def want_mappability(self, p: BuildPair) -> bool:
    #     return (
    #         self.refkey_to_strat(p.ref).strat_inputs.mappability is not None
    #         and len(self.buildkey_to_include(p).mappability) > 0
    #     )

    # def want_segdup_and_map(self, p: BuildPair) -> bool:
    #     return (
    #         self.buildkey_to_include(p).union
    #         and self.want_segdups(p)
    #         and self.want_mappability(p)
    #     )

    # def want_alldifficult(self, p: BuildPair) -> bool:
    #     return (
    #         self.want_segdup_and_map(p)
    #         and self.want_low_complexity(p)
    #         and self.want_gc(p)
    #     )

    # def want_benchmark(self, p: BuildPair) -> bool:
    #     return self.buildkey_to_build(p).bench is not None

    # def want_gaps(self, rk: RefKey) -> bool:
    #     return self.refkey_to_strat(rk).strat_inputs.gap is not None

    # def want_vdj(self, p: BuildPair) -> bool:
    #     cis = self.buildkey_to_chr_indices(p)
    #     vdj_chrs = {ChrIndex(i) for i in [2, 7, 14, 22]}
    #     return (
    #         self.buildkey_to_include(p).vdj
    #         and self.refkey_to_functional_ftbl_src(p.ref) is not None
    #         and self.refkey_to_functional_gff_src(p.ref) is not None
    #         and len(cis & vdj_chrs) > 0
    #     )

    # key lists for downloading resources

    @property
    def _all_haploid_builds(self) -> list[HaploidBuildPair]:
        return [
            HaploidBuildPair(rk, bk)
            for rk, s in self.haploid_stratifications.items()
            for bk in s.builds
        ]

    @property
    def _all_diploid1_builds(self) -> list[Diploid1BuildPair]:
        return [
            Diploid1BuildPair(rk, bk)
            for rk, s in self.diploid1_stratifications.items()
            for bk in s.builds
        ]

    @property
    def _all_diploid2_builds(self) -> list[Diploid2BuildPair]:
        return [
            Diploid2BuildPair(rk, bk)
            for rk, s in self.diploid2_stratifications.items()
            for bk in s.builds
        ]

    @property
    def _all_builds(self) -> list[BuildPair]:
        return (
            self._all_haploid_builds
            + self._all_diploid1_builds
            + self._all_diploid2_builds
        )

    @property
    def all_refkeys(self) -> list[RefKey]:
        return (
            [*self.haploid_stratifications]
            + [*self.diploid1_stratifications]
            + [*self.diploid2_stratifications]
        )

    # def strat_dict(
    #     self, rk: str, bk: str
    # ) -> (
    #     tuple[HapRefKey, HapBuildKey, HaploidStratification]
    #     | tuple[Dip1RefKey, Dip1BuildKey, Dip1Strat]
    #     | tuple[Dip2RefKey, Dip2BuildKey, Dip2Strat]
    # ):
    #     if rk in self.haploid_stratifications:
    #         _rk0 = HapRefKey(rk)
    #         return (
    #             _rk0,
    #             HapBuildKey(bk),
    #             self.haploid_stratifications[_rk0],
    #         )
    #     elif rk in self.diploid1_stratifications:
    #         _rk1 = Dip1RefKey(rk)
    #         return (
    #             _rk1,
    #             Dip1BuildKey(bk),
    #             self.diploid1_stratifications[_rk1],
    #         )
    #     elif rk in self.diploid2_stratifications:
    #         _rk2 = Dip2RefKey(rk)
    #         return (
    #             _rk2,
    #             Dip2BuildKey(bk),
    #             self.diploid2_stratifications[_rk2],
    #         )
    #     else:
    #         # TODO this seems sloppy, not sure if I want it here
    #         assert False, "this should not happen"

    def parse_refkey(self, rk: str) -> AnyRefKey:
        if HapRefKey_(rk) in self.haploid_stratifications:
            return RefKey_(HapRefKey_(rk))
        elif Dip1RefKey_(rk) in self.diploid1_stratifications:
            return RefKey_(Dip1RefKey_(rk))
        elif Dip2RefKey_(rk) in self.diploid2_stratifications:
            return RefKey_(Dip2RefKey_(rk))
        else:
            raise DesignError(f"could not get stratification data for ref key '{rk}'")

    def to_ref_data(self, rk: str) -> AnyStratification:
        try:
            return self.haploid_stratifications[HapRefKey_(rk)]
        except KeyError:
            pass

        try:
            return self.diploid1_stratifications[Dip1RefKey_(rk)]
        except KeyError:
            pass

        try:
            return self.diploid2_stratifications[Dip2RefKey_(rk)]
        except KeyError:
            pass

        raise DesignError(f"could not get stratification data for ref key '{rk}'")

    def to_build_data(self, rk: str, bk: str) -> AnyBuildData:
        try:
            return HaploidBuildData(
                *astuple(
                    self.haploid_stratifications.to_build_data_unsafe(
                        HapRefKey(rk), HapBuildKey(bk)
                    )
                )
            )
        except KeyError:
            pass

        try:
            return Diploid1BuildData(
                *astuple(
                    self.diploid1_stratifications.to_build_data_unsafe(
                        Dip1RefKey(rk), Dip1BuildKey(bk)
                    )
                )
            )
        except KeyError:
            pass

        try:
            return Diploid2BuildData(
                *astuple(
                    self.diploid2_stratifications.to_build_data_unsafe(
                        Dip2RefKey(rk), Dip2BuildKey(bk)
                    )
                )
            )
        except KeyError:
            pass

        raise DesignError(
            f"could not get build data for ref key '{rk}' and build key '{bk}'"
        )

    def with_ref_data_unsafe(
        self,
        rfk: str,
        hap_f: Callable[[HaploidStratification], X],
        dip1_f: Callable[[Dip1Strat], X],
        dip2_f: Callable[[Haplotype, Dip2Strat], X],
    ) -> X:
        rk_, hap = parse_refkey(rfk)
        rd = self.to_ref_data(rk_)
        if isinstance(rd.ref, HapChrSource) and hap is None:
            return hap_f(rd)
        elif isinstance(rd.ref, DipChrSource1) and hap is None:
            return dip1_f(rd)
        elif isinstance(rd.ref, DipChrSource2) and hap is not None:
            return dip2_f(hap, rd)
        else:
            raise DesignError(
                f"Invalid ref data with type '{type(rd)}' for key '{rfk}'"
            )

    # def with_ref_data_0_unsafe(
    #     self,
    #     rk: str,
    #     hap_f: Callable[[HaploidStratification], X],
    #     dip1_f: Callable[[Dip1Strat], X],
    #     hap1_f: Callable[[Haplotype, Dip1Strat], X],
    #     dip2_f: Callable[[Dip2Strat], X],
    #     hap2_f: Callable[[Haplotype, Dip2Strat], X],
    # ) -> X:
    #     rk_, hap = parse_refkey(rk)
    #     rd = self.to_ref_data(rk_)
    #     if isinstance(rd.ref, HapChrSource) and hap is None:
    #         return hap_f(rd)
    #     elif isinstance(rd.ref, DipChrSource1) and hap is None:
    #         return dip1_f(rd)
    #     elif isinstance(rd.ref, DipChrSource2) and hap is None:
    #         return dip2_f(rd)
    #     elif isinstance(rd.ref, DipChrSource1) and hap is not None:
    #         return hap1_f(hap, rd)
    #     elif isinstance(rd.ref, DipChrSource2) and hap is not None:
    #         return hap2_f(hap, rd)
    #     else:
    #         raise DesignError(f"Invalid ref data with type '{type(rd)}' for key '{rk}'")

    def with_build_data_ref_unsafe(
        self,
        rsk: str,
        bk: str,
        hap_f: Callable[[HaploidBuildData], X],
        dip1_f: Callable[[Diploid1BuildData], X],
        dip2_f: Callable[[Haplotype, Diploid2BuildData], X],
    ) -> X:
        return self.with_ref_data_unsafe(
            rsk,
            lambda s: hap_f(
                HaploidBuildData(
                    *astuple(
                        s.to_build_data_unsafe(HapBuildKey(bk)),
                    )
                )
            ),
            lambda s: dip1_f(
                Diploid1BuildData(
                    *astuple(
                        s.to_build_data_unsafe(Dip1BuildKey(bk)),
                    )
                )
            ),
            lambda hap, s: dip2_f(
                hap,
                Diploid2BuildData(
                    *astuple(
                        s.to_build_data_unsafe(Dip2BuildKey(bk)),
                    )
                ),
            ),
        )

    def with_build_data(
        self,
        rk: str,
        bk: str,
        hap_f: Callable[[HaploidBuildData], X],
        dip1_f: Callable[[Diploid1BuildData], X],
        dip2_f: Callable[[Haplotype, Diploid2BuildData], X],
    ) -> X:
        rk_, hap = parse_refkey(rk)
        bd = self.to_build_data(rk_, bk)
        if isinstance(bd, HaploidBuildData) and hap is None:
            return hap_f(bd)
        elif isinstance(bd, Diploid1BuildData) and hap is None:
            return dip1_f(bd)
        elif isinstance(bd, Diploid2BuildData) and hap is not None:
            return dip2_f(hap, bd)
        else:
            raise DesignError(
                f"Invalid build data with type '{type(bd)}' for keys '{rk}' and '{bk}'"
            )

    # TODO feed refkey to the HO-functions so that there is no confusion as to
    # when the refkey has been parsed
    def with_ref_data_and_bed(
        self,
        rk: RefKey_[RefKeyT_],
        get_bed_f: StratInputToBed,
        hap_f: Callable[[HapRefKey, HaploidStratification, HapBedFile], Z],
        dip_1to1_f: Callable[[Dip1RefKey, Dip1Strat, Dip1BedFile], Z],
        dip_1to2_f: Callable[[Dip2RefKey, Dip2Strat, Dip1BedFile], Z],
        dip_2to1_f: Callable[[Dip1RefKey, Dip1Strat, Dip2BedFile], Z],
        dip_2to2_f: Callable[[Dip2RefKey, Dip2Strat, Dip2BedFile], Z],
    ) -> Z:
        rd = self.to_ref_data(rk)
        # ref and bed are both haploid
        if (
            isinstance(rd.ref, HapChrSource)
            and (bf0 := get_bed_f(rd.strat_inputs)) is not None
        ):
            return hap_f(HapRefKey(rk), rd, bf0)
        # ref and bed are both diploid1
        elif (
            isinstance(rd.ref, DipChrSource1)
            and (bf1 := get_bed_f(rd.strat_inputs)) is not None
            and is_dip1_bed(bf1)
        ):
            return dip_1to1_f(Dip1RefKey(rk), rd, bf1)
        # ref is diploid2, bed is diploid1 (combine)
        elif (
            isinstance(rd.ref, DipChrSource2)
            and (bf3 := get_bed_f(rd.strat_inputs)) is not None
            and is_dip1_bed(bf3)
        ):
            return dip_1to2_f(Dip2RefKey(rk), rd, bf3)
        # ref is diploid1, bed is diploid2 (split)
        elif (
            isinstance(rd.ref, DipChrSource1)
            and (bf4 := get_bed_f(rd.strat_inputs)) is not None
            and is_dip2_bed(bf4)
        ):
            return dip_2to1_f(Dip1RefKey(rk), rd, bf4)
        # ref and bed are both diploid2
        elif (
            isinstance(rd.ref, DipChrSource2)
            and (bf2 := get_bed_f(rd.strat_inputs)) is not None
            and is_dip2_bed(bf2)
        ):
            return dip_2to2_f(Dip2RefKey(rk), rd, bf2)
        else:
            # TODO this errror seems wonky
            raise DesignError(f"Invalid ref data with type '{type(rd)}' for key '{rk}'")

    def with_ref_data_and_bed_hap(
        self,
        rfk: str,
        get_bed_f: StratInputToBed,
        hap_f: Callable[[HapRefKey, HaploidStratification, HapBedFile], Z],
        dip_1to1_f: Callable[[Dip1RefKey, Dip1Strat, Dip1BedFile], Z],
        dip_1to2_f: Callable[[Haplotype, Dip2RefKey, Dip2Strat, Dip1BedFile], Z],
        dip_2to1_f: Callable[[Dip1RefKey, Dip1Strat, Dip2BedFile], Z],
        dip_2to2_f: Callable[[Haplotype, Dip2RefKey, Dip2Strat, Dip2BedFile], Z],
    ) -> Z:
        rk, hap = parse_refkey(rfk)
        return self.with_ref_data_and_bed(
            rk,
            get_bed_f,
            hap_f,
            dip_1to1_f,
            lambda rk, rd, bf: not_none_unsafe(
                hap, lambda h: dip_1to2_f(h, rk, rd, bf)
            ),
            dip_2to1_f,
            lambda rk, rd, bf: not_none_unsafe(
                hap, lambda h: dip_2to2_f(h, rk, rd, bf)
            ),
        )

    def with_ref_data_bed_unsafe(
        self,
        rsk: str,
        get_bed_f: StratInputToBed,
        hap_f: Callable[[HapRefKey, HaploidStratification, HapBedFile], X],
        dip_1to1_f: Callable[[Dip1RefKey, Dip1Strat, Dip1BedFile], X],
        dip_1to2_f: Callable[[Dip2RefKey, Dip2Strat, Dip1BedFile], X],
        dip_2to1_f: Callable[[Haplotype, Dip1RefKey, Dip1Strat, Dip2BedFile], X],
        dip_2to2_f: Callable[[Haplotype, Dip2RefKey, Dip2Strat, Dip2BedFile], X],
    ) -> X:
        rk, hap = parse_refkey(rsk)
        return self.with_ref_data_and_bed(
            rsk,
            get_bed_f,
            lambda rk, rd, bd: none_unsafe(hap, hap_f(rk, rd, bd)),
            lambda rk, rd, bd: none_unsafe(hap, dip_1to1_f(rk, rd, bd)),
            lambda rk, rd, bd: none_unsafe(hap, dip_1to2_f(rk, rd, bd)),
            lambda rk, rd, bf: not_none_unsafe(
                hap, lambda h: dip_2to1_f(h, rk, rd, bf)
            ),
            lambda rk, rd, bf: not_none_unsafe(
                hap, lambda h: dip_2to2_f(h, rk, rd, bf)
            ),
        )

    def with_build_data_and_bed(
        self,
        rk: str,
        bk: str,
        get_bed_f: BuildDataToBed,
        hap_f: Callable[[HapRefKey, HaploidBuildData, HapBedFile], Z],
        dip_1to1_f: Callable[[Dip1RefKey, Diploid1BuildData, Dip1BedFile], Z],
        dip_1to2_f: Callable[[Dip2RefKey, Diploid2BuildData, Dip1BedFile], Z],
        dip_2to1_f: Callable[[Dip1RefKey, Diploid1BuildData, Dip2BedFile], Z],
        dip_2to2_f: Callable[[Dip2RefKey, Diploid2BuildData, Dip2BedFile], Z],
    ) -> Z:
        return self.with_ref_data_and_bed(
            rk,
            lambda bd: get_bed_f(bd.strat_inputs),
            lambda rk, rd, bf: hap_f(
                rk,
                HaploidBuildData(*astuple(rd.to_build_data_unsafe(HapBuildKey(bk)))),
                bf,
            ),
            lambda rk, rd, bf: dip_1to1_f(
                rk,
                Diploid1BuildData(*astuple(rd.to_build_data_unsafe(Dip1BuildKey(bk)))),
                bf,
            ),
            lambda rk, rd, bf: dip_1to2_f(
                rk,
                Diploid2BuildData(*astuple(rd.to_build_data_unsafe(Dip2BuildKey(bk)))),
                bf,
            ),
            lambda rk, rd, bf: dip_2to1_f(
                rk,
                Diploid1BuildData(*astuple(rd.to_build_data_unsafe(Dip1BuildKey(bk)))),
                bf,
            ),
            lambda rk, rd, bf: dip_2to2_f(
                rk,
                Diploid2BuildData(*astuple(rd.to_build_data_unsafe(Dip2BuildKey(bk)))),
                bf,
            ),
        )

    def with_build_data_and_bed_i(
        self,
        rk: str,
        bk: str,
        inputs: list[X],
        get_bed_f: BuildDataToBed,
        hap_f: Callable[[X, HapRefKey, HaploidBuildData, HapBedFile], Z],
        dip_1to1_f: Callable[[X, Dip1RefKey, Diploid1BuildData, Dip1BedFile], Z],
        dip_1to2_f: Callable[[X, Dip2RefKey, Diploid2BuildData, Dip1BedFile], Z],
        dip_2to1_f: Callable[
            [tuple[X, X], Dip1RefKey, Diploid1BuildData, Dip2BedFile], Z
        ],
        dip_2to2_f: Callable[
            [tuple[X, X], Dip2RefKey, Diploid2BuildData, Dip2BedFile], Z
        ],
    ) -> Z:
        return self.with_build_data_and_bed(
            rk,
            bk,
            get_bed_f,
            lambda rk, bd, bf: match1_unsafe(inputs, lambda i: hap_f(i, rk, bd, bf)),
            lambda rk, bd, bf: match1_unsafe(
                inputs, lambda i: dip_1to1_f(i, rk, bd, bf)
            ),
            lambda rk, bd, bf: match1_unsafe(
                inputs, lambda i: dip_1to2_f(i, rk, bd, bf)
            ),
            lambda rk, bd, bf: match2_unsafe(
                inputs, lambda i0, i1: dip_2to1_f((i0, i1), rk, bd, bf)
            ),
            lambda rk, bd, bf: match2_unsafe(
                inputs, lambda i0, i1: dip_2to2_f((i0, i1), rk, bd, bf)
            ),
        )

    def with_build_data_and_bed_io(
        self,
        rk: str,
        bk: str,
        inputs: list[X],
        output_f: Callable[[str], Y],
        get_bed_f: BuildDataToBed,
        hap_f: Callable[[X, Y, HaploidBuildData, HapBedFile], Z],
        dip_1to1_f: Callable[[X, Y, Diploid1BuildData, Dip1BedFile], Z],
        dip_1to2_f: Callable[[X, tuple[Y, Y], Diploid2BuildData, Dip1BedFile], Z],
        dip_2to1_f: Callable[[tuple[X, X], Y, Diploid1BuildData, Dip2BedFile], Z],
        dip_2to2_f: Callable[
            [tuple[X, X], tuple[Y, Y], Diploid2BuildData, Dip2BedFile], Z
        ],
    ) -> Z:
        return self.with_build_data_and_bed_i(
            rk,
            bk,
            inputs,
            get_bed_f,
            lambda i, rk, bd, bf: hap_f(i, output_f(rk), bd, bf),
            lambda i, rk, bd, bf: dip_1to1_f(i, output_f(rk), bd, bf),
            lambda i, rk, bd, bf: dip_1to2_f(
                i,
                (
                    output_f(f"{rk}_{Haplotype.HAP1.name}"),
                    output_f(f"{rk}_{Haplotype.HAP2.name}"),
                ),
                bd,
                bf,
            ),
            lambda i, rk, bd, bf: dip_2to1_f(i, output_f(rk), bd, bf),
            lambda i, rk, bd, bf: dip_2to2_f(
                i,
                (
                    output_f(f"{rk}_{Haplotype.HAP1.name}"),
                    output_f(f"{rk}_{Haplotype.HAP2.name}"),
                ),
                bd,
                bf,
            ),
        )

    def with_build_data_and_bed_io_(
        self,
        rk: str,
        bk: str,
        inputs: list[X],
        output_f: Callable[[str], Y],
        get_bed_f: BuildDataToBed,
        hap_f: Callable[[X, Y, HaploidBuildData, HapBedFile], None],
        dip_1to1_f: Callable[[X, Y, Diploid1BuildData, Dip1BedFile], None],
        dip_1to2_f: Callable[[X, tuple[Y, Y], Diploid2BuildData, Dip1BedFile], None],
        dip_2to1_f: Callable[[tuple[X, X], Y, Diploid1BuildData, Dip2BedFile], None],
        dip_2to2_f: Callable[[X, Y, Haplotype, Diploid2BuildData, Dip2BedFile], None],
    ) -> None:
        def _dip_2to2_f(
            i: tuple[X, X],
            o: tuple[Y, Y],
            bd: Diploid2BuildData,
            bf: Dip2BedFile,
        ) -> None:
            dip_2to2_f(i[0], o[0], Haplotype.HAP1, bd, bf)
            dip_2to2_f(i[1], o[1], Haplotype.HAP2, bd, bf)

        self.with_build_data_and_bed_io(
            rk,
            bk,
            inputs,
            output_f,
            get_bed_f,
            hap_f,
            dip_1to1_f,
            dip_1to2_f,
            dip_2to1_f,
            _dip_2to2_f,
        )

    # def with_build_src_data_unsafe(
    #     self,
    #     rk: str,
    #     bk: str,
    #     inputs: list[X],
    #     outputs: list[Y],
    #     hap: Haplotype | None,
    #     hap_f: Callable[[X, Y, HaploidBuildData], Z],
    #     dip1_f: Callable[[X, Y, Diploid1BuildData], Z],
    #     hap2_f: Callable[[X, Y, Haplotype, Diploid2BuildData], Z],
    #     dip_1to2_f: Callable[[X, tuple[Y, Y], Diploid2BuildData], Z],
    #     dip_1to2_f: Callable[[tuple[X, X], Y, Diploid1BuildData], Z],
    # ) -> None:
    #     # TODO make these errors more meaningful
    #     bd = self.to_build_data(rk, bk)
    #     match (inputs, outputs):
    #         case ([i], [o]):
    #             # one haplotype for both bed and ref (no combine)
    #             if isinstance(bd, HaploidBuildData) and hap is None:
    #                 hap_f(i, o, bd)
    #             # one bed with both haps in it; one reference with both haps (no combine)
    #             elif isinstance(bd, Diploid1BuildData) and hap is None:
    #                 dip1_f(i, o, bd)
    #             # one bed and one ref for a single haplotype in a diploid reference
    #             elif isinstance(bd, Diploid2BuildData) and hap is not None:
    #                 hap2_f(i, o, hap, bd)
    #         case ([i], [o0, o1]):
    #             # one bed with both haps in it; two references for both haps (split)
    #             if isinstance(bd, Diploid2BuildData) and hap is None:
    #                 dip_1to2_f(i, (o0, o1), bd)
    #             else:
    #                 assert False, "this should not happen"
    #         case ([i0, i1], [o]):
    #             # two beds for both haps; one reference with both haps (combine)
    #             if isinstance(bd, Diploid1BuildData) and hap is None:
    #                 dip_1to2_f((i0, i1), o, bd)
    #             else:
    #                 assert False, "this should not happen"
    #         case _:
    #             assert False, "this should not happen"

    # def to_build_pair(self, rk: str, bk: str) -> BuildPair:
    #     if rk in self.haploid_stratifications:
    #         return BuildPair_(HapRefKey(rk), HapBuildKey(bk))
    #     elif rk in self.diploid1_stratifications:
    #         return BuildPair_(Dip1RefKey(rk), Dip1BuildKey(bk))
    #     elif rk in self.diploid2_stratifications:
    #         return BuildPair_(Dip2RefKey(rk), Dip2BuildKey(bk))
    #     else:
    #         # TODO this seems sloppy, not sure if I want it here
    #         assert False, "this should not happen"

    def _all_refkey_from_want(
        self,
        f: Callable[[BuildPair], bool],
    ) -> list[RefKey]:
        return [*unique_everseen(x[0] for x in self._all_builds if f(x))]

    # TODO probably need to redo these
    # @property
    # def all_refkey_gap(self) -> list[RefKey]:
    #     return [k for k in self.all_refkeys if self.refkey_to_gap_src(k) is not None]

    # @property
    # def all_refkey_rmsk(self) -> list[RefKey]:
    #     return self._all_refkey_from_want(
    #         lambda p: self.want_low_complexity(p)
    #         and self.has_low_complexity_simreps(p.ref)
    #     )

    # @property
    # def all_refkey_trf(self) -> list[RefKey]:
    #     return self._all_refkey_from_want(
    #         lambda p: self.want_low_complexity(p)
    #         and self.has_low_complexity_rmsk(p.ref)
    #     )

    # @property
    # def all_refkey_censat(self) -> list[RefKey]:
    #     return self._all_refkey_from_want(
    #         lambda p: self.want_low_complexity(p)
    #         and self.has_low_complexity_censat(p.ref)
    #     )

    # @property
    # def all_refkey_functional(self) -> list[RefKey]:
    #     return self._all_refkey_from_want(self.want_functional)

    # @property
    # def all_refkey_segdups(self) -> list[RefKey]:
    #     return self._all_refkey_from_want(self.want_segdups)


# def narrow_bed_union(
#     x: BedFile[AnyBedT],
# ) -> HapBedFile | Dip1BedFile | Dip2BedFile | None:
#     if x is None:
#         return None
#     elif isinstance(x.data, HapChrSource):
#         return x
#     elif is_dip1_bed(x):
#         return x
#     elif is_dip2_bed(x):
#         return x
#     else:
#         # ASSUME this is unreachable
#         return None
#         # assert_never(y)
