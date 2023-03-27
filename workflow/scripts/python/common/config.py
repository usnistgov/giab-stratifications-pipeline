import gzip
from pathlib import Path
from pydantic import BaseModel as BaseModel_
from pydantic import validator, HttpUrl, FilePath
from enum import Enum, unique
from typing import NewType, Any, Callable, TypeVar
from typing_extensions import Self
from Bio import bgzf  # type: ignore

X = TypeVar("X")
Y = TypeVar("Y")


BuildKey = NewType("BuildKey", str)
RefKey = NewType("RefKey", str)


def fmap_maybe(f: Callable[[X], Y], x: X | None) -> None | Y:
    return None if x is None else f(x)


class RefFmt(Enum):
    NOZIP = "nozip"
    GZIP = "gzip"
    BGZIP = "bgzip"


class BedFmt(Enum):
    NOZIP = "nozip"
    GZIP = "gzip"


class XYFeature(Enum):
    XTR = "XTR"
    Ampliconic = "Ampliconic"


@unique
class ChrIndex(Enum):
    _ignore_ = "ChrIndex i"
    ChrIndex = vars()
    for i in range(1, 23):
        ChrIndex[f"CHR{i}"] = i
    CHRX = 23
    CHRY = 24

    @classmethod
    def from_name(cls, n: str) -> Self:
        try:
            return next(i for i in cls if i.chr_name == n)
        except StopIteration:
            raise ValueError(f"could make chr index from name '{n}'")

    def __init__(self, i: int) -> None:
        self.chr_name: str = "X" if i == 23 else ("Y" if i == 24 else str(i))

    def chr_name_full(self, prefix: str) -> str:
        return f"{prefix}{self.chr_name}"


class BaseModel(BaseModel_):
    class Config:
        frozen = True
        extra = "forbid"


class Paths(BaseModel):
    resources: Path
    results: Path


class Tools(BaseModel):
    repseq: HttpUrl
    gemlib: HttpUrl


# TODO non-negative ints which cannot equal each other
class BedColumns(BaseModel):
    chr: int = 0
    start: int = 1
    end: int = 2

    @property
    def columns(self) -> list[int]:
        return [self.chr, self.start, self.end]


class FileSrc_(BaseModel):
    filepath: FilePath


class BedFileSrc(FileSrc_):
    filepath: FilePath

    @validator("filepath")
    def is_gzip(cls, v: FilePath) -> FilePath:
        # test if gzip by trying to read first byte
        with gzip.open(v, "r") as f:
            try:
                f.read(1)
            except gzip.BadGzipFile:
                assert False, "not in gzip format"
        return v


def is_bgzip(p: Path) -> bool:
    # since bgzip is in blocks (vs gzip), determine if in bgzip by
    # attempting to seek first block
    with open(p, "rb") as f:
        try:
            next(bgzf.BgzfBlocks(f), None)
            return True
        except ValueError:
            return False


class RefFileSrc(FileSrc_):
    filepath: FilePath

    @validator("filepath")
    def path_is_bgzip(cls, v: FilePath) -> FilePath:
        assert is_bgzip(v), "not in bgzip format"
        return v


class HttpSrc_(BaseModel):
    url: HttpUrl


class BedHttpSrc(HttpSrc_):
    fmt: BedFmt = BedFmt.GZIP


class RefHttpSrc(HttpSrc_):
    fmt: RefFmt = RefFmt.BGZIP


RefSrc = RefFileSrc | RefHttpSrc

BedSrc = BedFileSrc | BedHttpSrc


class FuncFile(BaseModel):
    src: BedSrc
    chr_prefix: str = ""


class BedFile(BaseModel):
    src: BedSrc
    chr_prefix: str = "chr"
    bed_cols: BedColumns = BedColumns()
    skip_lines: int = 0


class RMSKFile(BedFile):
    class_col: int


class LowComplexity(BaseModel):
    rmsk: RMSKFile
    simreps: BedFile
    satellites: BedFile | None


class XYFeatures(BaseModel):
    x_bed: BedFile
    y_bed: BedFile
    regions: set[XYFeature]


class XYPar(BaseModel):
    start: tuple[int, int]
    end: tuple[int, int]

    def fmt(self, i: ChrIndex, prefix: str) -> str:
        # TODO this smells like something I'll be doing alot
        c = i.chr_name_full(prefix)
        return "\n".join(
            [
                f"{c}\t{self.start[0]}\t{self.start[1]}",
                f"{c}\t{self.end[0]}\t{self.end[1]}",
            ]
        )


class XY(BaseModel):
    features: XYFeatures
    x_par: XYPar
    y_par: XYPar

    def fmt_x_par(self, prefix: str) -> str:
        return self.x_par.fmt(ChrIndex.CHRX, prefix)

    def fmt_y_par(self, prefix: str) -> str:
        return self.y_par.fmt(ChrIndex.CHRY, prefix)


class SegDups(BaseModel):
    superdups: BedFile | None


class Include(BaseModel):
    low_complexity: bool
    xy: bool
    map: bool
    gc: bool
    functional: bool
    segdups: bool


class Build(BaseModel):
    chr_filter: set[ChrIndex]
    include: Include


class RefFile(BaseModel):
    src: RefSrc
    chr_prefix: str


class Functional(BaseModel):
    ftbl: FuncFile
    gff: FuncFile


class Stratification(BaseModel):
    ref: RefFile
    gap: BedFile | None
    low_complexity: LowComplexity
    xy: XY
    segdups: SegDups
    functional: Functional
    builds: dict[BuildKey, Build]


class GiabStrats(BaseModel):
    paths: Paths
    tools: Tools
    stratifications: dict[str, Stratification]

    # hack to make rmd scripts work with this (note this will totally kill
    # the config as it passes into an rmd script)
    def items(self) -> Any:
        return {}.items()

    def refkey_to_strat(self, k: RefKey) -> Stratification:
        return self.stratifications[k]

    # def refkey_to_src(
    #     self, k: RefKey, f: Callable[[Stratification], AnySrc | None]
    # ) -> AnySrc | None:
    #     return f(self.refkey_to_strat(k))

    def refkey_to_ref_src(self, k: RefKey) -> RefSrc:
        return self.refkey_to_strat(k).ref.src

    def refkey_to_gap_src(self, k: RefKey) -> BedSrc | None:
        return fmap_maybe(lambda x: x.src, self.stratifications[k].gap)

    def refkey_to_x_features_src(self, k: RefKey) -> BedSrc | None:
        return fmap_maybe(lambda x: x.features.x_bed.src, self.stratifications[k].xy)

    def refkey_to_y_features_src(self, k: RefKey) -> BedSrc | None:
        return fmap_maybe(lambda x: x.features.x_bed.src, self.stratifications[k].xy)

    # TODO not DRY
    def refkey_to_x_par_bed(self, k: RefKey) -> str:
        prefix = self.refkey_to_final_chr_prefix(k)
        return self.stratifications[k].xy.x_par.fmt(ChrIndex.CHRX, prefix)

    def refkey_to_y_par_bed(self, k: RefKey) -> str:
        prefix = self.refkey_to_final_chr_prefix(k)
        return self.stratifications[k].xy.y_par.fmt(ChrIndex.CHRY, prefix)

    def refkey_to_simreps_src(self, k: RefKey) -> BedSrc | None:
        return fmap_maybe(
            lambda x: x.simreps.src, self.stratifications[k].low_complexity
        )

    def refkey_to_rmsk_src(self, k: RefKey) -> BedSrc | None:
        return fmap_maybe(lambda x: x.rmsk.src, self.stratifications[k].low_complexity)

    def refkey_to_satellite_src(self, k: RefKey) -> BedSrc | None:
        return fmap_maybe(
            lambda x: fmap_maybe(lambda x: x.src, x.satellites),
            self.stratifications[k].low_complexity,
        )

    def refkey_to_superdups_src(self, k: RefKey) -> BedSrc | None:
        return fmap_maybe(
            lambda x: fmap_maybe(lambda x: x.src, x.superdups),
            self.stratifications[k].segdups,
        )

    def refkey_to_functional_ftbl_src(self, k: RefKey) -> BedSrc | None:
        return self.stratifications[k].functional.ftbl.src

    def refkey_to_functional_gff_src(self, k: RefKey) -> BedSrc | None:
        return self.stratifications[k].functional.gff.src

    def refkey_to_final_chr_prefix(self, k: RefKey) -> str:
        return self.stratifications[k].ref.chr_prefix

    def refkey_to_input_chr_prefix(
        self,
        f: Callable[[Stratification], BedFile],
        k: RefKey,
    ) -> str:
        return f(self.stratifications[k]).chr_prefix

    def buildkey_to_init_chr_mapping(
        self,
        f: Callable[[Stratification], BedFile],
        rk: RefKey,
        bk: BuildKey,
    ) -> dict[str, int]:
        p = self.refkey_to_input_chr_prefix(f, rk)
        cs = self.buildkey_to_chr_indices(rk, bk)
        return {c.chr_name_full(p): c.value for c in cs}

    def buildkey_to_final_chr_mapping(self, rk: RefKey, bk: BuildKey) -> dict[int, str]:
        p = self.refkey_to_final_chr_prefix(rk)
        cs = self.buildkey_to_chr_indices(rk, bk)
        return {c.value: c.chr_name_full(p) for c in cs}

    def buildkey_to_build(self, rk: RefKey, bk: BuildKey) -> Build:
        return self.stratifications[rk].builds[bk]

    def buildkey_to_chr_indices(self, rk: RefKey, bk: BuildKey) -> set[ChrIndex]:
        cs = self.stratifications[rk].builds[bk].chr_filter
        return set([x for x in ChrIndex]) if len(cs) == 0 else cs

    def buildkey_to_chr_names(self, rk: RefKey, bk: BuildKey) -> list[str]:
        # TODO don't hardcode this in the future
        prefix = "chr"
        return [i.chr_name_full(prefix) for i in self.buildkey_to_chr_indices(rk, bk)]

    def buildkey_to_chr_pattern(self, rk: RefKey, bk: BuildKey) -> str:
        return "\\|".join(self.buildkey_to_chr_names(rk, bk))

    def want_low_complexity_censat(self, rk: RefKey) -> bool:
        return self.stratifications[rk].low_complexity.satellites is not None

    def _want_chr_index(self, rk: RefKey, bk: BuildKey, i: ChrIndex) -> bool:
        cis = self.buildkey_to_chr_indices(rk, bk)
        return i in cis

    def want_xy_x(self, rk: RefKey, bk: BuildKey) -> bool:
        return (
            self._want_chr_index(rk, bk, ChrIndex.CHRX)
            and self.buildkey_to_build(rk, bk).include.xy
        )

    def want_xy_y(self, rk: RefKey, bk: BuildKey) -> bool:
        return (
            self._want_chr_index(rk, bk, ChrIndex.CHRY)
            and self.buildkey_to_build(rk, bk).include.xy
        )

    def wanted_xy_chr_names(self, rk: RefKey, bk: BuildKey) -> list[str]:
        return [
            i.chr_name
            for i in [ChrIndex.CHRX, ChrIndex.CHRY]
            if self._want_chr_index(rk, bk, i)
        ]

    def want_xy_sex(self, rk: RefKey, bk: BuildKey) -> bool:
        return self.want_xy_x(rk, bk) and self.want_xy_y(rk, bk)

    def want_xy_auto(self, rk: RefKey, bk: BuildKey) -> bool:
        cis = self.buildkey_to_chr_indices(rk, bk)
        return len(cis - set([ChrIndex.CHRX, ChrIndex.CHRY])) > 0

    def want_low_complexity(self, rk: RefKey, bk: BuildKey) -> bool:
        return self.buildkey_to_build(rk, bk).include.low_complexity

    def want_map(self, rk: RefKey, bk: BuildKey) -> bool:
        return self.buildkey_to_build(rk, bk).include.map

    def want_gc(self, rk: RefKey, bk: BuildKey) -> bool:
        return self.buildkey_to_build(rk, bk).include.gc

    def want_functional(self, rk: RefKey, bk: BuildKey) -> bool:
        return self.buildkey_to_build(rk, bk).include.functional

    def want_segdups(self, rk: RefKey, bk: BuildKey) -> bool:
        return self.buildkey_to_build(rk, bk).include.segdups