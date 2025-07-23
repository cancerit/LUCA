# LUCA
#
# Copyright (C) 2024, 2025 Genome Research Ltd.
#
# Author: Luca Barbon
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import abc
from collections import Counter
from contextlib import contextmanager, nullcontext
from dataclasses import dataclass, field
import gc
import json
import logging
import os
from typing import Any, Generator, Generic, IO

import numpy as np

from .app_info import AppInfo
from .experiment import CombinationInfo, CombinationRegionIndices, CombinationStrategy, Experiment, \
    get_unique_library_ids_from_read_templates, \
    LibraryInfo, Options, \
    ReadGroupId as ReadGroupIdEnum, ReadRegion, ReadTemplate, SequencingType
from .fs import get_combination_counts_file_name, get_library_counts_file_name, get_library_stats_file_name, \
    get_match_info_file_name, get_mm_match_info_file_name, get_swap_match_info_file_name, open_opt_output_file, \
    open_output_file
from .library import DynamicMultiTargetLibrary, fill_library_stats, LibraryBuilder, RawLibrary, TargetEncoder, \
    TargetLibrary, write_count
from .library_counts import get_library_counts_array
from .library_independent import write_read_counts
from .library_stats import LibraryStats
from .manifest import OutputFileManifest
from .path_bundle import PathBundle
from .readers.tsv import parse_tsv
from .readers.types import SendT, YieldT
from .stats import LibraryIndependentStats
from .utils import clamp_non_negative, sort_counter_desc

INVALID = -1

# How many values describe the match state of a region
MATCH_INFO_REGION_SIZE = 3

# Match state layout (offsets)
MATCH_INFO_LIBRARY_INDEX = 0
MATCH_INFO_TARGET_INDEX = 1
MATCH_INFO_STATUS = 2


ReadGroupId = int
READ_GROUP_IDS = {x: i for i, x in enumerate(ReadGroupIdEnum)}
READ_GROUP_DEFAULT = READ_GROUP_IDS[ReadGroupIdEnum.DEFAULT]
READ_GROUP_1 = READ_GROUP_IDS[ReadGroupIdEnum.READ_1]
READ_GROUP_2 = READ_GROUP_IDS[ReadGroupIdEnum.READ_2]


CombinationRegionOffsets = dict[ReadGroupId, list[dict[int, int]]]


# BEWARE: the mismatch must be the default value

ReadState = int
READ_STATE_MATCH = 1
READ_STATE_MISMATCH = INVALID
READ_STATE_SWAP = 2


def load_library(pb: PathBundle, library_info: LibraryInfo, library_index: int) -> RawLibrary:
    if library_info.values:

        # Use embedded values
        assert len(library_info.values) > 0
        return RawLibrary(library_index, library_info.reverse_on, library_info.values)

    # Load from file
    fp = pb.get_library_file_path(library_info.id)
    if not os.path.isfile(fp):
        raise FileNotFoundError(fp)
    return RawLibrary.load(fp, library_index, library_info.reverse_on)


def _write_match_info(matches: np.ndarray, fh: IO, row_index: int) -> None:
    fh.write(str(row_index))
    fh.write('\t')
    matches.tofile(fh, sep='\t')
    fh.write('\n')


def write_match_info(matches: np.ndarray, fh: IO | None, mm_fh: IO | None, swap_fh: IO | None, row_index: int, match_state: ReadState) -> None:
    """
    Write match information to file

    Reads where at least a swap event has been identified are always included
    in the mismatch file, and optionally copied in a separate file, for
    consistency with respect to the no swap test scenario (in which such reads
    would have been considered mismatching).
    """

    if fh is not None:

        # Write to matches file
        _write_match_info(matches, fh, row_index)

    if match_state == READ_STATE_MATCH:
        return

    if swap_fh and match_state == READ_STATE_SWAP:

        # Write to swap file
        _write_match_info(matches, swap_fh, row_index)

    # Write to mismatch file
    # TODO: count for logging purposes, and drop the file if empty?
    if mm_fh is not None:
        _write_match_info(matches, mm_fh, row_index)


@dataclass
class LibraryIndependentCounter:
    library_index: int
    library: DynamicMultiTargetLibrary
    counter: Counter[int]

    @classmethod
    def empty(cls, library_index: int, target_index_offset: int):
        return cls(
            library_index,
            DynamicMultiTargetLibrary(target_index_offset),
            Counter())

    @property
    def total_dynamic_targets(self) -> int:
        return len(self.counter)

    @property
    def is_empty(self) -> bool:
        return len(self.library) == 0

    def match_full(self, seq: str) -> tuple[int, int]:
        target_index, match_length = self.library.match_full(seq)
        self.counter[target_index] += 1
        return target_index, match_length

    def match(self, seq: str, offset: int, length: int, anchor_left: bool = True) -> tuple[int, int]:
        target_index, match_length = self.library.match(
            seq, offset=offset, length=length, anchor_left=anchor_left)
        self.counter[target_index] += 1
        return target_index, match_length

    def write_counts(self, fp: str) -> None:
        with open_output_file(fp) as f:
            fh = f.fh
            for target_index, count in self.counter.items():
                target = self.library.decode(target_index)
                write_count(fh, target_index, target, count)


@dataclass(slots=True)
class DataBundle:
    # TODO: add a constructor per setup?
    library_counts: np.ndarray
    library_stats: dict[int, LibraryStats]
    library_indep_counter: LibraryIndependentCounter
    combination_counts: list[Counter]
    mm_read_counts: Counter


@dataclass(slots=True, frozen=True)
class BaseRegion(abc.ABC):
    skip: int = 0

    @abc.abstractmethod
    def match(
        self,
        libraries: dict[int, TargetLibrary],
        bundle: DataBundle,
        matches: np.ndarray,
        seq: str,
        i: int,
        region_offset: int,
        anchor_left: bool,
        swap: bool = True
    ) -> tuple[int, int, bool]:
        pass

    def get_start(self, region_offset: int) -> int:
        return clamp_non_negative(region_offset + self.skip)


@dataclass(slots=True, frozen=True)
class Region(BaseRegion):
    libraries: list[int] = field(default_factory=list)
    max_offset: int = 0
    swap_libraries: list[int] = field(default_factory=list)

    def __post_init__(self) -> None:
        assert self.libraries

    @classmethod
    def from_read_region(cls, libraries: dict[str, int], t: ReadRegion):
        return cls(
            libraries=[
                libraries[x]
                for x in t.libraries
            ],
            skip=t.skip,
            max_offset=t.max_offset,
            swap_libraries=[
                libraries[x]
                for x in t.swap_libraries
            ])

    def match(
        self,
        libraries: dict[int, TargetLibrary],
        bundle: DataBundle,
        matches: np.ndarray,
        seq: str,
        i: int,
        region_offset: int,
        anchor_left: bool,
        swap: bool = True,
    ) -> tuple[int, int, bool]:
        match_length: int = 0
        is_swap: bool = False
        match_region_offset = region_offset

        # Try libraries
        for library_index in self.libraries:
            target_index, match_length, match_region_offset = match_read_region(
                seq, libraries[library_index], self, region_offset, anchor_left)
            if match_length > 0:
                update_lib_counts_match_info(
                    bundle.library_counts,
                    matches,
                    i,
                    library_index,
                    target_index,
                    READ_STATE_MATCH)
                break

        if match_length == 0 and swap:

            # Try alternate libraries (swap detection)
            for library_index in self.swap_libraries:
                target_index, match_length, match_region_offset = match_read_region(
                    seq, libraries[library_index], self, region_offset, anchor_left)
                if match_length > 0:
                    # NOTE: out-of-target matches contribute to the same count
                    update_lib_counts_match_info(
                        bundle.library_counts,
                        matches,
                        i,
                        library_index,
                        target_index,
                        READ_STATE_SWAP)
                    bundle.library_stats[library_index].swap_matching_reads += 1
                    is_swap = True
                    break

        return match_region_offset, match_length, is_swap


def match_lib_indep(
    bundle: DataBundle,
    matches: np.ndarray,
    seq: str,
    i: int,
    region_offset: int,
    region_length: int,
    anchor_left: bool
) -> tuple[int, bool]:
    # TODO: verify behaviour on out-of-bounds
    target_index, match_length = bundle.library_indep_counter.match(
        seq, region_offset, region_length, anchor_left=anchor_left)

    if match_length > 0:
        update_match_info(
            matches,
            i,
            bundle.library_indep_counter.library_index,
            target_index)

    return match_length, False


@dataclass(slots=True, frozen=True)
class LibraryIndependentRegion(BaseRegion):
    def match(
        self,
        libraries: dict[int, TargetLibrary],
        bundle: DataBundle,
        matches: np.ndarray,
        seq: str,
        i: int,
        region_offset: int,
        anchor_left: bool,
        swap: bool = True
    ) -> tuple[int, int, bool]:
        start = self.get_start(region_offset)
        k = len(seq) - start
        return start, *match_lib_indep(
            bundle, matches, seq, i, start, k, anchor_left)


@dataclass(slots=True, frozen=True)
class FixedLengthRegion(BaseRegion):
    length: int = 0

    def match(
        self,
        libraries: dict[int, TargetLibrary],
        bundle: DataBundle,
        matches: np.ndarray,
        seq: str,
        i: int,
        region_offset: int,
        anchor_left: bool,
        swap: bool = True
    ) -> tuple[int, int, bool]:
        k = self.length
        start = self.get_start(region_offset)
        return start, *match_lib_indep(
            bundle, matches, seq, i, start, k, anchor_left)


def build_region(t: ReadRegion, libraries: dict[str, int] | None = None) -> BaseRegion:
    if t.libraries:
        assert libraries
        return Region.from_read_region(libraries, t)
    elif t.length is not None:
        assert t.length > 0
        return FixedLengthRegion(length=t.length, skip=t.skip)
    else:
        return LibraryIndependentRegion(skip=t.skip)


def match_read_region(seq: str, library: TargetLibrary, region: Region, region_offset: int, anchor_left: bool) -> tuple[int, int, int]:
    start_region_offset = clamp_non_negative(region_offset + region.skip)
    for offset in range(region.max_offset + 1):
        matching_region_offset = start_region_offset + offset
        target_index, match_length = library.match(
            seq, offset=matching_region_offset, anchor_left=anchor_left)
        if match_length > 0:
            return target_index, match_length, matching_region_offset
    return 0, 0, region_offset


@dataclass(slots=True, frozen=True)
class Template:
    anchor_left: bool
    regions: list[BaseRegion]

    def __post_init__(self) -> None:
        assert self.regions
        assert sum([
            1 if isinstance(r, LibraryIndependentRegion) else 0
            for r in self.regions
        ]) <= 1

    @property
    def qc_info_slots(self) -> int:
        return len(self.regions)

    @classmethod
    def from_read_template(cls, libraries: dict[str, int], is_reverse: bool, t: ReadTemplate):
        # anchor_left = (t.anchor == Anchor.LEFT) or (t.anchor == Anchor.AUTO and not is_reverse)
        anchor_left = t.anchor.get_anchor_left(is_reverse)
        return cls(anchor_left, [
            build_region(region, libraries=libraries)
            for region in t.regions
        ])

    def match(self, libraries: dict[int, TargetLibrary], bundle: DataBundle, matches: np.ndarray, seq: str, swap: bool = True) -> tuple[int, bool]:
        region_offset = 0
        match_length = 0
        has_swap = False

        for i, region in enumerate(self.regions):
            match_region_offset, match_length, is_swap = region.match(
                libraries,
                bundle,
                matches,
                seq,
                i,
                region_offset,
                self.anchor_left,
                swap=swap)

            if is_swap:
                has_swap = True

            if match_length == 0:
                # Do not try to match anything past a mismatch
                break

            region_offset = match_region_offset + match_length

        return match_length, has_swap


def build_template(libraries: dict[str, int], is_reverse: bool, t: ReadTemplate) -> Template:
    return Template.from_read_template(libraries, is_reverse, t)


@dataclass(slots=True)
class ExperimentStats:
    read_counter: Counter[ReadState] = field(init=False, default_factory=Counter)

    # TODO: consider moving this to a paired-end-specific subclass
    pair_one_end_match: int = field(init=False, default=0)

    def to_dict(self, app_info: AppInfo) -> dict[str, Any]:
        return {
            **app_info.to_dict(),
            'read_match': self.read_counter.get(READ_STATE_MATCH, 0),
            'read_swap': self.read_counter.get(READ_STATE_SWAP, 0),
            'read_mismatch': self.read_counter.get(READ_STATE_MISMATCH, 0),
            'pair_one_end_match': self.pair_one_end_match
        }


def update_match_info(matches: np.ndarray, region_index: int, library_index: int, target_index: int) -> None:
    j = region_index * MATCH_INFO_REGION_SIZE
    matches[j] = library_index
    matches[j + 1] = target_index
    matches[j + 2] = READ_STATE_MATCH


def update_lib_counts_match_info(library_counts: np.ndarray, matches: np.ndarray, region_index: int, library_index: int, target_index: int, read_state: ReadState) -> None:
    library_counts[target_index] += 1
    j = region_index * MATCH_INFO_REGION_SIZE
    matches[j] = library_index
    matches[j + 1] = target_index
    matches[j + 2] = read_state


@dataclass(slots=True, frozen=True)
class BaseMatcher(abc.ABC):

    @property
    @abc.abstractmethod
    def match_info_field_count(self) -> int:
        pass

    @abc.abstractmethod
    def match_read_seq(self, bundle: DataBundle, matches: np.ndarray, seq: str, swap: bool = True) -> tuple[ReadState, int]:
        pass

    def get_match_info_array(self) -> np.ndarray:
        return np.empty(self.match_info_field_count * MATCH_INFO_REGION_SIZE, dtype=np.int64)


# @dataclass(slots=True, frozen=True)
# class SimpleMatcher(BaseMatcher):
#     library_index: int
#     library: TargetLibrary
#     region: Region
#     anchor_left: bool

#     @property
#     def match_info_field_count(self) -> int:
#         return 1

#     def __post_init__(self) -> None:
#         assert not self.region.swap_libraries

#     def match_read_seq(self, bundle: DataBundle, matches: np.ndarray, seq: str, swap: bool = True) -> tuple[ReadState, int]:
#         target_index, match_length = match_read_region(seq, self.library, self.region, 0, self.anchor_left)
#         if match_length == 0:
#             return READ_STATE_MISMATCH, 0
#         else:
#             update_lib_counts_match_info(
#                 bundle.library_counts, matches, 0, self.library_index, target_index, READ_STATE_MATCH)
#             return READ_STATE_MATCH, 0


@dataclass(slots=True, frozen=True)
class Matcher(BaseMatcher):
    libraries: dict[int, TargetLibrary]
    read_templates: list[Template]

    @property
    def match_info_field_count(self) -> int:
        return sum([
            read_template.qc_info_slots
            for read_template in self.read_templates
        ])

    def match_read_seq(self, bundle: DataBundle, matches: np.ndarray, seq: str, swap: bool = True) -> tuple[ReadState, int]:
        # Assumption: the matches array was reset to the invalid value!
        for template_index, read_template in enumerate(self.read_templates):

            match_length, has_swap = read_template.match(
                self.libraries, bundle, matches, seq, swap=swap)
            if has_swap:
                return READ_STATE_SWAP, template_index
            if match_length > 0:
                return READ_STATE_MATCH, template_index

        return READ_STATE_MISMATCH, 0


def build_matcher(exp: Experiment, read_group_id: ReadGroupIdEnum, library_builder: LibraryBuilder) -> BaseMatcher:
    if read_group_id not in exp.read_group_templates:
        logging.error("Read group ID not found: '%s'!" % read_group_id)
        raise ValueError

    read_templates = exp.get_read_group_templates(read_group_id)
    library_indices = exp.get_library_indices(
        get_unique_library_ids_from_read_templates(read_templates))

    # Build the strictly required libraries
    # Note: because being reverse complement or not is currently a property of the
    #  read group (vs. the template), the libraries may be pre-processed accordingly
    #  once (vs. having to potentially keep both forms of some libraries).
    read_group = exp.read_group_infos[read_group_id]
    # TODO: resolve per template vs. per read anchoring... may have to build both libraries is some cases
    is_reverse = read_group.is_reverse
    libraries = {
        library_index: library_builder.build(library_index, is_reverse=is_reverse)
        for library_index in library_indices.values()
    }

    # Encode library ID's into library indices
    templates = [
        build_template(library_indices, read_group.is_reverse, read_template)
        for read_template in read_templates
    ]

    # if len(templates) == 1:
    #     tpl = templates[0]
    #     if len(tpl.regions) == 1:
    #         region = tpl.regions[0]
    #         if (
    #             isinstance(region, Region) and
    #             len(region.libraries) == 1 and
    #             not region.swap_libraries
    #         ):
    #             library_index = region.libraries[0]
    #             logging.info("Building simple matcher")
    #             return SimpleMatcher(is_rc, library_index, libraries[library_index], region)

    return Matcher(libraries, templates)


def get_combination_region_offsets(
    combination_region_indices: CombinationRegionIndices,
    config: CombinationInfo
) -> CombinationRegionOffsets:
    """
    For each read group, for each template (by list index), map region indices
    to match info offsets

    The offset is relative to the in-memory representation, not the one written
    to disk (that may include the read index).
    """

    # map: read group ID -> region ID's
    read_group_regions = config.grouped_regions

    return {
        READ_GROUP_IDS[read_group_id]: [
            {
                # Convert the index of the region in the template into the offset
                #  of the matching target index in the match information array
                tpl_index: (tpl_regions[region_id] * MATCH_INFO_REGION_SIZE) + MATCH_INFO_TARGET_INDEX
                for tpl_index, tpl_regions in enumerate(combination_region_indices[read_group_id])
            }
            for region_id in read_group_regions[read_group_id]
        ]
        for read_group_id in combination_region_indices.keys()
    }


@dataclass
class BaseCombinationRules(abc.ABC):
    region_offsets: CombinationRegionOffsets

    @property
    def region_number(self) -> int:
        return len(self.region_offsets)

    def get_region_offsets(self, read_group_id: ReadGroupId, template_index: int) -> list[int]:
        return [r[template_index] for r in self.region_offsets[read_group_id]]

    @abc.abstractmethod
    def get_counter(self) -> Counter:
        pass

    @abc.abstractmethod
    def should_count_combination(self, combination: tuple) -> bool:
        pass


@dataclass
class CombinationRules(BaseCombinationRules):

    @classmethod
    def from_config(
        cls,
        config: CombinationInfo,
        combination_region_indices: CombinationRegionIndices
    ):
        region_offsets = get_combination_region_offsets(
            combination_region_indices, config)
        return cls(region_offsets)

    def should_count_combination(self, combination: tuple) -> bool:
        return True

    def get_counter(self) -> Counter:
        return Counter()


def build_combination_filter(
    config: CombinationInfo,
    pb: PathBundle,
    encoder: TargetEncoder
) -> set[tuple]:
    assert config.filters
    combinations: set[tuple] = set()

    # Parse valid combination files
    field_names = config.filter_field_names
    assert field_names
    for ft in config.filters:
        fp = pb.get_combination_filter_file_path(ft)
        assert os.path.isfile(fp)
        # NOTE: this is favouring the configuration order vs the TSV order
        for sequences in parse_tsv(fp, field_names):
            combinations.add(tuple([
                encoder.encode(sequence)
                for sequence in sequences
            ]))

    return combinations


@dataclass
class NCCombinationRules(BaseCombinationRules):
    combinations: set[tuple]

    @classmethod
    def from_config(
        cls,
        config: CombinationInfo,
        pb: PathBundle,
        combination_region_indices: CombinationRegionIndices,
        encoder: TargetEncoder
    ):
        combinations = build_combination_filter(
            config, pb, encoder)
        region_offsets = get_combination_region_offsets(
            combination_region_indices, config)

        return cls(region_offsets, combinations)

    def should_count_combination(self, combination: tuple) -> bool:
        return combination in self.combinations

    def get_counter(self) -> Counter:
        return Counter({k: 0 for k in self.combinations})


@dataclass
class PartialNCCombinationRules(BaseCombinationRules):
    filter_combinations: set[tuple]
    filter_region_indices: list[int]

    @classmethod
    def from_config(
        cls,
        config: CombinationInfo,
        pb: PathBundle,
        combination_region_indices: CombinationRegionIndices,
        encoder: TargetEncoder
    ):
        combinations = build_combination_filter(
            config, pb, encoder)
        region_offsets = get_combination_region_offsets(
            combination_region_indices, config)

        return cls(region_offsets, combinations, config.filtered_region_indices)

    def should_count_combination(self, combination: tuple) -> bool:
        return tuple([
            combination[i]
            for i in self.filter_region_indices
        ]) in self.filter_combinations

    def get_counter(self) -> Counter:
        return Counter()


def build_combination_rules(
    exp: Experiment, pb: PathBundle,
    library_builder: LibraryBuilder | None = None
) -> tuple[list[BaseCombinationRules], TargetEncoder | None]:

    if not exp.combinations:
        return [], TargetEncoder()

    # Map regions to indices to map them to match arrays
    combination_region_indices = exp.get_combination_region_indices()

    def get_combinatorial(ft: CombinationInfo) -> CombinationRules:
        return CombinationRules.from_config(
            ft, combination_region_indices)

    # A. No filtering
    if not exp.has_combination_filters:
        return [
            get_combinatorial(ft)
            for ft in exp.combinations
        ], None

    # B. Partial or full filtering
    assert library_builder
    encoder = library_builder.build_encoder()

    def get_combination_rules(ft: CombinationInfo) -> BaseCombinationRules:
        match ft.combination_strategy:

            case CombinationStrategy.COMBINATORIAL:
                return get_combinatorial(ft)

            case CombinationStrategy.NON_COMBINATORIAL:
                return NCCombinationRules.from_config(
                    ft, pb, combination_region_indices, encoder)

            case CombinationStrategy.PARTIAL_NON_COMBINATORIAL:
                return PartialNCCombinationRules.from_config(
                    ft, pb, combination_region_indices, encoder)

    return list(map(get_combination_rules, exp.combinations)), encoder


@dataclass
class MultiMatcher(Generic[YieldT, SendT], abc.ABC):
    libraries: LibraryBuilder
    matchers: dict[ReadGroupId, BaseMatcher]
    combination_rules: list[BaseCombinationRules] = field(default_factory=list)
    encoder: TargetEncoder | None = None

    @classmethod
    def from_experiment(cls, exp: Experiment, pb: PathBundle):

        library_indices = exp.assigned_library_indices
        if library_indices:

            # Load libraries assigned to at least one read template
            library_builder = LibraryBuilder.from_raw_libraries([
                load_library(pb, exp.libraries[library_index], library_index)
                for library_index in library_indices.values()
            ])

        else:
            library_builder = LibraryBuilder.empty()

        # Optionally load the valid combinations
        combination_rules, encoder = build_combination_rules(
            exp, pb, library_builder=library_builder)

        return cls(library_builder, {
            READ_GROUP_IDS[read_group_id]: build_matcher(exp, read_group_id, library_builder)
            for read_group_id in exp.read_group_templates
        }, combination_rules=combination_rules, encoder=encoder)

    def get_library_size(self, library_index: int) -> int:
        return self.libraries.get_library_size(library_index)

    @property
    def library_sizes(self) -> dict[int, int]:
        return self.libraries.library_sizes

    @abc.abstractmethod
    def _match(
        self,
        opt: Options,
        seqs: Generator[str, bool, Any],
        out_dir: str,
        exp_stats: ExperimentStats,
        bundle: DataBundle
    ) -> None:
        pass

    def get_data_bundle(self) -> DataBundle:
        library_counts = get_library_counts_array(sum(self.library_sizes.values()))

        library_indices = self.libraries.library_indices
        library_stats = {
            library_index: LibraryStats.empty()
            for library_index in library_indices
        }

        mm_read_counts = Counter()
        combination_counts = [cr.get_counter() for cr in self.combination_rules]
        target_n = library_counts.shape[0]

        lic_index = (max(library_indices) + 1) if library_indices else 0
        library_indep_counter = LibraryIndependentCounter.empty(
            lic_index,
            target_n)

        return DataBundle(
            library_counts,
            library_stats,
            library_indep_counter,
            combination_counts,
            mm_read_counts)

    def match_seqs(self, opt: Options, seqs: Generator[YieldT, SendT, Any], out_dir: str, profile: bool = False) -> tuple[LibraryIndependentStats, ExperimentStats]:
        # TODO: consider whether the library counts should be per mate, per region, or global...
        stats: LibraryIndependentStats | None = None
        exp_stats = ExperimentStats()
        bundle = self.get_data_bundle()

        if profile:
            from cProfile import Profile
            ctx = Profile()
        else:
            ctx = nullcontext()

        with ctx:
            try:
                self._match(opt, seqs, out_dir, exp_stats, bundle)
            except StopIteration as ex:
                stats = ex.value

            if profile:
                from pstats import SortKey, Stats
                Stats(ctx).strip_dirs().sort_stats(SortKey.CALLS).dump_stats('profile.dmp')

        assert stats
        manifest = OutputFileManifest.empty(sample_name=stats.sample_name)
        manifest.total_library_templates = bundle.library_counts.shape[0]

        # NOTE: the target sequences are expressed in whichever orientation
        #  they appeared in the original (raw) library file.
        for library_index, library_stats in bundle.library_stats.items():
            library = self.libraries.libraries[library_index]

            start, end = self.libraries.get_library_target_index_range(library_index)
            counts = bundle.library_counts[start:end]

            # Write library counts
            # TODO: let the user specify which libraries should have their counts dumped
            fnc = get_library_counts_file_name(library_index)
            fpc = os.path.join(out_dir, fnc)
            logging.info(f"Writing counts file: {fpc}")
            library.write_counts(fpc, counts, target_index_offset=start)

            # Write library stats
            fns = get_library_stats_file_name(library_index)
            fps = os.path.join(out_dir, fns)
            # BEWARE: the following function SORTS the counts IN PLACE
            fill_library_stats(library_stats, counts)
            logging.info(f"Writing statistics file: {fps}")
            with open(fps, 'w') as fh:
                json.dump(library_stats.to_dict(), fh)

            manifest.push_library_file_paths(library_index, fns, fnc)

        del bundle.library_counts

        # Handle library-independent targets (if any)
        lic = bundle.library_indep_counter
        if not lic.is_empty:
            logging.info("%d library-independent targets collected" % len(lic.library))

            # Write library-independent region counts
            fn = get_library_counts_file_name(lic.library_index)
            fp = os.path.join(out_dir, fn)
            logging.info(f"Writing library-independent region counts: {fp}")
            lic.write_counts(fp)

            # TODO: should the stats the computed?
            manifest.set_dynamic_target_file_paths(stats="", counts=fn)
            manifest.total_dynamic_targets = lic.total_dynamic_targets

        if self.combination_rules:

            if self.encoder:
                decoder = self.encoder.to_decoder()
                del self.encoder
                gc.collect()
            else:
                decoder = self.libraries.build_decoder()

            if not lic.is_empty:
                # Update decoder with library-independent targets
                # NOTE: duplicate targets are possible
                lic.library.update_decoder(decoder)

            for ci, cr in enumerate(self.combination_rules):

                # Write combination counts
                fn = get_combination_counts_file_name(ci)
                fp = os.path.join(out_dir, fn)
                logging.info(f"Writing combination counts file: {fp}")
                with open(fp, 'w') as fh:
                    for combination, count in bundle.combination_counts[ci].items():
                        for item in combination:
                            fh.write(decoder[item])
                            fh.write('\t')
                        fh.write(str(count))
                        fh.write('\n')

                manifest.push_combination_file_paths(cr.region_number, fn)

            del bundle.combination_counts, decoder

        if opt.count_mm_reads and bundle.mm_read_counts:
            if opt.sort_mm_read_counts:
                gc.collect()
                counts = sort_counter_desc(bundle.mm_read_counts)
            else:
                counts = bundle.mm_read_counts

            # Write mismatching read counts
            fp = write_read_counts(counts, out_dir, prefix='mm', compress=opt.compress_mm_read_counts)

            manifest.set_library_independent_mm_file_paths(fp, opt.sort_mm_read_counts, opt.compress_mm_read_counts)

        manifest.write(out_dir)
        return stats, exp_stats

    @contextmanager
    def open_match_info_file(self, opt: Options, out_dir: str, matcher_index: int = 0):
        with open_opt_output_file(
            get_match_info_file_name(matcher_index),
            out_dir=out_dir,
            no_op=not opt.full_match_info
        ) as f:
            yield f.fh if f else None

    @contextmanager
    def open_swap_match_info_file(self, opt: Options, out_dir: str, matcher_index: int = 0):
        with open_opt_output_file(
            get_swap_match_info_file_name(matcher_index),
            out_dir=out_dir,
            no_op=not opt.test_swaps
        ) as f:
            yield f.fh if f else None

    @contextmanager
    def open_mm_match_info_file(self, opt: Options, out_dir: str, matcher_index: int = 0):
        with open_opt_output_file(
            get_mm_match_info_file_name(matcher_index),
            out_dir=out_dir,
            no_op=not opt.out_mm_match_info
        ) as f:
            yield f.fh if f else None


def get_combination_target_indices(
    cr: BaseCombinationRules,
    read_group_id: ReadGroupId,
    template_index: int,
    matches: np.ndarray
) -> list[int]:
    return [
        int(matches[offset])
        for offset in cr.get_region_offsets(read_group_id, template_index)
    ]


@dataclass
class SingleEndMultiMatcher(MultiMatcher[str, bool]):

    def _match(
        self,
        opt: Options,
        seqs: Generator[str, bool, Any],
        out_dir: str,
        exp_stats: ExperimentStats,
        bundle: DataBundle
    ) -> None:
        m = self.matchers[READ_GROUP_DEFAULT]
        matches = m.get_match_info_array()
        test_swaps = opt.test_swaps
        count_mm_reads = opt.count_mm_reads
        out_match_info = opt.should_write_any_match_info

        with (
            self.open_match_info_file(opt, out_dir) as fh,
            self.open_mm_match_info_file(opt, out_dir) as mm_fh,
            self.open_swap_match_info_file(opt, out_dir) as swap_fh
        ):
            i = 0
            seq = next(seqs)
            while seq:
                matches.fill(INVALID)
                match_state, template_index = m.match_read_seq(bundle, matches, seq, swap=test_swaps)

                # NOTE: swaps are not counted as mismatches here (to consider)
                exp_stats.read_counter[match_state] += 1

                if out_match_info:
                    write_match_info(matches, fh, mm_fh, swap_fh, i, match_state)

                if self.combination_rules and match_state == READ_STATE_MATCH:
                    for ci, cr in enumerate(self.combination_rules):
                        combination = tuple(get_combination_target_indices(
                            cr, READ_GROUP_DEFAULT, template_index, matches))
                        if cr.should_count_combination(combination):
                            bundle.combination_counts[ci][combination] += 1

                is_not_match = match_state != READ_STATE_MATCH

                # Count mismatching read
                if count_mm_reads and is_not_match:
                    bundle.mm_read_counts[seq] += 1

                # Feed the mismatch status back to the generator and step forward
                # DISCUSS: have a separate output BAM file for swaps?
                seq = seqs.send(is_not_match)

                i += 1


@dataclass
class PairedEndMultiMatcher(MultiMatcher[tuple[str, str], tuple[bool, bool]]):
    def _match(
        self,
        opt: Options,
        seqs: Generator[str, tuple[bool, bool], Any],
        out_dir: str,
        exp_stats: ExperimentStats,
        bundle: DataBundle
    ) -> None:
        m1 = self.matchers[READ_GROUP_1]
        m2 = self.matchers[READ_GROUP_2]
        matches1 = m1.get_match_info_array()
        matches2 = m2.get_match_info_array()
        test_swaps = opt.test_swaps
        count_mm_reads = opt.count_mm_reads
        out_match_info = opt.should_write_any_match_info

        with (
            self.open_match_info_file(opt, out_dir, 1) as fh1,
            self.open_mm_match_info_file(opt, out_dir, 1) as mm_fh1,
            self.open_match_info_file(opt, out_dir, 2) as fh2,
            self.open_mm_match_info_file(opt, out_dir, 2) as mm_fh2,
            self.open_swap_match_info_file(opt, out_dir, 1) as swap_fh1,
            self.open_swap_match_info_file(opt, out_dir, 2) as swap_fh2
        ):
            i = 0
            pair = next(seqs)
            while pair:
                matches1.fill(INVALID)
                matches2.fill(INVALID)
                r1, r2 = pair
                match_state1, template_index1 = m1.match_read_seq(bundle, matches1, r1, swap=test_swaps)
                match_state2, template_index2 = m2.match_read_seq(bundle, matches2, r2, swap=test_swaps)

                # NOTE: swaps are not counted as mismatches here (to consider)
                exp_stats.read_counter[match_state1] += 1
                exp_stats.read_counter[match_state2] += 1

                if out_match_info:
                    write_match_info(matches1, fh1, mm_fh1, swap_fh1, i, match_state1)
                    write_match_info(matches2, fh2, mm_fh2, swap_fh2, i, match_state2)

                # TODO: in theory, the combinations could reflect only one of the groups
                if (
                    self.combination_rules and
                    match_state1 == READ_STATE_MATCH and
                    match_state2 == READ_STATE_MATCH
                ):
                    for ci, cr in enumerate(self.combination_rules):
                        combination_values = []

                        # NOTE: this is respecting the order in which the two read groups
                        #  are presented in the configuration.
                        # TODO: allow interleaving of regions from different read gropus
                        for read_group_id in cr.region_offsets.keys():
                            combination_values += (
                                get_combination_target_indices(cr, read_group_id, template_index1, matches1) if read_group_id == READ_GROUP_1 else
                                get_combination_target_indices(cr, read_group_id, template_index2, matches2)
                            )

                        combination = tuple(combination_values)
                        # TODO: respect the region order from the filter table?
                        if cr.should_count_combination(combination):
                            bundle.combination_counts[ci][combination] += 1

                is_not_match1 = match_state1 != READ_STATE_MATCH
                is_not_match2 = match_state2 != READ_STATE_MATCH

                # Test for one-end matches
                if is_not_match1 != is_not_match2:
                    exp_stats.pair_one_end_match += 1

                # Count mismatching reads
                if count_mm_reads:
                    if is_not_match1:
                        bundle.mm_read_counts[r1] += 1
                    if is_not_match2:
                        bundle.mm_read_counts[r2] += 1

                # Feed the mismatch status back to the generator and step forward
                pair = seqs.send((is_not_match1, is_not_match2))

                i += 1


def get_multi_matcher(exp: Experiment, pb: PathBundle) -> MultiMatcher:
    match exp.sequencing_type:
        case SequencingType.SINGLE_END:
            cls = SingleEndMultiMatcher
        case SequencingType.PAIRED_END:
            cls = PairedEndMultiMatcher
        case _:
            raise NotImplementedError("Unsupported sequencing type!")

    return cls.from_experiment(exp, pb)
