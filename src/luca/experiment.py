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

import logging
from collections import defaultdict
from enum import Enum, IntEnum
from itertools import chain, groupby
from typing import Any

from pydantic import BaseModel, ConfigDict, Field, model_validator

from .errors import InvalidExperiment
from .utils import eval_fail_conditions, get_offsets, get_set, greedy_all, greedy_any, has_duplicates, union


class ReadGroupId(Enum):
    DEFAULT = 'default'
    READ_1 = 'read_1'
    READ_2 = 'read_2'


def validate_read_group_ids(ids: set[ReadGroupId]) -> None:
    # TODO: consider the empty set
    if ReadGroupId.DEFAULT in ids:
        if len(ids) > 1:
            raise ValueError(
                f"Read group '{ReadGroupId.DEFAULT.value}' is mutually exclusive with any other!")


class LibraryReverseCondition(Enum):
    REVERSE_GROUP = 'reverse_group'
    ALWAYS = 'always'
    NEVER = 'never'


class BaseConfig(BaseModel, validate_assignment=True):
    model_config = ConfigDict(extra='forbid')  # noqa: F841


class LibraryInfo(BaseConfig):
    id: str
    values: list[str] | None = None
    reverse_on: LibraryReverseCondition = LibraryReverseCondition.REVERSE_GROUP

    @property
    def from_file(self) -> bool:
        return not bool(self.values)


class ReadRegion(BaseConfig):
    id: str
    libraries: list[str] = list()
    skip: int = 0
    max_offset: int = 0
    swap_libraries: list[str] = list()
    length: int | None = None

    @property
    def is_library_indep(self) -> bool:
        return not self.libraries

    @property
    def is_library_indep_unbound(self) -> bool:
        return self.is_library_indep and self.length is None

    def _raise_invalid_region(self, tpl_id: str, msg: str) -> None:
        logging.error(
            "Invalid read region '%s' in read template '%s': %s!" %
            (self.id, tpl_id, msg))

    def validate(self, tpl_id: str) -> bool:
        return eval_fail_conditions([
            (
                (self.is_library_indep and bool(self.swap_libraries)),
                "swap libraries set for library-independent region"
            ),
            (
                (self.length is not None and self.length < 1),  # type: ignore
                "the length is less than one"
            ),
            (
                ((self.length is not None or self.is_library_indep) and self.max_offset != 0),
                "library-independent quantification does not support a nonzero maximum offset"
            )
        ], lambda msg: self._raise_invalid_region(tpl_id, msg))


class SequencingType(Enum):
    SINGLE_END = 'single_end'
    PAIRED_END = 'paired_end'


class CombinationStrategy(IntEnum):
    COMBINATORIAL = 0
    NON_COMBINATORIAL = 1
    PARTIAL_NON_COMBINATORIAL = 2


class ExecMode(IntEnum):
    LIBRARY_INDEPENDENT = 0
    LIBRARY_DEPENDENT = 1


class Anchor(Enum):
    AUTO = 'auto'
    LEFT = 'left'
    RIGHT = 'right'

    def get_anchor_left(self, is_reverse: bool) -> bool:
        return (self == Anchor.LEFT) or (self == Anchor.AUTO and not is_reverse)


CombinationRegionIndices = dict[ReadGroupId, list[dict[str, int]]]


class ReadTemplate(BaseConfig):
    #  NOTE: Have reverse complement be a property
    #   of the template, for maximum flexibility
    #  (when both orientation are to be evaluated);
    #  the state of the read group may represent the default
    #  in that regard, and the template may override it.
    id: str
    regions: list[ReadRegion]
    anchor: Anchor = Anchor.AUTO

    def model_post_init(self, __context: Any) -> None:
        # TODO: consider whether this validation logic should be moved
        n = len(self.regions)
        assert n > 0
        assert not has_duplicates([r.id for r in self.regions])
        if n > 1:
            for r in self.regions[:-1]:
                assert not r.is_library_indep_unbound
        return super().model_post_init(__context)

    @property
    def library_indep_region_ids(self) -> set[str]:
        return {
            region.id
            for region in self.regions
            if region.is_library_indep
        }

    @property
    def library_ids(self) -> set[str]:
        library_ids = [
            set(region.libraries) | set(region.swap_libraries)
            for region in self.regions
        ]
        return union(library_ids) if library_ids else set()

    def get_region_indices(self, offset: int = 0) -> dict[str, int]:
        return {r.id: offset + i for i, r in enumerate(self.regions)}

    def validate(self) -> bool:
        return greedy_all(lambda x: x.validate(self.id), self.regions)


class ReadGroupOptions(BaseConfig):
    is_reverse: bool


class ReadGroupInfo(BaseConfig):
    """
    Information to classify reads in order to match them to templates

    mate: read 1 or 2, paired-end only
    read_groups: identifiers of the read groups to include (detect conflicts!)
    """

    id: ReadGroupId
    is_reverse: bool = False

    def update(self, opt: ReadGroupOptions) -> None:
        self.is_reverse = opt.is_reverse


class CombinationRegion(BaseConfig):
    id: str
    read_group: ReadGroupId = ReadGroupId.DEFAULT
    filter: bool = False

    @property
    def field_name(self) -> str:
        return (
            self.id if self.read_group == ReadGroupId.DEFAULT else
            f"{self.read_group.value}.{self.id}"
        )


class CombinationInfo(BaseConfig):
    id: str
    regions: list[CombinationRegion]
    filters: list[str] | None = None

    @property
    def combination_strategy(self) -> CombinationStrategy:
        return (
            CombinationStrategy.COMBINATORIAL if not self.filters else
            CombinationStrategy.PARTIAL_NON_COMBINATORIAL if self.has_partial_filtering else
            CombinationStrategy.NON_COMBINATORIAL
        )

    @property
    def has_filter(self) -> bool:
        return bool(self.filters)

    @property
    def region_ids(self) -> list[str]:
        return [r.id for r in self.regions]

    @property
    def has_partial_filtering(self) -> bool:
        return 0 < len(self.filtered_regions) < len(self.regions)

    @property
    def filtered_regions(self) -> list[CombinationRegion]:
        return [
            r
            for r in self.regions
            if r.filter
        ] if self.regions else []

    @property
    def filtered_region_indices(self) -> list[int]:
        return [
            i
            for i, r in enumerate(self.regions)
            if r.filter
        ]

    @property
    def filter_field_names(self) -> list[str]:
        return [
            r.field_name
            for r in self.filtered_regions
        ]

    @property
    def grouped_regions(self) -> dict[ReadGroupId, list[str]]:
        return {
            read_group: [r.id for r in regions]
            for read_group, regions in groupby(
                sorted(
                    self.regions,
                    key=lambda r: r.read_group.value),
                key=lambda r: r.read_group)
        }

    @property
    def has_duplicate_read_groups(self) -> bool:
        def on_duplicate(t: tuple[ReadGroupId, str]) -> None:
            logging.error(
                "Region '%s' from read group '%s' added multiple times to combination '%s'!" %
                (t[1], t[0].value, self.id))

        _, had_duplicates = get_set([
            (region.read_group, region.id)
            for region in self.regions
        ], on_duplicate)

        return had_duplicates


READ_GROUP_INFOS = {
    read_group_info.id: read_group_info
    for read_group_info in [
        ReadGroupInfo(id=ReadGroupId.DEFAULT),
        ReadGroupInfo(id=ReadGroupId.READ_1),
        ReadGroupInfo(id=ReadGroupId.READ_2, is_reverse=True)
    ]
}


class Options(BaseConfig):
    """
    Runtime options

    Can be overridden.

    test_swaps: Test alternative libraries on mismatch
    """

    # TODO: distinguish between library-dependent and library-independent-specific options?
    test_swaps: bool = False
    full_match_info: bool = False
    out_mm_match_info: bool = False
    out_mm_reads: bool = False
    rc_reverse: bool = False

    # Mismatching read counts
    count_mm_reads: bool = False
    sort_mm_read_counts: bool = False
    compress_mm_read_counts: bool = False

    @property
    def should_write_any_match_info(self) -> bool:
        # TODO: upon validation, mention interdependencies
        return self.out_mm_match_info or self.full_match_info

    def override(self, **kwargs) -> None:
        for k, v in kwargs.items():
            assert hasattr(self, k)
            if v is not None:
                self.__setattr__(k, v)

    @property
    def library_dependent_options(self) -> dict[str, bool]:
        return {
            'test_swaps': self.test_swaps,
            'full_match_info': self.full_match_info,
            'out_mm_reads': self.out_mm_reads,
            'count_mm_reads': self.count_mm_reads,
            'sort_mm_read_counts': self.sort_mm_read_counts
        }


def get_read_group_info(id: ReadGroupId) -> ReadGroupInfo:
    return READ_GROUP_INFOS[id]


def get_unique_library_ids_from_read_templates(read_templates: list[ReadTemplate]) -> set[str]:
    library_ids = [
        read_template.library_ids
        for read_template in read_templates
    ]
    return union(library_ids) if library_ids else set()


class Experiment(BaseConfig):
    # TODO: include species and assembly?
    sequencing_type: SequencingType
    libraries: list[LibraryInfo] = Field(default_factory=list)
    read_templates: list[ReadTemplate] = Field(default_factory=list)
    read_groups: dict[ReadGroupId, ReadGroupOptions] | None = None
    read_group_templates: dict[ReadGroupId, list[str]] = Field(default_factory=dict)
    combinations: list[CombinationInfo] | None = None
    default_options: Options = Options()

    @property
    def has_assigned_templates(self) -> bool:
        return len(self.read_group_templates) > 0 and any(
            len(tpl_ids) > 0
            for tpl_ids in self.read_group_templates.values()
        )

    @property
    def mode(self) -> ExecMode:
        return (
            ExecMode.LIBRARY_INDEPENDENT if len(self.libraries) == 0 else
            ExecMode.LIBRARY_DEPENDENT
        )

    @property
    def template_ids(self) -> list[str]:
        return [tpl.id for tpl in self.read_templates]

    @model_validator(mode='after')
    def _validate(self):
        success = True
        mode_label: str

        # Validate read group option overrides
        _ = self.read_group_infos

        match self.mode:
            case ExecMode.LIBRARY_INDEPENDENT:
                # TODO: provide better feedback
                assert len(self.libraries) == 0
                mode_label = "library-independent"

            case ExecMode.LIBRARY_DEPENDENT:
                assert len(self.libraries) > 0

                # Validate read template ID's
                if len(self.read_templates) == 0:
                    logging.error("No read templates defined in library-dependent mode!")
                    success = False

                tpl_ids, had_duplicates = get_set(self.template_ids, lambda tpl_id: logging.error(
                    "Duplicate read template ID '%s'!" % tpl_id))
                if had_duplicates:
                    success = False

                if not greedy_all(lambda x: x.validate(), self.read_templates):
                    success = False

                if self.read_group_templates:
                    for read_group_id, read_group_tpl_ids in self.read_group_templates.items():
                        def on_duplicate(tpl_id: str) -> None:
                            logging.error(
                                "Read template '%s' assigned multiple times to read group '%s'!" %
                                (tpl_id, read_group_id.value))

                        rg_tpl_ids, had_duplicates = get_set(read_group_tpl_ids, on_duplicate)
                        if had_duplicates:
                            success = False
                        for tpl_id in rg_tpl_ids:
                            if tpl_id not in tpl_ids:
                                logging.error(f"Read template '{tpl_id}' not defined!")
                                success = False

                else:
                    logging.error("No read groups defined!")
                    success = False

                # Validate combinations
                if self.combinations:
                    tpl_ids, had_duplicates = get_set(
                        [x.id for x in self.combinations],
                        lambda comb_id: logging.error(
                            "Duplicate combination ID '%s'!" % comb_id))
                    if had_duplicates:
                        success = False

                    if greedy_any(lambda x: x.has_duplicate_read_groups, self.combinations):
                        success = False

                    tpl_lir_ids = {
                        tpl.id: tpl.library_indep_region_ids
                        for tpl in self.read_templates
                    }

                    for cr in self.combinations:
                        read_group_ids: set[ReadGroupId] = {
                            region.read_group
                            for region in cr.regions
                        }
                        for id in read_group_ids:
                            if id not in self.read_group_templates:
                                logging.error("Read group '%s' in combination '%s' not assigned to any read template!" % (id.value, cr.id))
                                success = False

                        filtered_region_ids = set([
                            region.id
                            for region in cr.filtered_regions
                        ])
                        if cr.has_filter:
                            if not filtered_region_ids:
                                logging.error(
                                    "Combination %s has inclusion filters assigned to it " +
                                    "but none of its regions are marked as filtering!" % cr.id)
                                success = False
                        else:
                            if filtered_region_ids:
                                logging.error(
                                    "Combination %s has filtering regions (%s) " %
                                    (cr.id, ', '.join(sorted(filtered_region_ids))) +
                                    "but no inclusion filters assigned to it!")
                                success = False
                        for region_id in filtered_region_ids:
                            for tpl_id, lir_ids in tpl_lir_ids.items():
                                if region_id in lir_ids:
                                    logging.error(
                                        "Region %s is library-independent in template %s " % (region_id, tpl_id) +
                                        "and therefore may not be filtering in a combination!")
                                    success = False

                mode_label = "library-dependent"

            case _:
                raise NotImplementedError("Unsupported execution mode!")

        if not success:
            raise InvalidExperiment

        logging.info("Running in %s mode" % mode_label)

        return self

    @property
    def library_indep_region_ids(self) -> set[str]:
        return union([
            tpl.library_indep_region_ids
            for tpl in self.read_templates
        ])

    @property
    def combination_region_ids_by_read_group(self) -> dict[ReadGroupId, set[str]]:
        if not self.combinations:
            return {}

        d = defaultdict(set)
        for g in self.combinations:
            for read_group, region_ids in g.grouped_regions.items():
                d[read_group] |= set(region_ids)

        return d

    @property
    def combination_filters(self) -> set[str]:
        return union([
            set(g.filters)
            for g in self.combinations
            if g.filters
        ]) if self.combinations else set()

    @property
    def has_combination_filters(self) -> bool:
        return any(
            g.has_filter
            for g in self.combinations
        ) if self.combinations else False

    @property
    def read_group_count(self) -> int:
        return len(self.read_group_templates)

    @property
    def has_default_read_group(self) -> bool:
        # TODO: validate only the default read group is mutually exclusive with the mate groups, for now
        return (
            self.read_group_count == 1 and
            ReadGroupId.DEFAULT in self.read_group_templates
        )

    @property
    def read_group_infos(self) -> dict[ReadGroupId, ReadGroupInfo]:
        d = {
            read_group_id: get_read_group_info(read_group_id)
            for read_group_id in self.read_group_templates
        }

        # Override read group options
        if self.read_groups:
            try:
                validate_read_group_ids(set(self.read_groups.keys()))
            except ValueError as ex:
                logging.error(ex.args[0])
                raise InvalidExperiment

            for read_group_id, read_group_options in self.read_groups.items():
                if read_group_id not in d:
                    d[read_group_id] = get_read_group_info(read_group_id)
                d[read_group_id].update(read_group_options)

        return d

    @property
    def library_indices(self) -> dict[str, int]:
        return {
            library.id: i
            for i, library
            in enumerate(self.libraries)
        }

    @property
    def assigned_library_indices(self) -> dict[str, int]:
        read_template_ids = set(chain.from_iterable(self.read_group_templates.values()))
        library_ids = get_unique_library_ids_from_read_templates([
            t
            for t in self.read_templates
            if t.id in read_template_ids
        ])
        assigned_library_indices = self.get_library_indices(library_ids)

        # TODO: this test should be performed at validation
        if len(assigned_library_indices) != len(self.library_indices):
            logging.warning("Some libraries were not assigned to any read template!")

        return assigned_library_indices

    def get_library_indices(self, library_ids: set[str]) -> dict[str, int]:
        return {
            library_id: libray_index
            for library_id, libray_index in self.library_indices.items()
            if library_id in library_ids
        }

    def get_read_group_templates(self, read_group_id: ReadGroupId) -> list[ReadTemplate]:
        # Retrieve the read templates for the read group
        read_template_ids = set(self.read_group_templates[read_group_id])

        read_templates = [
            t
            for t in self.read_templates
            if t.id in read_template_ids
        ]

        if len(read_templates) != len(read_template_ids):
            for read_template_id in read_template_ids - {t.id for t in read_templates}:
                logging.error("Read template ID not found: '%s'!" % read_template_id)
                raise ValueError

        return read_templates

    def _get_read_group_template_region_indices(self, read_group_id: ReadGroupId) -> list[dict[str, int]]:
        read_templates = self.get_read_group_templates(read_group_id)
        region_offsets = get_offsets([
            len(read_template.regions)
            for read_template in read_templates[:-1]
        ]) if len(read_templates) > 1 else [0]

        # map: read template ID -> region ID -> region index (relative to the read group)
        return [
            tpl.get_region_indices(offset=region_offsets[tpl_index])
            for tpl_index, tpl in enumerate(read_templates)
        ]

    def get_combination_region_indices(self) -> CombinationRegionIndices:
        if not self.combinations:
            return {}

        # map: read group ID -> read template ID -> region ID -> region index (relative to the read group)
        read_group_template_region_indices = {
            read_group_id: self._get_read_group_template_region_indices(read_group_id)
            for read_group_id in self.read_group_templates.keys()
        }

        # Check the combination regions exist in all templates of the relevant read groups
        for read_group_id, read_group_region_ids in self.combination_region_ids_by_read_group.items():
            for tpl_region_indices in read_group_template_region_indices[read_group_id]:
                missing_region_ids = read_group_region_ids - tpl_region_indices.keys()
                if missing_region_ids:
                    for region_id in missing_region_ids:
                        logging.error(
                            "Region '%s' from combination not found in one of the templates assigned to read group '%s'!" % (
                                region_id,
                                read_group_id.value
                            ))
                    raise ValueError("Invalid combination setup!")

        return read_group_template_region_indices
