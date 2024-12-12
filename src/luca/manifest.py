# LUCA
#
# Copyright (C) 2024 Genome Research Ltd.
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

import os
from abc import ABC

from pydantic import BaseModel, model_validator


class BaseOutputFilePaths(BaseModel, ABC):
    counts: str


class CombinationOutputFilePaths(BaseOutputFilePaths):
    n: int


class LibraryOutputFilePaths(BaseOutputFilePaths):
    stats: str


class LibraryIndependentOutputFilePaths(BaseOutputFilePaths):
    pass


class LibraryIndependentMismatchingOutputFilePaths(BaseOutputFilePaths):
    is_sorted: bool
    is_compressed: bool


class OutputFileManifest(BaseModel):
    sample_name: str | None

    # Library-dependent counts
    library_count: int
    total_library_templates: int
    library_file_paths: dict[int, LibraryOutputFilePaths]

    # Library-independent counts
    library_independent_count_file_paths: LibraryIndependentOutputFilePaths | None
    library_independent_mm_count_file_paths: LibraryIndependentMismatchingOutputFilePaths | None
    total_dynamic_targets: int
    dynamic_target_file_paths: LibraryOutputFilePaths | None

    # Combinations
    combination_file_paths: list[CombinationOutputFilePaths]

    @model_validator(mode='after')
    def _validate(self):
        for i in range(self.library_count):
            assert i in self.library_file_paths
        return self

    @property
    def library_independent_counts_file_path(self) -> str | None:
        x = self.library_independent_count_file_paths
        return x.counts if x else None

    @property
    def library_independent_mm_counts_file_path(self) -> str | None:
        x = self.library_independent_mm_count_file_paths
        return x.counts if x else None

    @property
    def dynamic_target_counts_file_path(self) -> str | None:
        x = self.dynamic_target_file_paths
        return x.counts if x else None

    def set_library_independent_file_paths(self, counts: str) -> None:
        self.library_independent_count_file_paths = LibraryIndependentOutputFilePaths(
            counts=counts)

    def set_library_independent_mm_file_paths(self, counts: str, is_sorted: bool, is_compressed: bool) -> None:
        self.library_independent_mm_count_file_paths = LibraryIndependentMismatchingOutputFilePaths(
            counts=counts, is_sorted=is_sorted, is_compressed=is_compressed)

    def set_dynamic_target_file_paths(self, stats: str, counts: str) -> None:
        self.dynamic_target_file_paths = LibraryOutputFilePaths(
            stats=stats, counts=counts)

    def push_library_file_paths(self, library_index: int, stats: str, counts: str) -> None:
        assert library_index not in self.library_file_paths
        self.library_count += 1
        self.library_file_paths[library_index] = LibraryOutputFilePaths(
            stats=stats, counts=counts)

    def push_combination_file_paths(self, n: int, counts: str) -> None:
        self.combination_file_paths.append(
            CombinationOutputFilePaths(n=n, counts=counts))

    def write(self, out_dir: str) -> None:
        with open(os.path.join(out_dir, "manifest.json"), 'w') as fh:
            fh.write(self.model_dump_json())

    @classmethod
    def empty(cls, sample_name: str | None = None):
        return cls(
            sample_name=sample_name,
            library_count=0,
            total_library_templates=0,
            library_file_paths={},
            total_dynamic_targets=0,
            library_independent_count_file_paths=None,
            library_independent_mm_count_file_paths=None,
            dynamic_target_file_paths=None,
            combination_file_paths=[])
