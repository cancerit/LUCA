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

import logging
import os
from abc import ABC
from collections import Counter
from dataclasses import dataclass, field
from typing import Generator

import numpy as np

from .codec import Codec
from .errors import CustomFileException, EmptyFileError
from .fs import get_combination_counts_file_name, get_library_counts_file_name, get_library_stats_file_name, open_output_file, get_read_counts_file_name
from .library import fill_library_stats, write_count
from .library_stats import LibraryStats
from .manifest import OutputFileManifest
from .readers.tsv import parse_tsv_at
from .utils import get_range_from, sort_counter_desc, write_json


TSV_LIBRARY_INDEX = 0
TSV_LIBRARY_SEQUENCE = 1
TSV_LIBRARY_COUNT = 2

LIBRARY_INDEX = 0
LIBRARY_COUNT = 1
LIBRARY_SEQUENCE = 2


@dataclass(slots=True)
class CountFieldIndices:
    index: int | None
    sequence: int
    count: int

    def __post_init__(self) -> None:
        assert self.index != self.sequence != self.count


lib_count_field_indices = CountFieldIndices(TSV_LIBRARY_INDEX, TSV_LIBRARY_SEQUENCE, TSV_LIBRARY_COUNT)
lib_indep_count_field_indices = CountFieldIndices(None, 0, 1)


class InconsistentLibraries(CustomFileException, ABC):
    pass


class UnexpectedTarget(InconsistentLibraries):
    @property
    def message(self) -> str:
        return f"Unexpected target sequence at '{self.fp}'!"


class TargetIndexOutOfBounds(InconsistentLibraries):
    @property
    def message(self) -> str:
        return f"Unexpected target index at '{self.fp}': out of bounds!"


class EmptyCountsFileError(EmptyFileError):
    pass


def get_library_counts_array(n: int) -> np.ndarray:
    return np.zeros(n, dtype=np.uint64)


def add_library_counts_from_file(fp: str, a: np.ndarray, targets: list[str] | None = None) -> tuple[int, int]:
    field_indices = [TSV_LIBRARY_INDEX, TSV_LIBRARY_COUNT]
    if targets:
        field_indices.append(TSV_LIBRARY_SEQUENCE)
    it = parse_tsv_at(fp, field_indices)

    try:

        # Get first target index
        t = next(it)

    except StopIteration:
        raise EmptyCountsFileError(fp)

    buffer_size = a.shape[0]

    def proc_entry(t_: list[str]) -> int:

        # Validate target index
        i_ = int(t_[LIBRARY_INDEX])
        if i_ < 0 or i_ >= buffer_size:
            raise TargetIndexOutOfBounds(fp)

        # Validate target sequence
        if targets and t_[LIBRARY_SEQUENCE] != targets[i_]:
            raise UnexpectedTarget(fp)

        # Add count
        a[i_] += np.uint64(t_[LIBRARY_COUNT])
        return i_

    i = proc_entry(t)
    start = i
    n = 1

    for t in it:
        i = proc_entry(t)
        n += 1

    # Get last target index
    end = i

    # Test against gaps in the sequence
    assert n == end - start + 1

    return start, end


@dataclass(slots=True)
class MultiLibraryCounts:
    manifest: OutputFileManifest
    targets: list[str]
    counts: np.ndarray
    library_target_ranges: dict[int, tuple[int, int]]

    @classmethod
    def from_manifest(cls, manifest: OutputFileManifest, root_dir: str):
        counts = get_library_counts_array(manifest.total_library_templates)
        targets = []

        a = sorted([
            (library_index, library_file_paths)
            for library_index, library_file_paths
            in manifest.library_file_paths.items()
        ], key=lambda t: t[0])

        for library_index, library_file_paths in a:
            fpc = os.path.join(root_dir, library_file_paths.counts)
            targets += [r[0] for r in parse_tsv_at(fpc, [1])]

        mlc = cls(
            manifest=manifest,
            targets=targets,
            counts=counts,
            library_target_ranges={})
        mlc.library_target_ranges = mlc._update_from_manifest(
            root_dir, manifest, validate_sequences=False)
        return mlc

    def _update_library(self, root_dir: str, library_index: int, fns: str, fnc: str, validate_sequences: bool = False) -> tuple[int, int]:
        # Load counts
        fpc = os.path.join(root_dir, fnc)
        assert os.path.isfile(fpc)
        expected_targets = self.targets if validate_sequences else None
        try:
            start, end = add_library_counts_from_file(
                fpc, self.counts, targets=expected_targets)
        except IndexError:
            # TODO: log error
            raise
        n = end - start + 1

        # Validate number of targets (if stats are available)
        fps = os.path.join(root_dir, fns)
        if os.path.isfile(fps):
            # TODO: handle parsing erros
            stats = LibraryStats.load(fps)
            assert n == stats.total_templates
        else:
            logging.warning(
                "Number of targets in library counts file could not be validated: " +
                "stats file not found!")

        return start, end

    def _update_from_manifest(self, root_dir: str, manifest: OutputFileManifest, validate_sequences: bool = False) -> dict[int, tuple[int, int]]:
        return {
            library_index: self._update_library(
                root_dir,
                library_index,
                library_file_paths.stats,
                library_file_paths.counts,
                validate_sequences=validate_sequences)
            for library_index, library_file_paths in manifest.library_file_paths.items()
        }

    def update_from_manifest(self, manifest: OutputFileManifest, root_dir: str, validate_sequences: bool = False) -> None:
        d = self._update_from_manifest(root_dir, manifest, validate_sequences=validate_sequences)
        assert d == self.library_target_ranges

    def get_library_stats(self, library_index: int, start: int, end: int) -> LibraryStats:
        library_stats = LibraryStats.empty()
        # BEWARE: the following function SORTS the counts IN PLACE
        fill_library_stats(library_stats, self.counts[start:end + 1])
        return library_stats

    def get_all_library_stats(self) -> dict[int, LibraryStats]:
        return {
            library_index: self.get_library_stats(library_index, start, end)
            for library_index, (start, end) in self.library_target_ranges.items()
        }

    def _write_library_stats(self, out_dir: str) -> None:
        for library_index, library_stats in self.get_all_library_stats().items():
            fns = get_library_stats_file_name(library_index)
            fps = os.path.join(out_dir, fns)

            logging.info(f"Writing statistics file: {fps}")
            write_json(library_stats.to_dict(), fps)

    def _write_library_counts(self, out_dir: str) -> None:
        library_sizes = [
            (library_index, start, end - start + 1)
            for library_index, (start, end) in self.library_target_ranges.items()
        ]
        library_sizes.sort(key=lambda t: t[0])
        for library_index, start, library_size in library_sizes:
            fp = os.path.join(out_dir, get_library_counts_file_name(library_index))
            with open_output_file(fp) as f:
                fh = f.fh
                for target_index in get_range_from(start, library_size):
                    write_count(
                        fh,
                        target_index,
                        self.targets[target_index],
                        self.counts[target_index])

    def write(self, out_dir: str) -> None:
        # NOTE: do not change the order (stats calculation sorts the counts!)
        self._write_library_counts(out_dir)
        self._write_library_stats(out_dir)


def iter_counts_file(fp: str, field_indices: list[int]) -> tuple[Generator[list[str], None, None], list[str]]:
    """Return first row and iterator of a TSV file"""

    it = parse_tsv_at(fp, field_indices)
    try:
        return it, next(it)
    except StopIteration:
        raise EmptyCountsFileError(fp)


def cast_counts_row(t: list[str]) -> tuple[int, int, str]:
    return int(t[LIBRARY_INDEX]), int(t[LIBRARY_COUNT]), t[LIBRARY_SEQUENCE]


def cast_counts_no_index_row(t: list[str]) -> tuple[int, str]:
    return int(t[1]), t[0]


def codec_from_counts(fp: str, fields: CountFieldIndices) -> tuple[Codec, Counter[int]]:
    codec: Codec
    counts: Counter[int] = Counter()

    if fields.index is not None:
        it, t = iter_counts_file(fp, [fields.index, fields.count, fields.sequence])
        i, count, target = cast_counts_row(t)

        start = i

        codec = Codec(start=start)

        # Add first target
        counts[start] = count
        assert codec.add(target) == i

        n = 1
        for t in it:
            i, count, target = cast_counts_row(t)
            assert target not in codec
            counts[i] = count
            assert codec.add(target) == i
            n += 1

        # Get last target index
        end = i

        # Test against gaps in the sequence
        assert n == end - start + 1

    else:
        it, t = iter_counts_file(fp, [fields.sequence, fields.count])
        count, target = cast_counts_no_index_row(t)
        codec = Codec(start=0)
        counts[codec.add(target)] = count
        for t in it:
            count, target = cast_counts_no_index_row(t)
            assert target not in codec
            counts[codec.add(target)] = count

    return codec, counts


def update_codec_from_counts(codec: Codec, fp: str, validate_index: bool = True) -> None:
    for i, target in parse_tsv_at(fp, [TSV_LIBRARY_INDEX, TSV_LIBRARY_SEQUENCE]):
        target_index = codec.add(target)
        if validate_index:
            assert target_index == int(i)


@dataclass(slots=True)
class MultiDynamicLibraryCounts:
    fields: CountFieldIndices
    codec: Codec
    counts: Counter[int]

    @classmethod
    def from_counts_file(cls, fields: CountFieldIndices, fp: str):
        codec, counts = codec_from_counts(fp, fields)
        return cls(fields, codec, counts)

    def _update(self, target: str, count: int) -> None:
        self.counts[self.codec.add(target)] += count

    def _update_raw(self, t: list[str]) -> None:
        self._update(t[self.fields.sequence], int(t[self.fields.count]))

    def _get_field_indices(self) -> list[int]:
        return (
            [] if self.fields.index is None else
            [self.fields.index]
        ) + [
            self.fields.sequence,
            self.fields.count
        ]

    def update_from_counts_file(self, fp: str) -> None:
        it, t = iter_counts_file(fp, self._get_field_indices())

        self._update_raw(t)
        for t in it:
            self._update_raw(t)

    def write_as_regional(self, library_index: int, out_dir: str) -> None:
        fp = os.path.join(out_dir, get_library_counts_file_name(library_index))
        with open_output_file(fp) as f:
            fh = f.fh
            for target, target_index in self.codec.encoder.items():
                write_count(fh, target_index, target, self.counts[target_index])

    def write_as_global(self, out_dir: str, sort: bool = False, prefix: str | None = None, compress: bool = False) -> None:
        fp = os.path.join(out_dir, get_read_counts_file_name(prefix=prefix))
        counts = dict(sort_counter_desc(self.counts)) if sort else self.counts
        with open_output_file(fp, compress=compress) as f:
            fh = f.fh
            for target, target_index in self.codec.encoder.items():
                write_count(fh, target_index, target, counts[target_index])


@dataclass(slots=True)
class MultiCombinationCounts:
    n: int
    codec: Codec = field(init=False)
    counts: Counter[tuple] = field(init=False)

    def __post_init__(self) -> None:
        self.codec = Codec()
        self.counts = Counter()

    def update_from_counts_file(self, fp: str) -> None:
        n = self.n
        field_indices = list(range(n + 1))
        target_sl = slice(0, n)

        def update_combination(t_: list[str]) -> None:
            self.counts[tuple([
                self.codec.add(target)
                for target in t_[target_sl]
            ])] += int(t_[n])

        it, t = iter_counts_file(fp, field_indices)
        update_combination(t)
        for t in it:
            update_combination(t)

    def write(self, combination_index: int, out_dir: str) -> None:
        fp = os.path.join(out_dir, get_combination_counts_file_name(combination_index))
        with open_output_file(fp) as f:
            fh = f.fh
            for combination, count in self.counts.items():
                for item in combination:
                    fh.write(self.codec.decode(item))
                    fh.write('\t')
                fh.write(str(count))
                fh.write('\n')
