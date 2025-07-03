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

# Includes logic from pyQUEST (library module)

import abc
import gc
import logging
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Any, Sized

import numpy as np

from .codec import Codec
from .errors import InvalidLibraryError, TargetNotFound
from .experiment import LibraryReverseCondition
from .fs import open_output_file
from .library_stats import LibraryStats
from .readers.tsv import TsvEmpty, parse_tsv, TsvFieldNotFound, TsvTruncatedRow
from .utils import flip_dict, get_offsets, get_range_from, get_stats, get_unique_length, has_duplicates, is_dna, reverse_complement


# TSV field names
LIBRARY_SEQUENCE_FIELD = 'sequence'
LIBRARY_FIELDS = [LIBRARY_SEQUENCE_FIELD]


@dataclass
class TargetEncoder:
    _encoder: dict[str, int] = field(default_factory=dict)

    def encode(self, sequence: str) -> int:
        try:
            return self._encoder[sequence]
        except KeyError:
            raise TargetNotFound(sequence)

    def to_decoder(self) -> dict[int, str]:
        return flip_dict(self._encoder)


def get_target_encoder(targets: list[str], target_index_offset: int = 0) -> dict[str, int]:
    return {
        target: i + target_index_offset
        for i, target in enumerate(targets)
    }


def get_target_decoder(targets: list[str], target_index_offset: int = 0) -> dict[int, str]:
    return {
        i + target_index_offset: target
        for i, target in enumerate(targets)
    }


def write_count(fh, i: int, target: str, count: int) -> None:
    fh.write(str(i))
    fh.write('\t')
    fh.write(target)
    fh.write('\t')
    fh.write(str(count))
    fh.write('\n')


class Frozen:
    def _set(self, attr: str, value: Any) -> None:
        object.__setattr__(self, attr, value)


# TODO: store library identifier (e.g., to name output files?)
@dataclass(slots=True, frozen=True)
class TargetLibrary(abc.ABC, Frozen, Sized):
    target_index_offset: int

    @property
    @abc.abstractmethod
    def sequences(self) -> list[str]:
        pass

    @abc.abstractmethod
    def match(self, seq: str, offset: int = 0, length: int = 0, anchor_left: bool = True) -> tuple[int, int]:
        """Try to match from offset and return the target index and length"""
        pass


@dataclass(slots=True, frozen=True)
class DynamicMultiTargetLibrary(TargetLibrary):
    _codec: Codec = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, '_codec', Codec(start=self.target_index_offset))

    def __len__(self) -> int:
        return len(self._codec)

    def update_decoder(self, decoder: dict[int, str]) -> None:
        for i, k in enumerate(self._codec.items, start=self._codec.start):
            decoder[i] = k

    @property
    def sequences(self) -> list[str]:
        return self._codec.items

    def match_full(self, target: str) -> tuple[int, int]:
        return self._codec.add(target), len(target)

    def match(self, seq: str, offset: int = 0, length: int = 0, anchor_left: bool = True) -> tuple[int, int]:
        seq_len = len(seq)
        if anchor_left:
            end = offset + length
            if end > seq_len:
                return 0, 0
            start = offset
        else:
            end = seq_len - offset
            start = end - length
            if start < 0:
                return 0, 0

        target = seq[start:end]
        return self._codec.add(target), length

    def decode(self, target_index: int) -> str:
        return self._codec.decode(target_index)


@dataclass(slots=True, frozen=True)
class MonoTargetLibrary(TargetLibrary):
    target: str
    target_index: int = field(default=0)
    _target_length: int = field(init=False)

    def __len__(self) -> int:
        return 1

    def __post_init__(self) -> None:
        assert self.target
        self._set('_target_length', len(self.target))

    @property
    def sequences(self) -> list[str]:
        return [self.target]

    def match(self, seq: str, offset: int = 0, length: int = 0, anchor_left: bool = True) -> tuple[int, int]:
        is_match: bool
        if anchor_left:
            is_match = seq[offset:].startswith(self.target)
        else:
            end = len(seq) - offset
            start = end - self._target_length
            is_match = (seq[start:end] == self.target) if start >= 0 else False

        return (
            (self.target_index, self._target_length) if is_match else
            (0, 0)
        )


@dataclass(slots=True, frozen=True)
class BaseMultiTargetLibrary(TargetLibrary, abc.ABC):
    _sequences: list[str]

    @classmethod
    @abc.abstractmethod
    def from_sequences(cls, sequences: list[str], target_index_offset: int = 0):
        pass

    @property
    def sequences(self) -> list[str]:
        return self._sequences

    def __len__(self) -> int:
        return len(self._sequences)


@dataclass(slots=True, frozen=True)
class MultiTargetUniformLibrary(BaseMultiTargetLibrary):
    """Multiple targets of uniform lengths"""

    targets: dict[str, int]
    _target_length: int = field(init=False)

    def __post_init__(self) -> None:
        target_length: int = 0
        if self.targets:
            targets = iter(self.targets.keys())
            target_length = len(next(targets))
            for k in targets:
                if len(k) != target_length:
                    raise InvalidLibraryError("Non uniform target lengths!")

        self._set('_target_length', target_length)

    @classmethod
    def from_sequences(cls, sequences: list[str], target_index_offset: int = 0):
        return cls(target_index_offset, sequences, get_target_encoder(
            sequences, target_index_offset=target_index_offset))

    def match(self, seq: str, offset: int = 0, length: int = 0, anchor_left: bool = True) -> tuple[int, int]:
        seq_len = len(seq)

        if anchor_left:
            end = offset + self._target_length
            if end > seq_len:
                return 0, 0
            start = offset
        else:
            end = seq_len - offset
            start = end - self._target_length
            if start < 0:
                return 0, 0

        target = self.targets.get(seq[start:end])
        return (
            (target, self._target_length) if target is not None else
            (0, 0)
        )


@dataclass(slots=True, frozen=True)
class MultiTargetLibrary(BaseMultiTargetLibrary):
    length_targets: dict[int, dict[str, int]]
    _lengths: list[int] = field(init=False)

    def __post_init__(self) -> None:
        if 0 in self.length_targets:
            raise ValueError("zero-length target sequence")
        self._set_lengths()

    def _set_lengths(self) -> None:
        lengths = sorted(self.length_targets.keys(), reverse=True)
        for length in lengths:
            assert length > 0
        self._set('_lengths', lengths)

    @classmethod
    def from_sequences(cls, sequences: list[str], target_index_offset: int = 0):
        d = defaultdict(dict)
        for i, s in enumerate(sequences):
            d[len(s)][s] = i + target_index_offset
        return cls(target_index_offset, sequences, d)

    def match(self, seq: str, offset: int = 0, length: int = 0, anchor_left: bool = True) -> tuple[int, int]:
        seq_len = len(seq)

        # Assumption: the lengths are in descending order
        for target_length in self._lengths:
            if anchor_left:
                end = offset + target_length
                if end > seq_len:
                    continue
                start = offset
            else:
                end = seq_len - offset
                start = end - target_length
                if start < 0:
                    continue

            target = self.length_targets[target_length].get(seq[start:end])
            if target is not None:
                return target, target_length

        return 0, 0


@dataclass(slots=True, frozen=True)
class RawLibrary(Sized):
    # id: str
    index: int
    reverse_on: LibraryReverseCondition
    sequences: list[str]
    is_uniform_length: bool = field(init=False)
    fp: str | None = None

    def __len__(self) -> int:
        return len(self.sequences)

    @property
    def src(self) -> str:
        return self.fp or str(self.index)

    def should_apply_reverse_complement(self, is_reverse: bool) -> bool:
        return (
            self.reverse_on == LibraryReverseCondition.ALWAYS or
            self.reverse_on == LibraryReverseCondition.REVERSE_GROUP and is_reverse
        )

    def __post_init__(self) -> None:
        n = len(self.sequences)
        assert n > 0

        # Check for invalid targets
        for target in self.sequences:
            if not is_dna(target):
                raise InvalidLibraryError(
                    f"Invalid target '{target}' in library {self.src}!")

        if n == 1:
            is_uniform_length = True
        else:

            # Check for redundant targets
            unique_count = get_unique_length(self.sequences)
            if unique_count < n:
                logging.warning(
                    "%d redundant targets found in library %s!" %
                    (n - unique_count, self.src))
                raise InvalidLibraryError(f"Duplicate targets in library {self.src}!")

            # Assess target length uniformity
            target_length = len(self.sequences[0])
            is_uniform_length = all(
                len(seq) == target_length
                for seq in self.sequences[1:]
            )

        object.__setattr__(self, 'is_uniform_length', is_uniform_length)

    @property
    def forward_unique_targets(self) -> set[str]:
        # TODO: consider...
        return set(
            map(reverse_complement, self.sequences) if self.reverse_on == LibraryReverseCondition.ALWAYS else
            self.sequences
        )

    @classmethod
    def load(cls, fp: str, index: int, reverse_on: LibraryReverseCondition):
        try:
            sequences = list(map(lambda r: r[0], parse_tsv(fp, LIBRARY_FIELDS)))

            if not sequences:
                raise TsvEmpty(fp)

        except TsvFieldNotFound as ex:
            raise InvalidLibraryError(
                f"Field '{ex.field_name}' not found in library header at {ex.fp}!")

        except TsvTruncatedRow as ex:
            raise InvalidLibraryError(f"Truncated row in library at {ex.fp}!")

        except TsvEmpty as ex:
            raise InvalidLibraryError(f"Empty library at {ex.fp}!")

        return cls(index, reverse_on, sequences, fp=fp)

    def build(self, target_index_offset: int = 0, is_reverse: bool = False) -> TargetLibrary:
        # BEWARE: anything other than nonambiguous uppercase DNA base symbols
        #  would be silently ignored by the reverse complement function!

        sequences = self.sequences if not self.should_apply_reverse_complement(is_reverse) else list(
            map(reverse_complement, self.sequences))

        if len(self) == 1:
            return MonoTargetLibrary(
                target_index_offset, sequences[0], target_index=target_index_offset)

        return (
            MultiTargetUniformLibrary if self.is_uniform_length else
            MultiTargetLibrary
        ).from_sequences(
            sequences, target_index_offset=target_index_offset)

    def to_encoder(self, target_index_offset: int = 0) -> dict[str, int]:
        return get_target_encoder(self.sequences, target_index_offset=target_index_offset)

    def to_decoder(self, target_index_offset: int = 0) -> dict[int, str]:
        return get_target_decoder(self.sequences, target_index_offset=target_index_offset)

    def write_counts(self, fp: str, counts: np.ndarray, target_index_offset: int = 0) -> None:
        assert counts.shape[0] == len(self)
        with open_output_file(fp) as f:
            fh = f.fh
            for i, target, count in zip(
                get_range_from(target_index_offset, len(self)),
                self.sequences,
                counts
            ):
                write_count(fh, i, target, count)


@dataclass(slots=True, frozen=True)
class LibraryBuilder(Frozen, Sized):
    libraries: dict[int, RawLibrary]
    library_offsets: dict[int, int]
    library_shared_targets: dict[int, dict[int, dict[int, int]]] = field(init=False)

    @classmethod
    def empty(cls):
        return cls({}, {})

    def __post_init__(self) -> None:
        library_shared_target_map: dict[int, dict[int, dict[int, int]]] = defaultdict(dict)

        if len(self) > 1:
            library_target_sets = {
                library_index: library.forward_unique_targets
                for library_index, library in self.libraries.items()
            }

            # Collect the shared targets for each pair of libraries (forward sequence)
            library_shared_targets: dict[int, dict[int, set[str]]] = defaultdict(dict)
            for i, a in library_target_sets.items():
                for j, b in library_target_sets.items():
                    if j != i and j not in library_shared_targets:
                        shared_targets = a & b
                        if shared_targets:
                            n = len(shared_targets)
                            library_shared_targets[i][j] = shared_targets
                            library_shared_targets[j][i] = shared_targets
                            logging.warning(
                                "Libraries %d and %d share %d target%s!" %
                                (i, j, n, 's' if n > 1 else ''))

            overlapping_library_indices = list(library_shared_targets.keys())
            del library_target_sets
            gc.collect()

            # Build the required target encoders
            library_encoders = {
                library_index: self.libraries[library_index].to_encoder(
                    self.library_offsets[library_index])
                for library_index in overlapping_library_indices
            }

            # Map the encoded targets from one library to the other
            for i in overlapping_library_indices:
                for j in overlapping_library_indices:
                    if j != i:
                        library_shared_target_map[i][j] = {
                            library_encoders[i][target]: library_encoders[j][target]
                            for target in library_shared_targets[i][j]
                        }

        object.__setattr__(self, 'library_shared_targets', library_shared_target_map)

    @classmethod
    def from_raw_libraries(cls, raw_libraries: list[RawLibrary]):
        n = len(raw_libraries)
        match n:

            case 0:
                return cls({}, {})

            case 1:
                library = raw_libraries[0]
                return cls({library.index: library}, {library.index: 0})

            case _:
                if has_duplicates([x.index for x in raw_libraries]):
                    raise ValueError("Duplicate library indices!")

                libraries = {
                    raw_library.index: raw_library
                    for raw_library in raw_libraries
                }

                # Compute target index offsets
                library_indices = sorted(libraries.keys())
                library_offsets = get_offsets([
                    len(libraries[library_index])
                    for library_index in library_indices[:-1]
                ])

                return cls(libraries, dict(zip(library_indices, library_offsets)))

    def __len__(self) -> int:
        return len(self.libraries)

    def build(self, library_index: int, is_reverse: bool = False) -> TargetLibrary:
        return self.libraries[library_index].build(
            target_index_offset=self.library_offsets[library_index],
            is_reverse=is_reverse)

    def build_decoder(self, library_indices: set[int] | None = None) -> dict[int, str]:
        if library_indices:
            assert len(library_indices) <= len(self.libraries)
        return {
            k: v
            for library_index in (library_indices or self.libraries)
            for k, v in self.libraries[library_index].to_decoder(
                target_index_offset=self.library_offsets[library_index]).items()
        }

    def build_encoder(self, library_indices: set[int] | None = None) -> TargetEncoder:
        return TargetEncoder({
            k: v
            for library_index in (library_indices or self.libraries)
            for k, v in self.libraries[library_index].to_encoder(
                target_index_offset=self.library_offsets[library_index]).items()
        })

    def get_library_size(self, library_index: int) -> int:
        return len(self.libraries[library_index])

    @property
    def library_indices(self) -> list[int]:
        return sorted(list(self.libraries.keys()))

    @property
    def library_sizes(self) -> dict[int, int]:
        return {
            library_index: self.get_library_size(library_index)
            for library_index in self.libraries
        }

    def get_library_target_index_range(self, library_index: int) -> tuple[int, int]:
        start = self.library_offsets[library_index]
        end = start + self.get_library_size(library_index)
        return start, end


def fill_library_stats(stats: LibraryStats, template_counts: np.ndarray) -> None:
    # Sort the counts (required by `get_stats`)
    template_counts.sort()

    low_counts: Counter[int] = Counter()

    # Count templates with low read counts
    # TODO: take better advantage of the sorting...?
    low_counts[0] = int(np.count_nonzero(template_counts == 0))
    count_thresholds: list[int] = [15, 30]
    # if opt.custom_count_threshold is not None:
    #     count_thresholds.append(opt.custom_count_threshold)
    for t in count_thresholds:
        low_counts[t] = int(np.count_nonzero(template_counts < t))

    # Generate all stats
    mapped_to_template_reads, mean_count_per_template, median_count_per_template, gini_coefficient = get_stats(
        template_counts, gini_corr=False)

    stats.mapped_to_template_reads = mapped_to_template_reads
    stats.total_templates = template_counts.shape[0]

    # Derived stats
    stats.mean_count_per_template = mean_count_per_template
    stats.median_count_per_template = median_count_per_template
    stats.gini_coefficient = gini_coefficient

    # Count thresholds
    stats.zero_count_templates = low_counts.get(0, 0)
    stats.low_count_templates_lt_15 = low_counts.get(15, 0)
    stats.low_count_templates_lt_30 = low_counts.get(30, 0)
