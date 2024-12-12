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

from collections import Counter
from typing import Any, Generator
import os

from .fs import get_read_counts_file_name, open_output_file
from .readers.types import SendT, YieldT
from .stats import LibraryIndependentStats


def write_read_counts(counts: Counter | list[tuple[str | int]], out_dir: str, compress: bool = False, prefix: str | None = None) -> str:
    fn = get_read_counts_file_name(prefix=prefix)
    with open_output_file(fn, out_dir=out_dir, compress=compress) as f:
        fh = f.fh
        for sequence, count in (
            counts.items() if isinstance(counts, Counter) else
            counts
        ):
            fh.write(sequence)
            fh.write('\t')
            fh.write(str(count))
            fh.write('\n')

        return f.fp


def count_reads(seqs: Generator[str, SendT, Any], out_dir: str) -> LibraryIndependentStats:
    counts: Counter[str] = Counter()
    stats: LibraryIndependentStats | None = None
    try:
        # TODO: no need to feed back anything (add simpler generator?)
        seq = next(seqs)
        while seq:
            counts[seq] += 1
            seq = next(seqs)

    except StopIteration as ex:
        stats = ex.value

    assert stats

    # Write counts
    write_read_counts(counts, out_dir)

    return stats
