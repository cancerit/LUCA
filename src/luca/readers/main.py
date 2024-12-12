# LUCA
#
# Copyright (C) 2023 Genome Research Ltd.
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

# Adapted from pyCROQUET (readparser module)
# Adapted from pyQUEST (readers.main module)

import logging
from time import time
from typing import Generator

from ..stats import LibraryIndependentStats
from .read_file_info import ReadFileFormat, ReadFileInfo
from .hts import iter_paired_end, iter_single_end, parse_htsfile
from .types import SendT, YieldT


def parse_reads(
    read_file_info: ReadFileInfo,
    sample: str | None,
    cpus: int,
    reference: str | None = None,
    exclude_qcfail: bool = False,
    rc_reverse: bool = False,
    ofp: str | None = None,
    is_paired_end: bool = False,
    limit: int | None = None
) -> Generator[YieldT, SendT, LibraryIndependentStats]:
    """
    This function is for the initial collation of unique read sequences in the original orientation only (hts will do revcomp).
    There is no chunking of data, so this relies on large memory lookups at present.

    Selecting correct underlying parser is via file extension:
    - cram/bam/sam -> htslib processing
    """

    start = time()

    match read_file_info.fmt:

        case ReadFileFormat.HTS:
            logging.info(f"Sequence input detected as *{read_file_info.ext}")
            hts_reads = parse_htsfile(
                read_file_info,
                sample,
                cpus,
                reference=reference,
                exclude_qcfail=exclude_qcfail,
                rc_reverse=rc_reverse,
                ofp=ofp,
                limit=(limit * 2) if (is_paired_end and limit is not None) else limit,
                iter_reads=iter_paired_end if is_paired_end else iter_single_end)

            stats = yield from hts_reads

        case _:
            raise RuntimeError("Invalid reads format!")

    assert stats

    logging.info(f"Processed {stats.input_reads} reads in {int(time() - start)} s")

    return stats
