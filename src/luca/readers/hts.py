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
# Adapted from pyQUEST (readers.hts module)

import logging
from contextlib import contextmanager, nullcontext
from dataclasses import dataclass
from typing import Callable, Final, Generator, Literal

import pysam
from pysam.libcalignedsegment import SAM_FLAGS

from ..errors import InputReadError, InvalidHTSError, MissingMetadataError
from ..stats import LibraryIndependentStats
from .constants import LOAD_INFO_THRESHOLD
from .read_file_info import ReadFileInfo
from .read_info import ReadInfo
from .read_qc import non_dna_re, non_non_amb_dna_re
from .types import IterReadT, SendReadT
from ..utils import reverse_complement

CIGAR_SOFT_CLIPPING: Final[Literal['S']] = 'S'
MAX_HTS_CPUS = 4
SKIP_READ_FLAG: Final[int] = (
    SAM_FLAGS.FSECONDARY |
    SAM_FLAGS.FSUPPLEMENTARY
)

SAM_READ_GROUP = 'RG'
SAM_SAMPLE = 'SM'


def _cap_hts_cpus(cpus: int) -> int:
    return cpus if cpus < MAX_HTS_CPUS else MAX_HTS_CPUS


def _get_hts_sample_name(sam: pysam.AlignmentFile) -> str | None:
    header = sam.header.to_dict()
    sample_names: set[str] = {
        sample_name
        for sample_name in {
            rg[SAM_SAMPLE]
            for rg in header.get(SAM_READ_GROUP, [])
            if SAM_SAMPLE in rg
        }
        if sample_name
    }

    match len(sample_names):
        case 0:
            return None
        case 1:
            return sample_names.pop()
        case _:
            raise InvalidHTSError("Multiple different sample names found in header")


@contextmanager
def hts_reader(seq_file: str, mode, cpus: int, reference: str | None = None) -> Generator[pysam.AlignmentFile, None, None]:
    verbosity = pysam.set_verbosity(0)
    try:
        sam = pysam.AlignmentFile(
            seq_file,
            mode=mode,
            reference_filename=reference,
            require_index=False,
            threads=_cap_hts_cpus(cpus),
            check_sq=False)
    except ValueError as ex:
        raise InvalidHTSError("HTS file error: " + ex.args[0])

    pysam.set_verbosity(verbosity)
    try:
        yield sam
    finally:
        sam.close()


def get_hts_read_info(skip_read_flag: int, rc_reverse: bool, read: pysam.AlignedSegment) -> ReadInfo:

    # Fetch sequence
    seq: str | None = read.get_forward_sequence() if rc_reverse else read.query_sequence
    is_flagged = read.flag & skip_read_flag != 0

    if seq:
        # Does the sequence contain anything other than unambiguous nucleotide symbols?
        if non_non_amb_dna_re.search(seq):
            # Does the sequence contain anything other than nucleotide symbols?
            if non_dna_re.search(seq):
                raise InputReadError(f"Invalid sequence: '{seq}'!")
            is_ambiguous = True
        else:
            is_ambiguous = False

        is_masked = (
            read.cigarstring is not None and
            CIGAR_SOFT_CLIPPING in read.cigarstring
        )

        return ReadInfo(
            sequence=seq,
            is_flagged=is_flagged,
            is_qc_fail=read.is_qcfail,
            is_ambiguous=is_ambiguous,
            is_masked=is_masked,
            is_short=False,
            is_empty=False)

    else:
        return ReadInfo.empty(read.is_qcfail, is_flagged)


def eval_read(exclude_qcfail: bool, skip_read_flag: int, rc_reverse: bool, stats: LibraryIndependentStats, read: pysam.AlignedSegment) -> str:
    read_info: ReadInfo = get_hts_read_info(skip_read_flag, rc_reverse, read)
    stats.eval_read_info(exclude_qcfail, read_info)
    return read_info.sequence


def iter_single_end(exclude_qcfail: bool, skip_read_flag: int, rc_reverse: bool, has_output: bool, stats: LibraryIndependentStats, sam: pysam.AlignmentFile, out_bam: pysam.AlignmentFile | None, n: int | None) -> Generator[str, bool, None]:
    for read in sam.fetch(until_eof=True):
        seq = eval_read(exclude_qcfail, skip_read_flag, rc_reverse, stats, read)

        should_copy = (yield seq) and has_output
        if should_copy:
            out_bam.write(read)  # type: ignore

        if n is not None and stats.total_reads >= n:
            logging.warning("Terminating after parsing %d reads as per user request (--limit option)." % stats.total_reads)
            return

        if stats.total_reads % LOAD_INFO_THRESHOLD == 0:  # pragma: no cover
            logging.debug(f"Parsed {stats.total_reads} reads...")


def iter_paired_end(exclude_qcfail: bool, skip_read_flag: int, rc_reverse: bool, has_output: bool, stats: LibraryIndependentStats, sam: pysam.AlignmentFile, out_bam: pysam.AlignmentFile | None, n: int | None) -> Generator[tuple[str, str], tuple[bool, bool], None]:
    if rc_reverse:
        raise NotImplementedError("Forcing reverse complement not yet defined for paired-end sequencing data!")

    read_1: str = ''
    read_2: str = ''

    for read in sam.fetch(until_eof=True):
        # TODO: add stats about pairs?
        seq = eval_read(exclude_qcfail, skip_read_flag, rc_reverse, stats, read)

        if not read.query_name:
            continue

        # Assumption: read 2 always follows read 1 (samtools collate should guarantee that)
        if read.is_read1:
            read_1 = seq
            continue
        elif read.is_read2:
            # TODO: test QNAME, existing read 1?
            read_2 = seq
            should_copy = (yield read_1, read_2) and has_output
        else:
            # TODO: handle orphans
            continue

        # TODO: the mismatch status should be per mate
        # TODO: should the output be also collated with read 1 and 2 in the same order as in the original?
        if should_copy:
            out_bam.write(read)  # type: ignore

        if n is not None and stats.total_reads >= n:
            logging.warning("Terminating after parsing %d read pairs as per user request (--limit option)." % (stats.total_reads // 2))
            return

        if stats.total_reads % LOAD_INFO_THRESHOLD == 0:  # pragma: no cover
            logging.debug(f"Parsed {stats.total_reads} reads...")


IterReadSig = Callable[
    [bool, int, bool, bool, LibraryIndependentStats, pysam.AlignmentFile, pysam.AlignmentFile | None, int | None],
    Generator[IterReadT, SendReadT, None]
]


def _get_sample_name(sam: pysam.AlignmentFile, user_sample_name: str | None) -> str:
    read_group_sample_name = _get_hts_sample_name(sam)

    if user_sample_name:

        # User-provided sample name
        if read_group_sample_name and user_sample_name != read_group_sample_name:
            logging.warning(
                "Conflicting sample names: expected '%s', found '%s' in header!" %
                (user_sample_name, read_group_sample_name))

        return user_sample_name

    elif read_group_sample_name:

        # Read group sample name
        return read_group_sample_name

    else:
        raise MissingMetadataError(
            "No sample name found in input file header, please provide via '--sample'!")


def parse_htsfile(
    read_file_info: ReadFileInfo,
    sample: str | None,
    cpus: int,
    reference: str | None = None,
    exclude_qcfail: bool = False,
    rc_reverse: bool = False,
    ofp: str | None = None,
    iter_reads: IterReadSig = iter_single_end,
    limit: int | None = None
) -> Generator[IterReadT, bool, LibraryIndependentStats]:
    stats: LibraryIndependentStats
    skip_read_flag: int = (
        SKIP_READ_FLAG if not exclude_qcfail else
        (SKIP_READ_FLAG | SAM_FLAGS.FQCFAIL)
    ).value

    has_output = bool(ofp)

    with hts_reader(read_file_info.fp, read_file_info.read_mode, cpus, reference) as sam:
        sample_name = _get_sample_name(sam, sample)
        stats = LibraryIndependentStats.empty(sample_name)

        # TODO: check input CRAM can be used as template for an output BAM
        with pysam.AlignmentFile(ofp, mode='wb', template=sam) if ofp else nullcontext() as out_bam:

            try:
                yield from iter_reads(exclude_qcfail, skip_read_flag, rc_reverse, has_output, stats, sam, out_bam, limit)

            except OSError as ex:
                raise InvalidHTSError("HTS file error: " + ex.args[0])

    return stats
