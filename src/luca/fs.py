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

import gzip
import logging
import os
from contextlib import contextmanager, nullcontext
from dataclasses import dataclass
from typing import IO, Any, Iterable


PARTIAL_FILE_EXT = 'part'
PARTIAL_FILE_SUFFIX = '.' + PARTIAL_FILE_EXT

EXT_COMPRESSED = '.gz'


@dataclass(slots=True)
class OutputFile:
    fp: str
    fh: IO[Any] | gzip.GzipFile


@contextmanager
def open_input_file(fp: str):
    is_compressed = fp.endswith(EXT_COMPRESSED)
    with (gzip.open if is_compressed else open)(fp, 'rt') as fh:
        yield fh


def get_partial_file_path(fp: str) -> str:
    return fp + PARTIAL_FILE_SUFFIX


@contextmanager
def open_output_file(fp: str, out_dir: str | None = None, mode: str = 'wt', compress: bool = False):
    """
    Open file with partial extension, and delete it on close if empty

    Assumptions:
    - writing in text format when in compression mode
    """

    rfp = fp
    if compress:
        rfp += EXT_COMPRESSED
    fp = rfp if not out_dir else os.path.join(out_dir, rfp)
    fp_part = get_partial_file_path(fp)

    def rename_part():
        if os.path.getsize(fp_part) > 0:
            logging.debug("Renaming partial file '%s'...", fp_part)
            os.rename(fp_part, fp)

    try:
        with (
            (gzip.open if compress else open)(fp_part, mode)
        ) as fh:
            yield OutputFile(rfp, fh)

        # Rename partial file (if no iteration stop has been raised)
        rename_part()

    except StopIteration:

        # Rename partial file
        rename_part()

        # Allow the caller to handle the iteration stop
        raise

    finally:
        # TODO: should deletion be restricted to expected exceptions,
        #  for debugging purposes?
        if os.path.isfile(fp_part):
            logging.debug("Deleting partial file '%s'...", fp_part)
            os.remove(fp_part)


@contextmanager
def open_opt_output_file(fp: str, out_dir: str | None = None, mode: str = 'wt', compress: bool = False, no_op: bool = False):
    with (
        open_output_file(fp, out_dir=out_dir, mode=mode, compress=compress) if not no_op else
        nullcontext()
    ) as fh:
        yield fh


def get_match_info_file_name(i: int = 0) -> str:
    # E.g.: matches.0.tsv
    return f"matches.{i}.tsv"


def get_mm_match_info_file_name(i: int = 0) -> str:
    # E.g.: mm.matches.0.tsv
    return f"mm.{get_match_info_file_name(i=i)}"


def get_swap_match_info_file_name(i: int = 0) -> str:
    # E.g.: swap.matches.0.tsv
    return f"swap.{get_match_info_file_name(i=i)}"


def get_library_stats_file_name(library_index: int) -> str:
    return f"lib.{library_index}.stats.json"


def get_library_counts_file_name(library_index: int) -> str:
    return f"lib.{library_index}.counts.tsv"


def get_combination_counts_file_name(combination_index: int) -> str:
    # TODO: consider using the combination identifier instead
    return f"combination.{combination_index}.counts.tsv"


def get_read_counts_file_name(prefix: str | None = None) -> str:
    base = "read.counts.tsv"
    return f"{prefix}.{base}" if prefix else base


def get_abs_file_path(fp: str, root_dir: str | None = None) -> str:
    return (
        fp if os.path.isabs(fp) or not root_dir else
        os.path.join(root_dir, fp)
    )


def get_unique_abs_file_paths(fps: Iterable[str], root_dir: str | None = None) -> set[str]:
    return {get_abs_file_path(fp, root_dir=root_dir) for fp in fps}


def does_file_exist(fp: str, root_dir: str | None = None) -> bool:
    afp = get_abs_file_path(fp, root_dir=root_dir)
    if os.path.isfile(afp):
        return True
    else:
        logging.error("File not found: '%s'!", afp)
        return False
