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

from collections import Counter
from itertools import groupby
import json
import logging
import re
from typing import Callable, Iterable, TypeVar

import numpy as np
from pydantic import ValidationError

from .errors import UnsupportedData

T = TypeVar('T')
K = TypeVar('K')
V = TypeVar('V')

dna_re = re.compile('^[ACGT]*$')
dna_complement_tr_table = str.maketrans('ACGT', 'TGCA')
reverse_sl = slice(None, None, -1)


def load_text_from_file(fp: str) -> str:
    with open(fp, 'rt') as fh:
        return fh.read()


def clamp_non_negative(n: int) -> int:
    return n if n >= 0 else 0


def is_dna(s: str) -> bool:
    return dna_re.match(s) is not None


def reverse_complement(seq: str) -> str:
    return seq[reverse_sl].translate(dna_complement_tr_table)


def get_unique_length(a: list) -> int:
    return len(set(a))


def has_duplicates(a: list) -> bool:
    return get_unique_length(a) < len(a)


def get_range_from(start: int, length: int) -> range:
    return range(start, start + length)


def get_offsets(lengths: list[int]) -> list[int]:
    n = len(lengths)
    assert n > 0
    offsets = [0] + np.asarray(lengths, dtype=np.int64).cumsum(dtype=np.int64).tolist()
    return offsets


def greedy_any(f: Callable[[T], bool], a: Iterable[T]) -> bool:
    success = False
    for item in a:
        if f(item):
            success = True
    return success


def greedy_all(f: Callable[[T], bool], a: Iterable[T]) -> bool:
    return not greedy_any(lambda x: not f(x), a)


def sort_dict_by_value(d: dict[K, V], reverse: bool = False) -> list[tuple[K, V]]:
    items = [(k, v) for k, v in d.items()]
    items.sort(key=lambda t: t[1], reverse=reverse)
    return items


def sort_counter_desc(c: Counter[K]) -> list[tuple[K, int]]:
    return sort_dict_by_value(c, reverse=True)


def flip_dict(d: dict[K, V]) -> dict[V, K]:
    return {v: k for k, v in d.items()}


def get_median_from_sorted(sa: np.ndarray) -> float:
    """
    Calculate the median, assuming the input is sorted
    """

    n: int = sa.shape[0]

    if n == 0:
        raise UnsupportedData("Median: empty array!")

    if n == 1:
        return float(sa[0])

    # Even length
    if n % 2 == 0:
        i: int = n // 2
        return float(sa[[i - 1, i]].mean())

    # Odd length
    return float(sa[(n - 1) // 2])


def get_stats(sa: np.ndarray, gini_corr: bool = False) -> tuple[int, float, float, float]:
    """
    Calculate the following statistics, assuming the array is sorted:
    - total
    - mean
    - median
    - Gini coefficient
    """

    n: int = sa.shape[0]

    if n == 0:
        raise UnsupportedData("Stats calculation: empty array!")

    m: np.uint64 = sa.sum()

    if m == 0:
        logging.warning("No library matches!")
        return 0, 0.0, 0.0, 0.0

    # Mean
    mean: float = int(m) / n

    # Median
    median: float = get_median_from_sorted(sa)

    # Gini coefficient
    gini_coeff: float = float(
        (2 * np.sum(sa * np.arange(1, n + 1)) / m - (n + 1)) /
        ((n - 1) if gini_corr else n)
    )

    return int(m), mean, median, gini_coeff


def get_set(ls: list[T], f: Callable[[T], None]) -> tuple[set[T], bool]:
    s = set(ls)
    had_duplicates = False
    if len(s) < len(ls):
        had_duplicates = True
        ls.sort()
        for x, xs in groupby(ls):
            if len(list(xs)) > 1:
                f(x)

    return s, had_duplicates


def write_json(d: dict, fp: str) -> None:
    with open(fp, 'w') as fh:
        json.dump(d, fh)


def log_validation_error(ex: ValidationError) -> None:
    for err in ex.errors(include_url=False):
        logging.error(f"{err['msg']} ({err['type']}, at: {json.dumps(err['input'])})")
    if ex.args:
        logging.error(ex.args)


def eval_fail_conditions(a: list[tuple[bool, str]], f: Callable[[str], None]) -> bool:
    success = True
    for c, msg in a:
        if c:
            f(msg)
            success = False
    return success


def union(sets: list[set]) -> set:
    return set.union(*sets) if sets else set()
