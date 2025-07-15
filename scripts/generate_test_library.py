#!/usr/bin/env python3

# LUCA
#
# Copyright (C) 2025 Genome Research Ltd.
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

"""
Generate test library with randomised sequences
"""

import os
import random
import sys
from argparse import ArgumentParser
from dataclasses import dataclass, field


@dataclass(slots=True)
class DataGenOpt:
    sequence_number: int
    sequence_length: int


def get_random_sequence(k: int) -> str:
    return ''.join(random.choices('ACGT', k=k))


def main(opt: DataGenOpt, out_fp: str, seed: int | None = None):
    if seed is not None:
        random.seed(seed)

    if opt.sequence_number > 4 ** opt.sequence_length:
        sys.exit(
            f"Number of sequences ({opt.sequence_number}) is too large"
            " to generate the required number of unique sequences "
            f"of the required length ({opt.sequence_length})!")

    sequences: set[str] = set()
    while len(sequences) < opt.sequence_number:
        sequences.add(get_random_sequence(opt.sequence_length))

    with open(out_fp, 'w') as fh:
        fh.write('id\tsequence\n')
        for i, sequence in enumerate(sorted(sequences)):
            fh.write(str(i))
            fh.write('\t')
            fh.write(sequence)
            fh.write('\n')


if __name__ == '__main__':
    p = ArgumentParser(description="Generate synthetic CRISPR library (LUCA format)")
    p.add_argument('-o', '--output', required=True, help="Output library file path")
    p.add_argument('-l', '--length', required=True, type=int, help="Sequence length")
    p.add_argument('-n', '--number', required=True, type=int, help="Number of sequences")
    p.add_argument('-s', '--seed', type=int, help="Random seed for reproducible generation")
    args = p.parse_args()

    if os.path.isdir(args.output):
        sys.exit("Output path is a directory!")

    opt = DataGenOpt(sequence_number=args.number, sequence_length=args.length)
    if opt.sequence_number <= 0:
        sys.exit("Invalid number of sequences!")
    if opt.sequence_length <= 0:
        sys.exit("Invalid sequence length!")

    main(opt, args.output, seed=args.seed)
