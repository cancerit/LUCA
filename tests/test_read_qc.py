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

import pytest

from luca.readers.read_qc import non_dna_re, non_non_amb_dna_re


NON_DNA = 0
NON_AMB_DNA = 1
AMB_DNA = 2


def classify_dna(seq: str) -> int:
    # NOTE: marginally faster than a single expression
    #  in the non-ambiguous DNA case.

    # Non-ambiguous DNA (most common case)
    if not non_non_amb_dna_re.search(seq):
        return NON_AMB_DNA
    return (
        NON_DNA if non_dna_re.search(seq) else
        AMB_DNA
    )


@pytest.mark.parametrize('seq,cls', [
    ('ACGT', NON_AMB_DNA),
    ('NACGT', AMB_DNA),
    ('ANCGT', AMB_DNA),
    ('ACXGT', NON_DNA),
    ('XACGT', NON_DNA)
])
def test_classify_dna(seq, cls):
    assert classify_dna(seq) == cls
