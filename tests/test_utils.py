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

from luca.utils import clamp_non_negative, reverse_complement


@pytest.mark.parametrize('seq,rc', [
    ('ACGT', 'ACGT'),
    ('ACCT', 'AGGT'),
    ('GAT', 'ATC')
])
def test_reverse_complement(seq, rc):
    assert reverse_complement(seq) == rc


@pytest.mark.parametrize('n,m', [
    (0, 0),
    (-1, 0),
    (1, 1)
])
def test_clamp_non_negative(n, m):
    assert clamp_non_negative(n) == m
