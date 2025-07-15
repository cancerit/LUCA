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

from tempfile import NamedTemporaryFile

import pytest

from luca.readers.tsv import parse_tsv, TsvFieldNotFound, TsvTruncatedRow, TsvEmpty


def test_tsv_empty():
    with NamedTemporaryFile() as tmp:
        with pytest.raises(TsvEmpty):
            list(parse_tsv(tmp.name, ['a']))


def test_tsv_field_not_found():
    with NamedTemporaryFile() as tmp:
        with open(tmp.name, 'w') as fh:
            fh.write("a\tb\tc\n")
        with pytest.raises(TsvFieldNotFound):
            list(parse_tsv(tmp.name, ['MISSING_FIELD']))


def test_tsv_header_only():
    with NamedTemporaryFile() as tmp:
        with open(tmp.name, 'w') as fh:
            fh.write("a\tb\tc\n")
        data = list(parse_tsv(tmp.name, ['a']))
        assert data == []


def test_tsv_truncated():
    with NamedTemporaryFile() as tmp:
        with open(tmp.name, 'w') as fh:
            fh.write("a\tb\tc\n")
            fh.write("x\ty\tz\n")
            fh.write("x\n")
        with pytest.raises(TsvTruncatedRow):
            list(parse_tsv(tmp.name, ['b']))
