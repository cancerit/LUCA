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

import csv
from typing import Generator

from ..fs import open_input_file
from ..utils import has_duplicates


class TsvError(Exception):
    def __init__(self, fp: str, *args) -> None:
        self.fp = fp
        super().__init__(*args)


class TsvFieldNotFound(TsvError):
    def __init__(self, fp: str, field_name: str, *args) -> None:
        self.field_name = field_name
        super().__init__(fp, *args)


class TsvTruncatedRow(TsvError):
    pass


class TsvEmpty(TsvError):
    pass


def _parse_tsv(fp: str, reader, indices: list[int]) -> Generator[list[str], None, None]:
    try:

        # Emit rows
        for line in reader:
            yield [line[i] for i in indices]

    except IndexError:
        raise TsvTruncatedRow(fp)


def parse_tsv_at(fp: str, indices: list[int]) -> Generator[list[str], None, None]:
    with open_input_file(fp) as fh:
        reader = csv.reader(fh, delimiter='\t')
        yield from _parse_tsv(fp, reader, indices)


def parse_tsv(fp: str, fields: list[str]) -> Generator[list[str], None, None]:
    assert not has_duplicates(fields)
    with open(fp) as fh:
        reader = csv.reader(fh, delimiter='\t')

        try:

            # Parse header
            header = next(reader)

        except StopIteration:
            raise TsvEmpty(fp)

        column_headers = {s: i for i, s in enumerate(header)}
        try:
            indices = [column_headers[f] for f in fields]
        except KeyError as ex:
            raise TsvFieldNotFound(fp, ex.args[0])

        yield from _parse_tsv(fp, reader, indices)
