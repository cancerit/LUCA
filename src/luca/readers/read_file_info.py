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

from dataclasses import dataclass
from enum import Enum
import os
from typing import Final

from ..errors import InvalidReadFileExtension

EXT_TO_HTS: Final[dict[str, str]] = {
    '.bam': 'rb',
    '.cram': 'rc',
    '.sam': 'r',
}


class ReadFileFormat(Enum):
    HTS = 1


@dataclass(slots=True)
class ReadFileInfo:
    fp: str
    ext: str
    fmt: ReadFileFormat
    read_mode: str

    @classmethod
    def from_path(cls, fp: str):
        _, ext = os.path.splitext(fp)

        if ext not in EXT_TO_HTS:
            raise InvalidReadFileExtension(f"Unsupported file extension *{ext}!")

        return cls(fp, ext, ReadFileFormat.HTS, EXT_TO_HTS.get(ext, 'r'))
