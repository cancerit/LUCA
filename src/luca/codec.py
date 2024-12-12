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

from dataclasses import dataclass, field
from typing import Sized


@dataclass(slots=True)
class Codec(Sized):
    start: int = 0
    _encoder: dict[str, int] = field(init=False, default_factory=dict)
    _decoder: list[str] = field(init=False, default_factory=list)

    def __len__(self) -> int:
        return len(self._encoder)

    def __contains__(self, item: str) -> bool:
        return item in self._encoder

    @property
    def encoder(self) -> dict[str, int]:
        return self._encoder

    @property
    def items(self) -> list[str]:
        return self._decoder

    def encode(self, item: str) -> int:
        return self._encoder[item]

    def decode(self, i: int) -> str:
        return self._decoder[i - self.start]

    def add(self, item: str) -> int:
        if item in self._decoder:
            return self.encode(item)
        else:
            i = self.start + len(self._decoder)
            self._encoder[item] = i
            self._decoder.append(item)
            return i
