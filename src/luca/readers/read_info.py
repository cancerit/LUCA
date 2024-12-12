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


@dataclass(slots=True)
class ReadInfo:
    sequence: str
    is_masked: bool
    is_ambiguous: bool
    is_qc_fail: bool
    is_short: bool
    is_empty: bool
    is_flagged: bool

    @classmethod
    def empty(cls, is_qc_fail: bool, is_flagged: bool):
        return cls('', False, False, is_qc_fail, False, True, is_flagged)

    def to_discard(self, discard_qc: bool) -> bool:
        return (
            self.is_masked or
            self.is_short or
            self.is_empty or
            self.is_flagged or (
                discard_qc and self.is_qc_fail
            )
        )

    @property
    def is_sequence_invalid(self) -> bool:
        return self.is_ambiguous or self.is_masked
