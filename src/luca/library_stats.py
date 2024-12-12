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

from typing import Any, Iterable

from pydantic import BaseModel

from .utils import load_text_from_file


def _round_float(x: float) -> float:
    return round(x, ndigits=2)


class LibraryStats(BaseModel):

    # Reads
    mapped_to_template_reads: int
    mean_count_per_template: float
    median_count_per_template: float

    # Swap match
    swap_matching_reads: int

    # Templates
    total_templates: int
    zero_count_templates: int
    low_count_templates_lt_15: int
    low_count_templates_lt_30: int

    gini_coefficient: float

    @classmethod
    def empty(cls):
        return cls(
            mapped_to_template_reads=0,
            mean_count_per_template=0.0,
            median_count_per_template=0.0,
            swap_matching_reads=0,
            total_templates=0,
            zero_count_templates=0,
            low_count_templates_lt_15=0,
            low_count_templates_lt_30=0,
            gini_coefficient=0.0)

    @classmethod
    def load(cls, fp: str):
        return cls.model_validate_json(
            load_text_from_file(fp))

    def to_dict(self) -> dict[str, Any]:
        return {
            'mapped_to_template_reads': self.mapped_to_template_reads,
            'mean_count_per_template': _round_float(self.mean_count_per_template),
            'median_count_per_template': _round_float(self.median_count_per_template),
            'total_templates': self.total_templates,
            'swap_matching_reads': self.swap_matching_reads,
            'zero_count_templates': self.zero_count_templates,
            'low_count_templates_lt_15': self.low_count_templates_lt_15,
            'low_count_templates_lt_30': self.low_count_templates_lt_30,
            'gini_coefficient': _round_float(self.gini_coefficient)
        }
