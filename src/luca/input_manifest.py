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

import logging
import os

from pydantic import BaseModel, Field

from .fs import does_file_exist, get_unique_abs_file_paths
from .utils import greedy_all


def _get_file_path(d: dict[str, str], k: str, root_dir: str | None) -> str:
    # TODO: handle absolute paths
    fp = d[k]
    return os.path.join(root_dir, fp) if root_dir else fp


class InputManifest(BaseModel):
    libraries: dict[str, str] = Field(default_factory=dict)
    combination_filters: dict[str, str] = Field(default_factory=dict)

    @property
    def file_paths(self) -> set[str]:
        return (
            set(self.libraries.values()) |
            set(self.combination_filters.values())
        )

    @property
    def are_all_library_file_paths_absolute(self) -> bool:
        return all(os.path.isabs(x) for x in self.libraries.values())

    def do_all_files_exist(self, root_dir: str | None = None) -> bool:
        return greedy_all(
            lambda x: does_file_exist(x, root_dir=root_dir), self.file_paths)

    def get_library_file_path(self, id: str, root_dir: str | None = None) -> str:
        return _get_file_path(self.libraries, id, root_dir)

    def get_combination_filter_file_path(self, id: str, root_dir: str | None = None) -> str:
        return _get_file_path(self.combination_filters, id, root_dir)

    def is_valid(self, root_dir: str | None = None) -> bool:
        if not self.libraries and not self.combination_filters:
            logging.warning(
                "Input manifest lists neither libraries (`libraries`) " +
                "nor combination filters (`combination_filters`)!")

        def get_unique_paths(d: dict[str, str]) -> set[str]:
            return get_unique_abs_file_paths(d.values(), root_dir=root_dir)

        shared_file_paths = (
            get_unique_paths(self.libraries) &
            get_unique_paths(self.combination_filters)
        )

        for fp in shared_file_paths:
            logging.error(
                "Identifiers for both a library and a combination filter " +
                "map to the same file path: '%s'!" % fp)

        return not bool(shared_file_paths)
