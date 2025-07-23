# LUCA
#
# Copyright (C) 2024, 2025 Genome Research Ltd.
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
import os

from .input_manifest import InputManifest


@dataclass(slots=True, frozen=True)
class PathBundle:
    output_dir: str
    library_dir: str | None
    input_manifest: InputManifest | None

    def get_library_file_path(self, id: str) -> str:
        if self.input_manifest and id in self.input_manifest.libraries:
            return self.input_manifest.get_library_file_path(id, root_dir=self.library_dir)
        elif self.library_dir:
            return os.path.join(self.library_dir, id)
        else:
            return id

    def get_combination_filter_file_path(self, id: str) -> str:
        if self.input_manifest and id in self.input_manifest.combination_filters:
            return self.input_manifest.get_combination_filter_file_path(id, root_dir=self.library_dir)
        elif self.library_dir:
            return os.path.join(self.library_dir, id)
        else:
            return id

    @property
    def is_valid(self) -> bool:
        return (
            not self.input_manifest or
            self.input_manifest.is_valid(root_dir=self.library_dir)
        )
