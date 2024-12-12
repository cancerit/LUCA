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

import json
import logging

from .app_info import AppInfo
from .stats import LibraryIndependentStats


def write_stats(app_info: AppInfo, stats: LibraryIndependentStats, fp: str) -> None:
    logging.info(f"Writing statistics file: {fp}")
    with open(fp, 'w') as fh:
        json.dump(stats.to_dict(app_info), fh)
