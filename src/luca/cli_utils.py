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
import sys
from typing import NoReturn

import click


existing_file_path = click.Path(
    exists=True,
    file_okay=True,
    dir_okay=False,
    readable=True,
    resolve_path=True)


existing_directory_path = click.Path(
    exists=True,
    file_okay=False,
    dir_okay=True,
    readable=True,
    resolve_path=True)


existing_writable_directory_path = click.Path(
    exists=True,
    file_okay=False,
    dir_okay=True,
    readable=True,
    writable=True,
    resolve_path=True)


def abort(message: str, *args, level: int = logging.CRITICAL) -> NoReturn:
    logging.log(level, message, *args)
    sys.exit(1)
