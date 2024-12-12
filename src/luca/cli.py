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

from functools import wraps

import click
from click_option_group import OptionGroup

from .cli_utils import existing_writable_directory_path


HELP_OUTPUT = "Output directory"


LOG_LEVELS = [
    'WARNING',
    'INFO',
    'DEBUG'
]

debug_opts = OptionGroup("\nDebug", help="Options specific to troubleshooting, testing and debugging")


def common_cli_options(f):
    @click.option(
        '-o',
        '--output',
        required=True,
        type=existing_writable_directory_path,
        help=HELP_OUTPUT
    )
    @wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def common_cli_debug_options(f):
    @debug_opts.option(
        '--loglevel',
        default='INFO',
        show_default=True,
        type=click.Choice(LOG_LEVELS, case_sensitive=False),
        help="Set logging verbosity"
    )
    @wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper
