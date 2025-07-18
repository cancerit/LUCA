# LUCA
#
# Copyright (C) 2025 Genome Research Ltd.
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

import charset_normalizer
import logging


MAX_BYTES = 10000


def detect_encoding(fp: str) -> str | None:
    with open(fp, 'rb') as rfh:
        data = rfh.read(MAX_BYTES)

    if not data:
        return None

    encoding = charset_normalizer.detect(data)['encoding']
    logging.debug("File '%s' encoding: %s." % (fp, encoding))
    return encoding
