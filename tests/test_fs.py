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

import gzip
import os
from tempfile import TemporaryDirectory

from luca.fs import open_output_file


def test_open_output_file():
    s = 'ABC\n'

    def write_to_file(compress):
        with open_output_file(fp, compress=compress) as f:
            f.fh.write(s)

    with TemporaryDirectory() as out_dir:
        fp = os.path.join(out_dir, 'a.txt')
        cfp = fp + '.gz'

        write_to_file(True)
        write_to_file(False)
        assert os.path.isfile(fp)
        assert os.path.isfile(cfp)
        with open(fp) as ufh, gzip.open(cfp, 'rt') as cfh:
            assert ufh.read() == cfh.read()
