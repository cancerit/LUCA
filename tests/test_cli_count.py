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

import json
from tempfile import NamedTemporaryFile

from click.testing import CliRunner
import pytest

from luca.count import count


def get_experiment(lib_fp):
    return {
        'sequencing_type': 'single_end',
        'libraries': [
            {'id': lib_fp}
        ],
        'read_templates': [
            {
                'id': 'default',
                'regions': [{
                    'id': 'guide',
                    'libraries': [lib_fp]
                }]
            }
        ],
        'read_group_templates': {
            'default': ['default']
        }
    }


def test_count_help():
    runner = CliRunner()
    result = runner.invoke(count, ['--help'])
    assert result.exit_code == 0


def test_count_tsv_error():
    runner = CliRunner()
    with (
        NamedTemporaryFile() as exp,
        NamedTemporaryFile() as tsv,
        NamedTemporaryFile(suffix='.bam') as bam
    ):
        # Write experiment to JSON
        with open(exp.name, 'w') as fh:
            e = get_experiment(tsv.name)
            json.dump(e, fh)

        # Write truncated TSV
        with open(tsv.name, 'w') as fh:
            fh.write('id\tsequence\n0\n')

        # Run counting entry point
        result = runner.invoke(count, [
            '-l', '.',
            '-o', '.',
            exp.name,
            bam.name
        ])

    assert result.exit_code == 1
