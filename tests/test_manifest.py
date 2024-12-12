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

import os
from tempfile import NamedTemporaryFile, TemporaryDirectory

import pytest

from luca.experiment import CombinationInfo, CombinationRegion, Experiment, LibraryInfo, ReadGroupId, ReadGroupOptions, ReadRegion, ReadTemplate, SequencingType
from luca.input_manifest import InputManifest
from luca.path_bundle import PathBundle
from luca.count import validate_path_bundle


def test_path_bundle_is_valid():
    fn = 'a.tsv'
    with TemporaryDirectory() as temp_dir:

        # Empty manifest
        b = PathBundle('.', temp_dir, InputManifest())
        assert b.input_manifest
        assert b.is_valid

        # Manifest with different entities mapping to the same path
        b.input_manifest.libraries['a'] = fn
        b.input_manifest.combination_filters['a'] = os.path.join(temp_dir, fn)
        assert not b.is_valid
        b.input_manifest.libraries['a'] += '.x'
        assert b.is_valid


def test_validate_path_bundle():
    m = InputManifest()
    b = PathBundle('.', None, m)
    assert b.input_manifest

    exp = Experiment(
        sequencing_type=SequencingType.SINGLE_END,
        libraries=[
            LibraryInfo(id='embedded', values=['AAA']),
            LibraryInfo(id='external', values=None)
        ],
        read_templates=[
            ReadTemplate(id='default', regions=[
                ReadRegion(id='default', libraries=['external'])
            ])
        ],
        read_group_templates={
            ReadGroupId.DEFAULT: ['default']
        },
        read_groups={
            ReadGroupId.DEFAULT: ReadGroupOptions(is_reverse=False)
        },
        combinations=[
            CombinationInfo(
                id='default',
                regions=[
                    CombinationRegion(id='default', filter=True)
                ],
                filters=['ft1'])
        ])

    # Non existing input files
    with pytest.raises(SystemExit):
        validate_path_bundle(b, exp)

    # Existing input files
    with NamedTemporaryFile() as f, NamedTemporaryFile() as g:

        # Library and combination filter pointing to the same file path
        b.input_manifest.libraries['external'] = f.name
        b.input_manifest.combination_filters['ft1'] = f.name
        with pytest.raises(SystemExit):
            validate_path_bundle(b, exp)

        # Distinct file paths
        b.input_manifest.combination_filters['ft1'] = g.name
        validate_path_bundle(b, exp)
