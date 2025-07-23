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

from tempfile import NamedTemporaryFile

import pytest

from luca.errors import InvalidLibraryError
from luca.experiment import LibraryReverseCondition
from luca.library import DynamicMultiTargetLibrary, LibraryBuilder, MonoTargetLibrary, \
    MultiTargetLibrary, \
    MultiTargetUniformLibrary, RawLibrary

from .utils import mock_reads


def get_test_library():
    return MultiTargetLibrary(
        0, ['A', 'BBB', 'C'], {1: {'A': 0, 'C': 2}, 3: {'BBB': 1}})


def test_library_len():
    library = get_test_library()
    assert len(library) == 3


def test_library_from_sequences():
    library = MultiTargetLibrary.from_sequences(['A', 'BBB', 'C'])
    assert len(library) == 3
    assert len(library.length_targets) == 2
    for length in library._lengths:
        for target in library.length_targets[length]:
            assert len(target) == length


@pytest.mark.parametrize('seq,offset,exp_target,exp_length', [
    ('A', 0, 0, 1),
    ('XXXA', 3, 0, 1),
    ('XXBBB', 2, 1, 3)
])
def test_library_match(seq, offset, exp_target, exp_length):
    library = get_test_library()
    target, length = library.match(seq, offset=offset, anchor_left=True)
    assert target == exp_target
    assert length == exp_length


@pytest.mark.parametrize('seq,offset,exp_target,exp_length', [
    ('A', 0, 0, 1),
    ('AXXX', 3, 0, 1),
    ('BBBXX', 2, 1, 3)
])
def test_library_match_reverse(seq, offset, exp_target, exp_length):
    library = get_test_library()
    target, length = library.match(seq, offset=offset, anchor_left=False)
    assert target == exp_target
    assert length == exp_length


@pytest.mark.parametrize('anchor_left,query', [
    (True, 'XXAAA'),
    (False, 'AAAXX')
])
def test_mono_library(anchor_left, query):
    target = 'AAA'
    target_id, match_length = MonoTargetLibrary(0, target).match(query, offset=2, anchor_left=anchor_left)
    assert target_id == 0
    assert match_length == len(target)


def test_single_end_matcher_multi_target():
    rl = RawLibrary(0, False, ['AAA', 'CCC', 'GGG'])
    library = rl.build()
    targets = set(library.sequences)
    mm_reads = []
    for read in mock_reads([
        'AAAC',
        'AAAG',
        'CCCG',
        'TTTA'
    ], mm_reads):
        a, b = library.match(read)
        if b == 0:
            assert a == 0
            assert read[:3] not in targets
        else:
            assert read[:3] == ('AAA' if a == 0 else 'CCC')


def test_multi_target_uniform_library_from_sequences():
    MultiTargetUniformLibrary.from_sequences(['AAA', 'CCC'])
    with pytest.raises(InvalidLibraryError):
        MultiTargetUniformLibrary.from_sequences(['AAA', 'C'])


def test_dynamic_library():
    seq = 'AAACCC'
    region_length = 3
    library = DynamicMultiTargetLibrary(10)
    assert len(library) == 0

    target_index, length = library.match(seq, offset=0, length=region_length, anchor_left=True)
    assert length == region_length
    assert target_index == 10 + 0

    for _ in range(2):
        target_index, length = library.match(seq, offset=3, length=region_length, anchor_left=True)
        assert length == region_length
        assert target_index == 10 + 1

    assert library.decode(10 + 0) == 'AAA'
    assert library.decode(10 + 1) == 'CCC'
    assert len(library) == 2
    assert library.sequences == ['AAA', 'CCC']


def test_library_builder_shared_targets():
    lb = LibraryBuilder.from_raw_libraries([
        RawLibrary(0, False, ['AAA', 'CCC']),
        RawLibrary(1, False, ['TTT', 'CCC'])
    ])
    assert len(lb.library_shared_targets) == 2
    assert lb.library_shared_targets == {
        0: {1: {1: 3}},
        1: {0: {3: 1}}
    }


@pytest.mark.parametrize('sequences', [
    ['AZT'],  # invalid sequence
    ['AAA', 'AAA']  # duplicate sequence
])
def test_raw_library(sequences):
    with pytest.raises(InvalidLibraryError):
        RawLibrary(0, False, sequences)


def test_raw_library_load_empty_sequence():
    with NamedTemporaryFile() as tmp:

        # Library with one empty sequence
        with open(tmp.name, 'w') as fh:
            fh.write("id\tsequence\n")
            fh.write("0\tAAA\n")
            fh.write("1\t\n")

        with pytest.raises(InvalidLibraryError):
            RawLibrary.load(tmp.name, 0, LibraryReverseCondition.NEVER)
