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

from collections import Counter
from contextlib import nullcontext
import os
from tempfile import TemporaryDirectory, NamedTemporaryFile

import numpy as np
import pytest

from luca.errors import InvalidExperiment
from luca.experiment import Anchor, CombinationInfo, CombinationRegion, Experiment, LibraryInfo, LibraryReverseCondition, Options, \
    ReadGroupId, ReadRegion, \
    ReadTemplate, SequencingType
from luca.library import LibraryBuilder, MonoTargetLibrary, RawLibrary
from luca.matcher import LibraryIndependentRegion, Matcher, PairedEndMultiMatcher, Region, SingleEndMultiMatcher, Template, READ_GROUP_DEFAULT, DataBundle, FixedLengthRegion, LibraryIndependentCounter, match_read_region
from luca.path_bundle import PathBundle
from luca.utils import clamp_non_negative

from .utils import mock_reads


READS = [
    'AAACCC',
    'AAACCC',
    'CCCGG',
    'CCCGG',
    'AAACCC',
    'G',
    'AAAG',
    'TTTA'
]

EXP_OUT_FILES = {
    'lib.0.counts.tsv',
    'lib.0.stats.json',
    'matches.0.tsv',
    'mm.matches.0.tsv',
    'mm.read.counts.tsv',
    'manifest.json'
}


def yield_counts(fp):
    assert os.path.isfile(fp)
    with open(fp) as fh:
        while line := fh.readline().rstrip():
            index, target, count = line.split('\t')
            yield int(index), target, int(count)


def yield_combination_counts(fp):
    assert os.path.isfile(fp)
    with open(fp) as fh:
        while line := fh.readline().rstrip():
            t = line.split('\t')
            yield *t[:-1], int(t[-1])


@pytest.mark.parametrize('offset', [0, 3])
def test_fixed_length(offset):
    seq = 'AAACCC'
    n = 3
    region = FixedLengthRegion(length=n)
    lic = LibraryIndependentCounter.empty(0, 0)
    bundle = DataBundle(np.ndarray(0), {}, lic, [], Counter())
    region_offset, match_length, is_swap = region.match(
        {}, bundle, np.asarray([-1] * 3), seq, 0, offset, anchor_left=True)

    assert region_offset == offset
    assert match_length == n
    assert not is_swap
    assert len(lic.counter) == 1
    assert lic.counter[lic.library._codec.encode(seq[offset:offset + n])] == 1


def test_single_end_matcher():
    lib = MonoTargetLibrary(0, 'AAA')
    tpl = Template(True, [Region(libraries=[0])])
    m = Matcher({0: lib}, [tpl])
    mm_reads = []

    mm = SingleEndMultiMatcher(
        LibraryBuilder.from_raw_libraries([RawLibrary(0, LibraryReverseCondition.REVERSE_GROUP, ['AAA'])]),
        {READ_GROUP_DEFAULT: m})

    count_exp = sum(1 for r in READS if r.startswith(lib.target))
    mm_count_exp = len(READS) - count_exp
    with TemporaryDirectory() as out_dir:
        mm.match_seqs(Options(), mock_reads(READS, mm_reads), out_dir)

        # Verify the expected files were generated
        for fn in os.listdir(out_dir):
            assert fn in EXP_OUT_FILES

        # Check library counts file
        fn = 'lib.0.counts.tsv'
        fp = os.path.join(out_dir, fn)
        assert os.path.isfile(fp)
        for index, target, count in yield_counts(fp):
            assert int(index) >= 0
            assert target == lib.target
            assert int(count) == count_exp

    assert len(mm_reads) == mm_count_exp


def test_single_end_matcher_multi_target():
    trl = RawLibrary(0, LibraryReverseCondition.REVERSE_GROUP, ['ACT'])
    rl = RawLibrary(1, LibraryReverseCondition.REVERSE_GROUP, ['AAA', 'CCC', 'TTT', 'GGG'])
    raw_libs = [trl, rl]
    lib = rl.build(target_index_offset=1)
    tpl = Template(anchor_left=True, regions=[Region(libraries=[1])])
    m = Matcher({1: lib}, [tpl])
    mm = SingleEndMultiMatcher(
        LibraryBuilder.from_raw_libraries(raw_libs),
        {READ_GROUP_DEFAULT: m})
    mm_reads = []
    exp_counts = [0, 4, 2, 1, 0]
    with TemporaryDirectory() as out_dir:
        mm.match_seqs(Options(), mock_reads(READS, mm_reads), out_dir)
        fn = 'lib.1.counts.tsv'
        fp = os.path.join(out_dir, fn)
        assert os.path.isfile(fp)
        for index, target, count in yield_counts(fp):
            assert count == exp_counts[index]
            print(f"{index}\t{target}\t{count}")

    assert len(mm_reads) == 1
    assert mm_reads[0] == 'G'


def test_paired_end_matcher():
    exp = Experiment(
        sequencing_type=SequencingType.PAIRED_END,
        libraries=[
            LibraryInfo(id='a', values=['AAA']),
            LibraryInfo(id='b', values=['CCC'])
        ],
        read_templates=[
            ReadTemplate(id='forward', regions=[ReadRegion(id='x', libraries=['a'])]),
            ReadTemplate(id='reverse', regions=[ReadRegion(id='x', libraries=['b'])])
        ],
        read_group_templates={
            ReadGroupId.READ_1: ['forward'],
            ReadGroupId.READ_2: ['reverse']
        })
    PairedEndMultiMatcher.from_experiment(exp, PathBundle('.', '.', None))


@pytest.mark.parametrize('filtered', [True, False])
def test_paired_end_fixed_length_matcher(filtered):
    combination_info = CombinationInfo(id='default', regions=[
        CombinationRegion(id='x', read_group=ReadGroupId.READ_1),
        CombinationRegion(id='y', read_group=ReadGroupId.READ_2)
    ])

    if filtered:
        combination_info.regions[1].filter = True
        combination_info.filters = ['ft']
        ctx = pytest.raises(InvalidExperiment)
    else:
        ctx = nullcontext()

    with ctx:
        exp = Experiment(
            sequencing_type=SequencingType.PAIRED_END,
            libraries=[
                LibraryInfo(id='a', values=['AAA', 'TTT'])
            ],
            read_templates=[
                ReadTemplate(id='forward', regions=[ReadRegion(id='x', libraries=['a'])]),
                ReadTemplate(id='reverse', regions=[ReadRegion(id='y', length=3)], anchor=Anchor.LEFT)
            ],
            read_group_templates={
                ReadGroupId.READ_1: ['forward'],
                ReadGroupId.READ_2: ['reverse']
            },
            combinations=[combination_info])

    if filtered:
        return

    reads = [('AAAC', 'CCCA'), ('AAAC', 'CCCA')]
    mm_reads = []
    with TemporaryDirectory() as out_dir:
        with NamedTemporaryFile(prefix='ft', dir=out_dir) as ft:
            assert os.path.isfile(ft.name)
            mm = PairedEndMultiMatcher.from_experiment(exp, PathBundle('.', '.', None))
            print(mm)
            mm.match_seqs(Options(), mock_reads(reads, mm_reads), out_dir)
            fn = 'combination.0.counts.tsv'
            fp = os.path.join(out_dir, fn)
            for t in yield_combination_counts(fp):
                assert t == ('AAA', 'CCC', 2)


def test_negative_skip():
    tpl = Template(anchor_left=True, regions=[
        FixedLengthRegion(length=5),
        FixedLengthRegion(length=5, skip=-5),
        FixedLengthRegion(length=5)
    ])
    lic = LibraryIndependentCounter.empty(0, 0)
    bundle = DataBundle(np.ndarray(0), {}, lic, None, None)  # type: ignore
    matches = np.asarray([-1] * 3 * len(tpl.regions))
    tpl.match({}, bundle, matches, 'AAACCTTTTT')
    print(lic)
    assert len(lic.counter) == 2
    assert lic.counter[lic.library._codec.encode('AAACC')] == 2
    assert lic.counter[lic.library._codec.encode('TTTTT')] == 1
    assert matches[2] == 1


def test_match_read_region():
    library = MonoTargetLibrary(0, 'AAA')
    for skip in [-5, 0]:
        target_index, match_length, region_offset = match_read_region(
            'AAACCC', library, Region(skip=skip, libraries=[0]), 0, True)
        assert target_index == 0
        assert match_length == 3
        assert clamp_non_negative(region_offset) == 0


def test_match_read_region_fixed():
    s = 'AAAAAGGG'
    r = FixedLengthRegion(skip=5, length=3)
    matches = np.asarray([-1] * 3 * 1)
    bundle = DataBundle(None, None, LibraryIndependentCounter.empty(0, 0), None, None)  # type: ignore
    match_region_offset, match_length, is_swap = r.match({}, bundle, matches, s, 0, 0, True)
    assert match_region_offset == 5
    assert match_length == 3
    assert is_swap is False
    assert bundle.library_indep_counter.counter[0] == 1


def get_regions(skip):
    return [
        Region(skip=skip, libraries=[0]),
        LibraryIndependentRegion(skip=skip),
        FixedLengthRegion(skip=skip, length=3)
    ]


@pytest.mark.parametrize('region', get_regions(2))
def test_match_region_types(region):
    skip = 2

    rest = LibraryIndependentRegion(skip=0)

    libs = {
        0: MonoTargetLibrary(0, 'AAA')
    }
    s = 'T' * skip + 'AAA'

    bundle = DataBundle(
        np.zeros(1, dtype=np.uint64), None, LibraryIndependentCounter.empty(0, 0), None, None)  # type: ignore

    if isinstance(region, LibraryIndependentRegion):
        with pytest.raises(AssertionError):
            tpl = Template(True, [region, rest])
        tpl = Template(True, [region])
    else:
        s += 'GGG'
        tpl = Template(True, [region, rest])

    matches = np.asarray([-1] * 3 * len(tpl.regions))

    match_length, is_swap = tpl.match(libs, bundle, matches, s)
    assert match_length == 3
    assert is_swap is False
