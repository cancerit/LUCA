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

from luca.experiment import CombinationInfo, CombinationRegion, Experiment, LibraryInfo, ReadGroupId, ReadRegion, ReadTemplate, SequencingType
from luca.library import LibraryBuilder, RawLibrary
from luca.matcher import CombinationRules, NCCombinationRules, build_combination_rules, \
    PartialNCCombinationRules, READ_GROUP_1, READ_GROUP_2

exp = Experiment(
    sequencing_type=SequencingType.PAIRED_END,
    libraries=[
        LibraryInfo(id='a', values=['AAA']),
        LibraryInfo(id='b', values=['CCC']),
        LibraryInfo(id='c', values=['TTT']),
    ],
    read_templates=[
        ReadTemplate(id='forward', regions=[
            ReadRegion(id='y', libraries=['c']),
            ReadRegion(id='x', libraries=['a'])
        ]),
        ReadTemplate(id='reverse', regions=[ReadRegion(id='x', libraries=['b'])])
    ],
    read_group_templates={
        ReadGroupId.READ_1: ['forward'],
        ReadGroupId.READ_2: ['reverse']
    },
    combinations=[
        CombinationInfo(
            id='x',
            regions=[
                CombinationRegion(id='x', read_group=ReadGroupId.READ_1),
                CombinationRegion(id='x', read_group=ReadGroupId.READ_2)
            ])
    ])


def test_combination_rules_get_region_offsets():
    crs, _ = build_combination_rules(
        exp, '/missing', LibraryBuilder.from_raw_libraries([
            RawLibrary(0, False, ['AAA']),
            RawLibrary(1, False, ['CCC']),
            RawLibrary(2, False, ['TTT'])
        ]))
    assert len(crs) == 1
    cr = crs[0]
    assert cr.get_region_offsets(READ_GROUP_1, 0) == [1 * 3 + 1]
    assert cr.get_region_offsets(READ_GROUP_2, 0) == [0 * 3 + 1]


def test_combinatorial_combination_rules():
    cr = CombinationRules({})
    assert len(cr.get_counter()) == 0
    assert cr.should_count_combination((0, 1))


def test_non_combinatorial_combination_rules():
    combinations = {(0, 1), (0, 2)}
    cr = NCCombinationRules({}, combinations)
    c = cr.get_counter()
    assert len(c) == 2
    for combination in combinations:
        assert combination in c
        assert cr.should_count_combination(combination)
    assert not cr.should_count_combination((0, 3))
    assert not cr.should_count_combination((3, 3))


def test_partial_non_combinatorial_combination_rules():
    combinations = {(0, 1), (0, 2)}
    cr = PartialNCCombinationRules({}, combinations, [0, 3])
    assert len(cr.get_counter()) == 0
    assert cr.should_count_combination((0, 100, 100, 1))
    assert not cr.should_count_combination((0, 100, 100, 3))
    assert not cr.should_count_combination((3, 100, 100, 3))
