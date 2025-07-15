#!/usr/bin/env python3

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

"""
Generate test sequencing data from a configuration bundle
"""

from argparse import ArgumentParser
from dataclasses import dataclass, field
import os
import random
import sys

from pydantic import ValidationError
from pysam import AlignedSegment, AlignmentFile
from pysam.libcalignedsegment import SAM_FLAGS
import yaml
from yaml.scanner import ScannerError

from luca.cli_utils import abort
from luca.count import load_input_manifest, validate_path_bundle
from luca.errors import InvalidExperiment
from luca.experiment import (
    Experiment,
    ReadGroupId,
    ReadRegion,
    SequencingType,
)
from luca.library import LibraryBuilder, RawLibrary
from luca.matcher import load_library
from luca.path_bundle import PathBundle
from luca.readers.tsv import TsvError
from luca.utils import log_validation_error


@dataclass(slots=True)
class DataGenOpt:
    read_number: int
    read_length: int


# FROM MultiMatcher.from_experiment
def library_builder_from_experiment(exp: Experiment, pb: PathBundle) -> LibraryBuilder:
    library_indices = exp.assigned_library_indices
    if library_indices:

        # Load libraries assigned to at least one read template
        return LibraryBuilder.from_raw_libraries([
            load_library(pb, exp.libraries[library_index], library_index)
            for library_index in library_indices.values()
        ])

    else:
        return LibraryBuilder.empty()


# FROM luca.utils
def load_experiment(fp: str) -> Experiment:
    with open(fp) as fh:
        try:
            return Experiment.model_validate(yaml.safe_load(fh))
        except ScannerError as ex:
            m = ex.problem_mark
            abort(f"YAML parsing error: {ex.problem} at line {m.line} column {m.column} in {m.name}!")  # type: ignore
        except ValidationError as ex:
            log_validation_error(ex)
            abort("Failed to load experiment configuration!")
        except InvalidExperiment:
            abort("Failed to load experiment configuration!")


def load_path_bundle(exp: Experiment, output: str, library_dir: str | None, input_manifest: str | None) -> PathBundle:
    pb = PathBundle(
        output_dir=output,
        library_dir=library_dir,
        input_manifest=load_input_manifest(library_dir, input_manifest))
    validate_path_bundle(pb, exp)
    return pb


@dataclass(slots=True)
class GenRegion:
    prefix: str
    length: int  # if zero, assume library-dependent
    libraries: list[RawLibrary]
    _range: list[int] = field(init=False)

    def __post_init__(self) -> None:
        self._range = list(range(len(self.libraries)))

    def generate(self):
        if self.libraries:
            i = random.choice(range(len(self.libraries)))
            return self.prefix + random.choice(self.libraries[i].sequences)
        else:
            return self.prefix + ''.join(random.choices('ACGT', k=self.length))


def get_region_max_length(region: ReadRegion, library_indices: dict[str, int], lb: LibraryBuilder) -> int:
    # TODO: consider region.max_offset as well?
    if region.libraries:
        return region.skip + max(
            lb.libraries[library_indices[library_id]].max_sequence_length
            for library_id in region.libraries
        )

    elif region.length is not None:
        return region.skip + region.length

    else:
        raise Exception("Invalid region!")


def get_regions(exp: Experiment, lb: LibraryBuilder, rid: ReadGroupId) -> list[GenRegion]:
    library_indices = exp.library_indices

    tpl = exp.get_read_group_templates(rid)[0]
    regions = [
        GenRegion(
            prefix=region.skip * 'N',
            length=region.length or 0,
            libraries=[
                lb.libraries[library_indices[library_id]]
                for library_id in region.libraries
            ])
        for region in tpl.regions
    ]

    # Handle unbound last region
    # TODO: consider whether this is worth it
    if not regions[-1].length and not regions[-1].libraries:
        max_length: int = 0
        for region in tpl.regions:
            max_length += get_region_max_length(region, library_indices, lb)
        regions[-1].length = opt.read_length - max_length

    return regions


def get_read(flag: SAM_FLAGS | None = None) -> AlignedSegment:
    read: AlignedSegment = AlignedSegment()
    read.reference_id = -1
    read.flag = flag.value if flag is not None else 0
    read.query_name = 'r'
    return read


def get_read_pair() -> tuple[AlignedSegment, AlignedSegment]:
    return get_read(SAM_FLAGS.FREAD1), get_read(SAM_FLAGS.FREAD2)


def get_read_sequence(regions: list[GenRegion]) -> str:
    return ''.join([
        region.generate()
        for region in regions
    ])


def generate_single_end(exp: Experiment, opt: DataGenOpt, lb: LibraryBuilder, sam: AlignmentFile):
    regions = get_regions(exp, lb, ReadGroupId.DEFAULT)

    read = get_read()
    # read.template_length = opt.read_length
    for i in range(opt.read_number):
        read.query_sequence = get_read_sequence(regions)
        read.template_length = len(read.query_sequence)
        sam.write(read)


def generate_paired_end(exp: Experiment, opt: DataGenOpt, lb: LibraryBuilder, sam: AlignmentFile):
    regions_1 = get_regions(exp, lb, ReadGroupId.READ_1)
    regions_2 = get_regions(exp, lb, ReadGroupId.READ_2)

    read_1, read_2 = get_read_pair()
    for i in range(opt.read_number):
        read_1.query_sequence = get_read_sequence(regions_1)
        read_1.template_length = len(read_1.query_sequence)
        read_2.query_sequence = get_read_sequence(regions_2)
        read_2.template_length = len(read_2.query_sequence)
        sam.write(read_1)
        sam.write(read_2)


def main(opt: DataGenOpt, exp: Experiment, out_fp: str, pb: PathBundle, seed: int | None = None):
    if seed is not None:
        random.seed(seed)

    try:
        lb = library_builder_from_experiment(exp, pb)
    except TsvError as ex:
        abort("Invalid TSV file: '%s'!" % ex.fp)

    if not lb.libraries:
        raise NotImplementedError("Library-independent!")

    for templates in exp.read_group_templates.values():
        if len(templates) > 1:
            raise NotImplementedError("Multi-template!")

    # Create output BAM file
    header = {
        'HD': {'VN': '1.0'}
    }
    try:
        sam = AlignmentFile(out_fp, header=header, mode='wb')
    except ValueError as ex:
        sys.exit("HTS file error: " + ex.args[0])

    try:

        # Generate & write reads to BAM file
        match exp.sequencing_type:
            case SequencingType.SINGLE_END:
                generate_single_end(exp, opt, lb, sam)
            case SequencingType.PAIRED_END:
                generate_paired_end(exp, opt, lb, sam)

    finally:
        sam.close()


if __name__ == '__main__':
    p = ArgumentParser(description="Generate synthetic CRISPR data from a LUCA configuration")
    p.add_argument('-o', '--output', required=True, help="Output BAM file path")
    p.add_argument('-n', '--number', required=True, type=int, help="Number of reads or read pairs")
    p.add_argument('-s', '--seed', type=int, help="Random seed for reproducible generation")
    p.add_argument('-l', '--library-dir', help="LUCA library directory path")
    p.add_argument('-m', '--manifest', help="LUCA input manifest file path")
    p.add_argument('experiment', help="LUCA experiment configuration file path")
    args = p.parse_args()

    opt = DataGenOpt(read_number=args.number, read_length=60)
    if opt.read_number <= 0:
        sys.exit("Invalid read number!")

    if not os.path.isfile(args.experiment):
        sys.exit(f"File not found: '{args.experiment}'!")

    exp = load_experiment(args.experiment)
    pb = load_path_bundle(exp, '/tmp', args.library_dir, args.manifest)

    main(opt, exp, args.output, pb, seed=args.seed)
