# LUCA
#
# Copyright (C) 2023, 2024 Genome Research Ltd.
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

from json import JSONDecodeError
import logging
import os
import platform
import sys

import click
import yaml
from click_option_group import OptionGroup
from pydantic import ValidationError
from yaml.parser import ParserError
from yaml.scanner import ScannerError
from yaml.error import Mark

from . import __version__
from .app_info import AppInfo
from .cli import common_cli_debug_options, common_cli_options, debug_opts
from .cli_utils import abort, existing_directory_path, existing_file_path
from .errors import InvalidExperiment, InvalidHTSError, InvalidLibraryError, InvalidReadFileExtension, MissingMetadataError, TargetNotFound
from .experiment import Experiment, Options, ReadGroupId, SequencingType
from .fs import get_read_counts_file_name
from .input_manifest import InputManifest
from .library_independent import count_reads
from .manifest import OutputFileManifest
from .matcher import ExperimentStats, MultiMatcher, get_multi_matcher
from .path_bundle import PathBundle
from .readers.main import parse_reads
from .readers.read_file_info import ReadFileInfo
from .readers.tsv import TsvError
from .stats import LibraryIndependentStats
from .utils import load_text_from_file, log_validation_error
from .writer import write_stats


HELP_LIBRARY = "Expanded library definition TSV file with optional headers (common format for single/dual/other)"
HELP_SAMPLE = "Sample name, interrogate header for others when not defined"
HELP_OUTPUT = "Output directory"
HELP_REFERENCE = "Required for CRAM"
HELP_CPUS = "CPUs to use (0 to detect)"
HELP_COUNT_MM_READS = "Whether to count the mismatching reads (as in library-independent counting)"
HELP_SORT_MM_READ_COUNTS = "Whether to sort mismatching read counts (in descending order), if calculated"
HELP_COMPRESS_MM_READ_COUNTS = "Whether to compress the counts of the mismatching reads, if calculated"
HELP_RC_REVERSE = "Whether to apply the reverse complement to reverse reads before matching (currently for single-end only)"
HELP_FULL_MATCH_INFO = "Whether to generate the full match QC table (vs. only for mismatching reads)"
HELP_OUT_MM_MATCH_INFO = "Whether to generate the match QC table for mismatching reads"
HELP_TEST_SWAPS = "Whether to test swap libraries, if provided"
HELP_OUT_MM_READS = "Whether to write the mismatching reads to an output BAM file"
HELP_LIMIT = "Maximum number of reads (or read pairs) to process"
HELP_INPUT_MANIFEST = "Manifest mapping identifiers used in the experiment configuration to input file paths"


class PlatformError(Exception):
    pass


def get_cpus(cpus: int) -> int:
    if cpus > 0:
        return cpus
    if platform.system() == 'Darwin':
        logging.error("Process affinity detection not available on macOS: please set the CPU count!")
        raise PlatformError
    return len(os.sched_getaffinity(0))  # type: ignore[attr]


def setup_fs(output_dir: str) -> None:
    # TODO: verify file vs. directory as output
    outfolder = os.path.dirname(os.path.abspath(output_dir))
    if not os.path.exists(outfolder):
        os.makedirs(outfolder, exist_ok=True)


def load_matcher(exp: Experiment, pb: PathBundle) -> MultiMatcher:
    try:
        return get_multi_matcher(exp, pb)
    except FileNotFoundError as ex:
        abort("File not found: '%s'!", ex.filename or ex.args[0])
    except InvalidLibraryError as ex:
        abort(ex.message)
    except TsvError as ex:
        abort("Invalid TSV file: '%s'!" % ex.fp)


def library_dependent_counting(exp: Experiment, opt: Options, pb: PathBundle, iter_reads, profile: bool = False) -> tuple[LibraryIndependentStats, ExperimentStats]:
    assert exp.read_group_count > 0

    try:
        if exp.has_default_read_group:

            # Match all reads to the same templates
            mm = load_matcher(exp, pb)
            stats, exp_stats = mm.match_seqs(opt, iter_reads, pb.output_dir, profile=profile)

        else:

            match exp.sequencing_type:

                case SequencingType.SINGLE_END:
                    abort("Single-end experiments currently only support the default read group!")

                case SequencingType.PAIRED_END:
                    assert ReadGroupId.DEFAULT not in exp.read_group_templates

                    if (
                        ReadGroupId.READ_1 not in exp.read_group_templates or
                        ReadGroupId.READ_2 not in exp.read_group_templates
                    ):
                        abort("Paired-end experiments currently only support processing both mates!")

                    mm = load_matcher(exp, pb)
                    stats, exp_stats = mm.match_seqs(opt, iter_reads, pb.output_dir)

                case _:
                    abort(f"Unsupported sequencing type '{exp.sequencing_type}'!")

    except TargetNotFound as ex:
        abort(ex.message)

    return stats, exp_stats


def library_independent_counting(sample_name: str | None, output_dir: str, iter_reads) -> LibraryIndependentStats:
    stats = count_reads(iter_reads, output_dir)
    manifest = OutputFileManifest.empty(sample_name)
    manifest.set_library_independent_file_paths(get_read_counts_file_name())
    manifest.write(output_dir)
    return stats


def load_input_manifest(library_dir: str | None, input_manifest: str | None) -> InputManifest | None:
    if not input_manifest and library_dir:

        # Look for a manifest file in the library directory
        fp = os.path.join(library_dir, 'manifest.json')
        if os.path.isfile(fp):
            input_manifest = fp

    if input_manifest:
        try:
            m = InputManifest.model_validate_json(
                load_text_from_file(input_manifest))
        except JSONDecodeError:
            abort("Failed to load input manifest: invalid JSON!")
        except ValidationError as ex:
            log_validation_error(ex)
            abort("Failed to load input manifest!")

        if not m.do_all_files_exist(root_dir=library_dir):
            abort("Missing input files in manifest!")

        return m

    return None


def validate_path_bundle(pb: PathBundle, exp: Experiment) -> None:
    if not pb.is_valid:
        abort("Invalid input manifest!")
    for library in exp.libraries:
        if library.from_file:
            fp = pb.get_library_file_path(library.id)
            if not os.path.isfile(fp):
                abort("Library '%s' not found at '%s'!" % (library.id, fp))
    for ft in exp.combination_filters:
        fp = pb.get_combination_filter_file_path(ft)
        if not os.path.isfile(fp):
            abort("Combination filter '%s' not found at '%s'!" % (ft, fp))


lib_dep_opts = OptionGroup("\nLibrary-dependent", help="Options specific to library-dependent counting")
sample_opts = OptionGroup("\nInput sample metadata", help="Options adding information to the input")
perf_opts = OptionGroup("\nPerformance", help="Options to tune the performance")


@click.command()
@click.argument('experiment', required=True, type=existing_file_path)
@click.argument('queries', required=True, type=existing_file_path)
@common_cli_options
@click.option('--rc-reverse', is_flag=True, default=None, help=HELP_RC_REVERSE)
@click.option('-m', '--input-manifest', type=existing_file_path, help=HELP_INPUT_MANIFEST)
@sample_opts.option('-s', '--sample', type=str, help=HELP_SAMPLE)
@sample_opts.option('-r', '--reference', type=existing_file_path, help=HELP_REFERENCE)
@lib_dep_opts.option('-l', '--library-dir', type=existing_directory_path, help=HELP_LIBRARY)
@lib_dep_opts.option('--count-mm-reads', is_flag=True, default=None, help=HELP_COUNT_MM_READS)
@lib_dep_opts.option('--sort-mm-read-counts', is_flag=True, default=None, help=HELP_SORT_MM_READ_COUNTS)
@lib_dep_opts.option('--compress-mm-read-counts', is_flag=True, default=None, help=HELP_COMPRESS_MM_READ_COUNTS)
@lib_dep_opts.option('--full-match-info', is_flag=True, default=None, help=HELP_FULL_MATCH_INFO)
@lib_dep_opts.option('--out-mm-match-info', is_flag=True, default=None, help=HELP_OUT_MM_MATCH_INFO)
@lib_dep_opts.option('--test-swaps', is_flag=True, default=None, help=HELP_TEST_SWAPS)
@lib_dep_opts.option('--out-mm-reads', is_flag=True, default=None, help=HELP_OUT_MM_READS)
@lib_dep_opts.option('--profile', is_flag=True)
@perf_opts.option('-c', '--cpus', type=click.IntRange(min=0), default=1, show_default=True, help=HELP_CPUS)
@debug_opts.option('--limit', type=int, default=None, help=HELP_LIMIT)
@common_cli_debug_options
def count(
    input_manifest: str | None,
    library_dir: str | None,
    sample: str | None,
    output: str,
    reference: str | None,
    experiment: str,
    queries: str,
    cpus: int,
    loglevel: str,
    limit: int | None = None,
    profile: bool = False,

    # Experiment default option overrides
    test_swaps: bool | None = None,
    full_match_info: bool | None = None,
    out_mm_match_info: bool | None = None,
    out_mm_reads: bool | None = None,
    count_mm_reads: bool | None = None,
    sort_mm_read_counts: bool | None = None,
    compress_mm_read_counts: bool | None = None,
    rc_reverse: bool | None = None
):
    """
    Count reads or regions, optionally matching libraries.

    \b
    EXPERIMENT: Experiment configuration (YAML or JSON)
    QUERIES: Query sequence file (BAM or CRAM)
    """

    # DEVELOPMENT NOTES
    # When adding a new option, ensure it is added as:
    # - click option (in the appropriate option group)
    # - main function argument
    # - Options class attribute
    # - Options.override call

    # Setup logger
    logging.basicConfig(
        level=logging._nameToLevel[loglevel.upper()],
        format="%(levelname)s: %(message)s")

    # Get usable CPU's
    try:
        usable_cpus: int = get_cpus(cpus)
    except PlatformError:
        abort("Failed to set usable CPU number!")

    # Load experiment configuration
    with open(experiment) as fh:
        try:
            exp = Experiment.model_validate(yaml.safe_load(fh))
        except ScannerError as ex:
            m = ex.problem_mark
            abort(f"YAML parsing error: {ex.problem} at line {m.line} column {m.column} in {m.name}!")  # type: ignore
        except ParserError as ex:
            logging.error(str(ex).replace('\n', ''))
            abort("Failed to load experiment configuration!")
        except ValidationError as ex:
            log_validation_error(ex)
            abort("Failed to load experiment configuration!")
        except InvalidExperiment:
            abort("Failed to load experiment configuration!")

    # Validate sample name
    try:
        read_file_info = ReadFileInfo.from_path(queries)
    except InvalidReadFileExtension as ex:
        abort(ex.message)

    # Setup output directory
    setup_fs(output)

    app_info = AppInfo(
        __version__,
        ' '.join([os.path.basename(sys.argv[0]), *sys.argv[1:]])
    )

    mm_reads_fn = 'mm.reads.bam'
    mm_reads_fp = os.path.join(output, mm_reads_fn)

    # TODO: allow to override options from the command line interface
    opt = exp.default_options
    opt.override(
        test_swaps=test_swaps,
        full_match_info=full_match_info,
        out_mm_match_info=out_mm_match_info,
        out_mm_reads=out_mm_reads,
        rc_reverse=rc_reverse,
        count_mm_reads=count_mm_reads,
        sort_mm_read_counts=sort_mm_read_counts,
        compress_mm_read_counts=compress_mm_read_counts)

    if opt.sort_mm_read_counts and not opt.count_mm_reads:
        logging.warning(
            "The 'sort_mm_read_counts' option is set, but it has no effect " +
            "because 'count_mm_reads' is not")

    is_library_dependent = len(exp.libraries) > 0
    pb = PathBundle(
        output_dir=output,
        library_dir=library_dir,
        input_manifest=load_input_manifest(library_dir, input_manifest))
    validate_path_bundle(pb, exp)

    if not library_dir:
        if (
            is_library_dependent and
            not (
                pb.input_manifest and
                pb.input_manifest.are_all_library_file_paths_absolute
            )
        ):
            abort("Library directory or input manifest with absolute paths required in library-dependent mode!")
        if exp.has_combination_filters:
            # TODO: consider having a separate root directory for filter files
            abort("Library directory required when expecting combination filter files!")

    if not is_library_dependent:
        for opt_name, opt_value in opt.library_dependent_options.items():
            if opt_value:
                logging.warning(
                    "The '%s' option is set, but it has no effect in library-independent mode." % opt_name)

    iter_reads = parse_reads(
        read_file_info,
        sample,
        usable_cpus,
        limit=limit,
        ofp=mm_reads_fp if is_library_dependent and opt.out_mm_reads else None,
        reference=reference,
        rc_reverse=opt.rc_reverse,
        is_paired_end=exp.sequencing_type == SequencingType.PAIRED_END)

    try:
        if is_library_dependent or exp.has_assigned_templates:
            stats, exp_stats = library_dependent_counting(exp, opt, pb, iter_reads, profile=profile)
            stats_fp = os.path.join(output, 'exp.stats.json')
            write_stats(app_info, exp_stats, stats_fp)  # type: ignore

        else:
            stats = library_independent_counting(sample, output, iter_reads)

    except MissingMetadataError as ex:
        logging.error(ex.message)
        sys.exit(1)

    except InvalidHTSError as ex:
        logging.error(ex.message)
        sys.exit(1)

    except TsvError as ex:
        abort("Invalid TSV file: '%s'!" % ex.fp)

    # Write stats to file
    assert stats
    logging.info("Writing stats report...")
    stats_fp = os.path.join(output, 'stats.json')
    write_stats(app_info, stats, stats_fp)

    # Write experiment to file (includes default values)
    with open(os.path.join(output, 'exp.json'), 'w') as fh:
        fh.write(exp.model_dump_json())
