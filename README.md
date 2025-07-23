# Locus-based Universal CRISPR Aligner (LUCA)

Single- and multi-guide exact matching (with gaps) for single- and paired-end sequencing data, supporting both combinatorial and non-combinatorial setups.

This tool is primarily concerned with supporting the widest possible range of experimental designs.

Goals:

- support for BAM and CRAM inputs
- memory requirements independent of input size (read streaming)

Non-goals:

- fuzzy matching: it may be attempted downstream to rescue the mismatching reads, which are optionally collected in a separate BAM file for this purpose
- count normalisation
- support for redundant or overlapping libraries
- support for alternate formats for the input files (*e.g.*, libraries)

Please refer to the [examples](examples/) directory for some usage scenarios.

## Table of contents

- [Locus-based Universal CRISPR Aligner (LUCA)](#locus-based-universal-crispr-aligner-luca)
  - [Table of contents](#table-of-contents)
  - [Usage](#usage)
    - [Quantification](#quantification)
    - [Merge](#merge)
  - [Installation](#installation)
    - [Python virtual environment](#python-virtual-environment)
    - [Docker image](#docker-image)
  - [File formats](#file-formats)
    - [Per experiment](#per-experiment)
      - [Experiment](#experiment)
        - [Library](#library)
        - [Read group](#read-group)
        - [Read template](#read-template)
        - [Region](#region)
        - [Read group templates](#read-group-templates)
        - [Combination](#combination)
        - [Combination region](#combination-region)
        - [Options](#options)
      - [Experiment statistics](#experiment-statistics)
      - [Library-independent statistics](#library-independent-statistics)
      - [Match QC table](#match-qc-table)
      - [Additional QC tables](#additional-qc-tables)
      - [Read counts](#read-counts)
      - [Combination counts](#combination-counts)
      - [Input manifest](#input-manifest)
      - [Manifest](#manifest)
    - [Per library](#per-library)
      - [Library](#library-1)
      - [Counts](#counts)
      - [Library-dependent statistics](#library-dependent-statistics)

## Usage

### Quantification

The tool requires at least a BAM or CRAM file containing the reads to process, and a description of the experimental design as configuration. The default options in the configuration may be overridden via the corresponding command line options.

Input files:

- reads (BAM/CRAM)
- [experiment metadata](#experiment) (YAML or JSON)
- [libraries](#library) (TSV)

Output files:

- [experiment metadata](#experiment) (JSON)
- [counts](#counts) (TSV)
- [experiment stats](#experiment-statistics) (JSON)
- [library-independent stats](#library-independent-statistics) (JSON)
- [library-dependent stats](#library-dependent-statistics) (JSON)
- [match QC table](#match-qc-table) (optional; TSV)
- [all/mismatching read counts](#read-counts) (optional; TSV)
- [mismatch QC table](#additional-qc-tables) (optional; TSV)
- [swap QC table](#additional-qc-tables) (optional; TSV)
- fully/partially unmatched reads (optional; BAM)
- [combination counts](#combination-counts) (optional; TSV)

To run it, *e.g.*:

```sh
luca count -o out_dir -l lib_dir exp.yaml reads.bam
```

The output directory is expected to exist and to be writeable.

The library directory is expected to exist and to contain the library files, each named as the identifiers listed in the `libraries` section of the experiment configuration.

```
Usage: luca count [OPTIONS] EXPERIMENT QUERIES

  Count reads or regions, optionally matching libraries.

  EXPERIMENT: Experiment configuration (YAML or JSON)
  QUERIES: Query sequence file (BAM or CRAM)

Options:
  -o, --output DIRECTORY          Output directory  [required]
  --rc-reverse                    Whether to apply the reverse complement to
                                  reverse reads before matching (currently for
                                  single-end only)
  -m, --input-manifest FILE       Manifest mapping identifiers used in the
                                  experiment configuration to input file paths

Input sample metadata:         Options adding information to the input
    -s, --sample TEXT             Sample name, interrogate header for others
                                  when not defined
    -r, --reference FILE          Required for CRAM

Library-dependent:             Options specific to library-dependent
                                  counting
    -l, --library-dir DIRECTORY   Expanded library definition TSV file with
                                  optional headers (common format for
                                  single/dual/other)
    --count-mm-reads              Whether to count the mismatching reads (as
                                  in library-independent counting)
    --sort-mm-read-counts         Whether to sort mismatching read counts (in
                                  descending order), if calculated
    --compress-mm-read-counts     Whether to compress the counts of the
                                  mismatching reads, if calculated
    --full-match-info             Whether to generate the full match QC table
                                  (vs. only for mismatching reads)
    --out-mm-match-info           Whether to generate the match QC table for
                                  mismatching reads
    --test-swaps                  Whether to test swap libraries, if provided
    --out-mm-reads                Whether to write the mismatching reads to an
                                  output BAM file
    --profile

Performance:                   Options to tune the performance
    -c, --cpus INTEGER RANGE      CPUs to use (0 to detect)  [default: 1;
                                  x>=0]

Debug:                         Options specific to troubleshooting, testing
                                  and debugging
    --limit INTEGER               Maximum number of reads (or read pairs) to
                                  process
    --loglevel [WARNING|INFO|DEBUG]
                                  Set logging verbosity  [default: INFO]
  --help                          Show this message and exit.
```

The following command line options can override the defaults set in the [experiment configuration](#experiment) (see the [options](#options) reference documentation).

|CLI option|Experiment default option|
|-|-|
|`full-match-info`|`full_match_info`|
|`count-mm-reads`|`count_mm_reads`|
|`sort-mm-read-counts`|`sort_mm_read_counts`|
|`compress-mm-read-counts`|`compress_mm_read_counts`|
|`out-mm-match-info`|`out_mm_match_info`|
|`test-swaps`|`test_swaps`|
|`out-mm-reads`|`out_mm_reads`|
|`rc-reverse`|`rc_reverse`|

### Merge

To merge the raw counts and generate aggregate statistics, list the [manifest files](#manifest) generated by the quantification in their respective output directories. The other output files are expected to be in the same directories as the manifests.

To merge the results for two different lanes of the same sample, *e.g.*:

```sh
luca merge -o out_dir lane_1/manifest.json lane_2/manifest.json
```

```
Usage: luca merge [OPTIONS] MANIFEST...

  Merge counts across different runs and regenerate the statistics.

  MANIFEST: JSON manifest generated by the `count` subcommand

Options:
  -o, --output DIRECTORY          Output directory  [required]
  --validate-targets              Validate the target sequences when merging
  --skip-mm-read-counts           Do not merge mismatching read library-
                                  independent counts

Debug:                         Options specific to troubleshooting, testing
                                  and debugging
    --loglevel [WARNING|INFO|DEBUG]
                                  Set logging verbosity  [default: INFO]
  --help                          Show this message and exit.
```

The following changes may be observed when comparing the original with the aggregate counts files:

- dynamic targets may be assigned different identifiers, and the order in which the manifests are provided affects the result (such identifiers are only required to decode the optional [match QC table](#match-qc-table));
- combinations may appear in a different order, also depending on the order of the manifests.

## Installation

The instructions that follow apply to Linux and macOS.

### Python virtual environment

Please take care to read errors during the dependency installation step carefully. HTSlib (pysam) has system dependencies and will highlight the packages that need to be installed.

Requirements:

- Python 3.11 or above

To install in a virtual environment:

```sh
# Initialise the virtual environment
python -m venv .venv

# Activate the virtual environment
source .venv/bin/activate

# Install the luca package
pip install .
```

### Docker image

To build the Docker container:

```sh
docker build -t luca .
```

## File formats

For all character-delimited input files:

- the separator is tab (TSV)
- the header line, containing only field names, is mandatory
- comments are not allowed
- the order of the fields is free
- additional fields are allowed and ignored

For a [library](#library-1), *e.g.*:

```tsv
id	sequence	extra
ID1	ACGTTATC	<IGNORED>
```

For all character-delimited output files (*e.g.*, library counts):

- the separator is tab (TSV)
- no header is present (to facilitate concatenation when merging results)

All statistics are reported in JSON format.

### Per experiment

#### Experiment

**Format**: YAML, JSON

Description of the experiment. If no libraries are listed, **library-independent counting** will be performed.

|Property|Type|Description|
|-|-|-|
|`sequencing_type`|`single_end`\|`paired_end`|Sequencing type.|
|`libraries`|[library](#library) array|Libraries.|
|`read_groups`|[read group](#read-group) map|Read group properties.|
|`read_templates`|[read template](#read-template) array|Read templates.|
|`read_group_templates`|[read group to template](#read-group-templates) map|Association between read groups and templates.|
|`combinations`|[combination](#combination) array|(Optional) Combinations to evaluate.|
|`default_options`|[options](#options) object|(Optional) Default options.|

The sequencing type determines which read groups may be used:

- `single_end`: `default`
- `paired_end`: `read_1`, `read_2`

The sequences of `read_2` reads are expected to be reverse complements, and the regions of the corresponding read templates will be tested in reverse order against the reverse complement sequences of the targets.

Each read group has a default orientation that may be overridden via `read_groups`:

|Read group|Default orientation|
|-|-|
|`default`|forward|
|`read_1`|forward|
|`read_2`|reverse|

Paired-end data is expected to be collated (*e.g.*, via `samtool collate`).

##### Library

Set of DNA sequences with unique identifiers.

If the list of targets is not provided, the target sequences are expected to be found in a TSV file (on a path based on the library identifier).

|Property|Type|Default|Description|
|-|-|-|-|
|`id`|string||Library ID.|
|`reverse_on`|`reverse_group`\|`always`\|`never`|`reverse_group`|Condition upon which to apply the reverse complement to the targets before matching.|
|`values`|DNA string array|`null`|(Optional) List of targets.|

When setting the reverse condition to `reverse_group`, the reverse complement is only applied to the targets when the library is being queried on a read group marked as reverse (by default, `read_2`).

##### Read group

|Property|Type|Default|Description|
|-|-|-|-|
|`is_reverse`|Boolean||Whether the reads in the group are expected to be in the reverse orientation.|

##### Read template

Sequence of regions with a DNA sequence that are expected to match a target sequence.

|Property|Type|Default|Description|
|-|-|-|-|
|`id`|string||Read template ID.|
|`anchor`|`auto`\|`left`\|`right`|`auto`|Direction in which the regions are expected to be found.|
|`regions`|[region](#region) array||Regions.|

When setting `anchor` to `auto`, the regions are matched from the left in reads belonging to a forward group, and from the right in reads belonging to a reverse group. By default, the read sequence is considered as stored in the BAM/CRAM file; for mapped data, the `rc_reverse` option can be set to ensure all sequences are converted to forward before matching.

##### Region

A region is defined as a portion of a DNA sequence that is expected to start with a known sequence (usually out of a set of candidate targets, the library). Several bases may be ignored at the start for the purposes of the match. The length of a region is not defined, and so is therefore the start position of any region after the first. *E.g.*, the second region in a read template will start after the last position of the first region that matched the target. If no target matches any given regions, none of the regions downstream of it will be tested, as their start positions would not be definable.

Alternate offsets are always tested in ascending order, from zero to the maximum.

If none of the expected targets (from the collections listed in `libraries`) matches a region, alternate targets may be tested (`swap_libraries`); if any of the latter is a match, the read is still considered a mismatch, but is also counted as a potential **swap event** (known sequence found in a read or region it was not expected to appear in).

|Property|Type|Default|Description|
|-|-|-|-|
|`id`|||Region ID.|
|`libraries`|Library ID array||Identifiers of the libraries to be tested.|
|`swap_libraries`|Library ID array|empty array|Identifiers of the libraries to be tested as a fallback to detect swaps.|
|`skip`|integer|0|Number of bases to skip before matching.|
|`max_offset`|integer|0|Maximum number of bases that may be skipped before matching (after the amount expressed in the `skip` property).|
|`length`|integer|`null`|Exact length, mutually exclusive with maximum offset.|

If skip is 2 and the maximum offset is 3, *e.g.*:

```
      AACGATCAGTAC
2     **
2 + 0   ->
2 + 1    ->
2 + 2     ->
2 + 3      ->
```

The target sequences from the expected libraries will be compared with the start of the following subsequences within the read, until one matches:

1. `CGATCAGTAC`
2. `GATCAGTAC`
3. `ATCAGTAC`
4. `TCAGTAC`

If none of the above is a match, the read is flagged as mismatching.

When the `length` property is set, but no libraries are provided, any sequence may be captured at the region start offset. Capturing *e.g.* six bases after skipping three:

```
      AACGATCAGTAC
3     ***
3 + 0    ----->
```

If neither libraries nor a length are provided, the region will capture all remaining bases; only one such region is allowed in any given template, and it must be the last. This may be used to model library-independent quantification of paired-end data (see [example configuration](examples/dual_guide_library_independent.yaml)).

##### Read group templates

Mapping of read group identifiers to arrays of read template identifiers.

The valid read group identifiers depend on the sequencing type:

- only `default` for single-end
- both `read_1` and `read_2` for paired-end

##### Combination

Regions whose combined matches should be counted. The valid combinations for a subset of such regions (marked with the `filter` flag) may be restricted via inclusion lists, whose identifiers are to be listed in the `filters` property.

In paired-end dual-guide experiments, these regions would typically be the guides from read 1 and 2.

|Property|Type|Default|Description|
|-|-|-|-|
|`id`|string||Combination ID.|
|`regions`|[combination region](#combination-region) array||Regions to evaluate.|
|`filters`|string array|`null`|(Optional) Identifiers of the combination inclusion list files.|

##### Combination region

|Property|Type|Default|Description|
|-|-|-|-|
|`id`|string||Region ID.|
|`read_group`|`default`\|`read_1`\|`read_2`|`default`|Read group ID.|
|`filter`|Boolean|`false`|Whether it is expected to be included in the filter (if any).|

Combinations may be filtered by passing inclusion lists via files (`filters` property). The expected field names should match the pattern `[<READ GROUP ID>.]<REGION ID>` (unless the read group is `default`, in which case it should match just the region ID), *e.g.*:

- `tag` (the read group is assumed to be `default`)
- `read_1.guide`
- `read_2.guide`

##### Options

|Property|Type|Default|Description|
|-|-|-|-|
|`test_swaps`|Boolean|`false`|Whether to test swap libraries, if provided.|
|`full_match_info`|Boolean|`false`|Whether to generate the full match QC table (*vs.* only for mismatching reads).|
|`out_mm_match_info`|Boolean|`false`|Whether to generate the match QC table for mismatching reads.|
|`out_mm_reads`|Boolean|`false`|Whether to write the mismatching reads to an output BAM file.|
|`count_mm_reads`|Boolean|`false`|Whether to count the mismatching reads (as in library-independent counting).|
|`sort_mm_read_counts`|Boolean|`false`|Whether to sort mismatching read counts (in descending order), if calculated.|
|`compress_mm_read_counts`|Boolean|`false`|Whether to compress the counts of the mismatching reads, if calculated.|
|`rc_reverse`|Boolean|`false`|Whether to apply the reverse complement to reverse reads before matching.|

#### Experiment statistics

**Format**: JSON

**File name**: `exp.stats.json`

|Property|Type|Description|
|-|-|-|
|`read_match`|integer|Total matching reads.|
|`read_swap`|integer|Total mismatching reads with at least one swap event.|
|`read_mismatch`|integer|Total reads with at least one mismatching region.|
|`pair_one_end_match`|integer|(Paired-end only) Total read pairs in which only one read is matching.|

#### Library-independent statistics

**Format**: JSON

**File name**: `stats.json`

|Field|Format|Description|
|-|-|-|
|`sample_name`|string|Name of the sample.|
|`input_reads`|integer|Total input reads.|
|`total_reads`|integer|Total reads passed on to counting.|
|`discarded_reads`|integer|Total reads discarded before counting.|
|`vendor_failed_reads`|integer|Total reads with the `QCFAIL` flag.|
|`length_excluded_reads`|integer|Total reads discarded because shorter than a user-defined threshold.|
|`ambiguous_nt_reads`|integer|Total reads with ambiguous nucleotides.|
|`masked_reads`|integer|Total soft-masked reads.|
|`zero_length_reads`|integer|Total zero-length reads.|

#### Match QC table

**Format**: TSV

**File name**: `matches.<INDEX>.tsv`

One file per read group, one row per read, in input order.

|Index|Type|Description|
|-|-|-|
|1|integer|Read index.|
|2|integer|Matching library index.|
|3|integer|Matching target index.|
|4|integer|Match status.|

Fields two to four appear as many times as the regions in the template.

The match status may represent one of the following:

- expected match (`1`)
- unexpected match (prospective swap, `2`)
- mismatch (`-1`)

If any region is a mismatch, those that would follow are not evaluated and are considered mismatches, because their starting positions are relative to the immediately preceding regions, which do not have an intrinsic length (the matching target, if any, would determine it).

#### Additional QC tables

Rows from the match QC table representing mismatch and/or swap events are duplicated into separate files for ease of retrieval by downstream processes (*e.g.*, alignment-based rescue of mismatching reads).

The swap category is a subset of the mismatch category, and therefore its reads also appear in the mismatch file.

|Category|File name|
|-|-|
|Mismatch|`mm.matches.<INDEX>.tsv`|
|Swap|`swap.matches.<INDEX>.tsv`|

#### Read counts

**Format**: TSV

**File name (library-independent)**: `read.counts.tsv`

**File name (library-dependent)**: `mm.read.counts.tsv`

Count of the unique read sequences; all of them in library-independent counting, or just the mismatching ones in library-dependent counting. The latter may optionally be sorted by count to facilitate the identification of the most common mismatching reads.

If the `rc_reverse` option is set, the reverse reads won't be listed, and their counts will be added to those of their reverse complements.

|Index|Type|Description|
|-|-|-|
|1|string|Unique DNA sequence of the read.|
|2|integer|Total number of occurrences.|

#### Combination counts

**Format**: TSV

**File name**: `combination.<INDEX>.counts.tsv`

One file per entry in the [combinations](#combination) section of the experiment configuration.

|Index|Type|Description|
|-|-|-|
|1 to n|string|Unique DNA sequence.|
|n + 1|integer|Total number of matches.|

#### Input manifest

**Format**: JSON

**File name**: `manifest.json` (default)

Optional manifest to map identifiers to input files. File paths may be relative or absolute. Relative paths will use the library directory as their root, if provided. If a manifest is not explicitly provided, it will be looked for in the library directory (based on its default name).

|Index|Type|Description|
|-|-|-|
|`libraries`|map|Map of library identifiers to file paths.|
|`combination_filters`|map|Map of combination filter identifiers to file paths.|

*E.g.*:

```json
{
  "libraries": {
    "lib1": "hgi_1.tsv",
    "lib2": "hgi_2.tsv"
  },
  "combination_filters": {
    "ft1": "/absolute/path/to/filter.tsv"
  }
}
```

#### Manifest

**Format**: JSON

**File name**: `manifest.json`

Manifest required to inform the merge process.

|Index|Type|Description|
|-|-|-|
|`sample_name`|string|Sample name.|
|`library_count`|integer|Number of libraries.|
|`total_library_templates`|integer|Total number of targets in the libraries.|
|`library_file_paths`||Per-library relative output file paths.|
|`library_independent_count_file_paths`||Library-independent quantification relative output file paths.|
|`total_dynamic_targets`|integer|Total number of unique regions captured as they were.|
|`dynamic_target_file_paths`||Unique region relative output file paths.|
|`combination_file_paths`||Per-combination relative output file paths.|

For a library-depedent quantification, *e.g.*:

```json
{
  "sample_name": "SAMPLE1",
  "library_count": 4,
  "total_library_templates": 10,
  "library_file_paths": {
    "0": {
      "counts": "lib.0.counts.tsv",
      "stats": "lib.0.stats.json"
    },
    "1": {
      "counts": "lib.1.counts.tsv",
      "stats": "lib.1.stats.json"
    },
    "2": {
      "counts": "lib.2.counts.tsv",
      "stats": "lib.2.stats.json"
    },
    "3": {
      "counts": "lib.3.counts.tsv",
      "stats": "lib.3.stats.json"
    }
  },
  "library_independent_count_file_paths": null,
  "total_dynamic_targets": 0,
  "dynamic_target_file_paths": null,
  "combination_file_paths": [
    {
      "counts": "combination.0.counts.tsv",
      "n": 2
    }
  ]
}
```

### Per library

One file per library.

#### Library

**Format**: TSV

**File name**: same as the ID in the experiment metadata

|Field|Type|Description|
|-|-|-|
|`id`|integer|Unique identifier.|
|`sequence`|DNA string|Unambiguous unique DNA string.|

#### Counts

**Format**: TSV

**File name**: `lib.<INDEX>.counts.tsv`

Targets are reported in the same order they appear in the library of origin.

|Index|Type|Description|
|-|-|-|
|1|integer|Experiment-wide target index.|
|2|string|Unique DNA sequence.|
|3|integer|Total number of exact matches.|

The indices may be used to decode the matching targets in the [match QC tables](#match-qc-table).

The number of matches includes any deriving from fallback attempts at matching a region, *e.g.* to determine whether a swap has occurred (`swap_libraries`).

#### Library-dependent statistics

**Format**: JSON

**File name**: `lib.<INDEX>.stats.json`

|Field|Format|Description|
|-|-|-|
|`mapped_to_template_reads`|integer|Total reads mapping to the library.|
|`mean_count_per_template`|decimal|Mean reads per template.|
|`median_count_per_template`|decimal|Median reads per template.|
|`total_templates`|integer|Total number of templates.|
|`swap_matching_reads`|integer|Total number of swap matches.|
|`zero_count_templates`|integer|Total number of unique templates with no reads mapping to them.|
|`low_count_templates_lt_15`|integer|Total number of unique templates with less than 15 reads mapping to them.|
|`low_count_templates_lt_30`|integer|Total number of unique templates with less than 30 reads mapping to them.|
|`gini_coefficient`|decimal|Gini coefficient of the mapping read counts.|
