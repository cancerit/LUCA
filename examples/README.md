# Examples

To run all examples:

```shell
./run.sh
```

If no data is available for any given example, it will be generated from the corresponding experiment configuration.

Examples:

- single-end:
  - [library-independent](single_end/library_independent/experiment.json)
  - dual-guide:
    - [unfiltered](single_end/dual_guide/unfiltered/experiment.json)
    - [filtered](single_end/dual_guide/filtered/experiment.json)
    - [filtered & unfiltered](single_end/dual_guide/filtered_and_unfiltered/experiment.json)
- paired-end:
  - [library-independent](paired_end/library_independent/experiment.json)
  - dual-guide:
    - [unfiltered](paired_end/dual_guide/unfiltered/experiment.json)
    - [filtered](paired_end/dual_guide/filtered/experiment.json)
    - [filtered & unfiltered](paired_end/dual_guide/filtered_and_unfiltered/experiment.json)

## Adding a new example

Create a new subdirectory and a new experiment configuration (in JSON or YAML format), *e.g.*:

```shell
mkdir paired_end/triple_guide
touch paired_end/triple_guide/experiment.json
```

After filling in the experiment configuration as appropriate, generate test data, quantify it, and calculate the result checksums with the same script used for testing:

```shell
./run.sh
```

Verify the results are valid before commmitting the new directory.
