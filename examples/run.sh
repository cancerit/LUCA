#!/bin/bash

LIBRARY_DIR=libraries
SEED=32167
READ_LENGTH=60
READS=1000
MISMATCH_FRACTION=0.05  # [0.0, 1.0]

set -e

for fp in $(find . -name "experiment.json"); do
	dir=$(dirname $fp)
	data_fp="${dir}/data.bam"

	if [ ! -f $data_fp ]; then
		is_new=1

		printf "Generating test data in ${dir}...\n"
		../scripts/generate_test_data.py \
			-o ${data_fp} \
			-s ${SEED} \
			-n ${READS} \
			-r ${READ_LENGTH} \
			-M ${MISMATCH_FRACTION} \
			-l ${LIBRARY_DIR} \
			-m ${LIBRARY_DIR}/manifest.json \
			${fp}

	else
		is_new=0
	fi

	# Prepare output directory
	out_dir="${dir}/output"
	mkdir -p $out_dir

	printf "Testing ${dir}...\n"
	luca count \
		-l ${LIBRARY_DIR} \
		-o ${out_dir} \
		-s SAMPLE \
		--loglevel WARNING \
		${fp} \
		${data_fp}

	if [ "${is_new}" -eq "1" ]; then
		md5sum ${out_dir}/* > "${dir}/output.md5"
		mv -v "${out_dir}" "${out_dir}_exp"
	else
		md5sum -c "${dir}/output.md5" ${out_dir}
	fi

done

set +e
