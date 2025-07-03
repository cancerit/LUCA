#!/bin/bash

if [ -z ${TEST_DIR} ]; then
    echo "TEST_DIR not set!"
    exit 1
fi
PKG_DIR=$(python -c "import os;import luca;import inspect;print(os.path.dirname(inspect.getfile(luca)))")

echo "$(python --version)"
echo "Package source directory: ${PKG_DIR}"

pip install \
    pytest==8.4.1 \
    pytest-cov==6.2.1 && \
pytest --cov="${PKG_DIR}" "${TEST_DIR}"
