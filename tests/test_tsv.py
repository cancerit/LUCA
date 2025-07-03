from tempfile import NamedTemporaryFile

import pytest

from luca.readers.tsv import parse_tsv, TsvFieldNotFound, TsvTruncatedRow, TsvEmpty


def test_tsv_empty():
    with NamedTemporaryFile() as tmp:
        with pytest.raises(TsvEmpty):
            list(parse_tsv(tmp.name, ['a']))


def test_tsv_field_not_found():
    with NamedTemporaryFile() as tmp:
        with open(tmp.name, 'w') as fh:
            fh.write("a\tb\tc\n")
        with pytest.raises(TsvFieldNotFound):
            list(parse_tsv(tmp.name, ['MISSING_FIELD']))


def test_tsv_header_only():
    with NamedTemporaryFile() as tmp:
        with open(tmp.name, 'w') as fh:
            fh.write("a\tb\tc\n")
        data = list(parse_tsv(tmp.name, ['a']))
        assert data == []


def test_tsv_truncated():
    with NamedTemporaryFile() as tmp:
        with open(tmp.name, 'w') as fh:
            fh.write("a\tb\tc\n")
            fh.write("x\ty\tz\n")
            fh.write("x\n")
        with pytest.raises(TsvTruncatedRow):
            list(parse_tsv(tmp.name, ['b']))
