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

from tempfile import NamedTemporaryFile

import pytest

from luca.readers.utils import detect_encoding


@pytest.mark.parametrize('encoding,text', [
    ('ascii', 'Hello world!'),
    ('utf-8', 'Xin chào thế giới!'),
    ('utf-16', b'\xff\xfe'.decode('utf-16'))
])
def test_detect_encoding(encoding, text):
    with NamedTemporaryFile() as tmp:
        with open(tmp.name, 'w', encoding=encoding) as f:
            f.write(text)

        # Detect encoding
        detected_encoding = detect_encoding(tmp.name)
        assert detected_encoding is not None
        assert detected_encoding.lower() == encoding

        # Verify the original text can be read back with the detected encoding
        with open(tmp.name, encoding=detected_encoding, mode='r') as fh:
            assert fh.read() == text
