# LUCA
#
# Copyright (C) 2023 Genome Research Ltd.
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

import re

NT_UC = 'ACGT'
NT_LC = NT_UC.lower()
AMB_NT_UC = 'NRYKMSW'
AMB_NT_LC = AMB_NT_UC.lower()


def get_alphabet_re(symbols: str, ignore_case: bool = False) -> re.Pattern:
    return re.compile(
        f"[{symbols}]+",
        flags=re.IGNORECASE if ignore_case else re.NOFLAG)


dna_re = get_alphabet_re(NT_UC + AMB_NT_UC, ignore_case=True)
dna_masked_re = get_alphabet_re(NT_LC + AMB_NT_LC)
dna_ambiguous_re = get_alphabet_re(AMB_NT_UC + AMB_NT_LC)


non_non_amb_dna_re = re.compile(rf"[^{NT_UC}]", flags=re.IGNORECASE)
non_dna_re = re.compile(rf"[^{NT_UC}{AMB_NT_UC}]", flags=re.IGNORECASE)


def is_dna(seq: str) -> bool:
    return dna_re.fullmatch(seq) is not None


def is_masked(seq: str) -> bool:
    return dna_masked_re.search(seq) is not None


def is_ambiguous(seq: str) -> bool:
    return dna_ambiguous_re.search(seq) is not None
