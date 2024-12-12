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

import abc


class CustomException(Exception, abc.ABC):
    def __init__(self, message: str, *args: object) -> None:
        self.message = message
        super().__init__(*args)


class InvalidReadFileExtension(CustomException):
    pass


class InvalidInputFormatError(CustomException, abc.ABC):
    pass


class InvalidHTSError(InvalidInputFormatError):
    pass


class InputReadError(InvalidInputFormatError):
    pass


class MissingMetadataError(CustomException):
    pass


class InvalidLibraryError(InvalidInputFormatError):
    pass


class InvalidExperiment(Exception):
    pass


class UnsupportedData(CustomException):
    pass


class TargetNotFound(Exception):
    def __init__(self, target: str, *args: object) -> None:
        self.target = target
        super().__init__(*args)

    @property
    def message(self) -> str:
        return f"Target sequence '{self.target}' not found!"


class CustomFileException(Exception, abc.ABC):
    def __init__(self, fp: str, *args: object) -> None:
        self.fp = fp
        super().__init__(*args)

    @property
    @abc.abstractmethod
    def message(self) -> str:
        pass


class EmptyFileError(CustomFileException):
    @property
    def message(self) -> str:
        return f"Empty file at '{self.fp}'!"
