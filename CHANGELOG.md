# Changelog

## [1.0.2] - 2025-07-15

### Fixed

- Fix library stats serialisation for NumPy 2.3.*
- Fix ambiguous (non-canonical) base pattern
- Fix error message formatting on invalid combination configuration
- Handle input manifest failed look-ups
- Handle invalid DNA sequence errors

### Changed

- Forbid extra attributes in the experiment configuration

## [1.0.1] - 2025-04-07

### Fixed

- Fix combination counter initialisation: when a list of expected combinations was provided, all combination counts were incorrectly initialised to one

## [1.0.0] - 2024-12-12

First public release.
