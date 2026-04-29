# Changelog

## 0.7.1

- Fixed build isolation issues in `uv`/HPC environments by declaring explicit build requirements (`setuptools`, `wheel`, `packaging`) in `build-system.requires`.
- Switched `project.license` to SPDX string form (`"MIT"`) for compatibility with current setuptools behavior.

## 0.7

- Expanded Python support to 3.12 by updating package metadata to `>=3.10,<3.13`, improving compatibility with modern `uv` projects.
- Kept `PyFoam` as a required dependency and clarified installation guidance for `uv` users in the README.
- Simplified packaging metadata by removing unused `packaging` dependencies and using PEP 621 license table format.
