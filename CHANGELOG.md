# Changelog

All notable changes to this project will be documented in this file.

This project does not adhere to [Semantic Versioning](http://semver.org/spec/v2.0.0.html), but it follows some aspects inspired from it.
In particular (even though this is not a strict, always respected rule), given a version number `X.Y.Z`,
- `Z` is incremented for minor changes (e.g. bug fixes),
- `Y` is incremented for substantial refactoring and/or the introduction of minor new functionality and
- `X` for the introduction of substantial new features.

Refer also to the [TODO](TODO.md) file to get more information of the changes occurred since the last release.

### Legend

 * :heavy_plus_sign: New feature
 * :heavy_check_mark: Enhancement
 * :recycle: Refactoring
 * :boom: Substantial change
 * :sos: Bug fix
 * :heavy_minus_sign: Removed feature

---

## [Unreleased]

 * :boom: Version `1.60` of Boost is now required, so that unit tests command line options are better handled (user needs to separate boost and user options via `--`).
 * :heavy_check_mark: The codebase has been revised on modern architectures and code has been improved following compiler warnings and hints.

---

## [Version 1.0] &nbsp;&nbsp; <sub><sup>26 September 2018</sub></sup>

This has been the first release of the repository.

[Unreleased]: https://github.com/AG-Philipsen/cl2qcd/compare/v1.0...HEAD
[Version 1.0]: https://github.com/AG-Philipsen/cl2qcd/releases/tag/v1.0
