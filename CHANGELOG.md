# Changelog



## [0.1.1] -- 2018-05-14

### Added

- Mex files now compile also with Matlab R2018a

### Changed

- Fixed: Several potential memory violations in fnft_nsev, poly_chirpz and nsev_testcases
- Fixed: Some return codes in fnft_nsep had not been checked
- Fixed: misc_merge did not work correctly for empty input vectors
- Fixed: The continuous spectrum in mex_fnft_nsev_example.m was plotted over t, not xi
- Fixed: A superfluous parameter kappa was given in the documentation of mex_fnft_kdvv