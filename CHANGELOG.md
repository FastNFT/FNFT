# Changelog

## [0.5.0]

### Added

- The routine fnft_kdvv now can now also compute the discrete spectrum. To locate the bound states, either Newton's method or a grid search with additional Newton refinements are available.
- A NFT routine fnft_manakovv for the Manakov equation with vanishing boundaries was added (continuous spectrum only).
- The slow scattering methods for AKNS-type systems now include a normalization procedure to deal with numerical overflow (enabled by default).
- The periodic NFT routine fnft_nsep now also supports pure Newton refinement.
- The vanishing NFT routine fnft_nsev now also supports manual filtering.
- The routine fnft__kdv_finvscatter has been added. The plan is to later use it for a fast inverse KdV NFT.

### Changed

- New criteria for stopping Newton iterations in fnft_nsep and fnft_nsev.
- The tolerance for the Newton refinements in fnft_nsev can now be set by the user.
- The default number of iterations for the Newton refinements in fnft_nsev is now 100, and the user is warned if the number was too small.
- Reduced some tests to reduce run times.
- The code for AKNS scattering has been overhauled.

### Fixed

- Several memory leaks have been fixed.

## [0.4.1] -- 2020-07-13

### Changed

- Number of samples for fnft_nsep again has to be a power of two.
- misc_resample no longer issues a warning when the signal appears to be undersampled. This gave the wrong impression that CFx_y discrizations suffer more in such scenarios than the other ones, which do not use this routine.

### Fixed

- misc_downsample could return incorrect values for first_last_index[1].
- Some errors in fnft_nsev become meaningless when bound_states==NULL and should not be risen in that case.

## [0.4.0] -- 2020-07-08

### Added

- The routine fnft_nsev for the NFT of the vanishing nonlinear Schroedinger equation now supports higher-order discretizations. The new fourth-order discretizations 4SPLIT4A and 4SPLIT4B support fast computation. Furthermore, higher-order methods for the classical non-fast approach were added. See the documentation of nse_discretizations_t for more information.
- The Matlab routine mex_fnft_nsev has been updated accordingly.
- The routine fnft_nsep for the NFT of the periodic nonlinear Schroedinger equation can now find spectra of quasi-periodic signals. It also has support for 4SPLIT4A and 4SPLIT4B discretizations.
- The Matlab routine mex_fnft_nsep has been updated accordingly.

### Changed

- The routine fnft_nsep for the NFT of the periodic nonlinear Schroedinger equation now only requires number of samples to be even (instead of a power of two). It has an extra input that specifies the phase shift over one quasi-period of the signal (see doc of fnft_nsep).
- Improved computation of norming constants in fnft__nse_scatter_bound_states (see doc of fnft__nse_scatter_bound_states).

## [0.3.0] -- 2020-03-06

### Added

- The routine fnft_nsep for the NFT of the periodic nonlinear Schroedinger equation can now visualize spines (see the new points_per_spine option)
- The subsampling and refinement processes in fnft_nsep can now be better controlled using the new Dsub and max_evals options
- The Matlab routine mex_fnft_nsep has been updated accordingly

### Fixed

- The Matlab interface now builds also with the most recent versions of Matlab
- The refinement of the auxiliary spectrum in fnft_nsep was refining the complex conjugate of mu_k instead of mu_k

### Changed

- The refinement of the main spectrum in fnft_nsep has been improved
- The polynomial rootfinder in fnft__poly_roots_fasteigen now always uses the QR algorithm

## [0.2.2] -- 2018-12-13

### Added

- The routine fnft_nsev_inverse_XI now raises an error if D<2 to avoid division by zero
- A .gitattribute file to correct incorrect language detection on GitHub

### Fixed

- CMake failed under Windows due to an incomplete install command in CMakeLists.txt
- Removed incorrect "-SHELL=cmd" parameter for CMake from Windows build instructions
- Documentation for fnft_nsev_inverse now states that D should be a positive power of two
- Documentation for fnft_nsev and fnft_nsep now states that passing *K_ptr==0 and *M_ptr==0 is not sufficient in order to completely skip the computation of the corresponding spectra
- Reformatted changelog for 0.1.1

## [0.2.1] -- 2018-09-28

### Fixed

- Computation of norming constants from residues in fnft_nsev_inverse was incorrect when continuous spectrum was present.

## [0.2.0] -- 2018-09-21

### Added

- New public function fnft_nsev_inverse implementing a fast inverse NFT for the vanishing nonlinear Schroedinger equation (the inversion of the continuous spectrum is fast, bound states are added using a standard Darboux transform)
- New public function fnft_version returning current version of FNFT
- The version of FNFT can now be supplemented with a suffix
- The FFTW library can be used to compute internal FFTs instead of the standard Kiss FFT
- Building of the tests can now be turned off
- The number of samples D in fnft_nsev and fnft_kdvv no longer has to be a power of two
- The subsampling factor in the SUBSAMPLE_AND_REFINE method for locating the discrete spectrum in fnft_nsev can now be controlled by the user
- System-wide installation using "make install" is now possible
- New private function fnft__nse_finvscatter implementing fast inverse scattering
- New private function fnft__poly_specfact implementing Kolmogorov's method for spectral factorization
- New private function fnft__poly_eval implementing Horner's method
- New private functions fnft__akns_fscatter and fnft__akns_scatter_matrix implementing forward scattering for general AKNS systems
- New private functions fnft__akns_lambda_to_z and fnft__akns_z_to_lambda that centralize the continuous-time/discrete-time coordinate transforms

### Changed

- Forward scattering functions fnft__nse_fscatter, fnft__nse_scatter_matrix, fnft__kdv_fscatter and fnft__kdv_scatter_matrix now make use of their more general AKNS counter parts

### Fixed

- M could not be much larger than D in fnft_nsev
- dists[1] was not set correctly under some circumstances in nsev_compare_nfs
- Some public headers were including private ones
- Compiler warning about __Thread_local

## [0.1.1] -- 2018-05-14

### Fixed

- Mex files now compile also with Matlab R2018a
- Several potential memory violations in fnft_nsev, poly_chirpz and nsev_testcases
- Some return codes in fnft_nsep had not been checked
- misc_merge did not work correctly for empty input vectors
- The continuous spectrum in mex_fnft_nsev_example.m was plotted over t, not xi
- A superfluous parameter kappa was given in the documentation of mex_fnft_kdvv
