# FNFT: Fast Nonlinear Fourier Transforms

[![Build Status](https://travis-ci.org/FastNFT/FNFT.svg?branch=master)](https://travis-ci.org/FastNFT/FNFT) [![DOI](http://joss.theoj.org/papers/10.21105/joss.00597/status.svg)](https://doi.org/10.21105/joss.00597)

FNFT is a software library for the fast numerical computation of (inverse) nonlinear Fourier transforms, which are also known as (inverse) scattering transforms. FNFT is written in C and comes with a MATLAB interface. A [Python interface](https://github.com/xmhk/FNFTpy) is available seperately.

## Currently Implemented Cases

### Forward Transforms

* Nonlinear Schroedinger equation

    * Vanishing boundary conditions
      * Reflection coefficient and/or scattering coefficients (a and b)
      * Bound states (eigenvalues)
      * Norming constants and/or residues

    * (Quasi-)Periodic boundary conditions
      * Main spectrum
      * Auxiliary spectrum

* Korteweg-de Vries equation
    * Vanishing boundary conditions (reflection coefficient only)

### Inverse Transforms

* Nonlinear Schroedinger equation

    * Vanishing boundary conditions
      * Inversion of reflection coefficients, b-scattering coefficients or the inverse Fourier transform of the b-coefficient
      * Bound states (eigenvalues) can be added with arbitrary norming constants/residuals

## Mailing List

Please join the FNFT mailing list if you want to be notified about new releases of FNFT. You can subscribe either using the [web interface](https://listserv.tudelft.nl/mailman/listinfo/fnft-announcements), or by sending an email with the subject "subscribe" to <fnft-announcements-request@lists.tudelft.nl>.

## Citation

If you use FNFT for your academic work, please cite the accompanying [software paper](https://doi.org/10.21105/joss.00597).  Latex users can use the following BibTex entry.

```
@article{FNFT2018,
    author  = {S. Wahls and S. Chimmalgi and P.J. Prins},
    title   = {{FNFT: A Software Library for Computing Nonlinear Fourier Transforms}},
    journal = {{The Journal of Open Source Software}},
    year    = {2018},
    volume  = {3},
    issue   = {23},
    pages   = {597},
    doi     = {10.21105/joss.00597},
    url     = {https://doi.org/10.21105/joss.00597},
    issn    = {2475-9066}
}
```

## Installation

Please follow the instructions in the file [INSTALL.md](INSTALL.md).

## Getting started

Please read the file [Getting-Started.md](Getting-Started.md).

## Community Guidelines

Please use the [issue tracker](https://github.com/FastNFT/FNFT/issues) to report any problems with the software. If you want to contribute to the development of FNFT, please [email](mailto:s.wahls##at##tudelft.nl) Sander Wahls.

## Contributors

* Sander Wahls, TU Delft
* Shrinivas Chimmalgi, TU Delft
* Peter J. Prins, TU Delft
* Marius Brehler, TU Dortmund

## License

FNFT is provided under the terms of the [GNU General Public License, version 2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html).

## Acknowledgements

* This project has received funding from the European Research Council (ERC)
  under the European Union’s Horizon 2020 research and innovation programme
  (grant agreement No 716669).

* FNFT incorporates code from the Fortran library [eiscor](https://github.com/eiscor/eiscor).

* FNFT incorporates code from the C library [Kiss FFT](http://kissfft.sourceforge.net/).

## References

The algorithms in FNFT utilize ideas from the following references. More information can be found in the documentation of the individual routines.

- S. Wahls and H. V. Poor, ["Introducing the fast nonlinear Fourier transform"](http://dx.doi.org/10.1109/ICASSP.2013.6638772), Proc. IEEE International Conference on Acoustics, Speech, and Signal Processing (ICASSP), Vancouver, Canada, May 2013.
- S. Wahls and H. V. Poor, ["Fast Numerical Nonlinear Fourier Transforms"](http://dx.doi.org/10.1109/TIT.2015.2485944), IEEE Transactions on Information Theory, vol. 61, no. 12, pp. 6957-6974, Dec. 2015.
- P. J. Prins and S. Wahls, ["Higher order exponential splittings for the fast non-linear Fourier transform of the Korteweg-de Vries equation"](https://doi.org/10.1109/ICASSP.2018.8461708), Proc.  2018 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP). Piscataway, NJ, USA: IEEE. 2018. pp. 4524-4528.
- G. Boffetta and A. R. Osborne, ["Computation of the direct scattering transform for the nonlinear Schroedinger equation"](https://doi.org/10.1016/0021-9991(92)90370-E), Journal of Computational Physics, vol. 102, no. 2, pp. 252-264, Oct. 1992.
- V. Aref, ["Control and Detection of Discrete Spectral Amplitudes in Nonlinear Fourier Spectrum"](https://arxiv.org/abs/1605.06328), Preprint, arXiv:1605.06328v1, May 2016.
- S. Hari and F. R. Kschischang, ["Bi-Directional Algorithm for Computing Discrete Spectral Amplitudes in the NFT"](https://doi.org/10.1109/JLT.2016.2577702), Journal of Lightwave Technology, vol. 34, no. 15, pp. 3529-3537, Aug. 2016.
- J. L. Aurentz, T. Mach, L. Robol, R. Vandebril and D. S. Watkins, ["Fast and backward stable computation of roots of polynomials, Part IIa: general backward error analysis"](http://www.cs.kuleuven.be/publicaties/rapporten/tw/TW683.abs.html), Technical Report no. TW 683, KU Leuven, Oct. 2017.
- V. Aref, S. T. Le and H. Buelow, ["Modulation over Nonlinear Fourier Spectrum: Continuous and Discrete Spectrum"](https://dx.doi.org/10.1109/JLT.2018.2794475), Journal of Lightwave Technology, vol. 36, no. 6, pp. 1289--1295, Mar. 2018.
- W. K. McClary, ["Fast seismic inversion"](https://doi.org/10.1190/1.1441417), Geophysics, vol. 48, no. 10, pp. 1371--1372, Oct. 1983.
- S. Wahls and H. V. Poor, ["Fast Inverse Nonlinear Fourier Transform For Generating Multi-Solitons In Optical Fiber"](http://dx.doi.org/10.1109/ISIT.2015.7282741), Proc. IEEE International Symposium on Information Theory (ISIT’15), pp. 1676–1680, Hong Kong, China, Jun. 2015.
- S. Wahls and V. Vaibhav, ["Fast Inverse Nonlinear Fourier Transforms for Continuous Spectra of Zakharov-Shabat Type"](http://arxiv.org/abs/1607.01305v2), Withdrawn Preprint, Dec. 2016. arXiv:1607.01305v2 [cs.IT]
- S. Wahls, ["Generation of Time-Limited Signals in the Nonlinear Fourier Domain via b-Modulation"](https://doi.org/10.1109/ECOC.2017.8346231), Proc. European Conference on Optical Communcation (ECOC), Gothenburg, Sweden, Sep. 2017.
- J. Skaar, L. Wang and T. Erdogan, ["On the synthesis of fiber Bragg gratings by layer peeling"](https://doi.org/10.1109/3.903065), IEEE Journal of Quantum Electronics, vol. 37, no. 2, pp. 165--173, Feb. 2001.
- S. Chimmalgi, P. J. Prins and S. Wahls, ["Fast Nonlinear Fourier Transform Algorithms Using Higher Order Exponential Integrators"](https://doi.org/10.1109/ACCESS.2019.2945480), IEEE Access, vol. 7, pp. 145161--145176, Oct. 2019.
- P. J. Prins and S. Wahls, ["Soliton Phase Shift Calculation for the Korteweg–De Vries Equation"](https://doi.org/10.1109/ACCESS.2019.2932256), IEEE Access, vol. 7, pp. 122914--122930, July 2019.
- S. Medvedev, I. Vaseva, I. Chekhovskoy and M. Fedoruk, ["Exponential fourth order schemes for direct Zakharov-Shabat problem"](https://doi.org/10.1364/OE.377140), Optics Express, vol. 28, pp. 20--39, 2020.
- J. Mertsching, ["Quasiperiodie Solutions of the Nonlinear Schroedinger Equation"](https://doi.org/10.1002/prop.2190350704), Fortschritte der Physik, vol. 35, pp. 519--536, 1987.