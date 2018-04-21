---
title: 'FNFT: A Software Library for Computing Nonlinear Fourier Transforms'
tags:
- Nonlinear Fourier transform
- Scattering transform
authors:
- name: Sander Wahls
  orcid: 0000-0001-5159-0918
  affiliation: 1
- name: Shrinivas Chimmalgi
  orcid: 0000-0002-5868-3102
  affiliation: 1
- name: Peter J Prins
  orcid: 0000-0003-4692-5538
  affiliation: 1
affiliations:
 - name: Delft Center for Systems and Control, TU Delft, Delft, The Netherlands
   index: 1
date: 13 February 2018
bibliography: paper.bib
---

# Summary

The conventional Fourier transform was originally developed in order to solve the heat equation,
which is a standard example for a linear evolution equation. *Nonlinear Fourier transforms (NFTs)*^[NFTs are also known as *direct scattering transforms* in the literature.] are generalizations of the conventional Fourier transform that can be used to solve certain nonlinear evolution equations in a similar way [@Ablowitz1974]. An important difference to the conventional Fourier transform is that NFTs are equation-specific. The *Korteweg-de Vries (KdV) equation* [@Gardner1967] and the *nonlinear Schroedinger equation (NSE)* [@shabat1972exact] are two popular examples for nonlinear evolution equations that can be solved using appropriate NFTs.

NFTs are well-established theoretical tools in physics and mathematics, but in several areas of engineering they have only recently begun to draw significant attention. An idealized fiber-optic communication channel can be described by the NSE, which is solvable with a NFT. In the last few years, there has been much interest in
using NFTs for data transmission in optical fibers. We refer to @Turitsyn2017 for a recent review. Several invited papers and tutorials on the topic have been presented at major conferences such as the Optical Fiber Communication Conference and Exhibition (OFC) and the European Conference and Exhibition on Optical Communication (ECOC) in the last few years. Several new research projects, such as the [MSCA-ITN COIN](http://www.coinproject.eu/) or the ERC Starting Grant from which this project is funded (see below), have been initiated. Another area in which NFTs have found practical application is the analysis of waves in shallow water [@Osborne2010;@Bruehl2016]. Here, the KdV-NFT is typically used.

The implementation of a numerical algorithm for computing a NFT is however not as simple as for the conventional Fourier transform. While quite a few algorithms have been presented in the literature, *there is not a single publically available software library that implements numerical NFTs.* We believe that the lack of a reliable and efficient software library for computing NFTs is currently hindering progress in the field. For this reason, we have published *FNFT* on [GitHub](https://github.com/FastNFT/FNFT). FNFT, which is short for "Fast Nonlinear Fourier Transforms", is a software library that provides implementations of the fast NFT algorithms that were developed by some of the authors [@Wahls2013;@Wahls2015;@Prins2018]. Our goal was to develop an efficient, easy to use and reliable library. Therefore, FNFT was written in C and ships with a MATLAB interface as well as currently more than 60 unit and integration tests. To simplify the build process as much as possible, FNFT uses [CMake](http://www.cmake.org). FNFT has no external dependencies, but ships with some 3rd party code from the open source projects [EISCOR](https://github.com/eiscor/eiscor) (also see [@Aurentz2017]) and [KissFFT](http://kissfft.sourceforge.net/).    

# Features

FNFT currently provides algorithms for the numerical computation of the following NFTs: 

* NSE with vanishing boundary conditions,
* NSE with periodic conditions (main and auxiliary spectrum),
* KdV with vanishing boundary conditions (continuous spectrum only).

# Acknowledgements

This project has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (grant agreement No 716669).

# References
