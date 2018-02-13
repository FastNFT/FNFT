# FNFT: Fast Nonlinear Fourier Transforms

FNFT is a software library for the fast numerical computation of nonlinear Fourier transforms, which are also known as direct scattering transforms. FNFT is written in C and comes with a MATLAB interface. It also contains some third-party Fortran code.

The algorithms are mainly based on ideas in [1]-[3], but also use ideas from other papers such as [4]-[8]. The following nonlinear Fourier transforms are currently implemented in FNFT:

* Nonlinear Schroedinger equation
    * Focusing and defocusing case
    * Vanishing boundary conditions
    * Periodic boundary conditions (main and auxiliary spectrum)

* Korteweg-de Vries equation
    * Vanishing boundary conditions (reflection coefficient only)

It is planned to extend this list further and also add inverse transforms in future releases. Please join the [FNFT mailling list](https://listserv.tudelft.nl/mailman/listinfo/fnft-announcements) if you want to be notified about new releases of FNFT.

## Installation

### Required tools

The following tools are needed in order to build FNFT:

* [CMake](https://cmake.org/)
* A C compiler.
* A Fortran compiler.

For the MATLAB interface, additionally the following is needed:

* A C++ compiler.
* [MATLAB](https://www.mathworks.com/products/matlab.html).

We use the [GNU compiler collection](https://gcc.gnu.org/) to build FNFT. The documentation is built with [Doxygen](www.doxygen.org).

### Building under Linux

FNFT is mostly developed under Linux. We describe the build process for
a fresh installation of [Ubuntu 16.04 Desktop](https://www.ubuntu.com/download/desktop).
First, install the missing tools by running the command

    sudo apt install git cmake gfortran doxygen

in a terminal. Then, extract the source code of FNFT into a directory and change into that directory. To install (for example) in '~/FNFT', run the commands

	cd ~
	git clone https://github.com/FastNFT/FNFT.git

We are now ready to build FNFT:

    cd ~/FNFT/
    mkdir build
    cd build
    cmake ..
    make -j4

To test that everything works as expected, run the command

	make -j4 test

If MATLAB is installed, the MATLAB interface should have been built
automatically. It can be found in the 'matlab' folder.

### Building under Windows

The following instructions have been tested under Windows 7/8. A simple way to build FNFT is via the [Scoop](http://scoop.sh/) packet manager. Make sure that [PowerShell 3](https://docs.microsoft.com/en-us/powershell/scripting/setup/installing-windows-powershell?view=powershell-6) is available. Run the command

	powershell.exe -ExecutionPolicy RemoteSigned

in the PowerShell to change the execution policy for the current Powershell session to allow installation of Scoop (Read more about [execution policies](https://docs.microsoft.com/en-us/powershell/module/microsoft.powershell.core/about/about_execution_policies?view=powershell-6)). Run the command

	iex (new-object net.webclient).downloadstring('https://get.scoop.sh')

in the PowerShell to install Scoop. Then, run the commands

	scoop install git
	scoop install cmake
	scoop install gcc

in the PowerShell to install the required tools. Change (for example) to your home directory,
download FNFT, and change into the new dir:

	git clone https://github.com/FastNFT/FNFT.git
	cd FNFT

Create a build directory and run cmake there

	mkdir build
	cd build
	cmake .. -G"MinGW Makefiles"

Build the library:

	mingw32-make SHELL=cmd -j4

Run the tests:

	mingw32-make SHELL=cmd -j4 test

_Note:_ If MATLAB is installed, the MATLAB interface should have been built
automatically. It can be found in the 'matlab' folder. Please copy the libfnft... dll file and the mex_... mexw... files manually from the 'build' into the 'matlab' folder.

### Building under MacOS

The following instructions have been tested under MacOS 10.12 (Sierra). The
simplest way to build FNFT is with the help of the [Homebrew](https://brew.sh) package manager. To install it, open a terminal and enter the command

	/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

Then, run the commands

	brew install cmake
	brew install gcc@7

Now, download the library in (for example) your home directory by running the command

    cd ~
    git clone https://github.com/FastNFT/FNFT.git

in the terminal. Finally, change into the directory and build using the commands

	cd FNFT
	mkdir build
	cd build
	CC=gcc-7 FC=gfortran-7 cmake ..
	make -j4

To test if the build was successful, run the command

	make -j4 test

## Getting started

### Examples

FNFT comes with several examples, both in C and MATLAB, which are
located in the 'examples' directory. For example, run the following
commands in a terminal.

	cd ~/FNFT/examples/
	./fnft_nsev_example
	./fnft_nsep_example
	./fnft_kdvv_example

The sources for the examples can be found in the same directory. The example binaries have already been built together with the library. If you want to see the exact commands that were used to compile the examples, replace the 'make -j4' command in the build process with

    VERBOSE=1 make -j4

If the MATLAB interface has been built, try the following in MATLAB

	addpath ~/FNFT/matlab/
	cd ~/FNFT/examples
	mex_fnft_nsev_example
	mex_fnft_nsep_example
	mex_fnft_kdvv_example

### Documentation

The C interface is separated in a public ('fnft_' prefix) and a private part ('fnft__' prefix). To get started with the public part, read the documentation in the public header files in the 'include' folder. It is also possible to build a html version of the documentation. To build it, run doxygen in the main folder of the library. It can then be found in the
doc folder. Simply open the file '~/FNFT/doc/html/index.html' in your web browser.

The MATLAB interface is documented in the usual way. Run the commands

	help mex_fnft_nsev
	help mex_fnft_nsep
	help mex_fnft_kdvv

in MATLAB to get more information.

## Contributors

* Sander Wahls, TU Delft
* Shrinivas Chimmalgi, TU Delft
* Peter J. Prins, TU Delft

## License

FNFT is provided under the terms of the [GNU General Public License, version 2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html). 

## Acknowledgements

* This project has received funding from the European Research Council (ERC)
  under the European Union’s Horizon 2020 research and innovation programme
  (grant agreement No 716669).

* FNFT incorporates code from the Fortran library 
[eiscor](https://github.com/eiscor/eiscor).

* FNFT incorporates code from the C library
[Kiss FFT](http://kissfft.sourceforge.net/).

## References

1. S. Wahls and H. V. Poor, ["Introducing the fast nonlinear Fourier transform"](http://dx.doi.org/10.1109/ICASSP.2013.6638772), Proc. IEEE International Conference on Acoustics, Speech, and Signal Processing (ICASSP), Vancouver, Canada, May 2013.
2. S. Wahls and H. V. Poor, [Fast Numerical Nonlinear Fourier Transforms](http://dx.doi.org/10.1109/TIT.2015.2485944), IEEE Transactions on Information Theory, vol. 61, no. 12, pp. 6957-6974, Dec. 2015.
3. P. J. Prins and S. Wahls, "Higher order exponential splittings for the fast non-linear Fourier transform of the Korteweg-de Vries equation", Accepted for presentation at the IEEE International Conference on Acoustics, Speech, and Signal Processing (ICASSP), Calgary, Canada, April 2018.
4. G. Boffetta and A. R. Osborne, ["Computation of the direct scattering transform for the nonlinear Schroedinger equation"](https://doi.org/10.1016/0021-9991(92)90370-E), Journal of Computational Physics, vol. 102, no. 2, pp. 252-264, Oct. 1992.
5. V. Aref, ["Control and Detection of Discrete Spectral Amplitudes in Nonlinear Fourier Spectrum"](https://arxiv.org/abs/1605.06328), Preprint, arXiv:1605.06328v1, May 2016.
6. S. Hari and F. R. Kschischang, ["Bi-Directional Algorithm for Computing Discrete Spectral Amplitudes in the NFT"](https://doi.org/10.1109/JLT.2016.2577702), Journal of Lightwave Technology, vol. 34, no. 15, pp. 3529-3537, Aug. 2016.
7. J. L. Aurentz, T. Mach, L. Robol, R. Vandebril and D. S. Watkins, ["Fast and backward stable computation of roots of polynomials, Part IIa: general backward error analysis"](http://www.cs.kuleuven.be/publicaties/rapporten/tw/TW683.abs.html), Technical Report no. TW 683, KU Leuven, Oct. 2017.
8. V. Aref, S. T. Le and H. Buelow, ["Modulation over Nonlinear Fourier Spectrum: Continuous and Discrete Spectrum"](https://doi-org.tudelft.idm.oclc.org/10.1109/JLT.2018.2794475), Journal of Lightwave Technology, to appear.
