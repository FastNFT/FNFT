# Installation

## Required tools

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

## Building under Windows

The following instructions have been tested under Windows 7/8/10. A simple way to build FNFT is via the [Scoop](http://scoop.sh/) packet manager. Make sure that [PowerShell 3](https://docs.microsoft.com/en-us/powershell/scripting/setup/installing-windows-powershell?view=powershell-6) is available. Run the command

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
    cmake .. -SHELL=cmd -G"MinGW Makefiles"

Build the library:

    mingw32-make -j4

Run the tests:

    mingw32-make -j4 test

_Note:_ If MATLAB is installed, the MATLAB interface should have been built
automatically. It can be found in the 'matlab' folder. Please copy the libfnft... dll file and the mex_... mexw... files manually from the 'build' into the 'matlab' folder.

## Building under MacOS

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
