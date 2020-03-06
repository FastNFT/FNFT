# Getting started

## Examples

FNFT comes with several examples, both in C and MATLAB, which are
located in the 'examples' directory. For example, run the following
commands in a terminal.

    cd ~/FNFT/examples/
    ./fnft_nsev_example
    ./fnft_nsep_example
    ./fnft_kdvv_example
    ./fnft_nsev_inverse_example

The sources for the examples can be found in the same directory. The example binaries have already been built together with the library. If you want to see the exact commands that were used to compile the examples, replace the 'make -j4' command in the build process with

    VERBOSE=1 make -j4

If the MATLAB interface has been built, try the following in MATLAB

    addpath ~/FNFT/matlab/
    cd ~/FNFT/examples
    mex_fnft_nsev_example
    mex_fnft_nsep_example
    mex_fnft_kdvv_example
    mex_fnft_inverse_example_1
    mex_fnft_inverse_example_2
    mex_fnft_inverse_example_3

## Documentation

The documentation of the C interface is available online at

    [https://fastnft.github.io/FNFT/](https://fastnft.github.io/FNFT/)

The C interface is separated into a public part ('fnft_' prefix) and a private part ('fnft__' prefix). The public part corresponds to the header files directly inside the 'include' folder. The header files in the subdirectories of the 'include' folder are considered private. To build the documentation yourself, run doxygen in the main folder of the library. It can then be found in the doc folder. Simply open the file '~/FNFT/doc/html/index.html' in your web browser.

The MATLAB interface is documented in the usual way. Run the commands

    help mex_fnft_nsev
    help mex_fnft_nsep
    help mex_fnft_kdvv
    help mex_fnft_nsev_inverse

in MATLAB to get more information.
