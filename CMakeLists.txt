# This file is part of FNFT.
#
# FNFT is free software; you can redistribute it and/or
# modify it under the terms of the version 2 of the GNU General
# Public License as published by the Free Software Foundation.
#
# FNFT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# Contributors:
# Sander Wahls (TU Delft) 2017-2018.
# Marius Brehler (TU Dortmund) 2018.

cmake_minimum_required(VERSION 3.5)
project(fnft C)

set(FNFT_VERSION_MAJOR 0)
set(FNFT_VERSION_MINOR 1)
set(FNFT_VERSION_PATCH 1)

include(CheckIncludeFiles)
include(CheckFunctionExists)
include(CheckCSourceCompiles)
include(CheckLanguage)

# options
option(DEBUG "Compile with debugging symbols" OFF)
option(WITH_MATLAB "Build the Matlab interface" ON)
option(MACHINE_SPECIFIC_OPTIMIZATION "Activate optimizations specific for this machine" ON)
option(ADDRESS_SANITIZER "Enable address sanitzer for known compilers" OFF)
option(ENABLE_FFTW "Use FFTW if it is available" OFF)
option(BUILD_TESTS "Build tests" ON)

# check for complex.h
check_include_files(complex.h HAVE_COMPLEX_H)
if (NOT HAVE_COMPLEX_H)
    message(FATAL_ERROR "Compiler does not provide complex.h")
endif()

# check if math code requires us to link libm
check_function_exists(log2 HAVE_LOG2)
if (NOT HAVE_LOG2)
	message("++ Adding libm to the required libraries, repeating test.")
	list(APPEND CMAKE_REQUIRED_LIBRARIES m)
	check_function_exists(log2 HAVE_LOG2_2)
	if (NOT HAVE_LOG2_2)
        message(FATAL_ERROR "no log2")
	else()
		set(LIBM m)	
	endif()
else()
	set(LIBM )
endif()

# check if thread local storage is available
check_c_source_compiles("_Thread_local int g; int main() { g = 1; return g; }" HAVE__THREAD_LOCAL)
if (NOT HAVE__THREAD_LOCAL)
	check_c_source_compiles("__thread int g; int main() { g = 1; return g; }" HAVE___THREAD)
	if (NOT HAVE___THREAD)
		message(WARNING "Thread local storage is not available. fnft_errwarn_setprintf will not be thread-safe.")
	endif()
endif()

# check if FFTW3 is available
find_library(FFTW3_LIB fftw3)
find_path(FFTW3_INCLUDE fftw3.h)
if (FFTW3_LIB AND FFTW3_INCLUDE)
    if (ENABLE_FFTW)
        message("++ FFTW3 found and enabled. Run cmake with \"-DENABLE_FFTW=OFF\" to disable.")
        set(HAVE_FFTW3 1) # for updating fnft_config.h
    else()
        message("++ FFTW3 found but NOT enabled. Run cmake with \"-DENABLE_FFTW=ON\" to enable.")
    endif()
else()
    if (ENABLE_FFTW)
        message(FATAL_ERROR "FFTW3 NOT found but set to enabled the by user. Run cmake WITHOUT \"-DENABLE_FFTW=ON\" to disable.")
    else()
        set(FFTW3_LIB "")
        set(FFTW3_INCLUDE "")
        message("++ FFTW3 NOT found. Using Kiss FFT instead.")
    endif()
endif()

# header files
include_directories(include)
include_directories(include/3rd_party/eiscor)
include(CheckCCompilerFlag)

# configure C compiler
if (CMAKE_COMPILER_IS_GNUCC) # gcc
	if (DEBUG)
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
        message("++ Enabling debugging symbols in the C compiler")
	endif()
    if (ADDRESS_SANITIZER)
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
        message("++ Enabling address sanitizer")
    endif()
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
	if (MACHINE_SPECIFIC_OPTIMIZATION)
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
        message("++ Enabling machine specific optimization in the C compiler")
	endif()
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wextra -pedantic -Werror=implicit-function-declaration")
else()
	message(WARNING "++ Compiler is not gcc. Will try to set flags anyway.")
	if (DEBUG)    
		check_c_compiler_flag("-g" HAS_C_G_FLAG)
		if (HAS_C_G_FLAG)
			set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS}$ -g")
            message("++ Enabled -g flag in the C compiler")
		endif()
	else()
		check_c_compiler_flag("-O3" HAS_C_O3_FLAG)
       	if (HAS_C_O3_FLAG)
			set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
            message("++ Enabled -O3 flag in the C compiler")
		endif()
		if (MACHINE_SPECIFIC_OPTIMIZATION)
       		check_c_compiler_flag("-march=native" HAS_C_MARCH_NATIVE_FLAG)
       		if (HAS_C_MARCH_NATIVE_FLAG)
   				set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
                message("++ Enabled -march=native flag in the C compiler")
       		endif()
		endif()
	endif()
endif()

# generate config header file
configure_file(include/fnft_config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/include/fnft_config.h)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/private)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/3rd_party/kiss_fft)
include_directories(${FFTW3_INCLUDE})

# configure Fortran support
enable_language(Fortran)
if (CMAKE_Fortran_COMPILER_ID MATCHES GNU) # gfortran
    set (CMAKE_Fortran_FLAGS " -O3 -cpp -ffree-line-length-none")
    if (DEBUG)
        set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g")
        message("++ Enabling debugging symbols in the Fortran compiler")
	endif()
	if (MACHINE_SPECIFIC_OPTIMIZATION)
        set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native")
        message("++ Enabling machine specific optimization in the Fortran compiler")
	endif()
endif ()

if (BUILD_TESTS)
	enable_testing()
endif()

# sources
file(GLOB SOURCES "src/*.c")
file(GLOB PRIVATE_SOURCES "src/private/*.c")
file(GLOB KISS_FFT_SOURCES "src/3rd_party/kiss_fft/*.c")
file(GLOB EISCOR_SOURCES "src/3rd_party/eiscor/*.f90")
file(GLOB EXAMPLE_SOURCES "examples/*.c")
file(GLOB_RECURSE TEST_SOURCES "test/*_test*.c")

# generate shared library
add_library(fnft SHARED ${SOURCES} ${PRIVATE_SOURCES} ${KISS_FFT_SOURCES} ${EISCOR_SOURCES})
target_link_libraries(fnft ${FFTW3_LIB})
set_target_properties(fnft PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")

# installation under Linux
install(TARGETS fnft DESTINATION lib)

# generate examples
foreach (srcfile ${EXAMPLE_SOURCES})
	get_filename_component(example ${srcfile} NAME_WE)
	add_executable(${example} ${srcfile})
    target_link_libraries(${example} fnft ${LIBM} ${FFTW3_LIB})
	set_target_properties(${example} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/examples")
endforeach()

# generate tests
if (BUILD_TESTS)
	foreach (srcfile ${TEST_SOURCES})
		get_filename_component(dir ${srcfile} DIRECTORY)
		get_filename_component(test ${srcfile} NAME_WE)
		add_executable(${test} ${srcfile})
	    target_link_libraries(${test} fnft ${LIBM} ${FFTW3_LIB})
		add_test(NAME ${test} COMMAND ${test})
		set_target_properties(${test} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${dir}")
	endforeach()
endif()

# Try to build the Matlab interface if requested
if (WITH_MATLAB)
    check_language(CXX)
    if (CMAKE_CXX_COMPILER)
    	enable_language(CXX) # sometimes required by FindMatlab below
    else()
        message("++ Could not find a C++ compiler, which is sometimes required to localize Matlab.")
        message("++ If the build process fails, deactivate the localization with the -DWITH_MATLAB=OFF switch.")
    endif()
	find_package(Matlab COMPONENTS MX_LIBRARY)
	if (Matlab_FOUND)
		file(GLOB MEX_SOURCES "matlab/mex_*.c")
		foreach (srcfile ${MEX_SOURCES})
			get_filename_component(mexfile ${srcfile} NAME_WE)
      matlab_add_mex(NAME ${mexfile} SRC ${srcfile} LINK_TO fnft ${LIBM} ${FFTW3_LIB})
			set_target_properties(${mexfile} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/matlab")
      target_compile_options(${mexfile} PRIVATE -DMEX_DOUBLE_HANDLE)
      # -DMEX_DOUBLE_HANDLE avoids the new complex interleaved API      
		endforeach()
	endif()
endif()
