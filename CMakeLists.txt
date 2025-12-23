cmake_minimum_required ( VERSION 3.21 )
set ( CMAKE_Fortran_COMPILER gfortran )
set ( CMAKE_C_COMPILER gcc )
project ( hello LANGUAGES Fortran VERSION 1.0.0 )
enable_language ( Fortran C )
enable_testing ()

add_subdirectory ( test )
