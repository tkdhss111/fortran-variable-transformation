set ( EXE unit_test )
set ( SRCS ./unit_test.f90
           ../src/variable_transformation_mo.f90
    )

set ( CMAKE_Fortran_FLAGS
  "-ffree-line-length-none -Wno-uninitialized -Wno-maybe-uninitialized -fcoarray=single")
set ( CMAKE_Fortran_FLAGS_DEBUG
  "-Ddebug -O0 -fbounds-check -fbacktrace -Wall" )

add_executable ( ${EXE} ${SRCS} )
add_test ( NAME ${EXE} COMMAND ${EXE} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )
