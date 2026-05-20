cmake_minimum_required (VERSION 3.0)

project(main Fortran)

# Lapack
find_package( LAPACK REQUIRED )
if ( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" )
    set (BLA_VENDOR OpenBLAS)
elseif ( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" )
    set (BLA_VENDOR Intel)
endif()
# set (BLA_VENDOR OpenBLAS)
# set (BLA_VENDOR Intel)
if ( LAPACK_FOUND STREQUAL false )
    set (BLA_VENDOR All)
endif()
if ( BLA_VENDOR STREQUAL "OpenBLAS" )
    set( LAPACK_LINKER_FLAGS "-lopenblas" )
elseif ( BLA_VENDOR STREQUAL "Intel")
    set( LAPACK_LINKER_FLAGS "-lmkl_rt -lmkl_core -lmkl_sequential" )
elseif( BLA_VENDOR STREQUAL "All" )
    set( LAPACK_LINKER_FLAGS "-llapack" )
endif()

# Source
file(GLOB SRC_LIST "*.f90")


# compiler and compile flag
set (CMAKE_Fortran_FLAGS "-O3 ")
#set (CMAKE_Fortran_FLAGS "-O2 -xAVX -ip ")
#set (CMAKE_Fortran_FLAGS "-O2 -xSSE4.2 -static ")
#set (CMAKE_Fortran_FLAGS "-O2 -xSSE4.2 -unroll ")

# debug
set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O0 -g -C" )
#set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wunused-variable" )

# add Lapack
set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${LAPACK_LINKER_FLAGS}" )

# main
if ( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" )
    set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none" )
elseif ( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" )
    set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -free" )
endif()

# mesg
MESSAGE( "" )
MESSAGE( "Fortran COMPILER      = " ${CMAKE_Fortran_COMPILER_ID} )
MESSAGE( "Fortran COMPILE FLAGS = " ${CMAKE_Fortran_FLAGS} )
MESSAGE( "Lapack Vendor         = " ${BLA_VENDOR} )
MESSAGE( "" )

# compile
add_executable (main ${SRC_LIST})
