cmake_minimum_required(VERSION 3.2)

project(NEVPT2dmrg Fortran C)
# WARNING: Currently this does NOT work without ENABLE_MOLCAS
# HOWTO compile the sources from the nevpt2 source code directory:
#
# $ mkdir build && cd build
# $ ccmake ..
# $ make -j4
#
# done! you may want to set the compile type to Debug or Release in the ccmake GUI
#
# in order to set the fortran/c/c++ compiler use in the build directory:
# $ FC=ifort CC=icc CXX=icpc ccmake -DENABLE_DMRG=ON -DENABLE_64BIT_INTEGERS=ON/OFF ..

# The MKL search adapted from the Molcas CMakeLists.txt file

option(ENABLE_64BIT_INTEGERS  "Enable 64-bit integers"              ON)
option(ENABLE_DMRG            "Compile with DMRG interface support" OFF)
option(ENABLE_DEBUG_DMRG      "Compile with DEBUG printing for the DMRG interface development" OFF)
option(ENABLE_OPENMP          "Compile with OpenMP parallelization" OFF)
option(ENABLE_BOUNDS_CHECK    "Enable bounds check"                 OFF)
option(ENABLE_INTEL_MKL       "Enable Intel MKL library, requires $MKLROOT to be set" OFF)
option(ENABLE_MOLCAS          "Enable compilation within OpenMOLCAS framework" OFF)

set(MOLCAS_MATHLIBS "" CACHE STRING "Linear algebra libraries from OpenMOLCAS")
mark_as_advanced(MOLCAS_MATHLIBS)


# Replace the improvised separators with proper semicolons to recover the lists
if (MOLCAS_MATHLIBS)
  string(REPLACE "<->" ";" MOLCAS_MATHLIBS "${MOLCAS_MATHLIBS}")
endif(MOLCAS_MATHLIBS)

if (GA_LIBS)
  string(REPLACE "<->" ";" GA_LIBS "${GA_LIBS}")
endif(GA_LIBS)


set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
   )

if (NOT ENABLE_MOLCAS)
   include(ConfigCompilerFlags)
endif()

if (ENABLE_DMRG)
  project(NEVPT2dmrg Fortran C CXX)
endif()
######################################################################
# Version information
######################################################################
set(NEVPT_YEAR 2018)
set(NEVPT_VERSION_MAJOR 2)
set(NEVPT_VERSION_MINOR 1)
set(NEVPT_VERSION_BUILD "SCINE-QCMaquis")

if(EXISTS ${PROJECT_SOURCE_DIR}/.git)
    find_package(Git)
    if(GIT_FOUND)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-list --max-count=1 HEAD
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_REVISION
            ERROR_QUIET
            )
        if(NOT ${GIT_REVISION} STREQUAL "")
            string(STRIP ${GIT_REVISION} GIT_REVISION)
        endif()

        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_BRANCH
            ERROR_QUIET
            )
        if(NOT ${GIT_BRANCH} STREQUAL "")
            string(STRIP ${GIT_BRANCH} GIT_BRANCH)
        endif()

        set(NEVPT_VERSION_BUILD "${GIT_REVISION} (${GIT_BRANCH})")
    endif()
endif(EXISTS ${PROJECT_SOURCE_DIR}/.git)

set(NEVPT_VERSION "${NEVPT_VERSION_MAJOR}.${NEVPT_VERSION_MINOR}")
if(NEVPT_VERSION_BUILD)
  set(NEVPT_GIT_VERSION "${NEVPT_VERSION_BUILD}")
else(NEVPT_VERSION_BUILD)
  set(NEVPT_GIT_VERSION "")
endif(NEVPT_VERSION_BUILD)
set(NEVPT_VERSION_STRING "NEVPT version ${NEVPT_VERSION}")
MESSAGE(STATUS "NEVPT version: ${NEVPT_VERSION}")

######################################################################
# Configure files
######################################################################

configure_file(${PROJECT_SOURCE_DIR}/src/version.h.in ${PROJECT_SOURCE_DIR}/src/version.h)
######################################################################
# Install
######################################################################

install(FILES ${PROJECT_SOURCE_DIR}/src/version.h
        DESTINATION include COMPONENT headers)
install(PROGRAMS ${PROJECT_SOURCE_DIR}/tools-distributed-rdm/jobmanager.py
              ${PROJECT_SOURCE_DIR}/tools-distributed-rdm/prepare_rdm_template.sh
              ${PROJECT_SOURCE_DIR}/tools-distributed-rdm/submit-3rdm.sh
              ${PROJECT_SOURCE_DIR}/tools-distributed-rdm/submit-4rdm.sh
              DESTINATION Tools/distributed-4rdm COMPONENT runtime)



set (AUXLIB_F_SOURCE_FILES
     ${PROJECT_SOURCE_DIR}/src/info_state_energy.F90
     ${PROJECT_SOURCE_DIR}/src/nevpt2_cfg.F90
     ${PROJECT_SOURCE_DIR}/src/input_reader.F90
     ${PROJECT_SOURCE_DIR}/src/input_reader_sections.F90
     ${PROJECT_SOURCE_DIR}/src/parse_input.F90
     ${PROJECT_SOURCE_DIR}/src/inforb.F90
     ${PROJECT_SOURCE_DIR}/src/infsym.F90
     ${PROJECT_SOURCE_DIR}/src/indices_norm.F90
     ${PROJECT_SOURCE_DIR}/src/swap_util.F90
     ${PROJECT_SOURCE_DIR}/src/nevpt_header.F90
     ${PROJECT_SOURCE_DIR}/src/ord_utils.F90
     ${PROJECT_SOURCE_DIR}/src/eexx_functions.F90
     ${PROJECT_SOURCE_DIR}/src/ijkl_utils.F90
     ${PROJECT_SOURCE_DIR}/src/rdm_utils.F90
     ${PROJECT_SOURCE_DIR}/src/namel.f
     ${PROJECT_SOURCE_DIR}/src/cput.linux.f
     ${PROJECT_SOURCE_DIR}/src/filesplit.f
    )

set (AUXLIB_C_SOURCE_FILES
     ${PROJECT_SOURCE_DIR}/src/gctime.c
    )

set (KOOPMANSQD_SOURCE_FILES
     ${PROJECT_SOURCE_DIR}/src/koop_matrices.F90
     ${PROJECT_SOURCE_DIR}/src/koopro4QD.F90
    )
set (KOOPMANSQD_EXE_SOURCE_FILES
     ${PROJECT_SOURCE_DIR}/src/koopro4QD_main.F90
    )
set (QDNEVPT2_SOURCE_FILES
     ${PROJECT_SOURCE_DIR}/src/qdnevpt.F
    )
set (QDNEVPT2_EXE_SOURCE_FILES
     ${PROJECT_SOURCE_DIR}/src/qdnevpt_main.F90
    )


if(ENABLE_DEBUG_DMRG)
  add_definitions(-DDEBUG_DMRG_NEVPT)
endif(ENABLE_DEBUG_DMRG)

if(ENABLE_OPENMP)
  find_package (OpenMP)
  if (OPENMP_FOUND)
    add_definitions(-D_OPENMP_)
    add_definitions(${OpenMP_Fortran_FLAGS})
  endif (OPENMP_FOUND)
endif(ENABLE_OPENMP)

if (ENABLE_MOLCAS)
# Ignore MKL flag and obtain BLAS libraries from outside
  if (MOLCAS_MATHLIBS)
    list (APPEND LINK_LIBS ${MOLCAS_MATHLIBS})
    message(STATUS "LAPACK libraries from OpenMOLCAS: ${MOLCAS_MATHLIBS}")
  else()
    message(FATAL_ERROR "Please set MOLCAS_MATHLIBS when compiling with OpenMOLCAS.")
  endif(MOLCAS_MATHLIBS)
else()
if(ENABLE_INTEL_MKL)
  set (MKLROOT $ENV{MKLROOT} CACHE PATH "MKL root directory" FORCE)
  message ("-- using Intel Math Kernel Library (MKL)")
# 	mark_as_advanced (CLEAR MKLROOT)
        if (MKLROOT STREQUAL "")
                        message (FATAL_ERROR
                                "You must set the environment variable MKLROOT, "
                                "or specify -DMKLROOT=/path/to/mkl_root_dir     "
                                "when running cmake."
                                )
        endif ()
        # at this point, MKLROOT should be defined and not empty
        message ("-- MKLROOT = ${MKLROOT}")

        # here we check if MKLROOT has changed after it was last used
        # succesfully, if so then we need to reprocess the checks here.
        if (NOT MKLROOT STREQUAL MKLROOT_LAST)
                # reset MKL paths
                set (MKL_INCLUDE_PATH "${MKLROOT}/include" CACHE PATH
                        "location of MKL include files" FORCE)
                if (ENABLE_64BIT_INTEGERS)
                        set (libpath "${MKLROOT}/lib/intel64")
                else ()
                        set (libpath "${MKLROOT}/lib/intel64")
                endif ()
                set (MKL_LIBRARY_PATH ${libpath} CACHE PATH
                        "location of MKL libraries" FORCE)
                # uncache variables
                unset (LIBMKL_CORE CACHE)
                unset (LIBMKL_INTERFACE CACHE)
                unset (LIBMKL_SEQUENTIAL CACHE)
                unset (LIBMKL_THREADING CACHE)
                # cache last used MKLROOT
                set (MKLROOT_LAST ${MKLROOT} CACHE INTERNAL "last value" FORCE)
        endif ()
        message ("-- MKL_INCLUDE_PATH = ${MKL_INCLUDE_PATH}")
        message ("-- MKL_LIBRARY_PATH = ${MKL_LIBRARY_PATH}")

        # core library
        find_library (LIBMKL_CORE NAMES "mkl_core"
                PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
        # compiler-specific library interface
        if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
                find_library (LIBMKL_INTERFACE NAMES "mkl_gf_ilp64"
                        PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
        elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
                find_library (LIBMKL_INTERFACE NAMES "mkl_intel_ilp64"
                        PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
        endif ()
        # sequential/compiler-specific threading interface
        find_library (LIBMKL_SEQUENTIAL NAMES "mkl_sequential"
                PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
        if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
                find_library (LIBMKL_THREADING NAMES "mkl_gnu_thread"
                        PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
        elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
                find_library (LIBMKL_THREADING NAMES "mkl_intel_thread"
                        PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
        endif ()
        list (APPEND MKL_LIBRARIES ${LIBMKL_INTERFACE})
        list (APPEND MKL_LIBRARIES ${LIBMKL_CORE})
        if (ENABLE_OPENMP)
                list (APPEND MKL_LIBRARIES ${LIBMKL_THREADING})
        else ()
                list (APPEND MKL_LIBRARIES ${LIBMKL_SEQUENTIAL})
        endif ()

        include_directories (${MKL_INCLUDE_PATH})
        list (APPEND LINK_LIBS ${MKL_LIBRARIES})

else ()
	find_package(LAPACK REQUIRED)
	find_package(BLAS REQUIRED)
	message (STATUS "LAPACK libraries for NEVPT2: ${LAPACK_LIBRARIES}")
	message (STATUS "${LAPACK_LIBRARIES}")
endif(ENABLE_INTEL_MKL)
endif(ENABLE_MOLCAS)

if(ENABLE_MOLCAS)
  message (STATUS "Enabling compilation within OpenMOLCAS")
  mark_as_advanced(MOLCAS_BUILD_DIR)
  if ("${MOLCAS_BUILD_DIR}" STREQUAL "")
    message (FATAL_ERROR "Please set MOLCAS_BUILD_DIR before compiling with MOLCAS support!")
  endif ()
else ()
  message (FATAL_ERROR "This NEVPT2 version can only be compiled from within OpenMOLCAS.")
endif(ENABLE_MOLCAS)



  FIND_PACKAGE (HDF5) # find non-cmake built HDF5
  if(NOT HDF5_FOUND)
    #then search first the cmake-built HDF5
    FIND_PACKAGE (HDF5 NAMES hdf5)
  endif()

# HDF5 support
if(HDF5_FOUND)
    INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})
    list (APPEND LINK_LIBS qcmaquis-hdf5-interface ${HDF5_C_LIBRARIES})
else ()
  message (STATUS "Could not find HDF5 support" FATAL_ERROR)
endif(HDF5_FOUND)

if(ENABLE_DMRG)
  add_definitions(-DDMRG_NEVPT)

  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")

  set(ExternalProjectCMakeArgs
           -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
           -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
           )

    if(ENABLE_MOLCAS)
       set (proj_bin ${MOLCAS_BUILD_DIR}/External/qcmaquis/src/qcmaquis-driver-build/)
       if (ENABLE_DMRG)
          if (NOT MAQUIS_DMRG_DIR)
            set(MAQUIS_DMRG_DIR ${MOLCAS_BUILD_DIR}/qcmaquis/share/cmake)
          endif()

          find_package(MAQUIS_DMRG)
          if (NOT MAQUIS_DMRG_FOUND)
            message(FATAL_ERROR "Cannot found QCMaquis DMRG installation.")
          endif()
       endif(ENABLE_DMRG)

       link_directories (${proj_bin})
       include_directories (${DMRG_INCLUDE})

       message ("-- DMRG_INCLUDE_PATH: ${DMRG_INCLUDE}")
       message ("-- DMRG_LIBRARIES: ${MAQUIS_DMRG_LIBRARIES}")

    endif (ENABLE_MOLCAS)

    if(DEFINED GA_LIBS AND NOT "${GA_LIBS}" STREQUAL "")
      list (APPEND LINK_LIBS ${GA_LIBS})
      add_definitions(-D_MOLCAS_MPP_)
    endif()

  list (APPEND LINK_LIBS ${MAQUIS_DMRG_LIBRARIES} qcmaquis-driver)
  message (STATUS "DMRG interface succesfully added to target libs")
endif(ENABLE_DMRG)

add_library(AUXLIB_F ${AUXLIB_F_SOURCE_FILES})
add_library(AUXLIB_C ${AUXLIB_C_SOURCE_FILES})

add_library(koopro OBJECT  ${KOOPMANSQD_SOURCE_FILES})
add_library(qdnevpt OBJECT ${QDNEVPT2_SOURCE_FILES})

add_dependencies(koopro AUXLIB_F AUXLIB_C)
add_dependencies(qdnevpt AUXLIB_F AUXLIB_C)

macro(compile_standalone _executable _source LIST_OF_EXECUTABLES)
    add_executable(
        ${_executable}
        ${_source}
    )
    set_property(TARGET ${_executable} PROPERTY LINKER_LANGUAGE Fortran)
    set(LIST_OF_EXECUTABLES
        ${LIST_OF_EXECUTABLES}
        ${_executable}
    )
endmacro()

# set standalone executables

set(KOOPMANSQD_SRC $<TARGET_OBJECTS:koopro> ${KOOPMANSQD_EXE_SOURCE_FILES})
set(QDNEVPT_SRC $<TARGET_OBJECTS:qdnevpt> ${QDNEVPT2_EXE_SOURCE_FILES})

compile_standalone(koopro4QD.x "${KOOPMANSQD_SRC}" "${LIST_OF_EXECUTABLES}")
compile_standalone(QDnevpt2.x  "${QDNEVPT_SRC}" "${LIST_OF_EXECUTABLES}")

if (ENABLE_MOLCAS)
  add_library(qdnevpt2 $<TARGET_OBJECTS:koopro> $<TARGET_OBJECTS:qdnevpt>)
  add_dependencies(qdnevpt2 AUXLIB_F AUXLIB_C)
  target_link_libraries(qdnevpt2 ${LINK_LIBS} AUXLIB_F AUXLIB_C)
endif()

foreach(
    EXECUTABLE
    ${LIST_OF_EXECUTABLES}
    )
    if(${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
        SET_TARGET_PROPERTIES(${EXECUTABLE} PROPERTIES LINK_FLAGS "-Wl,-bbigtoc")
    endif()
    if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
        SET_TARGET_PROPERTIES(${EXECUTABLE} PROPERTIES LINK_FLAGS "-Wl,-E")
    endif()
    target_link_libraries(
        ${EXECUTABLE}
        AUXLIB_F
        AUXLIB_C
        ${LAPACK_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${LINK_LIBS}
        )
endforeach()
