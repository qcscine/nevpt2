# This variable is not guaranteed to be defined for all compilers or languages
if(NOT DEFINED CMAKE_Fortran_COMPILER_ID)
	message(FATAL_ERROR "CMAKE_Fortran_COMPILER_ID variable is not defined! (CMake Error)")
endif()

if(NOT DEFINED DEFAULT_Fortran_FLAGS_SET OR RESET_FLAGS)

message("NOTE: 64-bit integer..." ${ENABLE_64BIT_INTEGERS})

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU) # this is gfortran
    set(CMAKE_Fortran_FLAGS         "-g -fno-align-commons -fbacktrace -DVAR_GFORTRAN -fno-range-check")
    if(ENABLE_STATIC_LINKING)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -static"
            )
    endif()
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fdefault-integer-8"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fbounds-check"
            )
    endif()
    if(ENABLE_CODE_COVERAGE)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fprofile-arcs -ftest-coverage"
            )
    endif()
    if(ENABLE_DMRG)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -cpp"
            )
    endif()
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops -w")
    set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE} -pg")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(CMAKE_Fortran_FLAGS         "-w -assume byterecl -DINTEL -g -traceback")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip")
    set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE} -g -pg")

    if(DEFINED MKL_FLAG)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MKL_FLAG}")
    endif()

    if(ENABLE_DMRG)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fpp"
            )
    endif()

    if(ENABLE_STATIC_LINKING)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -static-libgcc -static-intel"
            )
    endif()
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -check bounds -fpstkchk -check pointers -check uninit -check output_conversion"
            )
    endif()

    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        message("--Switch off warnings due to incompatibility XCode 4 and Intel 11 on OsX 10.6")
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -Qoption,ld,-w"
            )
    endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
    set(CMAKE_Fortran_FLAGS         "-DVAR_PGF90")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE} -g -pg")
    if(ENABLE_STATIC_LINKING)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -Bstatic"
            )
    endif()
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} "
            )
    endif()
    if(ENABLE_CODE_COVERAGE)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} "
            )
    endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set(CMAKE_Fortran_FLAGS         "-qzerosize -qextname -qsuppress=cmpmsg")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE} -g -pg")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -qintsize=8 -q64"
            )
    endif()

    set_source_files_properties(${FREE_FORTRAN_SOURCES}
        PROPERTIES COMPILE_FLAGS
        "-qfree=f90"
        )
    set_source_files_properties(${FIXED_FORTRAN_SOURCES}
        PROPERTIES COMPILE_FLAGS
        "-qfixed"
        )
else()
	message(FATAL_ERROR "Vendor of your Fortran compiler is not supported")
endif()

if(DEFINED EXTRA_Fortran_FLAGS)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTRA_Fortran_FLAGS}")
endif()

save_compiler_flags(Fortran)
endif()
