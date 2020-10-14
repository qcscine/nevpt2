# This variable is not guaranteed to be defined for all compilers or languages
if(NOT DEFINED CMAKE_C_COMPILER_ID)
	message(FATAL_ERROR "CMAKE_C_COMPILER_ID variable is not defined! (CMake Error)")
endif()

if(NOT DEFINED DEFAULT_C_FLAGS_SET OR RESET_FLAGS)

if(CMAKE_C_COMPILER_ID MATCHES GNU)
    set(CMAKE_C_FLAGS         "-g")
    set(CMAKE_C_FLAGS_DEBUG   "-O0")
    set(CMAKE_C_FLAGS_RELEASE "-O2 -Wno-unused")
    if(ENABLE_STATIC_LINKING)
        set(CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} -static -fpic"
            )
    endif()
elseif(CMAKE_C_COMPILER_ID MATCHES Intel)
    set(CMAKE_C_FLAGS         "-g -wd981 -wd279 -wd383 -vec-report0 -wd1572 -wd177")
    set(CMAKE_C_FLAGS_DEBUG   "-O0")
    set(CMAKE_C_FLAGS_RELEASE "-O2")
    set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -shared-intel")

    if(DEFINED MKL_FLAG)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MKL_FLAG}")
    endif()
elseif(CMAKE_C_COMPILER_ID MATCHES PGI)
    set(CMAKE_C_FLAGS         "-g")
    set(CMAKE_C_FLAGS_DEBUG   "-O0")
    set(CMAKE_C_FLAGS_RELEASE "-O3")
elseif(CMAKE_C_COMPILER_ID MATCHES XL)
    set(CMAKE_C_FLAGS         "-qcpluscmt")
    set(CMAKE_C_FLAGS_DEBUG   " ")
    set(CMAKE_C_FLAGS_RELEASE "-O3")
elseif(CMAKE_CXX_COMPILER_ID MATCHES Clang)
    set(CMAKE_C_FLAGS         "-g")
    set(CMAKE_C_FLAGS_DEBUG   "-O0")
    set(CMAKE_C_FLAGS_RELEASE "-O2 -Wno-unused")
    if(ENABLE_STATIC_LINKING)
        set(CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} -Bstatic -fpic"
            )
    endif()
else()
	message(FATAL_ERROR "Vendor of your C compiler is not supported")
endif()

if(DEFINED EXTRA_C_FLAGS)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${EXTRA_C_FLAGS}")
endif()

save_compiler_flags(C)
endif()
