include(SaveCompilerFlags)

if(CMAKE_C_COMPILER_WORKS AND CMAKE_Fortran_COMPILER_WORKS)
    include(CFlags)
#   include(CXXFlags)
    include(FortranFlags)
endif()
