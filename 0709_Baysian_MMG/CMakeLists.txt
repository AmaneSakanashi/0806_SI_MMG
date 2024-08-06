cmake_minimum_required(VERSION 3.5)
# indentify minimum required version. 

enable_language(Fortran)
# designate language. Cmake mainly designed for C++, so for fortran use, set language,
find_package(OpenMP REQUIRED)


project(build_test Fortran)

# Comment out the following if using OpenMP
set(CMAKE_Fortran_FLAGS         "-O3 -p")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -fbacktrace -g -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow -Wall")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")

set(EXECUTABLE main.out)
# name of executable file.

add_executable(${EXECUTABLE}
  cmaes.f90
  eigen.f90
  random.f90
  problem_SI.f90
  MMG_model.f90
  main.f90
  #  problem.f90
  )
# files to build

# add link to library 
# comment out hereafter for openmp compllie
find_package(OpenMP)
if(OpenMP_Fortran_FOUND)
    target_link_libraries(${EXECUTABLE}
      PRIVATE
      OpenMP::OpenMP_Fortran
      )
endif()