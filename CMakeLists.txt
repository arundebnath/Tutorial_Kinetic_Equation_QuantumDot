cmake_minimum_required(VERSION 3.18)
project(QuantumDotHEOM Fortran)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)

add_executable(qdot
  constants.f90
  bath_model.f90
  dynamics.f90
  ode_solver.f90
  main.f90
)

target_link_libraries(qdot PRIVATE OpenMP::OpenMP_Fortran)
target_compile_options(qdot PRIVATE -O3 -Wall -Wextra -pedantic)
