cmake_minimum_required(VERSION 3.25)
project(geo_tools VERSION 0.1.0 LANGUAGES Fortran CXX C)

# Se configuran las variables estandar de fortran y c++
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_Fortran_PREPROCESS ON)

# Se configuran las banderas de fortran
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fbacktrace -ffree-form -ffree-line-length-none")

# Se configuran los directorios para guardar los modulos y ejecutables
set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
set(PROJECT_MODULE_BINARY_DIR "${CMAKE_BINARY_DIR}/modules")

# Crear directorio de módulos
file(MAKE_DIRECTORY "${PROJECT_MODULE_BINARY_DIR}")

# ← CORRECCIÓN: Módulos en BUILD dir, no SOURCE dir
set(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_MODULE_BINARY_DIR}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# Busqueda de dependencias base
set(BLA_VENDOR OpenBLAS)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

if(CMAKE_BUILD_TYPE MATCHES Release)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -march=native")
elseif(CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O0 -Wall -Wextra -Wunderflow")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -Wall -Wextra")
endif()

add_subdirectory(src)