# Buscar todos los módulos en la carpeta modules/
file(GLOB FORTRAN_MODULE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/modules/*.f90")

# Crear librería con los módulos
add_library(geodesia_core ${FORTRAN_MODULE_SOURCES})

# Linkear con BLAS/LAPACK
target_link_libraries(geodesia_core PUBLIC
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
)

# Configurar propiedades
set_target_properties(geodesia_core PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
    LINKER_LANGUAGE Fortran
)

# Crear ejecutable desde main.f90
add_executable(geo_tool "${CMAKE_CURRENT_SOURCE_DIR}/main.f90")

# El ejecutable necesita los módulos compilados
target_include_directories(geo_tool PRIVATE
    "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

# Linkear el ejecutable con la librería y dependencias
target_link_libraries(geo_tool PRIVATE
    geodesia_core
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
)