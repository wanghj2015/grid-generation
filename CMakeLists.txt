cmake_minimum_required(VERSION 3.16)

project(grid_generation
        VERSION 1.0.0
        DESCRIPTION "Curvilinear field-line-following coordinate grid generation"
        LANGUAGES Fortran)

# ---------------------------------------------------------------------------
# Build type & flags
# ---------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo)
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-g -fcheck=all -fbacktrace")
endif()

# All binaries land in build/bin next to a copy of the IGRF coefficients.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ---------------------------------------------------------------------------
# Shared sources (used by both executables)
# ---------------------------------------------------------------------------
# mod_igrf's `igrf` subroutine calls the Schmidt Legendre routines, so those
# are part of the shared set even though the dipole program never invokes them.
set(COMMON_SRC
  src/mod_precision.f90
  src/mod_const.f90
  src/PlmSchmidt_d1.f95
  src/PlmIndex.f95
  src/mod_igrf.f90)

# ---------------------------------------------------------------------------
# Target 1: analytical dipole coordinates
# ---------------------------------------------------------------------------
add_executable(dipole_coordinates
  ${COMMON_SRC}
  src/dipole_mu.f90
  src/dipole_psi.f90
  src/dipole2spherical.f90
  src/dipole_coordinates.f90)

# ---------------------------------------------------------------------------
# Target 2: general IGRF field-line-following coordinates
# ---------------------------------------------------------------------------
add_executable(general_coordinates
  ${COMMON_SRC}
  src/mod_function.f90
  src/rksuite_90.f90
  src/zero_rc.f90
  src/mod_fieldline.f90
  src/general_coordinates.f90)

# Give each target a private module directory so the two builds of the shared
# modules never clobber one another's .mod files.
foreach(tgt dipole_coordinates general_coordinates)
  set_target_properties(${tgt} PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/${tgt})
endforeach()

# ---------------------------------------------------------------------------
# Runtime data: both programs read ./igrf13coeffs.txt from the working dir.
# Copy it next to the executables for convenient direct runs.
# ---------------------------------------------------------------------------
add_custom_command(
  OUTPUT  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/igrf13coeffs.txt
  COMMAND ${CMAKE_COMMAND} -E copy
          ${CMAKE_SOURCE_DIR}/data/igrf13coeffs.txt
          ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/igrf13coeffs.txt
  DEPENDS ${CMAKE_SOURCE_DIR}/data/igrf13coeffs.txt)
add_custom_target(copy_igrf_data ALL
  DEPENDS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/igrf13coeffs.txt)

# ---------------------------------------------------------------------------
# Install rules (optional): `cmake --install build --prefix <dir>`
# ---------------------------------------------------------------------------
install(TARGETS dipole_coordinates general_coordinates RUNTIME DESTINATION bin)
install(FILES data/igrf13coeffs.txt DESTINATION share/grid-generation)
