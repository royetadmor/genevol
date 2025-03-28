cmake_minimum_required (VERSION 3.0.0)
project (genevol CXX)

# Compile options
set (CMAKE_CXX_FLAGS "-std=c++14 -Wall -Wshadow -Wconversion")

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

IF(NOT NO_DEP_CHECK)
  SET(NO_DEP_CHECK FALSE CACHE BOOL
      "Disable dependencies check for building distribution only."
      FORCE)
ENDIF(NOT NO_DEP_CHECK)

IF(NO_DEP_CHECK)
	MESSAGE(STATUS "Dependencies checking disabled. Only distribution can be built.")
ELSE(NO_DEP_CHECK)

# include openmp for parallelization of likelihood computation
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    set (OMP_NUM_THREADS "${OMP_NUM_THREADS}")
    set (OMP_DYNAMIC "${OMP_DYNAMIC}")
endif()

#### Optional Debug

option(USE_DEBUG "Enter Debug Mode" OFF)

IF (USE_DEBUG)
  MESSAGE("-- Enter Debug mode")
  add_compile_options(-DDEBUG)
ENDIF()


option(USE_VERSION_EIGEN_PIJT "Use Eigen-based PIJT version" OFF)

IF(USE_VERSION_EIGEN_PIJT)
  MESSAGE("-- Compiling with Eigen-decomposition-based PIJT version")
  add_compile_options(-DUSE_VERSION_EIGEN_PIJT)
ENDIF()

# Libtool-like version number
# CURRENT:REVISION:AGE => file.so.(C-A).A.R
# current:  The most recent interface number that this library implements.
# revision: The implementation number of the current interface.
# age:      The difference between the newest and oldest interfaces that this
#           library implements.
# In other words, the library implements all the interface numbers in the
# range from number current - age to current.
SET(${PROJECT_NAME}_VERSION_CURRENT "1")
SET(${PROJECT_NAME}_VERSION_REVISION "0")
SET(${PROJECT_NAME}_VERSION_AGE "0")

#static linkage?
IF(NOT BUILD_STATIC)
	SET(BUILD_STATIC FALSE CACHE BOOL
	  "Compile and install static libraries."
      FORCE)
ENDIF()
IF(BUILD_STATIC)
  MESSAGE(STATUS "Static libraries requested.")
ENDIF()

# Effective version number computation
MATH(EXPR ${PROJECT_NAME}_VERSION_MAJOR "${${PROJECT_NAME}_VERSION_CURRENT} - ${${PROJECT_NAME}_VERSION_AGE}")
SET(${PROJECT_NAME}_VERSION_MINOR ${${PROJECT_NAME}_VERSION_AGE})
SET(${PROJECT_NAME}_VERSION_PATCH ${${PROJECT_NAME}_VERSION_REVISION})
SET(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}")

set (PROJECT_VERSION ${${PROJECT_NAME}_VERSION})

# Find dependencies (add install directory to search)
if (CMAKE_INSTALL_PREFIX)
  set (CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
endif (CMAKE_INSTALL_PREFIX)

include (GNUInstallDirs)
find_package (bpp-core3 1.0.0 REQUIRED)
find_package (bpp-seq3 1.0.0 REQUIRED)
find_package (bpp-phyl3 1.0.0 REQUIRED)

# CMake package
set (cmake-package-location ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
include (CMakePackageConfigHelpers)
configure_package_config_file (
  package.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/package.cmake
  INSTALL_DESTINATION ${cmake-package-location}
  )
write_basic_package_version_file (
  ${CMAKE_CURRENT_BINARY_DIR}/package-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
  )
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/package.cmake DESTINATION ${cmake-package-location}
  RENAME ${PROJECT_NAME}-config.cmake)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/package-version.cmake DESTINATION ${cmake-package-location}
  RENAME ${PROJECT_NAME}-config-version.cmake)

# Use eigen
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

# Define the libraries
add_subdirectory (src)



ENDIF(NO_DEP_CHECK)

# Packager
SET(CPACK_PACKAGE_NAME "genevol")
SET(CPACK_PACKAGE_VERSION "1.0.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "GenEvol")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING.txt")
SET(CPACK_SOURCE_GENERATOR "TGZ")

SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_DEBSOURCE_PACKAGE_FILE_NAME "lib${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.orig")
INCLUDE(CPack)

#This adds the 'dist' target
ADD_CUSTOM_TARGET(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)

IF(NOT NO_DEP_CHECK)
IF (UNIX)

#This creates rpm packages:
ADD_CUSTOM_TARGET(rpm rpmbuild -ta ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(rpm dist)

ENDIF()

SET(CTEST_UPDATE_TYPE git)
SET(UPDATE_COMMAND "git")
SET(UPDATE_OPTIONS "")

# add_executable ( genEvol ./GenEvol/genEvol)
# target_link_libraries (genEvol)


# Tests
enable_testing ()
include (CTest)
if (BUILD_TESTING)
  add_subdirectory ( GenEvol )
endif (BUILD_TESTING)

ENDIF(NOT NO_DEP_CHECK)