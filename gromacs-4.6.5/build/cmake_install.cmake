# Install script for directory: /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local/gromacs")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "data")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/gromacs/share/gromacs/COPYING")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/usr/local/gromacs/share/gromacs" TYPE FILE FILES "/Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/COPYING")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "data")

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/cmake_install.cmake")
  INCLUDE("/Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/share/cmake_install.cmake")
  INCLUDE("/Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/include/cmake_install.cmake")
  INCLUDE("/Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/cmake_install.cmake")
  INCLUDE("/Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/scripts/cmake_install.cmake")
  INCLUDE("/Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/tests/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

IF(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
ELSE(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
ENDIF(CMAKE_INSTALL_COMPONENT)

FILE(WRITE "/Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/${CMAKE_INSTALL_MANIFEST}" "")
FOREACH(file ${CMAKE_INSTALL_MANIFEST_FILES})
  FILE(APPEND "/Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
ENDFOREACH(file)
