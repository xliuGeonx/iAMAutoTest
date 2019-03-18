# iAMAutoTest 2019
# Developped By Xin
# Project CMakeLists
# -----------------------------------------------------------------------------
# ----------------------------- Project configuration -------------------------
# -----------------------------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

FUNCTION(COLLECT_SRC_FROM_DIRS DIRS)
  FOREACH(folder ${${DIRS}})
    FILE(
      GLOB
      source_files_one_dir
      ${folder}/[^.]*.h      
      ${folder}/[^.]*.hxx          
      ${folder}/[^.]*.ipp  
      ${folder}/[^.]*.hpp  
      ${folder}/[^.]*.hh
      ${folder}/[^.]*.cc      
      ${folder}/[^.]*.c
      ${folder}/[^.]*.cxx
      ${folder}/[^.]*.cpp
      ${folder}/[^.]*.occa
    )
    LIST(APPEND SOURCE_FILES_LIST ${source_files_one_dir})
  ENDFOREACH()
  SET(SOURCE_FILES ${SOURCE_FILES_LIST} PARENT_SCOPE)
ENDFUNCTION()

PROJECT(iAMAutoTest)
#-----------------------------------------------------------------------------
# ----------------------------- Project configuration -------------------------
# -----------------------------------------------------------------------------
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
MESSAGE(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
# find Boost
SET(Boost_NO_SYSTEM_PATHS true)
SET(Boost_USE_STATIC_LIBS ON CACHE BOOL "use static libraries from Boost")
FIND_PACKAGE(Boost REQUIRED COMPONENTS regex serialization iostreams system filesystem)
INCLUDE_DIRECTORIES(SYSTEM ${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

# find Barracuda
FIND_PACKAGE(Barracuda REQUIRED)

# find iAM
SET(USE_iAM_NIGHTBUILD OFF CACHE BOOL "AutoTest is configured for NightBuild")
SET(USE_iAM_RELEASE OFF CACHE BOOL "Using release version of iAM")
FIND_PACKAGE(iAM REQUIRED)

# vstest.console.exe
SET(VSTEST_EXE vstest.console.exe)

# ctest
SET(CTEST_USE_LAUNCHERS 1)
INCLUDE(CTest)
SET(CTEST_ENV_VARS )
LIST(APPEND CTEST_ENV_VARS "MSV_VERSION=${MSVC_VERSION}")

# add Check tool
ADD_SUBDIRECTORY(compare_results)

#add UnitTest
ADD_SUBDIRECTORY(UnitTest)