# Barracuda - Copyright (C) 2016 geonX
#
# Please report all bugs and problems to <support@geonx.com>.
#
# ParsingYaml test case CMakeLists

# -----------------------------------------------------------------------------
# ----------------------------- Project configuration -------------------------
# -----------------------------------------------------------------------------

# Environment
SET(VALIDATE_EXE "${PROJECT_NAME}.compare_vectors")
SET(VALIDATE_EXE ${VALIDATE_EXE} PARENT_SCOPE)
# -- Include directories --
LIST(APPEND VALIDATE_EXE_INCLUDE_DIRS  ${CMAKE_CURRENT_SOURCE_DIR}/src)
SET(VALIDATE_EXE_INCLUDE_DIR ${VALIDATE_EXE_INCLUDE_DIRS})
SET(VALIDATE_EXE_INCLUDE_DIR ${VALIDATE_EXE_INCLUDE_DIR} PARENT_SCOPE)

IF(TARGET ${VALIDATE_EXE})
  MESSAGE("${VALIDATE_EXE} is already defined.")
ELSE()
  # -- Source FILEs list ---
  collect_src_from_dirs(VALIDATE_EXE_INCLUDE_DIR)
  ADD_EXECUTABLE(${VALIDATE_EXE} ${SOURCE_FILES})
  # -- Target --
  SET_PROPERTY(TARGET ${VALIDATE_EXE} PROPERTY CXX_STANDARD 98)
  TARGET_LINK_LIBRARIES(${VALIDATE_EXE} LINK_PUBLIC  ${Boost_LIBRARIES})
ENDIF()

ADD_SUBDIRECTORY(TestCases)
