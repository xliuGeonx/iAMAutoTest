# BARRACUDA - Copyright (C) 2015-2016 geonX
#
# Please report all bugs and problems to <support@geonx.com>.
#
# Test Cases of Solvers library CMakeLists
# -----------------------------------------------------------------------------
# ----------------------------- Project configuration -------------------------
# -----------------------------------------------------------------------------

#SET(${Solvers}_TEST FALSE CACHE BOOL "Generate testcases for $${Solvers} library - default: OFF")
SET(RUN_INTETEST TRUE CACHE BOOL "Run iAM_InteTest - default: ON")
IF(RUN_INTETEST)
	LIST(APPEND PorjectVM ProjectVMTest)	
ENDIF()

LIST(APPEND TEST_DIRS ${PorjectVM} Canonical)
FOREACH(T_DIR ${TEST_DIRS})
	ADD_SUBDIRECTORY(${T_DIR})
ENDFOREACH()
