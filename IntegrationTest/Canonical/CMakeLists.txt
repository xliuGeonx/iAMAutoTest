# BARRACUDA - Copyright (C) 2015-2016 geonX
#
# Please report all bugs and problems to <support@geonx.com>.
#
# Test Cases of Solvers library CMakeLists
# -----------------------------------------------------------------------------
# ----------------------------- Project configuration -------------------------
# -----------------------------------------------------------------------------

#SET(${Solvers}_TEST FALSE CACHE BOOL "Generate testcases for $${Solvers} library - default: OFF")
LIST(APPEND CanonicalTest Thermal)
FOREACH(T_DIR ${CanonicalTest})
	ADD_SUBDIRECTORY(${T_DIR})
ENDFOREACH()

