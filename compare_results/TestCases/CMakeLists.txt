# BARRACUDA - Copyright (C) 2015-2016 geonX
#
# Please report all bugs and problems to <support@geonx.com>.
#
# Compare_results library CMakeLists
# -----------------------------------------------------------------------------
# ----------------------------- Project configuration -------------------------
# -----------------------------------------------------------------------------

LIST(APPEND testCases CompareTwoFiles CompareFourFiles FileCheck)


FOREACH(T_DIR ${testCases})
  ADD_SUBDIRECTORY(${T_DIR})
ENDFOREACH()

