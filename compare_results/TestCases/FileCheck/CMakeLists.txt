# Barracuda - Copyright (C) 2015 geonX
#
# Please report all bugs and problems to <support@geonx.com>.
#
# CompareTwoFiles test case CMakeLists

# -----------------------------------------------------------------------------
# ----------------------------- Project configuration -------------------------
# -----------------------------------------------------------------------------


SET(TESTNAME BARRACUDA_FileCheck)
SET(QUALITY_TESTNAME QA_${TESTNAME})
# -- Execution --
# --------------- ----------------------------------
SET(RUN_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(QUALITY_RESULTS_FILE ${RUN_DIR}/ref.txt)
#SET(QUALITY_REFERENCE_RESULTS_FILE ${RUN_DIR}/ref_values_check.txt)

# validate testcase ---------------------------------------------------
ADD_TEST(NAME ${QUALITY_TESTNAME} 
  COMMAND $<TARGET_FILE:${VALIDATE_EXE}> --filecheck true -f ${QUALITY_RESULTS_FILE}  WORKING_DIRECTORY  ${RUN_DIR}
 )
SET_PROPERTY(TEST ${QUALITY_TESTNAME} PROPERTY LABELS ${VALIDATE_EXE})