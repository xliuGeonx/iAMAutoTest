# Barracuda - Copyright (C) 2015 geonX
#
# Please report all bugs and problems to <support@geonx.com>.
#
# CompareTwoFiles test case CMakeLists

# -----------------------------------------------------------------------------
# ----------------------------- Project configuration -------------------------
# -----------------------------------------------------------------------------


SET(TESTNAME BARRACUDA_CompareFourFiles)
SET(QUALITY_TESTNAME QA_${TESTNAME})
# -- Execution --
# --------------- ----------------------------------
SET(RUN_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(QUALITY_RESULTS_FILE1 ${RUN_DIR}/values_check_01.txt)
SET(QUALITY_REFERENCE_RESULTS_FILE1 ${RUN_DIR}/ref_values_check_01.txt)
SET(QUALITY_RESULTS_FILE2 ${RUN_DIR}/values_check_02.txt)
SET(QUALITY_REFERENCE_RESULTS_FILE2 ${RUN_DIR}/ref_values_check_02.txt)
SET(QUALITY_RESULTS_FILE3 ${RUN_DIR}/values_check_03.txt)
SET(QUALITY_REFERENCE_RESULTS_FILE3 ${RUN_DIR}/ref_values_check_03.txt)
SET(QUALITY_RESULTS_FILE4 ${RUN_DIR}/values_check_04.txt)
SET(QUALITY_REFERENCE_RESULTS_FILE4 ${RUN_DIR}/ref_values_check_04.txt)

# validate testcase ---------------------------------------------------
ADD_TEST(NAME ${QUALITY_TESTNAME} 
  COMMAND $<TARGET_FILE:${VALIDATE_EXE}> -v --l2_norm -f ${QUALITY_RESULTS_FILE1}  ${QUALITY_RESULTS_FILE2}  ${QUALITY_RESULTS_FILE3}  ${QUALITY_RESULTS_FILE4} -r ${QUALITY_REFERENCE_RESULTS_FILE1} ${QUALITY_REFERENCE_RESULTS_FILE2} ${QUALITY_REFERENCE_RESULTS_FILE3} ${QUALITY_REFERENCE_RESULTS_FILE4} WORKING_DIRECTORY  ${RUN_DIR}
 )
SET_PROPERTY(TEST ${QUALITY_TESTNAME} PROPERTY LABELS ${VALIDATE_EXE})