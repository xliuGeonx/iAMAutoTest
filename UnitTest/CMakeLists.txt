# Barracuda - Copyright (C) 2015 geonX
#
# Please report all bugs and problems to <support@geonx.com>.
#
# UnitTest test case CMakeLists

# -----------------------------------------------------------------------------
# ----------------------------- Project configuration -------------------------
# -----------------------------------------------------------------------------


SET(TESTNAME iAM_UnitTest)
#SET(QUALITY_TESTNAME QA_${TESTNAME})
# -- Execution --
# --------------- ----------------------------------
#SET(RUN_DIR ${CMAKE_CURRENT_SOURCE_DIR})
#SET(QUALITY_RESULTS_FILE ${RUN_DIR}/values_check.txt)
#SET(QUALITY_REFERENCE_RESULTS_FILE ${RUN_DIR}/ref_values_check.txt)
SET(errorRegex "ERROR")
# validate testcase ---------------------------------------------------
ADD_TEST(NAME ${TESTNAME}_KernelType_Test
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:"TestCategory=KernelType" /Inisolation}
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_KernelType_Test PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_KernelType_Test PROPERTIES 
                         FAIL_REGULAR_EXPRESSION "${errorRegex}")
#SET_PROPERTY(TEST ${QUALITY_TESTNAME} PROPERTY LABELS ${VALIDATE_EXE})