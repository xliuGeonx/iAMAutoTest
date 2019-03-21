# Barracuda - Copyright (C) 2015 geonX
#
# Please report all bugs and problems to <support@geonx.com>.
#
# UnitTest test case CMakeLists

# -----------------------------------------------------------------------------
# ----------------------------- Project configuration -------------------------
# -----------------------------------------------------------------------------


SET(TESTNAME CanonicalThermalTester)
SET(TESTNAME_CreateMeshAndRun ${TESTNAME}_CreateMeshAndRun)
SET(errorRegex "Test Run Failed")
SET(warningRegex "No test is available")
SET(errorRegexBarracuda "ERROR")
SET(warningRegexBarracuda "WARNING")
SET(inputnotfound "The input data file is not found")

SET(TEST_RESULTS_DIR ${CMAKE_CURRENT_BINARY_DIR}/TestResults)
SET(QUALITY_RESULTS_FILE_CreateMeshAndRun ${TEST_RESULTS_DIR}/Barracuda/results/0/Temperature_proc0_000032.qa)
SET(QUALITY_REFERENCE_RESULTS_FILE_CreateMeshAndRun ${TEST_RESULTS_DIR}/Barracuda/reference/Ref.qa)

FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/reference DESTINATION ${TEST_RESULTS_DIR}/Barracuda )
# validate testcase ---------------------------------------------------
ADD_TEST(NAME ${TESTNAME}_Copy_CR
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${iAM_SB_INTETEST_TEMP}/${TESTNAME} ${TEST_RESULTS_DIR}
 )

# Add CanonicalThmeralTester_LoadMeshAndRun_Barracuda Test
ADD_TEST(NAME ${TESTNAME_CreateMeshAndRun}_Barracuda
  COMMAND ${BARRACUDA_EXE} -i ${TEST_RESULTS_DIR}/Barracuda/CreateMeshAndRun.yml -l CreateMeshAndRun.log WORKING_DIRECTORY ${TEST_RESULTS_DIR}/Barracuda
 )
SET_TESTS_PROPERTIES(${TESTNAME_CreateMeshAndRun}_Barracuda PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME_CreateMeshAndRun}_Barracuda PROPERTIES  FAIL_REGULAR_EXPRESSION "${warningRegexBarracuda};${errorRegexBarracuda};${inputnotfound}")

# Add CanonicalThmeralTester_LoadMeshAndRun_Barracuda_QA Test
ADD_TEST(NAME ${TESTNAME_CreateMeshAndRun}_Barracuda_QA
    COMMAND $<TARGET_FILE:${VALIDATE_EXE}> -v -t 1.e-4 --l2_norm  -f ${QUALITY_RESULTS_FILE_LoadMeshAndRun}  -r ${QUALITY_REFERENCE_RESULTS_FILE_LoadMeshAndRun} WORKING_DIRECTORY  ${TEST_RESULTS_DIR}/Barracuda
   )