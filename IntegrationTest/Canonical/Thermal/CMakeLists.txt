# Barracuda - Copyright (C) 2015 geonX
#
# Please report all bugs and problems to <support@geonx.com>.
#
# UnitTest test case CMakeLists

# -----------------------------------------------------------------------------
# ----------------------------- Project configuration -------------------------
# -----------------------------------------------------------------------------


SET(TESTNAME_CreatMesh CanonicalThmeralTesterCreateMesh)
SET(errorRegex "Test Run Failed")
SET(errorRegexBarracuda "Warning")
SET(warningRegex "WARNING")
SET(TEST_RESULTS_DIR ${CMAKE_CURRENT_BINARY_DIR}/TestResults)
SET(QUALITY_RESULTS_FILE ${TEST_RESULTS_DIR}/meshes/Mesh_Collection.xmf)
IF(EXISTS QUALITY_RESULTS_FILE)
  FILE(REMOVE ${QUALITY_RESULTS_FILE})
ENDIF()
# validate testcase ---------------------------------------------------

# Add CanonicalThmeralTesterCreateMesh Test
ADD_TEST(NAME ${TESTNAME_CreatMesh}
  COMMAND ${VSTEST_EXE} ${iAM_SB_INTETEST_DLL} /TestCaseFilter:TestCategory=${TESTNAME_CreatMesh} /Inisolation /Settings:${RUNSETTINGS_FILE}
 )
SET_TESTS_PROPERTIES(${TESTNAME_CreatMesh} PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME_CreatMesh} PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add CanonicalThmeralTesterCreateMesh_Barracuda Test
ADD_TEST(NAME ${TESTNAME_CreatMesh}_Barracuda
  COMMAND ${BARRACUDA_EXE} -i ${iAM_SB_INTETEST_TEMP}/${TESTNAME_CreatMesh}/Barracuda/CreateMesh.yml -l CreateMesh.log WORKING_DIRECTORY ${TEST_RESULTS_DIR}
 )
SET_TESTS_PROPERTIES(${TESTNAME_CreatMesh}_Barracuda PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME_CreatMesh}_Barracuda PROPERTIES  FAIL_REGULAR_EXPRESSION "${warningRegex};${errorRegexBarracuda}")


# Add CanonicalThmeralTesterCreateMesh_Barracuda_QA Test
ADD_TEST(NAME ${TESTNAME_CreatMesh}_Barracuda_QA 
    COMMAND $<TARGET_FILE:${VALIDATE_EXE}> --filecheck true -f ${QUALITY_RESULTS_FILE}
)
SET_TESTS_PROPERTIES(${TESTNAME_CreatMesh}_Barracuda_QA PROPERTIES DEPENDS ${TESTNAME_CreatMesh}_Barracuda)