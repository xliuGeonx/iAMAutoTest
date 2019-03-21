# Barracuda - Copyright (C) 2015 geonX
#
# Please report all bugs and problems to <support@geonx.com>.
#
# UnitTest test case CMakeLists

# -----------------------------------------------------------------------------
# ----------------------------- Project configuration -------------------------
# -----------------------------------------------------------------------------


SET(TESTNAME CanonicalThermalTester)
SET(TESTNAME_CreateMesh ${TESTNAME}_CreateMesh)
SET(TESTNAME_LoadMeshAndRun ${TESTNAME}_LoadMeshAndRun)
SET(errorRegex "Test Run Failed")
SET(warningRegex "No test is available")
SET(errorRegexBarracuda "ERROR")
SET(warningRegexBarracuda "WARNING")
SET(inputnotfound "The input data file is not found")

SET(TEST_RESULTS_DIR ${CMAKE_CURRENT_BINARY_DIR}/TestResults)
SET(QUALITY_RESULTS_FILE_CreateMesh ${TEST_RESULTS_DIR}/Barracuda/meshes/Mesh_Collection.xmf)
SET(QUALITY_RESULTS_FILE_LoadMeshAndRun ${TEST_RESULTS_DIR}/Barracuda/results/0/Temperature_proc0_000032.qa)
SET(QUALITY_REFERENCE_RESULTS_FILE_LoadMeshAndRun ${TEST_RESULTS_DIR}/Barracuda/reference/Ref.qa)
SET(QUALITY_RESULTS_FILE_LoadMeshAndRun_Postpro ${TEST_RESULTS_DIR}/Barracuda/results/TemperatureOnSkins_After_Heating_00000032.xmf)
IF(EXISTS QUALITY_REFERENCE_RESULTS_FILE_LoadMeshAndRun)
  FILE(REMOVE ${QUALITY_REFERENCE_RESULTS_FILE_LoadMeshAndRun})
ENDIF()
FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/reference DESTINATION ${TEST_RESULTS_DIR}/Barracuda )
# validate testcase ---------------------------------------------------

# Add CanonicalThermalTester Test
ADD_TEST(NAME ${TESTNAME}
  COMMAND ${VSTEST_EXE} ${iAM_SB_INTETEST_DLL} /TestCaseFilter:TestCategory=${TESTNAME} /Inisolation /Settings:${RUNSETTINGS_FILE}
 )
SET_TESTS_PROPERTIES(${TESTNAME} PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME} PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex};${warningRegex}")
ADD_TEST(NAME ${TESTNAME}_Copy
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${iAM_SB_INTETEST_TEMP}/${TESTNAME} ${TEST_RESULTS_DIR}
 )

# Add CanonicalThmeralTester_CreateMesh_Barracuda Test
ADD_TEST(NAME ${TESTNAME_CreateMesh}_Barracuda
  COMMAND ${BARRACUDA_EXE} -i ${TEST_RESULTS_DIR}/Barracuda/CreateMesh.yml -l CreateMesh.log WORKING_DIRECTORY ${TEST_RESULTS_DIR}/Barracuda
 )
SET_TESTS_PROPERTIES(${TESTNAME_CreateMesh}_Barracuda PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME_CreateMesh}_Barracuda PROPERTIES  FAIL_REGULAR_EXPRESSION "${warningRegexBarracuda};${errorRegexBarracuda};${inputnotfound}")


# Add CanonicalThmeralTester_CreateMesh_Barracuda_QA Test
ADD_TEST(NAME ${TESTNAME_CreateMesh}_Barracuda_QA 
    COMMAND $<TARGET_FILE:${VALIDATE_EXE}> --filecheck true -f ${QUALITY_RESULTS_FILE_CreateMesh}
)
SET_TESTS_PROPERTIES(${TESTNAME_CreateMesh}_Barracuda_QA PROPERTIES DEPENDS ${TESTNAME_CreateMesh}_Barracuda)


# Add CanonicalThmeralTester_LoadMeshAndRun_Barracuda Test
ADD_TEST(NAME ${TESTNAME_LoadMeshAndRun}_Barracuda
  COMMAND ${BARRACUDA_EXE} -i ${TEST_RESULTS_DIR}/Barracuda/LoadMeshAndRun.yml -l LoadMeshAndRun.log WORKING_DIRECTORY ${TEST_RESULTS_DIR}/Barracuda
 )
SET_TESTS_PROPERTIES(${TESTNAME_LoadMeshAndRun}_Barracuda PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME_LoadMeshAndRun}_Barracuda PROPERTIES  FAIL_REGULAR_EXPRESSION "${warningRegexBarracuda};${errorRegexBarracuda};${inputnotfound}")

# Add CanonicalThmeralTester_LoadMeshAndRun_Barracuda_QA Test
ADD_TEST(NAME ${TESTNAME_LoadMeshAndRun}_Barracuda_QA
    COMMAND $<TARGET_FILE:${VALIDATE_EXE}> -v -t 1.e-4 --l2_norm  -f ${QUALITY_RESULTS_FILE_LoadMeshAndRun}  -r ${QUALITY_REFERENCE_RESULTS_FILE_LoadMeshAndRun} WORKING_DIRECTORY  ${TEST_RESULTS_DIR}/Barracuda
   )
   
# Add CanonicalThmeralTester_LoadMeshAndRun_Postpro Test
ADD_TEST(NAME ${TESTNAME_LoadMeshAndRun}_Postpro
  COMMAND ${POSTPROCESSING_EXE} -i ${TEST_RESULTS_DIR}/Barracuda/LoadMeshAndRun.yml WORKING_DIRECTORY ${TEST_RESULTS_DIR}/Barracuda
 )
SET_TESTS_PROPERTIES(${TESTNAME_LoadMeshAndRun}_Postpro PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME_LoadMeshAndRun}_Postpro PROPERTIES  FAIL_REGULAR_EXPRESSION "${warningRegexBarracuda};${errorRegexBarracuda};${inputnotfound}")

# Add CanonicalThmeralTester_LoadMeshAndRun_Postpro_QA Test
ADD_TEST(NAME ${TESTNAME_LoadMeshAndRun}_Postpro_QA 
    COMMAND $<TARGET_FILE:${VALIDATE_EXE}> --filecheck true -f ${QUALITY_RESULTS_FILE_LoadMeshAndRun_Postpro}
)
SET_TESTS_PROPERTIES(${TESTNAME_LoadMeshAndRun}_Postpro_QA PROPERTIES DEPENDS ${TESTNAME_LoadMeshAndRun}_Postpro)
