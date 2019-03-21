# Barracuda - Copyright (C) 2015 geonX
#
# Please report all bugs and problems to <support@geonx.com>.
#
# UnitTest test case CMakeLists

# -----------------------------------------------------------------------------
# ----------------------------- Project configuration -------------------------
# -----------------------------------------------------------------------------


SET(TESTNAME iAM_InteTest)
SET(errorRegex "Test Run Failed")
# validate testcase ---------------------------------------------------

# Add InteTest Test
ADD_TEST(NAME ${TESTNAME}
  COMMAND ${VSTEST_EXE} ${iAM_INTETEST_DLL} /Inisolation /Settings:${RUNSETTINGS_FILE}
 )
 
SET_TESTS_PROPERTIES(${TESTNAME} PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME} PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add ProjectVMCretor Test
ADD_TEST(NAME ${TESTNAME}_ProjectVMCretor
  COMMAND ${VSTEST_EXE} ${iAM_SB_INTETEST_DLL} /TestCaseFilter:TestCategory=ProjectVMCretor /Inisolation /Settings:${RUNSETTINGS_FILE}
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_ProjectVMCretor PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_ProjectVMCretor PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")



# Add AddSTLTester Test
ADD_TEST(NAME ${TESTNAME}_AddSTLTester
  COMMAND ${VSTEST_EXE} ${iAM_SB_INTETEST_DLL} /TestCaseFilter:TestCategory=AddSTLTester /Inisolation /Settings:${RUNSETTINGS_FILE}
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_AddSTLTester PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_AddSTLTester PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add AddMaterialTester Test
ADD_TEST(NAME ${TESTNAME}_AddMaterialTester
  COMMAND ${VSTEST_EXE} ${iAM_SB_INTETEST_DLL} /TestCaseFilter:TestCategory=AddMaterialTester /Inisolation /Settings:${RUNSETTINGS_FILE}
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_AddMaterialTester PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_AddMaterialTester PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")


# Add AddMappingTester Test
ADD_TEST(NAME ${TESTNAME}_AddMappingTester
  COMMAND ${VSTEST_EXE} ${iAM_SB_INTETEST_DLL} /TestCaseFilter:TestCategory=AddMappingTester /Inisolation /Settings:${RUNSETTINGS_FILE}
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_AddMappingTester PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_AddMappingTester PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")


# Add AddLayersTester Test
ADD_TEST(NAME ${TESTNAME}_AddLayersTester
  COMMAND ${VSTEST_EXE} ${iAM_SB_INTETEST_DLL} /TestCaseFilter:TestCategory=AddLayersTester /Inisolation /Settings:${RUNSETTINGS_FILE}
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_AddLayersTester PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_AddLayersTester PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")


# Add AddMacroLayersTester Test
ADD_TEST(NAME ${TESTNAME}_AddMacroLayersTester
  COMMAND ${VSTEST_EXE} ${iAM_SB_INTETEST_DLL} /TestCaseFilter:TestCategory=AddMacroLayersTester /Inisolation /Settings:${RUNSETTINGS_FILE}
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_AddMacroLayersTester PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_AddMacroLayersTester PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")


# Add CallSolverBridgeTester Test
ADD_TEST(NAME ${TESTNAME}_CallSolverBridgeTester
  COMMAND ${VSTEST_EXE} ${iAM_SB_INTETEST_DLL} /TestCaseFilter:TestCategory=CallSolverBridgeTester /Inisolation /Settings:${RUNSETTINGS_FILE}
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_CallSolverBridgeTester PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_CallSolverBridgeTester PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

