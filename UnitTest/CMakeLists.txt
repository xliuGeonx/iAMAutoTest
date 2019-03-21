# Barracuda - Copyright (C) 2015 geonX
#
# Please report all bugs and problems to <support@geonx.com>.
#
# UnitTest test case CMakeLists

# -----------------------------------------------------------------------------
# ----------------------------- Project configuration -------------------------
# -----------------------------------------------------------------------------


SET(TESTNAME iAM_UnitTest)
SET(errorRegex "Test Run Failed")
# validate testcase ---------------------------------------------------

# Add KernelType Test
ADD_TEST(NAME ${TESTNAME}_KernelType
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=KernelType /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_KernelType PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_KernelType PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")


# Add NumberofThread Test
ADD_TEST(NAME ${TESTNAME}_NumberofThread
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=NumberofThread /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_NumberofThread PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_NumberofThread PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")
 
 
# Add VoxelizationMaximumSubdomainsLevel Test
ADD_TEST(NAME ${TESTNAME}_VoxelizationMaximumSubdomainsLevel
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=VoxelizationMaximumSubdomainsLevel /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_VoxelizationMaximumSubdomainsLevel PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_VoxelizationMaximumSubdomainsLevel PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add Octrees Test
ADD_TEST(NAME ${TESTNAME}_Octrees
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=Octrees /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_Octrees PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_Octrees PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add Trajectories Test
ADD_TEST(NAME ${TESTNAME}_Trajectories
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=Trajectories /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_Trajectories PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_Trajectories PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add TrajectoryDomains Test
ADD_TEST(NAME ${TESTNAME}_TrajectoryDomains
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=TrajectoryDomains /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_TrajectoryDomains PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_TrajectoryDomains PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add Skins Test
ADD_TEST(NAME ${TESTNAME}_Skins
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=Skins /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_Skins PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_Skins PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add ComputeVolumeDefinitions Test
ADD_TEST(NAME ${TESTNAME}_ComputeVolumeDefinitions
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=ComputeVolumeDefinitions /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_ComputeVolumeDefinitions PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_ComputeVolumeDefinitions PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add Materials Test
ADD_TEST(NAME ${TESTNAME}_Materials
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=Materials /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_Materials PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_Materials PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add SkinAssociations Test
ADD_TEST(NAME ${TESTNAME}_SkinAssociations
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=SkinAssociations /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_SkinAssociations PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_SkinAssociations PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add Layers Test
ADD_TEST(NAME ${TESTNAME}_Layers
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=Layers /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_Layers PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_Layers PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add Mappings Test
ADD_TEST(NAME ${TESTNAME}_Mappings
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=Mappings /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_Mappings PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_Mappings PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add ParametersSets Test
ADD_TEST(NAME ${TESTNAME}_ParametersSets
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=ParametersSets /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_ParametersSets PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_ParametersSets PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add MacroLayersThermalModelingParametersSets Test
ADD_TEST(NAME ${TESTNAME}_MacroLayersThermalModelingParametersSets
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=MacroLayersThermalModelingParametersSets /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_MacroLayersThermalModelingParametersSets PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_MacroLayersThermalModelingParametersSets PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add MacroLayers Test
ADD_TEST(NAME ${TESTNAME}_MacroLayers
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=MacroLayers /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_MacroLayers PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_MacroLayers PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add SourceEvaluators Test
ADD_TEST(NAME ${TESTNAME}_SourceEvaluators
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=SourceEvaluators /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_SourceEvaluators PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_SourceEvaluators PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add SourceTerms Test
ADD_TEST(NAME ${TESTNAME}_SourceTerms
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=SourceTerms /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_SourceTerms PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_SourceTerms PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add FieldProjectionParameters Test
ADD_TEST(NAME ${TESTNAME}_FieldProjectionParameters
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=FieldProjectionParameters /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_FieldProjectionParameters PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_FieldProjectionParameters PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add Points Test
ADD_TEST(NAME ${TESTNAME}_Points
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=Points /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_Points PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_Points PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")


# Add Domains Test
ADD_TEST(NAME ${TESTNAME}_Domains
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=Domains /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_Domains PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_Domains PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")


# Add Bodies Test
ADD_TEST(NAME ${TESTNAME}_Bodies
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=Bodies /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_Bodies PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_Bodies PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add Keys Test
ADD_TEST(NAME ${TESTNAME}_Keys
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=Keys /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_Keys PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_Keys PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")


# Add NonLinearSolverParameters Test
ADD_TEST(NAME ${TESTNAME}_NonLinearSolverParameters
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=NonLinearSolverParameters /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_NonLinearSolverParameters PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_NonLinearSolverParameters PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")


# Add LinearSolverParameters Test
ADD_TEST(NAME ${TESTNAME}_LinearSolverParameters
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=LinearSolverParameters /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_LinearSolverParameters PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_LinearSolverParameters PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add OctreeAdaptationParameters Test
ADD_TEST(NAME ${TESTNAME}_OctreeAdaptationParameters
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=OctreeAdaptationParameters /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_OctreeAdaptationParameters PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_OctreeAdaptationParameters PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")

# Add Fields Test
ADD_TEST(NAME ${TESTNAME}_Fields
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=Fields /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_Fields PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_Fields PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")


# Add BoundaryConditions Test
ADD_TEST(NAME ${TESTNAME}_BoundaryConditions
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=BoundaryConditions /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_BoundaryConditions PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_BoundaryConditions PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")


# Add BoundaryConditionsSets Test
ADD_TEST(NAME ${TESTNAME}_BoundaryConditionsSets
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=BoundaryConditionsSets /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_BoundaryConditionsSets PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_BoundaryConditionsSets PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")


# Add LayerMonitorings Test
ADD_TEST(NAME ${TESTNAME}_LayerMonitorings
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=LayerMonitorings /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_LayerMonitorings PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_LayerMonitorings PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")


# Add TimeStepperEndCriterias Test
ADD_TEST(NAME ${TESTNAME}_TimeStepperEndCriterias
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=TimeStepperEndCriterias /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_TimeStepperEndCriterias PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_TimeStepperEndCriterias PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")


# Add TimeStepperDeltaCriterias Test
ADD_TEST(NAME ${TESTNAME}_TimeStepperDeltaCriterias
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=TimeStepperDeltaCriterias /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_TimeStepperDeltaCriterias PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_TimeStepperDeltaCriterias PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")


# Add TimeSteppers Test
ADD_TEST(NAME ${TESTNAME}_TimeSteppers
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=TimeSteppers /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_TimeSteppers PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_TimeSteppers PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")


# Add CreateSteadyComputationCreateMesh Test
ADD_TEST(NAME ${TESTNAME}_CreateSteadyComputationCreateMesh
  COMMAND ${VSTEST_EXE} ${iAM_SB_UNITTEST_DLL} /TestCaseFilter:TestCategory=CreateSteadyComputationCreateMesh /Inisolation
 )
 
SET_TESTS_PROPERTIES(${TESTNAME}_CreateSteadyComputationCreateMesh PROPERTIES  ENVIRONMENT "${CTEST_ENV_VARS}")
SET_TESTS_PROPERTIES(${TESTNAME}_CreateSteadyComputationCreateMesh PROPERTIES  FAIL_REGULAR_EXPRESSION "${errorRegex}")