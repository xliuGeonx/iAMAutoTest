
# Try to find iAM UnitTest And IntegrationTest
#

unset(iAM_FOUND CACHE)
unset(iAM_SB_UNITTEST_DLL CACHE)
unset(iAM_SB_INTETEST_DLL CACHE)
unset(iAM_INTETEST_DLL CACHE)
unset(iAM_SB_INTETEST_TEMP CACHE)

SET(iAM_FOUND "NO")

# IF(NOT DEFINED BARRACUDA_DIR)
	# MESSAGE(FATAL_ERROR "Barracuda Package directory is missing, please define BARRACUDA_DIR")
# ENDIF


IF(iAM_DIR)
	IF(NOT USE_iAM_RELEASE)
		SET(iAM_COMPLIE_TAG Debug)
	ELSE()
		SET(iAM_COMPLIE_TAG Release)
	ENDIF()
	SET(iAM_SB_UNITTEST_NAME BarracudaSolverBridgeUnitTest)
	SET(iAM_SB_INTETEST_NAME BarracudaSolverBridgeInteTest)
	SET(iAM_INTETEST_NAME IntegrationTests)
	SET(iAM_SB_UNITTEST_BIN_DIR ${iAM_DIR}/${iAM_SB_UNITTEST_NAME}/bin/${iAM_COMPLIE_TAG})	
	SET(iAM_SB_INTETEST_BIN_DIR ${iAM_DIR}/${iAM_SB_INTETEST_NAME}/bin/${iAM_COMPLIE_TAG})
	SET(iAM_INTETEST_BIN_DIR ${iAM_DIR}/${iAM_INTETEST_NAME}/bin/${iAM_COMPLIE_TAG})
	SET(iAM_SB_INTETEST_TEMP ${iAM_SB_INTETEST_BIN_DIR}/Temp)
			
	#Find iAM SlvBridge UnitTest dll
	FIND_FILE(iAM_SB_UNITTEST_DLL NAME "${iAM_SB_UNITTEST_NAME}.dll" PATHS ${iAM_SB_UNITTEST_BIN_DIR} NO_DEFAULT_PATH)
	IF(NOT iAM_SB_UNITTEST_DLL)
		MESSAGE(WARNING "${iAM_SB_UNITTEST_NAME}.dll is not found in ${iAM_SB_UNITTEST_BIN_DIR}")
		SET(iAM_FOUND "NO")
	ELSE()
		MESSAGE(STATUS "${iAM_SB_UNITTEST_NAME}.dll is found: ${iAM_SB_UNITTEST_DLL}")
		SET(iAM_FOUND "YES")
	ENDIF()
	
	#Find iAM SlvBridge InteTest dll
	FIND_FILE(iAM_SB_INTETEST_DLL NAME "${iAM_SB_INTETEST_NAME}.dll" PATHS ${iAM_SB_INTETEST_BIN_DIR} NO_DEFAULT_PATH)
	IF(NOT iAM_SB_INTETEST_DLL)
		MESSAGE(WARNING "${iAM_SB_INTETEST_NAME}.dll is not found")
		SET(iAM_FOUND "NO")
	ELSE()
		MESSAGE(STATUS "iAM_SB_INTETEST_DLL is found: ${iAM_SB_INTETEST_DLL}")
		SET(iAM_FOUND "YES")
	ENDIF()
	#Find Runsettings file
	FIND_FILE(RUNSETTINGS_FILE NAME "TestConfigure.runsettings" PATHS ${iAM_SB_INTETEST_BIN_DIR}/DefaultAsset NO_DEFAULT_PATH)
	IF(NOT RUNSETTINGS_FILE)
		MESSAGE(WARNING "RUNSETTINGS_FILE is not found")
		SET(iAM_FOUND "NO")
	ELSE()
		MESSAGE(STATUS "RUNSETTINGS_FILE is found: ${RUNSETTINGS_FILE}")
		SET(iAM_FOUND "YES")
	ENDIF()
	
	#Find iAM InteTest dll
	FIND_FILE(iAM_INTETEST_DLL NAME "${iAM_INTETEST_NAME}.dll" PATHS ${iAM_INTETEST_BIN_DIR} NO_DEFAULT_PATH)
	IF(NOT iAM_INTETEST_DLL)
		MESSAGE(WARNING "${iAM_INTETEST_NAME}.dll is not found")
		SET(iAM_FOUND "NO")
	ELSE()
		MESSAGE(STATUS "iAM_INTETEST_DLL is found: ${iAM_INTETEST_DLL}")
		SET(iAM_FOUND "YES")
	ENDIF()
ENDIF()

IF(${iAM_FOUND} STREQUAL "NO")
	MESSAGE(FATAL_ERROR
		"Please set iAM_DIR to the root location of 
		iAM containing ${iAM_SB_UNITTEST_NAME}.dll and ${iAM_SB_INTETEST_NAME}.dll.\n
		Not found at iAM_DIR: ${iAM_DIR}\n"
	)
ENDIF()

SET(iAM_SB_UNITTEST_DLL ${iAM_SB_UNITTEST_DLL} CACHE PATH "iAM_SB_UNITTEST_DLL")
SET(iAM_SB_INTETEST_DLL ${iAM_SB_INTETEST_DLL} CACHE PATH "iAM_SB_INTETEST_DLL")
SET(iAM_INTETEST_DLL ${iAM_INTETEST_DLL} CACHE PATH "iAM_INTETEST_DLL")
SET(RUNSETTINGS_FILE ${RUNSETTINGS_FILE} CACHE PATH "RUNSETTINGS_FILE")
SET(iAM_SB_INTETEST_TEMP ${iAM_SB_INTETEST_TEMP} CACHE PATH "iAM_SB_INTETEST_TEMP")
mark_as_advanced( iAM_FOUND )

