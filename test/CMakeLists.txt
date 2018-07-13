#===============================================================================
# Usage:
#   BULLET_ADD_TEST("unit" test_UnitTestA) # assumed source is test_UnitTestA.cpp
#   BULLET_ADD_TEST("unit" test_UnitTestB test_SourceB1.cpp)
#   BULLET_ADD_TEST("unit" test_UnitTestA test_SourceC1.cpp test_SourceC2.cpp)
#===============================================================================
FUNCTION(BULLET_ADD_TEST target_name) # ARGN for source files
  IF(${ARGC} GREATER 2)
    SET(sources ${ARGN})
  ELSE()
    SET(sources "${target_name}.cpp")
  ENDIF()

  ADD_EXECUTABLE(${target_name} ${sources})

  ADD_TEST(${target_name} ${target_name})

  TARGET_LINK_LIBRARIES(${target_name}
      gtest
      LinearMath
      BulletCollision
      BulletDynamics
  )
  IF(NOT WIN32)
    TARGET_LINK_LIBRARIES(${target_name} pthread)
  ENDIF()

  IF(INTERNAL_ADD_POSTFIX_EXECUTABLE_NAMES)
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES DEBUG_POSTFIX "_Debug")
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES MINSIZEREL_POSTFIX "_MinsizeRel")
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES RELWITHDEBINFO_POSTFIX "_RelWithDebugInfo")
  ENDIF(INTERNAL_ADD_POSTFIX_EXECUTABLE_NAMES)
ENDFUNCTION()

INCLUDE_DIRECTORIES(${target_name} BEFORE "${CMAKE_SOURCE_DIR}/src")
INCLUDE_DIRECTORIES(${target_name} "${CMAKE_CURRENT_SOURCE_DIR}/gtest-1.7.0/include")

SUBDIRS(  gtest-1.7.0 collision BulletDynamics )

IF(BUILD_BULLET3)
	#SUBDIRS( TestBullet3OpenCL )
        SUBDIRS(LinearMath InverseDynamics SharedMemory)
ENDIF(BUILD_BULLET3)
