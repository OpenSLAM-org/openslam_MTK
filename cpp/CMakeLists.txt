CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# If SLoM is imported from somewhere else, don't change the FLAGS:
IF(NOT SLOM_BUILD_LIBONLY)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-function -Wno-attributes -fno-strict-aliasing -march=native -mfpmath=sse")
	SET(SLOM_BUILD_EXAMPLES ON)
ENDIF(NOT SLOM_BUILD_LIBONLY)



######################
# find include dirs:
######################

SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)

# FIND_PACKAGE(CXSparse REQUIRED) 
FIND_PACKAGE(Eigen3 3.2.0 REQUIRED)
FIND_PACKAGE(Boost 1.36.0 REQUIRED COMPONENTS program_options unit_test_framework)

INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR} )
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} )

SET(MTK_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

ADD_LIBRARY( slom STATIC 
			slom/src/Estimator.cpp
			slom/src/SparseFunction.cpp
			slom/src/SparseInterface.cpp
			slom/src/Optimizer.cpp
			slom/solvers/GolubKahanLanczos.cpp
			slom/solvers/CGNRSolver.cpp 
			slom/solvers/LSMRSolver.cpp 
#			slom/solvers/LSMRSolver_simplified.cpp 
			slom/solvers/LSQRSolver.cpp
#			slom/solvers/LSQRSolver_simplified.cpp 
			slom/solvers/LSQRSolverMatlab.cpp
			slom/algorithms/GaussNewton.cpp
			slom/algorithms/LevenbergMarquardt.cpp
			slom/algorithms/PowellsDogLeg.cpp
			slom/tools/ExportSparse.cpp
)

MESSAGE(STATUS ${EIGEN3_VERSION})

IF(${EIGEN3_VERSION} VERSION_LESS 3.1.90 AND NOT SLOM_DISABLE_CHOLMOD)
    MESSAGE("Warning: You need a newer Eigen version for Cholmod support")
    SET(SLOM_DISABLE_CHOLMOD true)
ENDIF(${EIGEN3_VERSION} VERSION_LESS 3.1.90 AND NOT SLOM_DISABLE_CHOLMOD)


IF(NOT SLOM_DISABLE_CHOLMOD)
    # only look for SuiteSparse, if Cholmod is not disabled
    FIND_PACKAGE(SuiteSparse)
ENDIF(NOT SLOM_DISABLE_CHOLMOD)



IF(CHOLMOD_FOUND AND NOT SLOM_DISABLE_CHOLMOD)
    MESSAGE("Found CholMod library")
    TARGET_LINK_LIBRARIES(slom
        ${CHOLMOD_LIBRARIES}
#   If the above line does not work, try this or some variation of it.
#   And send me a FindSuiteSparse that works everywhere ;)
#        ${CHOLMOD_LIBRARY}  metis lapack camd colamd ccolamd amd ${BLAS_LIBRARIES}
    )
    INCLUDE_DIRECTORIES(${CHOLMOD_INCLUDE_DIR})
    add_definitions(-DSLOM_CHOLMOD_SUPPORT) 
ENDIF(CHOLMOD_FOUND AND NOT SLOM_DISABLE_CHOLMOD)


IF(SLOM_BUILD_EXAMPLES)
	ENABLE_TESTING()
#	ADD_EXECUTABLE( min_example doc/example.cpp)
#	TARGET_LINK_LIBRARIES(min_example slom)
	
#	ADD_SUBDIRECTORY (test)
	ADD_SUBDIRECTORY (example)
	ADD_SUBDIRECTORY (unit_test)
ENDIF()

INSTALL(DIRECTORY mtk   DESTINATION include)
INSTALL(DIRECTORY ukfom DESTINATION include)
INSTALL(DIRECTORY slom  DESTINATION include)


install(TARGETS slom
    LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
    )

