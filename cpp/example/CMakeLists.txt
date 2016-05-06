include_directories( ${Eigen3_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR})
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/slom)

include_directories( ${CXSPARSE_INCLUDE_DIR} )

# This needs some clean-up before publication, sorry
# ADD_EXECUTABLE( dlr_data_set dlr/dlr_data_set.cpp)
# TARGET_LINK_LIBRARIES( dlr_data_set slom ${Boost_LIBRARIES})

ADD_EXECUTABLE( sba sba/slom_sba.cc)
TARGET_LINK_LIBRARIES( sba slom)

ADD_EXECUTABLE( relation relation/relation.cpp)
TARGET_LINK_LIBRARIES( relation slom ${Boost_LIBRARIES})

ADD_EXECUTABLE( relation_slameval relation/relation_slameval.cpp)
TARGET_LINK_LIBRARIES( relation_slameval slom)
