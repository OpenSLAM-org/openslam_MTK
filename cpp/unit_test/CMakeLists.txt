include_directories( ${EIGEN3_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR})

INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/slom)
#include_directories( ${CXSPARSE_INCLUDE_DIR} )

ADD_LIBRARY( test_base STATIC test_base.cpp)


MACRO(MAKE_TEST name scalar)
	ADD_TEST(${name}_${scalar} ${name}_${scalar})
	ADD_EXECUTABLE(${name}_${scalar} ${name}.cc)
	SET_TARGET_PROPERTIES(${name}_${scalar} PROPERTIES COMPILE_FLAGS -DSCALAR=${scalar})
	TARGET_LINK_LIBRARIES(${name}_${scalar} test_base ${Boost_LIBRARIES} slom)
ENDMACRO()


MAKE_TEST(vectview_test float)
MAKE_TEST(vectview_test double)

MAKE_TEST(IO_test float)
MAKE_TEST(IO_test double)

MAKE_TEST(mtk_test float)
MAKE_TEST(mtk_test double)
MAKE_TEST(slom_choleskytest float)
MAKE_TEST(slom_choleskytest double)

MAKE_TEST(slom_covariance double)
MAKE_TEST(slom_dynamic_measurement double)
