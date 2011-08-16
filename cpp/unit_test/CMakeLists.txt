include_directories( ${EIGEN2_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR})

INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/slom)
include_directories( ${CXSPARSE_INCLUDE_DIR} )


MACRO(MAKE_TEST name scalar)
	ADD_TEST(${name}_${scalar} ${name}_${scalar})
	ADD_EXECUTABLE(${name}_${scalar} ${name}.cc)
	SET_TARGET_PROPERTIES(${name}_${scalar} PROPERTIES COMPILE_FLAGS -DSCALAR=${scalar})
	TARGET_LINK_LIBRARIES(${name}_${scalar} boost_unit_test_framework-mt)
ENDMACRO()


MAKE_TEST(mtk_test float)
MAKE_TEST(mtk_test double)
MAKE_TEST(slom_choleskytest double)

