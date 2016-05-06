/**
 * @file /mtk-trunk/unit_test/vectview_test.cc
 * @brief Brief description
 * 
 */


#include "test_base.hpp"
typedef MTK::vectview<scalar, Eigen::Dynamic> vv_dyn;
typedef vv_dyn::matrix_type vect_dyn;

BOOST_AUTO_TEST_SUITE( mtk_test_suite )

BOOST_AUTO_TEST_CASE( dynamic_vectview )
{
	for(int i=1; i<100; ++i) {
		vect_dyn ran = vect_dyn::Random(i);
		vv_dyn  v1(ran), v2(ran.data(), ran.size());
		BOOST_CHECK_EQUAL(i, v1.size());
		BOOST_CHECK_EQUAL(i, v2.size());
		
		MTK_CHECK_CLOSE(ran, v1, 1e-9);
		MTK_CHECK_CLOSE(ran, v2, 1e-9);
	}
	
}




BOOST_AUTO_TEST_SUITE_END()

