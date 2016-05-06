/**
 * @file /mtk-trunk/unit_test/vectview_test.cc
 * @brief Brief description
 * 
 */


#include "test_base.hpp"
typedef MTK::vectview<scalar, Eigen::Dynamic> vv_dyn;
typedef MTK::vect<Eigen::Dynamic,scalar> vect_dyn;

BOOST_AUTO_TEST_SUITE( mtk_test_suite )

BOOST_AUTO_TEST_CASE( dynamic_IO )
{
	
	for(int i=1; i<100; ++i) {
		vect_dyn::base ran = vect_dyn::Random(i);
		vect_dyn  v1(ran), v2;
		std::ostringstream out;
		out.precision(16);
		out << v1.format(MTK::IO_no_spaces);
		
		std::istringstream in(out.str());
		in >> v2;
		BOOST_CHECK(in.good());
		
		BOOST_CHECK_EQUAL(i, v1.size());
		BOOST_CHECK_EQUAL(i, v2.size());
		
		
		MTK_CHECK_CLOSE(ran, v1, 1e-9);
		MTK_CHECK_CLOSE(ran, v2, 1e-9);
	}
	
}




BOOST_AUTO_TEST_SUITE_END()

