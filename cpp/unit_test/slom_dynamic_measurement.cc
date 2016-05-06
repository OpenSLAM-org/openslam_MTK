/**
 * @file /mtk-trunk/unit_test/slom_dynamic_measurement.cc
 * @brief Brief description
 * 
 */




#include "test_base.hpp"
#include <slom/Estimator.hpp>
#include <slom/BuildMeasurement.hpp>
#include <mtk/types/vect.hpp>
typedef MTK::vectview<scalar, Eigen::Dynamic> vv_dyn;
typedef vv_dyn::matrix_type vect_dyn;

typedef MTK::Scalar<double> Scalar;

BOOST_AUTO_TEST_SUITE( slom_test_suite )

SLOM_BUILD_MEASUREMENT_DYNAMIC(DynTest, ((Scalar, x)), )
SLOM_IMPLEMENT_MEASUREMENT(DynTest, ret) {
	ret.fill(*x);
}


BOOST_AUTO_TEST_CASE( increasing_size )
{
	
	for(int len = 1; len<10; ++len) {
		SLOM::Estimator est;
		
		SLOM::VarID<Scalar> x_id = est.insertRV(Scalar(1.0));
		BOOST_CHECK_EQUAL(1, est.getN());
		
		est.insertMeasurement(DynTest(len, x_id));
		
		BOOST_CHECK_EQUAL(len, est.getM());
		
		MTK_CHECK_CLOSE(est.getRes(), Eigen::VectorXd::Ones(len), 1e-6);
		
		MTK_CHECK_CLOSE(est.getRSS(), len*1.0, 1e-6);
		
		est.optimizeStep();
		
		BOOST_CHECK_SMALL(est.getRSS(),  1e-6);
		BOOST_CHECK_SMALL(double(*x_id), 1e-6);
		}
}

BOOST_AUTO_TEST_CASE( add_remove )
{
	int len1 = 10, len2=4;
	SLOM::Estimator est;
	
	SLOM::VarID<Scalar> x_id[2] = {est.insertRV(Scalar(1.0)), est.insertRV(Scalar(2.0))};
	BOOST_CHECK_EQUAL(2, est.getN());
	
	DynTest::id m1 = est.insertMeasurement(DynTest(len1, x_id[0]));
	
	BOOST_CHECK_EQUAL(len1, est.getM());
	
	MTK_CHECK_CLOSE(est.getRes(), Eigen::VectorXd::Ones(len1), 1e-6);
	
	MTK_CHECK_CLOSE(est.getRSS(), len1*1.0, 1e-6);
	
	DynTest::id m2 = est.insertMeasurement(DynTest(len2, x_id[1])); (void)m2;
	
	BOOST_CHECK_EQUAL(len1 + len2, est.getM());
	
	MTK_CHECK_CLOSE(est.getRes().head(len1).eval(), Eigen::VectorXd::Ones(len1), 1e-6);
	MTK_CHECK_CLOSE(est.getRes().tail(len2).eval(), Eigen::VectorXd::Ones(len2)*2.0, 1e-6);
	
	MTK_CHECK_CLOSE(est.getRSS(), len1*1.0 + len2*4.0, 1e-6);
	
	est.removeMeasurement(m1);
	
	BOOST_CHECK_EQUAL(len2, est.getM());
	MTK_CHECK_CLOSE(est.getRes(), Eigen::VectorXd::Ones(len2)*2.0, 1e-6);
	
	MTK_CHECK_CLOSE(est.getRSS(), len2*4.0, 1e-6);
	
	est.removeRV(x_id[0]);
	est.optimizeStep();
	
	BOOST_CHECK_SMALL(est.getRSS(),  1e-6);
	BOOST_CHECK_SMALL(double(*x_id[1]), 1e-6);
}




BOOST_AUTO_TEST_SUITE_END()


