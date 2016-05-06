/**
 * @file unit_test/slom_covariance.cc
 * @brief Brief description
 * 
 */




#include "test_base.hpp"

#include <slom/CholeskyCovariance.hpp>
#include <Eigen/Cholesky>
#include <Eigen/LU>

#include <Estimator.hpp>
#include <BuildMeasurement.hpp>


typedef MTK::vect<4> vec4;
typedef SLOM::VarID<vec4> vec4_id;
SLOM_BUILD_MEASUREMENT(Zero_Meas, 4, ((vec4, v)), )
SLOM_IMPLEMENT_MEASUREMENT(Zero_Meas, ret) {
	ret = *v;
}


// Implement previously required manual covariance normalization:

#define SLOM_ZERO_MEASUREMENT_WITH_STD(name, std_type, std_op) \
	SLOM_BUILD_MEASUREMENT(    Zero_with_ ## name, 4, ((vec4, v)), ((std_type, std))) \
	SLOM_IMPLEMENT_MEASUREMENT(Zero_with_ ## name, ret) { \
		ret = *v; \
		std_op; \
	}



SLOM_ZERO_MEASUREMENT_WITH_STD(std,     double, ret /= std)
SLOM_ZERO_MEASUREMENT_WITH_STD(inv_std, double, ret *= std)

SLOM_ZERO_MEASUREMENT_WITH_STD(std_vec,     vec4, ret.array() /= std.array())
SLOM_ZERO_MEASUREMENT_WITH_STD(inv_std_vec, vec4, ret.array() *= std.array())

SLOM_ZERO_MEASUREMENT_WITH_STD(std_mat,     SLOM::CholeskyCovariance<4>, std.invApply(ret))
// FIXME this does not work, because chol(X)^-1 != chol(X^-1)
SLOM_ZERO_MEASUREMENT_WITH_STD(inv_std_mat, SLOM::CholeskyCovariance<4>, std.apply(ret))




BOOST_AUTO_TEST_SUITE( slom_test_suite )




struct CovarianceTester {
	SLOM::Estimator e;
	vec4_id v;
	
	CovarianceTester() {
		v = e.insertRV(vec4());
	}
	
	void compare_residuum(){
		for(int i=0; i<100; ++i) {
			e.reinitRV(v, r.get_normal<4>());
			const Eigen::VectorXd & res = e.getRes();
			Eigen::Vector4d comp = res.head<4>();
			for(int k=4; k<res.rows(); k+=4) {
				BOOST_CHECK_CLOSE_FRACTION(comp, res.segment<4>(k), tol);
			}
		}
	}
	
	void insert_matrix_covariances(const Eigen::Matrix4d& cov) {
		Eigen::LLT<Eigen::Matrix4d> chol(cov);
//		Eigen::Matrix4d cov_inv = chol.solve(Eigen::Matrix4d::Identity());
		Eigen::Matrix4d cov_inv = cov.inverse();
		Eigen::Matrix4d dev = chol.matrixL();
		Eigen::Matrix4d dev_inv = chol.matrixL().solve(Eigen::Matrix4d::Identity());
		
		e.insertMeasurement(Zero_Meas(v), SLOM::Covariance(cov));
		e.insertMeasurement(Zero_Meas(v), SLOM::InvCovariance(cov_inv));
		e.insertMeasurement(Zero_Meas(v), SLOM::StandardDeviation(dev));
		e.insertMeasurement(Zero_Meas(v), SLOM::InvStandardDeviation(dev_inv));
		
		e.insertMeasurement(Zero_with_std_mat(    v, SLOM::CholeskyCovariance<4>(cov.data(), SLOM::CholeskyMode::CHOLESKY_FULL)));
		// FIXME this does not work (see above)
//		e.insertMeasurement(Zero_with_inv_std_mat(v, SLOM::CholeskyCovariance<4>(cov_inv.data(), SLOM::CholeskyMode::CHOLESKY_FULL)));
	}
	
	
	void insert_vector_covariances(const Eigen::Vector4d& cov) {
		Eigen::Vector4d dev = cov.cwiseSqrt();
		e.insertMeasurement(Zero_Meas(v), SLOM::Covariance(cov));
		e.insertMeasurement(Zero_Meas(v), SLOM::InvCovariance(cov.cwiseInverse()));
		e.insertMeasurement(Zero_Meas(v), SLOM::StandardDeviation(dev));
		e.insertMeasurement(Zero_Meas(v), SLOM::InvStandardDeviation(dev.cwiseInverse()));
		
		e.insertMeasurement(Zero_with_std_vec(v,dev));
		e.insertMeasurement(Zero_with_inv_std_vec(v, dev.cwiseInverse()));
		
		insert_matrix_covariances(cov.asDiagonal());
	}
	
	void insert_scalar_covariances(const double &cov) {
		double dev = sqrt(cov);
		e.insertMeasurement(Zero_Meas(v), SLOM::Covariance(cov));
		e.insertMeasurement(Zero_Meas(v), SLOM::InvCovariance(1./cov));
		e.insertMeasurement(Zero_Meas(v), SLOM::StandardDeviation(dev));
		e.insertMeasurement(Zero_Meas(v), SLOM::InvStandardDeviation(1./dev));
		
		e.insertMeasurement(Zero_with_std(v, dev));
		e.insertMeasurement(Zero_with_inv_std(v, 1./dev));
		
		
		insert_vector_covariances(Eigen::Vector4d::Constant(cov));
	}
	
};








BOOST_AUTO_TEST_CASE( ScalarCovariance )
{
	using std::exp;
	using std::sqrt;
	for(int i=0; i<100; ++i) {
		double cov = exp(r.uniform());
		CovarianceTester ct;
		ct.insert_scalar_covariances(cov);
		ct.compare_residuum();
	}
}


BOOST_AUTO_TEST_CASE( VectorCovariance )
{
	using std::exp;
	using std::sqrt;
	for(int i=0; i<100; ++i) {
		vec4 cov = Eigen::exp(r.get_uniform<4>().array());
		CovarianceTester ct;
		ct.insert_vector_covariances(cov);
		ct.compare_residuum();
	}
}

BOOST_AUTO_TEST_CASE( MatrixCovariance )
{
	using std::exp;
	using std::sqrt;
	for(int i=0; i<100; ++i) {
		Eigen::Matrix4d cov = generateCovariance<4>();
		CovarianceTester ct;
		ct.insert_matrix_covariances(cov);
		ct.compare_residuum();
	}
}

BOOST_AUTO_TEST_SUITE_END()
