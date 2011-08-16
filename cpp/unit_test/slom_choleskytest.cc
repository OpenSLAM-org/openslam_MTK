/**
 * @file slom_choleskytest.cc
 * @brief Test cases for SLOM::CholeskyCovariance
 */

#define BOOST_TEST_MODULE template_for_test_test
#include <boost/test/included/unit_test.hpp>

#include "boost_check_close_vector.hh"

#include "random_vector.hh"

#include <iostream>

#include <slom/CholeskyCovariance.hpp>

#include <Eigen/Cholesky>

BOOST_AUTO_TEST_SUITE( slom_test_suite )


// FIXME currently only SCALAR = double is supported
typedef SCALAR scalar;

random_vector<scalar> r;


template<int dim>
void testCholesky(const Eigen::Matrix<double, dim, dim> &sigma)
{
	typedef Eigen::Matrix<double, dim, dim> m_type;
	typedef Eigen::Matrix<double, dim, 1> c_type;
	typedef SLOM::CholeskyCovariance<dim> Cov;
	m_type L = sigma.llt().matrixL();
	
	Cov c1(sigma.data(), SLOM::CholeskyMode::CHOLESKY_FULL);
	
	for(int i=0; i<dim; ++i)
	{
		c_type col = L.col(i);
		c1.invApply(col);
		BOOST_CHECK_CLOSE_FRACTION(col, c_type::Unit(i).eval(), 1e-9);
	}
	
	for(int i=0; i<10; ++i)
	{
		c_type x = r.get_normal<dim>();
		c_type y = x;
		c1.invApply(x);
		c1.apply(x);
		BOOST_CHECK_CLOSE_FRACTION(x, y, 1e-9);
		c1.apply(y);
		c1.invApply(y);
		BOOST_CHECK_CLOSE_FRACTION(x, y, 1e-9);
	}
	
	// now check what happens if L was already decomposed:
#if MTK_EIGEN < 300
	L = sigma.template part<Eigen::LowerTriangular>();
#else
	L = sigma.template triangularView<Eigen::Lower>();
#endif
	Cov c2(L.data(), SLOM::CholeskyMode::COPY_UPPER_FULL);
	for(int i=0; i<dim; ++i)
	{
		c_type col = L.col(i);
		c2.invApply(col);
		BOOST_CHECK_CLOSE_FRACTION(col, c_type::Unit(i).eval(), 1e-9);
	}
}

template<int dim>
Eigen::Matrix<double, dim, dim> generateCovariance(){
	Eigen::Matrix<double, dim, dim> S;
	S.setIdentity();
	for(int i=0; i<10 * dim; ++i)
	{
		Eigen::Matrix<double, dim, 1> v = r.get_normal<dim>();
		S += v * v.transpose();
	}
	
	return S;
}


BOOST_AUTO_TEST_CASE( Cholesky )
{
	testCholesky(generateCovariance<1>());
	testCholesky(generateCovariance<2>());
	testCholesky(generateCovariance<3>());
	testCholesky(generateCovariance<4>());
	testCholesky(generateCovariance<5>());
	testCholesky(generateCovariance<10>());
	testCholesky(generateCovariance<31>());
	
//	BOOST_CHECK(false);

}

BOOST_AUTO_TEST_SUITE_END()


