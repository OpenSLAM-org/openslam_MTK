/**
 * @file slom_choleskytest.cc
 * @brief Test cases for SLOM::CholeskyCovariance
 */

#include "test_base.hpp"

#include <slom/CholeskyCovariance.hpp>

#include <Eigen/Cholesky>

BOOST_AUTO_TEST_SUITE( slom_test_suite )




template<int dim>
void testCholesky(const Eigen::Matrix<scalar, dim, dim> &sigma)
{
	typedef Eigen::Matrix<scalar, dim, dim> m_type;
	typedef Eigen::Matrix<scalar, dim, 1> c_type;
	typedef SLOM::CholeskyCovariance<dim, scalar> Cov;
	Eigen::LLT<m_type> chol(sigma);
	m_type inv_cov = chol.solve(m_type::Identity());
	m_type L = sigma.llt().matrixL();
	
	Cov c1(sigma.data(), SLOM::CholeskyMode::CHOLESKY_FULL);
	Cov c1_inv = c1; c1_inv.invert();
	// FIXME calculate Cholesky factor of inverse
//	Cov c2_inv(inv_cov.data(), SLOM::CholeskyMode::CHOLESKY_FULL);
//	Cov c2 = c2_inv; c2.invert();
	for(int i=0; i<dim; ++i)
	{
		c_type col = L.col(i);
		c1.invApply(col);
		BOOST_CHECK_CLOSE_FRACTION(col, c_type::Unit(i).eval(), tol);
		col = L.col(i);
		c1_inv.apply(col);
		BOOST_CHECK_CLOSE_FRACTION(col, c_type::Unit(i).eval(), tol);
		col = L.col(i);
//		c2.invApply(col);
//		BOOST_CHECK_CLOSE_FRACTION(col, c_type::Unit(i).eval(), tol);
//		col = L.col(i);
//		c2_inv.apply(col);
//		BOOST_CHECK_CLOSE_FRACTION(col, c_type::Unit(i).eval(), tol);
	}
	
	for(int i=0; i<10; ++i)
	{
		c_type x = r.get_normal<dim>();
		c_type y = x;
		c1.invApply(x);
		c1.apply(x);
		BOOST_CHECK_CLOSE_FRACTION(x, y, tol);
		c1.apply(y);
		c1.invApply(y);
		BOOST_CHECK_CLOSE_FRACTION(x, y, tol);
		c1_inv.apply(x);
		c1.invApply(y);
		BOOST_CHECK_CLOSE_FRACTION(x, y, tol);
		c1_inv.invApply(y);
		c1.apply(x);
		BOOST_CHECK_CLOSE_FRACTION(x, y, tol);
		c1_inv.invApply(y);
		c1.invApply(y);
		BOOST_CHECK_CLOSE_FRACTION(x, y, tol);
		c1_inv.apply(x);
		c1.apply(x);
		BOOST_CHECK_CLOSE_FRACTION(x, y, tol);
	}
	
	// now check what happens if L was already decomposed:
#if MTK_EIGEN < 300
	L = sigma.template part<Eigen::LowerTriangular>();
#else
	L = sigma.template triangularView<Eigen::Lower>();
#endif
	Cov c3(L.data(), SLOM::CholeskyMode::COPY_UPPER_FULL);
	for(int i=0; i<dim; ++i)
	{
		c_type col = L.col(i);
		c3.invApply(col);
		BOOST_CHECK_CLOSE_FRACTION(col, c_type::Unit(i).eval(), tol);
	}
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
	testCholesky(generateCovariance<32>());
	
//	BOOST_CHECK(false);
}

BOOST_AUTO_TEST_SUITE_END()


