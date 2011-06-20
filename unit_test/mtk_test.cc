#define BOOST_TEST_MODULE template_for_test_test
#include <boost/test/included/unit_test.hpp>

#include "boost_check_close_vector.hh"

#include "random_vector.hh"

#include <iostream>

#include <mtk/types/SOn.hpp>
#include <mtk/types/S2.hpp>
#include <mtk/types/pose.hpp>

#include "mtk_axioms.hh"

typedef SCALAR scalar;

BOOST_AUTO_TEST_SUITE( mtk_test_suite )


random_vector<scalar> r;
const scalar tol = MTK::tolerance<scalar>();
const scalar one = 1.0;


BOOST_AUTO_TEST_CASE( SO3 )
{
	typedef MTK::SO3<scalar> SO3;
	for(int i=0; i<16; ++i)
	{
		Eigen::Matrix<scalar, 3, 1> rand;
		do{
			rand = r.get_normal<3>();
		} while(rand.norm() < 1 || rand.norm() > 3);
		SO3 x, y = SO3::exp(rand);
		SO3 z = y.inverse();
		SO3 minus_one(-2*one, 0, 0, 0);
		MTK_CHECK_CLOSE(x.norm(), one, tol);
		MTK_CHECK_CLOSE(y.norm(), one, tol);
		MTK_CHECK_CLOSE(z.norm(), one, tol);
		
		check_equal(x, minus_one);
		MTK_CHECK_CLOSE(y * rand, (y*minus_one) * rand, tol);
		std::cerr << (y * minus_one).norm() << ", " << (y * minus_one).coeffs().transpose() << "; ";
		std::cerr  << (y * rand).transpose() << " "<< ((y*minus_one) * rand).transpose() << " " << ((minus_one * y) * rand).transpose() << " " << std::endl;
		check_equal(y, SO3(minus_one * y));
		check_equal(z, SO3(minus_one * z));
		
		check_axioms(x, y, z, rand, r.get_normal<3>(), MTK::pi);
		
		z = SO3::exp(r.get_normal<3>());
		check_group(x, y, z);
	}
}

BOOST_AUTO_TEST_CASE( SO2 )
{
	typedef MTK::SO2<scalar> SO2;
	for(int i=0; i<16; ++i)
	{
		Eigen::Matrix<scalar, 1, 1> rand;
		do{
			rand = r.get_normal<1>();
		} while(rand.norm() < 1 || rand.norm() > 3);
		SO2 x, y = SO2(rand(0));
		SO2 z = y.inverse();
		SO2 two_pi(2*MTK::pi);
		check_equal(x, two_pi);
		check_equal(y, SO2(two_pi * y));
		check_equal(z, SO2(two_pi * z));
		
		check_axioms(x, y, z, rand, r.get_normal<1>(), MTK::pi);
		
		z = SO2(r.get_normal<1>()(0));
		check_group(x, y, z);
	}
}

BOOST_AUTO_TEST_CASE( SE2 )
{
	typedef MTK::SO2<scalar> SO2;
	typedef MTK::trafo<SO2> SE2;
	for(int i=0; i<16; ++i)
	{
		SE2 x, y(r.get_normal<1>()(0),  r.get_normal<2>());
		SE2 z(r.get_normal<1>()(0), r.get_normal<2>());
		check_group(x, y, z);
	}
}

BOOST_AUTO_TEST_CASE( SE3 )
{
	typedef MTK::SO3<scalar> SO3;
	typedef MTK::trafo<SO3> SE3;
	for(int i=0; i<16; ++i)
	{
		SE3 x, y(SO3::exp(r.get_normal<3>()),  r.get_normal<3>());
		SE3 z(SO3::exp(r.get_normal<3>()),  r.get_normal<3>());
		check_group(x, y, z);
	}
}


BOOST_AUTO_TEST_CASE( vect3 )
{
	typedef MTK::vect<3, scalar> vect3;
	for(int i=0; i<16; ++i)
	{
		Eigen::Matrix<scalar, 3, 1> rand;
		do{
			rand = r.get_normal<3>();
		} while(rand.norm() < 1 );
		vect3 x, y = rand;
		vect3 z = -y;
		
		check_axioms(x, y, z, rand, r.get_normal<3>(), -1);
	}
}


BOOST_AUTO_TEST_CASE( S2 )
{
	typedef MTK::S2<scalar> S2;
	for(int i=0; i<16; ++i)
	{
		Eigen::Matrix<scalar, 2, 1> rand;
		do{
			rand = r.get_normal<2>();
		} while(rand.norm() < 1 || rand.norm() > 3);
		S2 x, y; y.boxplus(rand);
		S2 z(-2.0, 0,0);
		MTK_CHECK_CLOSE(x.get_vect().norm(), one, tol);
		MTK_CHECK_CLOSE(y.get_vect().norm(), one, tol);
		MTK_CHECK_CLOSE(z.get_vect().norm(), one, tol);
		
		check_axioms(x, y, z, rand, r.get_normal<2>(), scalar(MTK::pi));
		
	}
}



BOOST_AUTO_TEST_SUITE_END()
