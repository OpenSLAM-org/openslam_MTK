/**
 * @file /mtk-trunk/unit_test/test_base.hpp
 * @brief Brief description
 * 
 */


#ifndef TEST_BASE_HPP_
#define TEST_BASE_HPP_


#define BOOST_TEST_MODULE template_for_test_test
#include <boost/test/unit_test.hpp>

#include "boost_check_close_vector.hh"

#include "random_vector.hh"

#include <iostream>

#include <mtk/src/mtkmath.hpp>


#ifndef SCALAR
#define SCALAR double
#endif

typedef SCALAR scalar;

random_vector<scalar> r;
const scalar tol = MTK::tolerance<scalar>();
const scalar one = 1.0;


template<int dim>
Eigen::Matrix<scalar, dim, dim> generateCovariance(int count = 10*dim){
	Eigen::Matrix<scalar, dim, dim> S;
	S.setIdentity();
	for(int i=0; i<count; ++i)
	{
		Eigen::Matrix<scalar, dim, 1> v = r.get_normal<dim>();
		S += v * v.transpose();
	}
	
	return S;
}





#endif /* TEST_BASE_HPP_ */
