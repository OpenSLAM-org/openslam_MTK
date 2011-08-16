/*
 * random_vector.hh
 *
 *  Created on: 11.05.2010
 *      Author: chtz
 */

#ifndef RANDOM_VECTOR_HH_
#define RANDOM_VECTOR_HH_


#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
//#include <boost/random/uniform_on_sphere.hpp>

#include <mtk/src/eigen.hpp>

/**
 * Generate normal or uniform distributed random vectors of arbitrary dimension
 */

template<class scalar, class RandomGenerator = boost::mt19937>
struct random_vector
{
	random_vector(unsigned int seed) 
		: rng_(RandomGenerator(seed)), 
		  uniform(rng_, boost::uniform_real<scalar>(-1.0, 1.0)),
		  normal(rng_, boost::normal_distribution<scalar>())
	{}
	random_vector(const RandomGenerator &rng = RandomGenerator()) 
		: rng_(rng), 
		  uniform(rng_, boost::uniform_real<scalar>(-1.0, 1.0)), 
		  normal(rng_, boost::normal_distribution<scalar>())
	//, 		  sphere(rng_, boost::uniform_on_sphere<scalar>())
	{}
	
	template <int dim>
	Eigen::Matrix<scalar, dim, 1> get_uniform()
	{
		Eigen::Matrix<scalar, dim, 1> ret;
		for(int i=0; i<dim; ++i){
			ret[i] = uniform();
		}
		return ret;
	}

	template <int dim>
	Eigen::Matrix<scalar, dim, 1> get_normal()
	{
		Eigen::Matrix<scalar, dim, 1> ret;
		for(int i=0; i<dim; ++i){
			ret[i] = normal();
		}
		return ret;
	}
	
	template <int dim>
	Eigen::Matrix<scalar, dim, 1> get_uniform_sphere()
	{
		Eigen::Matrix<scalar, dim, 1> ret = get_normal<dim>();
		ret.normalize();
		return ret;
	}

	
	RandomGenerator rng_;
	boost::variate_generator<RandomGenerator, boost::uniform_real<scalar> > uniform;
	boost::variate_generator<RandomGenerator, boost::normal_distribution<scalar> > normal;
//	boost::variate_generator<RandomGenerator, boost::uniform_on_sphere<scalar> > sphere;
};



#endif /* RANDOM_VECTOR_HH_ */
