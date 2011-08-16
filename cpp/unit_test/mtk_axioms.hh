/**
 * @file mtk_axioms.hh
 * @brief Generic testcases for axioms
 */


#ifndef MTK_AXIOMS_HH_
#define MTK_AXIOMS_HH_

#include "boost_check_close_vector.hh"

#include "random_vector.hh"
#include <mtk/types/vect.hpp>


/**
 * Checks that x == y.
 */
template<class Manifold>
void check_equal(const Manifold &x, const Manifold &y)
{
	typedef typename Manifold::scalar scalar;
	typedef MTK::vect<Manifold::DOF, scalar> vect_type;
	scalar tol = MTK::tolerance<scalar>();
	vect_type delta;
	
	x.boxminus(delta, y);
	BOOST_CHECK_SMALL(delta.norm(), tol);
	y.boxminus(delta, x);
	BOOST_CHECK_SMALL(delta.norm(), tol);
}

/**
 * Checks that x != y.
 */
template<class Manifold>
void check_not_equal(const Manifold &x, const Manifold &y)
{
	typedef typename Manifold::scalar scalar;
	typedef MTK::vect<Manifold::DOF, scalar> vect_type;
	scalar tol = MTK::tolerance<scalar>();
	vect_type delta;
	
	x.boxminus(delta, y);
	BOOST_CHECK_GT(delta.norm(), tol);
	y.boxminus(delta, x);
	BOOST_CHECK_GT(delta.norm(), tol);
}



/**
 * Check x [+] 0 == x and x [-] x == 0.
 */
template<class Manifold>
void check_axiom_zero(const Manifold &x)
{
	typedef typename Manifold::scalar scalar;
	typedef MTK::vect<Manifold::DOF, scalar> vect_type;
	scalar tol = MTK::tolerance<scalar>();
	vect_type delta;
	x.boxminus(delta, x);
	BOOST_CHECK_SMALL(delta.norm(), tol);
	Manifold y = x;
	delta.setZero();
	y.boxplus(delta);
	check_equal(x, y);
}

/**
 * Check x [+] (y [-] x) == y.
 */
template<class Manifold>
void check_axiom_plus_minus(const Manifold &x, const Manifold &y)
{
	typedef typename Manifold::scalar scalar;
	typedef MTK::vect<Manifold::DOF, scalar> vect_type;
	vect_type delta;
	y.boxminus(delta, x);
	Manifold z = x;
	z.boxplus(delta);
	vect_type diff;
	check_equal(z, y);
}

/**
 * Check norm(y [-] x) < maxNorm 
 */
template<class Manifold, class vector>
void check_minus_le_maxnorm(const Manifold &x, const vector &v, const typename Manifold::scalar &maxNorm)
{
	typedef typename Manifold::scalar scalar;
	typedef MTK::vect<Manifold::DOF, scalar> vect_type;
	scalar tol = MTK::tolerance<scalar>();
	vect_type delta = v;
	delta.normalize();
	for(scalar i=0.0; i <= 100 * maxNorm; i+=maxNorm/16)
	{
		Manifold y = x;
		y.boxplus(delta, i);
		typename vect_type::base diff;
		y.boxminus(diff, x);
		BOOST_CHECK_LE(diff.norm(), maxNorm*(1.0 + tol));
	}
}

/**
 * Check (x [+] delta) [-] x == delta, for any delta.
 */
template<class Manifold, class vector>
void check_axiom_minus_plus(const Manifold &x, const vector &v)
{
	typedef typename Manifold::scalar scalar;
	typedef MTK::vect<Manifold::DOF, scalar> vect_type;
	scalar tol = MTK::tolerance<scalar>();
	vect_type delta = v;
	for(scalar i=0; i < 10.0; i+=1.0/16.0)
	{
		Manifold y = x;
		y.boxplus(delta, i);
		typename vect_type::base diff;
		y.boxminus(diff, x);
		MTK_CHECK_CLOSE(diff, (i*delta).eval(), tol);
	}
}


/**
 * Check (x [+] delta) [-] x == delta, for delta small enough.
 */
template<class Manifold, class vector>
void check_axiom_minus_plus(const Manifold &x, const vector &v, const double &maxNorm)
{
	typedef typename Manifold::scalar scalar;
	typedef MTK::vect<Manifold::DOF, scalar> vect_type;
	scalar tol = MTK::tolerance<scalar>();
	vect_type delta = v;
	delta.normalize();
	for(scalar i=0; i < maxNorm; i+=maxNorm/32.0)
	{
		Manifold y = x;
		y.boxplus(delta, i);
		typename vect_type::base diff;
		y.boxminus(diff, x);
		MTK_CHECK_CLOSE(diff, (i*delta).eval(), tol);
	}
}

template<class Manifold, class vector>
void check_axiom_triangle(const Manifold &x, const vector &d1, const vector &d2)
{
	typedef typename Manifold::scalar scalar;
	typedef MTK::vect<Manifold::DOF, scalar> vect_type;
	scalar tol = MTK::tolerance<scalar>();
	//vect_type delta = v;
	for(double i=0; i <= 10; i+=0.5)
	{
		Manifold y1 = x, y2 = x;
		y1.boxplus(d1, i);
		y2.boxplus(d2, i);
		vect_type diff;
		y1.boxminus(diff, y2);
		BOOST_CHECK_LE(diff.norm(), (i*d1 - i*d2).norm() + tol);
	}
}



template<class Manifold, class vector>
void check_axioms(const Manifold &x, const Manifold &y, const Manifold &z, 
                  const vector  &v1, const vector  &v2, const typename Manifold::scalar &maxNorm)
{
	check_not_equal(x, y);
	check_not_equal(x, z);
	check_not_equal(y, z);
	
	check_axiom_zero(x);
	check_axiom_zero(y);
	check_axiom_zero(z);
	
	check_axiom_plus_minus(x, x);
	check_axiom_plus_minus(y, y);
	check_axiom_plus_minus(z, z);
	
	check_axiom_plus_minus(x, y);
	check_axiom_plus_minus(x, z);
	check_axiom_plus_minus(y, z);
	
	
	if(maxNorm > 0)
	{
		check_minus_le_maxnorm(x, v1, maxNorm);
		check_minus_le_maxnorm(y, v1, maxNorm);
		check_minus_le_maxnorm(z, v1, maxNorm);
		
		check_axiom_minus_plus(x, v1, maxNorm);
		check_axiom_minus_plus(y, v1, maxNorm);
		check_axiom_minus_plus(z, v1, maxNorm);
		check_axiom_minus_plus(x, v2, maxNorm);
		check_axiom_minus_plus(y, v2, maxNorm);
		check_axiom_minus_plus(z, v2, maxNorm);
	}
	else
	{
		check_axiom_minus_plus(x, v1);
		check_axiom_minus_plus(y, v1);
		check_axiom_minus_plus(z, v1);
		
		check_axiom_minus_plus(x, v2);
		check_axiom_minus_plus(y, v2);
		check_axiom_minus_plus(z, v2);
	}
	
	check_axiom_triangle(x, v1, v2);
	check_axiom_triangle(y, v1, v2);
	check_axiom_triangle(z, v1, v2);

}

/**
 * Check correctness of multiplicative groups.
 */
template<class Manifold>
void check_group(const Manifold &x, const Manifold &y, const Manifold &z)
{
	Manifold one; // default constructor shall construct neutral element
	
	// Neutral element:
	check_equal<Manifold>(one * x, x);
	check_equal<Manifold>(x * one, x);
	check_equal<Manifold>(one * y, y);
	check_equal<Manifold>(y * one, y);
	check_equal<Manifold>(one * z, z);
	check_equal<Manifold>(z * one, z);
	
	// associativity:
	check_equal<Manifold>(x * (y * z), (x * y) * z);
	check_equal<Manifold>(y * (x * z), (y * x) * z);
	check_equal<Manifold>(z * (y * x), (z * y) * x);
	
	// inverse:
	check_equal<Manifold>(x * x.inverse(), one);
	check_equal<Manifold>(y * y.inverse(), one);
	check_equal<Manifold>(z * z.inverse(), one);
	check_equal<Manifold>(x.inverse() * x, one);
	check_equal<Manifold>(y.inverse() * y, one);
	check_equal<Manifold>(z.inverse() * z, one);
	
	check_equal<Manifold>(x.inverse().inverse(), x);
	check_equal<Manifold>(y.inverse().inverse(), y);
	check_equal<Manifold>(z.inverse().inverse(), z);
	
	// inverse of product:
	check_equal<Manifold>((x*y).inverse(), y.inverse() * x.inverse());
	check_equal<Manifold>((x*z).inverse(), z.inverse() * x.inverse());
	check_equal<Manifold>((z*y).inverse(), y.inverse() * z.inverse());
	
	// division operator:
	check_equal<Manifold>(x / y, x * y.inverse());
	check_equal<Manifold>(y / z, y * z.inverse());
	check_equal<Manifold>(z / x, z * x.inverse());
	
	// % operator:
	check_equal<Manifold>(x % y, x.inverse() * y);
	check_equal<Manifold>(y % z, y.inverse() * z);
	check_equal<Manifold>(z % x, z.inverse() * x);
}



#endif /* MTK_AXIOMS_HH_ */
