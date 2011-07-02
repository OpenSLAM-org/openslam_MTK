/*
 *  Copyright (c) 2008--2011, Universitaet Bremen
 *  All rights reserved.
 *
 *  Author: Christoph Hertzberg <chtz@informatik.uni-bremen.de>
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the Universitaet Bremen nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */
/**
 * @file mtk/src/mtkmath.hpp
 * @brief several math utility functions.
 */

#ifndef MTKMATH_H_
#define MTKMATH_H_

#include <cmath>

#include "vectview.hpp"

#ifndef M_PI
#define M_PI  3.1415926535897932384626433832795
#endif


namespace MTK {

namespace internal {

template<class Manifold>
struct traits {
	typedef typename Manifold::scalar scalar;
	enum {DOF = Manifold::DOF};
	typedef vect<DOF, scalar> vectorized_type;
	typedef Eigen::Matrix<scalar, DOF, DOF> matrix_type;
};

template<>
struct traits<float> : traits<Scalar<float> > {};
template<>
struct traits<double> : traits<Scalar<double> > {};

}  // namespace internal

/**
 * \defgroup MTKMath Mathematical helper functions
 */
//@{

//! constant @f$ \pi @f$
const double pi = M_PI;

template<class scalar> inline scalar tolerance();

template<> inline float  tolerance<float >() { return 1e-5f; }
template<> inline double tolerance<double>() { return 1e-11; }


/**
 * normalize @a x to @f$[-bound, bound] @f$.
 * 
 * result for @f$ x = bound + 2\cdot n\cdot bound @f$ is arbitrary @f$\pm bound @f$.
 */
template<class scalar>
inline scalar normalize(scalar x, scalar bound){
	if(std::fabs(x) <= bound) return x;
	int r = (int)(x / bound);
	return x - ((r + (r>>31) + 1) & ~1)*bound; 
}


template<class scalar, int n>
scalar exp(vectview<scalar, n> result, vectview<const scalar, n> vec, const scalar& scale = 1) {
	scalar norm = vec.norm();
	// protect against division by zero:
	norm = std::max(norm,  tolerance<scalar>());
	scalar alpha = scale * norm;
	scalar mult = std::sin(alpha) / norm;
	result = mult * vec;
	return std::cos(alpha);
}


/**
 * Inverse function to @c exp.
 * 
 * @param result @c vectview to the result
 * @param w      scalar part of input
 * @param vec    vector part of input
 * @param scale  scale result by this value
 * @param plus_minus_periodicity if true values @f$[w, vec]@f$ and @f$[-w, -vec]@f$ give the same result 
 */
template<class scalar, int n>
void log(vectview<scalar, n> result,
		const scalar &w, const vectview<const scalar, n> vec,
		const scalar &scale, bool plus_minus_periodicity)
{
	scalar nv = vec.norm();
	if(nv < tolerance<scalar>()) {
		if(!plus_minus_periodicity) {
			// find the maximal entry:
			int i;
			vec.maxCoeff(&i);
			result = scale * std::atan2(0, w) * vect<n, scalar>::Unit(i);
			return;
		}
		nv = tolerance<scalar>();
	}
	scalar s = scale / nv * (plus_minus_periodicity ? std::atan(nv / w) : std::atan2(nv, w) );
	
	result = s * vec;
}

//@}

} // namespace MTK


#endif /* MTKMATH_H_ */
