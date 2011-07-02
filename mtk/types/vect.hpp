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
 * @file mtk/types/vect.hpp
 * @brief Basic vectors interpreted as manifolds.
 * 
 * This file also implements a simple wrapper for matrices, for arbitrary scalars
 * and for positive scalars.
 */
#ifndef VECT_H_
#define VECT_H_

#include <iosfwd>

#include "../src/vectview.hpp"

namespace MTK {

/**
 * A simple vector class.
 * Implementation is basically a wrapper around Eigen::Matrix with manifold 
 * requirements added.
 */
template<int D, class _scalar = double>
struct vect : public Eigen::Matrix<_scalar, D, 1> {
	typedef Eigen::Matrix<_scalar, D, 1> base;
	enum {DOF = D};
	typedef _scalar scalar;
	
	using base::operator=;
	
	/** Standard constructor. Sets all values to zero. */
	vect() {
		base::setZero();
	}
	
	/** Constructor copying the value of the expression \a other */
	template<typename OtherDerived>
	EIGEN_STRONG_INLINE vect(const Eigen::MatrixBase<OtherDerived>& other) : base(other) {}
	
	/** Construct from memory. */
	vect(const scalar* src) : base(src) { }
	
	void boxplus(MTK::vectview<const scalar, DOF> vec, scalar scale=1) {
		*this += scale * vec;
	}
	void boxminus(MTK::vectview<scalar, DOF> res, const vect<D, scalar>& other) const {
		res = *this - other;
	}
	
	
	friend std::ostream& operator<<(std::ostream &os, const vect<D, scalar>& v){
		// Old versions of Eigen messed up with the streams flags, so do it manually:
		for(int i=0; i<DOF; ++i)
			os << v(i) << " ";
		return os;
	}
	friend std::istream& operator>>(std::istream &is, vect<D, scalar>& v){
		for(int i=0; i<DOF; ++i)
			is >> v(i);
		return is;
	}
	
	template<int dim>
	vectview<scalar, dim> tail(){
		BOOST_STATIC_ASSERT(0< dim && dim <= DOF);
		return base::data() + (DOF - dim);
	}
	template<int dim>
	vectview<const scalar, dim> tail() const{
		BOOST_STATIC_ASSERT(0< dim && dim <= DOF);
		return base::data() + (DOF - dim);
	}
	template<int dim>
	vectview<scalar, dim> head(){
		BOOST_STATIC_ASSERT(0< dim && dim <= DOF);
		return base::data();
	}
	template<int dim>
	vectview<const scalar, dim> head() const{
		BOOST_STATIC_ASSERT(0< dim && dim <= DOF);
		return base::data();
	}
};


/**
 * A simple matrix class.
 * Implementation is basically a wrapper around Eigen::Matrix with manifold 
 * requirements added, i.e., matrix is viewed as a plain vector for that.
 */
template<int M, int N, class _scalar = double, int _Options = Eigen::Matrix<_scalar, M, N>::Options>
struct matrix : public Eigen::Matrix<_scalar, M, N, _Options> {
	typedef Eigen::Matrix<_scalar, M, N, _Options> base;
	enum {DOF = M * N};
	typedef _scalar scalar;
	
	using base::operator=;
	
	/** Standard constructor. Sets all values to zero. */
	matrix() {
		base::setZero();
	}
	
	/** Constructor copying the value of the expression \a other */
	template<typename OtherDerived>
	EIGEN_STRONG_INLINE matrix(const Eigen::MatrixBase<OtherDerived>& other) : base(other) {}
	
	/** Construct from memory. */
	matrix(const scalar* src) : base(src) { }
	
	void boxplus(MTK::vectview<const scalar, DOF> vec, scalar scale=1) {
		*this += scale * base::Map(vec.data());
	}
	void boxminus(MTK::vectview<scalar, DOF> res, const matrix& other) const {
		base::Map(res.data()) = *this - other;
	}
	
	friend std::ostream& operator<<(std::ostream &os, const matrix<M, N, scalar, _Options>& mat){
		for(int i=0; i<DOF; ++i){
			os << mat.data()[i] << " ";
		}
		return os;
	}
	friend std::istream& operator>>(std::istream &is, matrix<M, N, scalar, _Options>& mat){
		for(int i=0; i<DOF; ++i){
			is >> mat.data()[i];
		}
		return is;
	}
};



/**
 * A simple scalar type.
 */
template<class _scalar = double>
struct Scalar {
	enum {DOF = 1};
	typedef _scalar scalar;
	
	scalar value;
	
	Scalar(const scalar& value = scalar(0)) : value(value) {}
	operator const scalar&() const { return value; }
	operator scalar&() { return value; }
	Scalar& operator=(const scalar& val) { value = val; return *this; }
	
	void boxplus(MTK::vectview<const scalar, DOF> vec, scalar scale=1) {
		value += scale * vec[0];
	}
	void boxminus(MTK::vectview<scalar, DOF> res, const Scalar& other) const {
		res[0] = *this - other;
	}
};

/**
 * Positive scalars.
 * Boxplus is implemented using multiplication by @f$x\boxplus\delta = x\cdot\exp(\delta) @f$.
 */
template<class _scalar = double>
struct PositiveScalar {
	enum {DOF = 1};
	typedef _scalar scalar;
	
	scalar value;
	
	PositiveScalar(const scalar& value = scalar(1)) : value(value) {
		assert(value > scalar(0));
	}
	operator const scalar&() const { return value; }
	PositiveScalar& operator=(const scalar& val) { assert(val>0); value = val; return *this; }
	
	void boxplus(MTK::vectview<const scalar, DOF> vec, scalar scale=1) {
		value *= std::exp(scale * vec[0]);
	}
	void boxminus(MTK::vectview<scalar, DOF> res, const PositiveScalar& other) const {
		res[0] = std::log(*this / other);
	}
	
	friend std::istream& operator>>(std::istream &is, PositiveScalar<scalar>& s){
		is >> s.value;
		assert(s.value > 0);
		return is;
	}
};

}  // namespace MTK


#endif /*VECT_H_*/
