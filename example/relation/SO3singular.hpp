/*
 *  Copyright (c) 2010--2011, Universitaet Bremen
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
 * @file SO3singular.hpp
 * @brief Representations of SO(3) which you should NOT use.
 * 
 * This was implemented mainly for experimental reasons.
 */


#ifndef SO3SCALEDAXIS_HPP_
#define SO3SCALEDAXIS_HPP_

namespace MTK {

template<class Storage>
struct SO3singular : public Storage{
	enum {DOF = Storage::DOF, DIM = 3};
	typedef typename Storage::scalar scalar;
	
	typedef Eigen::Quaternion<scalar> quat;
	typedef typename MTK::SO3<scalar> SO3;
	typedef typename SO3::vect_type vect_type;
	
	typedef Storage base;
	
	operator SO3() const {
		return toSO3();
	}
	
	SO3 toSO3() const {
		return base::operator SO3();
	}
	
	typename quat::Coefficients coeffs() const {
		return toSO3().coeffs();
	}
	SO3 operator%(const SO3 &r) const {
		return toSO3() % r;
	}
	//! Calculate @c this->inverse() * @c r
	SO3 operator*(const SO3 &r) const {
		return toSO3() * r;
	}
	
	template<class Derived>
	vect_type operator*(const Eigen::MatrixBase<Derived> &vec) const {
		return toSO3() * vec;
	}
	//! Calculate @c this->inverse() * @c r
	template<class Derived>
	vect_type operator%(const Eigen::MatrixBase<Derived> &vec) const {
		return toSO3() % vec;
	}
	
	//! Calculate @c this * @c r.conjugate()
	SO3 operator/(const SO3& r) const {
		return toSO3() * r.conjugate();
	}
	
	/**
	 * Construct from real part and three imaginary parts.
	 * Quaternion is normalized after construction.
	 */
	SO3singular(const scalar& w, const scalar& x, const scalar& y, const scalar& z) : base(quat(w, x, y, z)) { }
	
	/**
	 * Construct from Eigen::Quaternion.
	 * @note Non-normalized input may result result in spurious behavior.
	 */
	SO3singular(const quat& src = quat::Identity()) : base(src) {}
	
	/**
	 * Construct from rotation matrix.
	 * @note Invalid rotation matrices may lead to spurious behavior.
	 */
	template<class Derived>
	SO3singular(const Eigen::MatrixBase<Derived>& matrix) : base(quat(matrix)) {}
	
	//! @name Manifold requirements
	//{
	void boxplus(MTK::vectview<const scalar, DOF> vec, scalar scale=1) {
		base::boxplus(vec, scale);
	}
	void boxminus(MTK::vectview<scalar, DOF> res, const SO3singular<Storage>& other) const {
		base::boxminus(res, other);
	}
	//}
	
	quat inverse() const {
		return toSO3().conjugate();
	}
	
	friend std::ostream& operator<<(std::ostream& os, const SO3singular<Storage>& q){
		return os << q.toSO3();
	}
	friend std::istream& operator>>(std::istream& is, SO3singular<Storage>& res){
		SO3 q;
		is >> q;
		res = q;
		return is;
	}
};


template<class _scalar>
struct scaledAxis {
	enum {DOF = 3};
	
	typedef _scalar scalar;
	typedef Eigen::Quaternion<scalar> quat;
	typedef typename MTK::SO3<scalar> SO3;
	typedef Eigen::Matrix<scalar, DOF, 1> vect_type;
	vect_type _scaled_axis;
	
	operator SO3() const {
		return SO3::exp(_scaled_axis);
	}
	
	/**
	 * Construct from Eigen::Quaternion.
	 * @note Non-normalized input may result result in spurious behavior.
	 */
	scaledAxis(const quat& src = quat::Identity()) : _scaled_axis(SO3::log(src)) {}
	
	//! @name Manifold requirements
	//{
	void boxplus(MTK::vectview<const scalar, DOF> vec, scalar scale=1) {
		_scaled_axis += scale * vec;
	}
	void boxminus(MTK::vectview<scalar, DOF> res, const scaledAxis<scalar>& other) const {
		res = _scaled_axis - other._scaled_axis;
	}
	//}
};

template<class _scalar>
struct euler {
	enum {DOF = 3};
	
	typedef _scalar scalar;
	typedef Eigen::Quaternion<scalar> quat;
	typedef typename MTK::SO3<scalar> SO3;
	typedef Eigen::Matrix<scalar, DOF, 1> vect_type;
	vect_type _ang;
	
	operator SO3() const {
		return toQuat();
	}
	
	//! @name Manifold requirements
	//{
	void boxplus(MTK::vectview<const scalar, DOF> vec, scalar scale=1) {
		_ang += scale * vec;
	}
	void boxminus(MTK::vectview<scalar, DOF> res, const euler<scalar>& other) const {
		res = _ang - other._ang;
	}
	//}
	
	euler(const quat &q = quat::Identity()) {
		setFromQuat(q);
	}
	
private:
	/**
	 * Convert Euler-Angle to Quaternion.
	 * Code based on: 
	 * http://www.euclideanspace.com/maths/geometry/rotations/conversions/eulerToQuaternion/
	 */
	quat toQuat() const
//			double heading, double attitude, double bank)
	{
		scalar c1=std::cos(_ang[0]/2), c2=std::cos(_ang[1]/2), c3=std::cos(_ang[2]/2);
		scalar s1=std::sin(_ang[0]/2), s2=std::sin(_ang[1]/2), s3=std::sin(_ang[2]/2);
		
		scalar w = c1 * c2 * c3 - s1 * s2 * s3;
		scalar x = s1 * s2 * c3 + c1 * c2 * s3;
		scalar y = s1 * c2 * c3 + c1 * s2 * s3;
		scalar z = c1 * s2 * c3 - s1 * c2 * s3;
		
		return quat(w, x, y,z);
	}
	
	void setFromQuat(const quat &q) {
		scalar w=q.w(), x=q.x(), y=q.y(), z=q.z();
		scalar sqw = w*w;
		scalar sqx = x*x;
		scalar sqy = y*y;
		scalar sqz = z*z;
		scalar unit = sqx + sqy + sqz + sqw; // if normalized is one, otherwise is correction factor
		scalar test = x*y + z*w;
		if (test > 0.49999*unit) { // singularity at north pole
			_ang << 2 * std::atan2(x, w), M_PI/2, 0;
//			std::cerr << test/unit << std::endl << q.coeffs().transpose() << std::endl << toQuat().coeffs().transpose() << std::endl;
			return;
		}
		if (test < -0.49999*unit) { // singularity at south pole
			_ang << -2 * std::atan2(x, w), -M_PI/2, 0;
//			std::cerr << test/unit << std::endl << q.coeffs().transpose() << std::endl << toQuat().coeffs().transpose() << std::endl;
			return;
		}
		_ang << 
				std::atan2(2*y*w-2*x*z , sqx - sqy - sqz + sqw),
				std::asin (2*test/unit),
				std::atan2(2*x*w-2*y*z , -sqx + sqy - sqz + sqw);
	}
};

template<class _scalar>
struct vec4SO3 {
	enum {DOF = 4};
	typedef _scalar scalar;
	typedef Eigen::Quaternion<scalar> quat;
	typedef typename MTK::SO3<scalar> SO3;
	
	quat q;
	operator SO3() const {
		return q.normalized();
	}
	vec4SO3(const quat &src = quat::Identity()) : q(src) {};
	
	//! @name Manifold requirements
	//{
	void boxplus(MTK::vectview<const scalar, DOF> vec, scalar scale=1) {
		q.coeffs() += scale * vec;
	}
	void boxminus(MTK::vectview<scalar, DOF> res, const vec4SO3<scalar>& other) const {
		res = q.coeffs() - other.q.coeffs();
	}
	//}
	
	scalar squaredNorm() const {
		return q.squaredNorm();
	}
	scalar norm() const {
		return q.norm();
	}
	
};



}  // namespace MTK



#endif /* SO3SCALEDAXIS_HPP_ */
