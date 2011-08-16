/*
 *  Copyright (c) 2011, Universitaet Bremen
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
 * @file euler_stable.hpp
 * @brief Stable conversions from Quaterion to Euler angles and back
 * 
 * @author Christoph Hertzberg
 * 
 * Did I mention that I don't like Euler angles?
 */


#ifndef EULER_STABLE_HPP_
#define EULER_STABLE_HPP_

#include <cmath>
#include <Eigen/Geometry>

/**
 * Euler to Quaternion conversion.
 */
template<class Scalar>
static void euler_to_quat(Eigen::Quaternion<Scalar> &res, const Scalar& roll, const Scalar& pitch, const Scalar& yaw){
	using std::cos;
	using std::sin;
	
	Scalar sr = sin(roll*0.5);
	Scalar cr = cos(roll*0.5);
	Scalar sp = sin(pitch*0.5);
	Scalar cp = cos(pitch*0.5);
	Scalar sy = sin(yaw*0.5);
	Scalar cy = cos(yaw*0.5);
	res.w() = cr*cp*cy + sr*sp*sy;
	res.x() = sr*cp*cy - cr*sp*sy;
	res.y() = cr*sp*cy + sr*cp*sy;
	res.z() = cr*cp*sy - sr*sp*cy;
}

/**
 * Convert Quaternion quat to Euler angles, gracefully handling the singularity
 * for abs(pitch) = PI/2.
 * 
 * Passing the result of this method to the method above should always result 
 * in quat or -quat assuming quat was normalized.
 * This method also has no problems handling non-normalized Quaternions.
 * 
 * I REALLY DON'T LIKE EULER ANGLES!
 */
template<class Scalar>
static void quat_to_euler(Scalar& roll, Scalar& pitch, Scalar& yaw, const Eigen::Quaternion<Scalar>& quat){
	using std::cos;
	using std::sin;
	using std::atan2;
	
	// Get yaw angle:
	Scalar qx=quat.x(), qy=quat.y(), qz=quat.z(), qw=quat.w();
	Scalar qx2=qx*qx, qy2=qy*qy, qz2=qz*qz, qw2=qw*qw;
	// for abs(pitch) = PI/2 this will lead to atan2(0,0)
	// i.e. for noisy values, result will be arbitrary
	yaw = atan2(2*(qw*qz + qx*qy), qw2 + qx2 - qy2 - qz2);
	
	// Now rotate the original Quaternion backwards by yaw:
	Scalar c = cos(yaw/2), s=sin(yaw/2);
	Scalar px=c*qx+s*qy, py=c*qy-s*qx, pz=c*qz-s*qw, pw=c*qw+s*qz;
	Scalar px2=px*px, py2=py*py, pz2=pz*pz, pw2=pw*pw;
	
	// Now calculating pitch and roll does not have singularities anymore:
	pitch = atan2(2*(py*pw - px*pz), px2 + pw2 - py2 - pz2);
	roll  = atan2(2*(px*pw - py*pz), py2 + pw2 - px2 - pz2);
}






#endif /* EULER_STABLE_HPP_ */
