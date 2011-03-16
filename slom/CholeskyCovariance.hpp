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
 * @file slom/CholeskyCovariance.hpp
 * @brief Helper class to apply Covariance or Information matrix to measurements
 */
#ifndef CHOLESKYCOVARIANCE_H_
#define CHOLESKYCOVARIANCE_H_

#include <algorithm>
#include <numeric>
#include <cassert>

#include "../mtk/src/vectview.hpp"

namespace SLOM {



namespace CholeskyMode{
/**
 * How to initialize a CholeskyCovariance.
 * COPY_* means source already is a Cholesky factor,
 * CHOLESKY_* does the decomposition itself.
 * *_UPPER means only the upper half of the matrix is stored,
 * *_FULL means zeros or symmetric values from lower half are stored.
 */
enum CM {
	COPY_UPPER,
	COPY_UPPER_FULL,
	COPY_LOWER,
	COPY_LOWER_FULL,
	CHOLESKY_UPPER,
	CHOLESKY_FULL
};
};

/**
 * This class provides a possibility to normalize a measurement.
 * SLoM requires measurements to have standard normal distribution, so for a
 * measurement having distribution @f$f(x) \sim \Ndp\mu\Sigma @f$ the function
 * @f$ \tilde f(x) := L\inv(f(x) \mmnus \mu @f$ is distributed 
 * @f$ \tilde f(x)\sim\Ndp0\I @f$ if @f$ \Sigma = L\trans L @f$.
 * In some cases not the covariance but the information matrix 
 * @f$\Omega = \Sigma\inv @f$ is available, then for @f$\Omega = L\trans L @f$
 * one has to calculate @f$ \tilde f(x) := L(f(x) \mmnus \mu @f$.
 * 
 * Performing the necessary Cholesky decomposition and multiplying by 
 * @f$ L @f$ (@ref apply) or @f$ L\inv @f$ (@ref invApply) is done by this class
 * 
 * 
 * Notice that often @f$\Sigma @f$ or @f$\Omega @f$ have diagonal structure, 
 * in which case simple element-wise multiplication/division is more efficient.
 */
template<int dim>
struct CholeskyCovariance
{
	enum {DIM = dim, SIZE = (dim*(dim+1))/2};
	// chol is a lower triangular matrix, saved in row major order
	double chol[SIZE];
	
	
	CholeskyCovariance() {};
	
	CholeskyCovariance(const double *A, CholeskyMode::CM mode){
		bool skip = false;
		switch(mode){
		case CholeskyMode::COPY_LOWER:
			skip = true;
		case CholeskyMode::COPY_LOWER_FULL:
			copyCholeskyLower(A, skip);
			break;
		case CholeskyMode::COPY_UPPER:
			skip = true;
		case CholeskyMode::COPY_UPPER_FULL:
			copyCholeskyUpper(A, skip);
			break;
		case CholeskyMode::CHOLESKY_UPPER:
			skip = true;
		case CholeskyMode::CHOLESKY_FULL:
			calculateCholesky(A, skip);
			break;
		default:
			assert(false); //Unknown CholeskyMode
		}
	}
	/**
	 * Initialize by copying lower values from data.
	 */
	void copyCholeskyLower(const double *A, bool skip){
		if(skip){
			std::copy(A, A+SIZE, chol);
			return;
		}
		const double *a = A;
		double *rowi = chol;
		for(int i=1; i<=DIM; i++){
			for(int j=0; j<i; j++){
				*rowi++ = *a++;
			}
			a += DIM-i;
		}
		assert(a==A+DIM*DIM);
		assert(rowi == chol+SIZE);
	}
	/**
	 * Initialize by copying upper values from data.
	 */
	void copyCholeskyUpper(const double *A, bool skip){
		const double *a = A;
		for(int i=0; i<DIM; ++i)
		{ // copy row i from A to col i from chol
			double *coli = chol + ( (i*(i+3)) >> 1);
			if(!skip){ // skip the first elements of this row in A:
				a += i;
			}
			for(int j=i; j<DIM; ++j)
			{ // chol(j,i) = A(i,j);
				*coli = *a++;
				coli += j+1; 
			}
		}
	}
	
	/**
	 * Initialize by Cholesky decomposition of matrix A.
	 * A is stored row major and only the upper half of A is considered.
	 * If parameter skip == true, it is assumed, that only the upper half of A is stored.
	 * I.e: A = [a00 a01 .. a0n a11 a12 .. a1n a22 ...]
	 * (As A is symmetric, you can substitute "row major" by "column major"
	 * and every "upper" by "lower")
	 */
	void calculateCholesky(const double *A, bool skip){
		const double *a = A;
		double *coli = chol;
		for(int i=0; i<DIM; i++){
			if(!skip){ // skip the first elements of this row in A:
				a += i;
			}
			coli += i; // start of column i
			// C(i,i) = sqrt(A(i,i) - C(0:i,i)'*C(0:i,i))
			double sum = *a++ - std::inner_product(coli, coli+i, coli, 0.0 );
			assert(sum>0);
			double sqrtSum = coli[i] = sqrt(sum);
			sqrtSum = 1/sqrtSum;
			double *colj = coli;
			for(int j=i+1; j<DIM; j++){
				colj += j;
				// C(i,j) = (A(i,j) - C(0:i,i)'*C(0:i,j))/C(i,i)
				colj[i] = (*a++ - std::inner_product(coli, coli+i, colj, 0.0)) * sqrtSum;
			}
		}
		assert(skip ? (a==A+SIZE) : (a==A+DIM*DIM));
	}
	
	/**
	 * multiplies arr[0,dim) by inverse of this Cholesky factor.
	 * I.e. computes @f$ chol^{-1}\cdot arr@f$;
	 */
	void invApply(MTK::vectview<double, DIM> arr) const {
		const double *c = chol;
		for(double *a = arr.data(); a < arr.data() + DIM; ++a){
			double ai = *a;
			for(double *b=arr.data(); b<a; ++b){
				ai -= *c++ * *b;
			}
			*a = ai / *c++;
		}
	}
	
	/**
	 * multiplies arr[0,dim) by this Cholesky factor.
	 * I.e. computes chol*arr;
	 */
	void apply(MTK::vectview<double, DIM> arr) const {
		const double *c = chol + SIZE -1;
		for(double *a=arr.data() + DIM - 1; a >= arr.data(); --a){
			*a *= *c--;
			for(double *b = a-1; b>= arr.data(); --b){
				*a += *c-- * *b;
			}
		}
		
	}
};

}  // namespace SLOM

#endif /*CHOLESKYCOVARIANCE_H_*/
