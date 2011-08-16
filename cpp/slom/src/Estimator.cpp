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
 * @file slom/src/Estimator.cpp
 * @brief implements the Estimator class
 */

#include "../Estimator.hpp"

#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>

#include <iostream>



#include "../CallBack.hpp"


namespace SLOM {



void Estimator::freeWorkspace(){
//	JtJ      = cs_spfree(JtJ);
	symbolic = cs_sfree(symbolic);
	numeric  = cs_nfree(numeric);
}


Estimator::~Estimator(){
	freeWorkspace();
}


double Estimator::evaluate(double * result) const{
	CALLBACK_START(evaluate);
	return func.evaluate(result);
}


void Estimator::updateSparse(){
	{
		CALLBACK_START(calculateJacobian);
		func.calculateJacobian();
	}
	func.updateDiagonal(lamda, usedAlgorithm == LevenbergMarquardt ? func.cholCovariance.data() : 0);
}

bool Estimator::qrSolve(double* delta, const double *res){
	CALLBACK_START(qrSolve);
	const cs * jacobian = func.jacobian;
	assert(symbolic);
	cs_nfree(numeric);
	numeric = cs_qr(jacobian, symbolic);
	assert(numeric);
	int m=jacobian->m, n=jacobian->n;
	
	// The following code essentially does a cs_qrsol(3, matrix, workspace);
	// but it doesn't recalculate the symbolic decomposition
	
	cs_ipvec(symbolic->pinv, res, workspace.data(), m) ; /* x(0:m-1) = b(p(0:m-1) */
	for (int k = 0; k < n; k++) /* apply Householder refl. to x */
	{
		cs_happly(numeric->L, k, numeric->B [k], workspace.data()) ;
	}
	cs_usolve(numeric->U, workspace.data()) ; /* x = R\x */
	
	cs_ipvec(symbolic->q, workspace.data(), delta, n) ; /* b(q(0:n-1)) = x(0:n-1) */
	// end of qrsol
	return true; // qrSolve only fails, if out of memory (result could be non-finite, though)
}


bool Estimator::choleskySolve(double *delta, const double *res){
	CALLBACK_START(choleskySolve);
//	const cs * jacobian = func.jacobian;
	int m = func.getM(), n=func.getN(); (void)m;
	assert((delta - res >= m || res - delta >= n) 
	       && "delta and residuum shall not overlap!");
	
	std::fill(delta, delta + n, 0);
	func.add_Jtu(delta, res);
	
	const cs *JtJ = func.get_JtJ();
	
	
	if (true) { 
		CALLBACK_START(choleskyDecompose);
		csn *N = cs_di_chol (JtJ, symbolic) ;                    /* numeric Cholesky factorization */
		if(!N) return false;
		
		double *x = workspace.data();
		assert (workspace.size()>=n);
		assert (n==N->L->n);		
		cs_ipvec (symbolic->pinv, delta, x, n) ;   /* x = P*b */
		cs_lsolve (N->L, x) ;           /* x = L\x */
		cs_ltsolve (N->L, x) ;          /* x = L'\x */
		cs_pvec (symbolic->pinv, x, delta, n) ;    /* b = P'*x */
		cs_nfree (N) ;
		return true;
	}
	else return cs_di_cholsol(1, JtJ, delta);	
}


void Estimator::init(){
	if(!func.initialize(usedAlgorithm !=GaussNewton))
		return; // if SparseFunction has not changed, we don't need to re-initialize.
	
	// free Estimators workspace:
	freeWorkspace();
	
	std::pair<int, int> size = func.getSize();
	
	switch(usedSolver){
	case QR:
		assert(false && "QR not supported currently");
		symbolic = cs_sqr(3, func.jacobian, true);
		assert(symbolic);
		break;
	case Cholesky:
		symbolic = cs_schol(1, func.get_JtJ());
		assert(symbolic);
		break;
	}
	
	res.resize(size.first);
	newRes.resize(size.first);
	lastRSS = -1;
	workspace.resize(size.second);
	delta.resize(size.second);
}


void Estimator::calculate_delta()
{
	CALLBACK_START(calculate_delta);
	
	init();
	assert(func.jacobian);
	
	int n=func.jacobian->n;
	
	updateSparse();
	
	if(lastRSS < 0){
		lastRSS = evaluate(res.data());
	}
	
	bool success;
	// Solve the linear system and store the result in delta:
	switch(usedSolver){
	case QR:
		success = qrSolve(delta.data(), res.data());
		break;
	case Cholesky:
		success = choleskySolve(delta.data(), res.data());
		break;
	default:
		success = false; // should never happen!
	}
	if(!success) {
		CALLBACK_FINISH(optimizeStep);
		throw "Solving linear system failed";
	}
	
	//some information about the delta_vector:
	double normInf= delta.lpNorm<Eigen::Infinity>(), norm2=delta.squaredNorm();
	cb.output(3) << "Update vector: rms: " << std::sqrt(norm2/n) << " normInf: " << normInf; 
	cb.output(4) << "\ndelta:";
	for(int i = 0; i < std::min(20, n); ++i)
		cb.output(4) << " " << delta[i];
	cb.output(3) << " ...\n";
}

double Estimator::apply_delta(double scale){
	func.apply_delta(delta, scale);
	// calculate the new RSS:
	return evaluate(newRes.data());
}

double Estimator::store_or_restore(double newRSS, scoped_callback& cb)
{
	int m=func.jacobian->m;
	double gain = (lastRSS - newRSS)/newRSS;
	cb.output(2) << ", RSS: " << newRSS << ", RMS: " << std::sqrt(newRSS/m) << ", Gain: " << gain;
	if(gain > 0 || usedAlgorithm == GaussNewton){
		// positive gain or GaussNewton: 
		// Store modified variables permanently, current RSS to res, and reduce lamda.
		func.store_or_restore(true);
		
		res.swap(newRes);
		lastRSS = newRSS;
		lamda *= sqrt(0.1);
	} else {
		// Restore variables, increase lamda
		func.store_or_restore(false);
		lamda *= sqrt(10.0);
	}
	if(usedAlgorithm != GaussNewton){
		cb.output(2) << ", lamda = " << lamda;
	}
	return gain;
	
}

double Estimator::optimizeStep(bool linear_refine){
	CALLBACK_START(optimizeStep);
	// TODO better parameter control for LMA
	// TODO factor out into Algorithm class
	
	
	// calculate delta = (J^T * J + lambda*I) \ (J^T * res):
	calculate_delta();
	
	double newRSS;
	if(!linear_refine)
	{
		newRSS = apply_delta(-1);
	}
	else
	{
		// try linear refinement (highly experimental and not optimized)
		cb.output(5) << "RSS(scale): ";
		double bestRSS = lastRSS;
		double bestScale = 0;
		for(double scale = -2; scale<0; scale+=1.0/16)
		{
			newRSS = apply_delta(scale);
			cb.output(5) << newRSS << " ";
			if(newRSS < bestRSS)
			{
				bestRSS = newRSS;
				bestScale = scale;
			}
		}
		cb.output(5) << lastRSS << "\n";
		
		newRSS = apply_delta(bestScale);
	}
	
	cb.output(2) << "lastRSS: " << lastRSS; 
	double gain = store_or_restore(newRSS, cb);
	
	cb.output(2) << std::endl;

	return gain;
}

std::ostream& Estimator::dumpJacobian(std::ostream& out)
{
	const cs * jacobian = func.getJ();
	assert(jacobian);
	
	int n = jacobian->n;
	int *Ap = jacobian->p, *Ai = jacobian->i;
	double *Ax = jacobian->x;
	for (int j = 0 ; j < n ; ++j)
	{
		for (int p = Ap [j] ; p < Ap [j+1] ; ++p)
		{
			out << Ai[p] << " " << j << " " << Ax[p] << "\n";
		}
	}
	return out;
}


}  // namespace SLOM
