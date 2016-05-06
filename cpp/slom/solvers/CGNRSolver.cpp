/**
 * @file /mtk-trunk/slom/solvers/CGNRSolver.cpp
 * @brief Brief description
 * 
 */

#include "CGNRSolver.hpp"

#include <slom/CallBack.hpp>

namespace SLOM {

// FIXME this file needs a lot of clean-up


bool CGNRSolver::solve(VectorType& delta, const double & lambda){
	CALLBACK_START_M();
	
	int m = getM(), n = getN(); (void)(m);
	
	double tol = 1e-7; // FIXME add parameter
	
	delta.setZero(n);
	const SparseType &A = getJ();
	
	
	int zn = z.size();
	if(recycle && zn > 0 && n2s > 0) { // n2s <=0 means, last solution was exact to numeric accuracy, zn==0 means no previous iterations
		const VectorType& g = getGradiant();
		assert(zn == p.size() && n2s > 0 && "Something went wrong in previous iteration");
		double sTz = (zn==n) ? z.dot(g) : z.dot(g.head(zn));
		z = getGradiant();
		applyPreconditioner(z, true);
		double n2s_new = z.head(zn).squaredNorm();
		double beta = (n2s_new - sTz)/n2s;
		cb(5) << "CGNR> beta=" << beta << "\tn2s=" << n2s_new << "\ts^tz=" << sTz << "\tn2s_old=" << n2s << std::endl;
		n2s = n2s_new + z.tail(n - zn).squaredNorm();
		applyPreconditioner(z, false);
		
		if(true && zn >0) {
			// NOTE new approach to calculate recycled beta, ensures that p_new is orthogonal to p_old
			// FIXME initialize last elements of p by linearized initialization function of new elements
			VectorType w = A.middleCols(0, zn) * p;
			double n2w = w.squaredNorm();
			if(!(n2w > 0)){
				beta = 0;
			} else {
				VectorType v = A.transpose() * w;
				beta = -v.dot(z)/n2w;
			}
			
			cb(5) << "CGNR> better beta=" << beta << std::endl;
			if(beta < 0 || beta > 10) beta = 0; // FIXME what's a reasonable upper bound for beta?
			VectorType p_old;
			p.swap(p_old);
			p = z;
			p.head(zn) += beta * p_old;
//			p = z + beta * p;
		} else if(beta < 0.0 || beta > 20.0) {
			cb(5) << "CGNR> Reset search direction\n";
			p = z;
			beta = 0;
		} else if( zn == n) {
			p = z + beta * p;
		} else {
			VectorType p_old;
			p.swap(p_old);
			p = z;
			p.head(zn) += beta * p_old;
		}
	} else {
		z = getGradiant();
		applyPreconditioner(z, true);
		n2s = z.squaredNorm();
		applyPreconditioner(z, false);
		p = z;
	}
	
	double n2s0 = n2s; // initial norm(s)^2, for convergence testing
	
	VectorType v;
	VectorType w; // workspace
	for(int j=0; j<max_it; ++j){
		w = A * p;
		double n2w = w.squaredNorm() + lambda * p.squaredNorm();
		if(n2w < tol*tol) {
			n2s = 0; // reset search direction
			break;
		}
		double alpha = n2s / n2w;
		delta += alpha * p;
		
		if(j==max_it-1) break;
		
		v = A.transpose() * w + lambda * p;
		double zs = z.dot(v);
		applyPreconditioner(v, true);
		double n2s_new = n2s - 2*alpha * zs + alpha*alpha* v.squaredNorm();
		if(n2s_new<= tol * n2s0) {
			cb(2) <<  "CGNR> Early Break\tj=" << j << "\tn = " << n << "\tn2s_new = " << n2s_new << "\tres=" << z.norm() << "\n";
			n2s = -1;
			break;
		}
		applyPreconditioner(v, false);
		z -= alpha * v;
		double beta = n2s_new / n2s;
		n2s = n2s_new;
		if(!(j&0) && !true)
			cb(10) <<"CGNR>>>" << j << ":\ta: " << alpha << ",\tb: " << beta << ",\tn2s: " << n2s << ",\tzs: " << zs << std::endl; 
		p = z + beta * p;
	}
	
	
	return true;
}



} /* namespace SLOM */
