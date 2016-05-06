/**
 * @file LSQRSolver.cpp
 * @brief Refactored and simplified version of the LSQR solver from Paige and Saunders
 * 
 * LSQR uses an iterative (conjugate-gradient-like) method.
 * For further information, see 
 * 1. C. C. Paige and M. A. Saunders (1982a).
 *    LSQR: An algorithm for sparse linear equations and sparse least squares,
 *    ACM TOMS 8(1), 43-71.
 * 2. C. C. Paige and M. A. Saunders (1982b).
 *    Algorithm 583.  LSQR: Sparse linear equations and least squares problems,
 *    ACM TOMS 8(2), 195-209.
 * 3. M. A. Saunders (1995).  Solution of sparse rectangular systems using
 *    LSQR and CRAIG, BIT 35, 588-604.
 * 
 */

#include "LSQRSolver.hpp"

#include <cmath>
#include <cstdio>

namespace SLOM {

using std::min;
using std::sqrt;

bool LSQRSolver::solve(VectorType& x, const double &lambda2) {
	int m = getM(), n=getN();   // Determine dimensions m and n
	int itnlim = min(min(m, n), this->max_it); // iteration limit
	
	double lambda = sqrt(lambda2); // LSQR works with $\sqrt{\lambda}$:
	
	init();  // Initialize GKL process
	
	// These variables are updated at each GKL step:
	const double &alpha = getAlpha(), &beta  = getBeta();
	const VectorType &v = getZ();   // For LSQR v==z
	
	// Initialize variables for 1st iteration.
	double rhobar = alpha,  phibar = beta;
	VectorType w = v;  x = VectorType::Zero(n);
	
	for(int k = 0; k < itnlim; ++k){
		step(); // After this call, beta = $\beta_{k+1}$, alpha = $\alpha_{k+1}$.
		
		// Use a plane rotation to eliminate the damping parameter $\lambda$.
		// This alters the diagonal $\bar\rho$ of the lower-bidiagonal matrix.
		double cs1, sn1;
		double rhobar1 = rotation(cs1, sn1, rhobar, lambda);
		phibar         = cs1*phibar;
		
		// Use a plane rotation to eliminate the subdiagonal element $\beta$
		// of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
		double cs, sn;
		double rho   = rotation(cs, sn, rhobar1, beta);
		double theta =   sn*alpha;
		rhobar       = - cs*alpha;
		double phi   =   cs*phibar;
		phibar       =   sn*phibar;
		
		// Update $x$ and $w$.
		x = x + (phi  /rho) * w;
		w = v - (theta/rho) * w;
	}
	return true;
}

} /* namespace SLOM */
