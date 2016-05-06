/**
 * @file LSMRSolver.cpp
 * @brief LSMR Solver
 * 
 * 
 * Ported from BSD-licensed Matlab code at:
 * http://www.mathworks.com/matlabcentral/fileexchange/27183-lsmr-an-iterative-algorithm-for-least-squares-problems
 * 
 * Original code was written by:
 * 
 * David Chin-lung Fong            clfong@stanford.edu
 * Institute for Computational and Mathematical Engineering
 * Stanford University
 * 
 * Michael Saunders                saunders@stanford.edu
 * Systems Optimization Laboratory
 * Dept of MS&E, Stanford University.
 * 
 * 
 * Shared code with LSQR solver has been outfactored to GolubKahanLanczos base class
 * 
 * 
 */

#include "LSMRSolver.hpp"

#include <stdio.h>
#include "../TicToc.hpp"

namespace SLOM {

using std::min;
using std::sqrt;

bool LSMRSolver::solve(VectorType &x, const double &lambda2) {
	int m = getM(), n=getN();   // Determine dimensions m and n
	int itnlim = min(min(m, n), this->max_it); // iteration limit
	
	double lambda = sqrt(lambda2); // LSMR works with $\sqrt{\lambda}$:
	
	init();  // Initialize GKL process
	
	// These variables are updated at each GKL step:
	const double &alpha = getAlpha(), &beta  = getBeta();
	const VectorType &v = getZ();   // For LSMR v==z
	
	// Initialize variables for 1st iteration.
	double alphabar = alpha, zetabar = alpha*beta, rho = 1, rhobar = 1, cbar = 1, sbar = 0;
	VectorType h = v, hbar = VectorType::Zero(n);
	x = VectorType::Zero(n);
	
	for(int k = 0; k < itnlim; ++k){
		step(); // After this call, beta = $\beta_{k+1}$, alpha = $\alpha_{k+1}$.
		
		// Construct rotation $\hat Q_{k,2k+1}$ to eliminate $\lambda$.
		double chat, shat;
		double alphahat = rotation(chat, shat, alphabar, lambda);
		
		// Use a plane rotation $Q_k$ to turn $B_k$ to $R_k$.
		double rhoold   = rho;
		double c, s;
		       rho      = rotation(c, s, alphahat, beta);
		double thetanew = s*alpha;
		       alphabar = c*alpha;
		
		// Use a plane rotation $\bar Q_i$ to turn $R_i\trans$ to $\bar R_i$.
		double rhobarold = rhobar;
		double thetabar  = sbar*rho;
		       rhobar    = rotation(cbar, sbar, cbar*rho, thetanew);
		double zeta      =  cbar*zetabar;
		       zetabar   = -sbar*zetabar;
		
		// Update $\bar h$, $x$, $h$.
		hbar = h - (thetabar*rho/(rhoold*rhobarold))*hbar;
		x    = x + (zeta/(rho*rhobar))*hbar;
		h    = v - (thetanew/rho)*h;
	}
	return true;
}

} /* namespace SLOM */
