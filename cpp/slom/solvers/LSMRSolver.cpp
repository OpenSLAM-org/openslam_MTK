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


bool LSMRSolver::solve(VectorType &x, const double &lambda2) {
	//if nargin < 3 || isempty(lambda)   , lambda    = 0;          end
	double atol      = 1e-8;
	double btol      = 1e-6;
	double conlim    = 1e+8;
	
	// Set default parameters.
	int m = getM(), n=getN(); (void)m;
	
	int itnlim    = std::min(std::min(m, n), this->max_it);
	// LSMR works with $\sqrt{\lambda}$:
	double lambda = sqrt(lambda2);

	
	int show      = 0;
	
	const char hdg1[] = "   itn      x(1)       norm r    norm A'r";
	const char hdg2[] = " compatible   LS      norm A   cond A";
	int pfreq  = 200000;   // print frequency (for repeating the heading)
	int pcount = 0;    // print counter
	
	x.setZero(n);
	
	// Determine dimensions m and n, and
	// form the first vectors u and v.
	// These satisfy  beta*u = b,  alpha*v = A'u.
	
	init();
	
	// variables are updated at each GKL step
	const double &beta = getBeta();
	const double &alpha = getAlpha();
	const Eigen::VectorXd &z = getZ();
	
	if(show > 10) {
		printf("\n\nLSMR            Least-squares solution of  Ax = b");
		printf("\nVersion 1.11                          09 Jun 2010");
		printf("\nThe matrix A has %8d rows  and %8d cols", m,n);
		printf("\nlambda = %16.10e", lambda );
		printf("\natol   = %8.2e               conlim = %8.2e", atol,conlim);
		printf("\nbtol   = %8.2e               itnlim = %8d"  , btol,itnlim);
	}
	
	
	
	// Initialization for local reorthogonalization.
	
	
	// Initialize variables for 1st iteration.
	
	int itn      = 0;
	double zetabar  = alpha*beta;
	double alphabar = alpha;
	double rho      = 1;
	double rhobar   = 1;
	double cbar     = 1;
	double sbar     = 0;
	
	Eigen::VectorXd h    = z;
	Eigen::VectorXd hbar = Eigen::VectorXd::Zero(n); //zeros(n,1);
	
	// Initialize variables for estimation of ||r||.
	
	double betadd      = beta;
	double betad       = 0;
	double rhodold     = 1;
	double tautildeold = 0;
	double thetatilde  = 0;
	double zeta        = 0;
	double d           = 0;
	
	// Initialize variables for estimation of ||A|| and cond(A).
	
	double normA2  = alpha * alpha;
	double maxrbar = 0;
	double minrbar = 1e+100;
	
	double normA=0, condA=0, normx=0;
	
	// Items for use in stopping rules.
	double normb  = beta;
	int istop  = 0;
	double ctol   = (conlim > 0) ? 1/conlim : 0;
	double normr  = beta;
	
	// Exit if b=0 or A'b = 0.
	
	double normAr = alpha * beta;
	if(normAr == 0) {
		std::cerr << "The exact solution is  x = 0";
		return true;
	}
	
	// Heading for iteration log.
	
	if(show > 3){
		double test1 = 1;
		double test2 = alpha/beta;
		printf("\n\n%s%s"      , hdg1 , hdg2   );
		printf("\n%6d %12.5e"  , itn  , x(1)   );
		printf(" %10.3e %10.3e", normr, normAr );
		printf("  %8.1e %8.1e" , test1, test2  );
	}
	
	//------------------------------------------------------------------
	//     Main iteration loop.
	//------------------------------------------------------------------
	while(itn < itnlim){
		itn++;
		
		// Perform the next step of the bidiagonalization to obtain the
		// next beta, u, alpha, v.  These satisfy the relations
		//      beta*u  =  A*v  - alpha*u,
		//      alpha*v  =  A'*u - beta*v.
		
		step();
		
		
		// calculating A * z and A' * u could be done with one sweep over the measurements
		// At this point, beta = beta_{k+1}, alpha = alpha_{k+1}.
		
		// Construct rotation Qhat_{k,2k+1}.
		
		double chat, shat;
		double alphahat = rotation(chat, shat, alphabar, lambda);
		
		// Use a plane rotation (Q_i) to turn B_i to R_i.
		
		double rhoold   = rho;
		double c, s;
		rho      = rotation(c, s, alphahat, beta);
		double thetanew = s*alpha;
		alphabar = c*alpha;
		
		// Use a plane rotation (Qbar_i) to turn R_i^T to R_i^bar.
		
		double rhobarold = rhobar;
		double zetaold   = zeta;
		double thetabar  = sbar*rho;
		double rhotemp   = cbar*rho;
		rhobar    = rotation(cbar, sbar, cbar*rho, thetanew);
		zeta      =   cbar*zetabar;
		zetabar   = - sbar*zetabar;
		// Update h, h_hat, x.
		
		hbar      = h - (thetabar*rho/(rhoold*rhobarold))*hbar;
		x         = x + (zeta/(rho*rhobar))*hbar;
		h         = z - (thetanew/rho)*h;
		// Estimate of ||r||.
		
		// Apply rotation Qhat_{k,2k+1}.
		double betaacute =   chat* betadd;
		double betacheck = - shat* betadd;
		
		// Apply rotation Q_{k,k+1}.
		double betahat   =   c*betaacute;
		betadd    = - s*betaacute;
		
		// Apply rotation Qtilde_{k-1}.
		// betad = betad_{k-1} here.
		
		double thetatildeold = thetatilde;
		double ctildeold, stildeold;
		double rhotildeold   = rotation(ctildeold, stildeold, rhodold, thetabar);
		thetatilde    = stildeold* rhobar;
		rhodold       =   ctildeold* rhobar;
		betad         = - stildeold*betad + ctildeold*betahat;
		
		// betad   = betad_k here.
		// rhodold = rhod_k  here.
		
		tautildeold   = (zetaold - thetatildeold*tautildeold)/rhotildeold;
		double taud          = (zeta - thetatilde*tautildeold)/rhodold;
		d             = d + std::pow(betacheck,2);
		normr         = std::sqrt(d + std::pow(betad - taud,2) + std::pow(betadd,2));
		
		// Estimate ||A||.
		normA2        = normA2 + std::pow(beta,2);
		normA         = std::sqrt(normA2);
		normA2        = normA2 + std::pow(alpha,2);
		
		// Estimate cond(A).
		maxrbar       = std::max(maxrbar,rhobarold);
		if (itn>1) 
			minrbar     = std::min(minrbar,rhobarold);
		condA         = std::max(maxrbar,rhotemp)/std::min(minrbar,rhotemp);
		
		// Test for convergence.
		
		// Compute norms for convergence testing.
		normAr  = std::abs(zetabar);
		normx   = x.norm();
		
		// Now use these norms to estimate certain other quantities,
		// some of which will be small near a solution.
		
		double test1   = normr /normb;
		double test2   = normAr/(normA*normr);
		double test3   =      1/condA;
		double t1      =  test1/(1 + normA*normx/normb);
		double rtol    = btol + atol*normA*normx/normb;
		
		// The following tests guard against extremely small values of
		// atol, btol or ctol.  (The user may have set any or all of
		// the parameters atol, btol, conlim  to 0.)
		// The effect is equivalent to the normAl tests using
		// atol = eps,  btol = eps,  conlim = 1/eps.
		
		if(itn >= itnlim )  istop = 7;
		if(1 + test3  <= 1) istop = 6;
		if(1 + test2  <= 1) istop = 5;
		if(1 + t1     <= 1) istop = 4;
		
		// Allow for tolerances set by the user.
		
		if(test3 <= ctol)  istop = 3;
		if(test2 <= atol)  istop = 2;
		if(test1 <= rtol)  istop = 1;
		
		// See if it is time to print something.
		
		if(show > 3) {
			bool prnt = false;
			if(n     <= 40       ) prnt = true;
			if(itn   <= 10       ) prnt = true;
			if(itn   >= itnlim-10) prnt = true;
			if(itn % 20 == 0     ) prnt = true;
			if(test3 <= 1.1*ctol ) prnt = true;
			if(test2 <= 1.1*atol ) prnt = true;
			if(test1 <= 1.1*rtol ) prnt = true;
			if(istop !=  0       ) prnt = true;
			
			if(prnt){
				if(pcount >= pfreq){
					pcount = 0;
					printf("\n\n%s%s"    , hdg1 , hdg2  );
				}  
				pcount = pcount + 1;
				printf("\n%6d %12.5e"  , itn  , x(1)  );
				printf(" %10.3e %10.3e", normr, normAr);
				printf("  %8.1e %8.1e" , test1, test2 );
				printf(" %8.1e %8.1e"  , normA, condA );
			}
		}
		
		if(istop > 0) break;
	} // iteration loop
	
	// Print the stopping condition.
	
	if(show){
		printf("\nLSMR finished");
//		printf("\n%d", istop);
		printf("\nistop =%8d    normr =%8.1e"     , istop, normr );
		printf("    normA =%8.1e    normAr =%8.1e", normA, normAr);
		printf("\nitn   =%8d    condA =%8.1e"     , itn  , condA );
		printf("    normx =%8.1e\n\n", normx);
	}
	//atol = std::min(atol, normx*normx);
	
	// end function lsmr
	return true;
}

} /* namespace SLOM */
