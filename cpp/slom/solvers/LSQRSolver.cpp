/**
 * @file LSQRSolver.cpp
 * @brief Refactored version of the LSQR solver from Paige and Saunders
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



bool LSQRSolver::solve(VectorType& x, const double & dampsq) {
	
	// Default parameters:
	int show = 0;
	double atol      = 1e-10;
	double btol      = 1e-7;
	double conlim    = 1e+8;
	
	int m = getM(), n=getN(); (void)m;
	
	int itnlim    = std::min(std::min(m, n), this->max_it);
	
	
	// Initialize.
	
	
	double damp = std::sqrt(dampsq);
	
//	bool wantvar   = false;
	//	if wantvar, var = zeros(n,1); end //TODO optional;
	
	const char *msg[]={
			"The exact solution is  x = 0                              ",
			"Ax - b is small enough, given atol, btol                  ",
			"The least-squares solution is good enough, given atol     ",
			"The estimate of cond(Abar) has exceeded conlim            ",
			"Ax - b is small enough for this machine                   ",
			"The least-squares solution is good enough for this machine",
			"Cond(Abar) seems to be too large for this machine         ",
			"The iteration limit has been reached                      ",};
	
	if(show > 10) {
		printf("\n\nLSQR            Least-squares solution of  Ax = b");
		printf("\nThe matrix A has %8d rows  and %8d cols", m,n);
		printf("\ndamp = %20.14e", damp);
		printf("\natol   = %8.2e               conlim = %8.2e", atol,conlim);
		printf("\nbtol   = %8.2e               itnlim = %8d"  , btol,itnlim);
	}
	
	int itn    = 0;            int istop  = 0;
	double ctol   = conlim > 0 ? 1/conlim : 0;
	double Anorm  = 0,             Acond  = 0;
	double ddnorm = 0,             res2   = 0;
	double xnorm  = 0,             xxnorm = 0;
	double cs2    = -1,            sn2    = 0;
	
	x.setZero(n);
	
	// Set up the first vectors u and v for the bidiagonalization.
	// These satisfy  beta*u = b,  alfa*v = A'u.
	
	VectorType oldV = getZ();
	
	init();
	const VectorType & v=getZ();
	if(show > 0 && oldV.size() == v.size()) {
		std::cout << "LSQR> Diff between old and new v: " << (v - oldV).norm() << std::endl;
	}
	const double &alfa = getAlpha();
	const double &beta = getBeta();
	double z      = 0;
	
	
	double Arnorm = alfa*beta;
	if(Arnorm == 0) {
		if(show > 2)
			std::cout << msg[0];
		return true;
	}
	
	if(recycle && v.size() == w.size() && nextGamma>0) {
		const SparseType &J = getJ();
		VectorType Jw = J*w;
		double gamma = v.dot(J.transpose() * Jw) / Jw.squaredNorm();
		if(show > 0) {
			std::cout << "LSQR> Recycling with gamma = min(" << gamma << ", " << nextGamma << ")\n";
		}
		gamma = std::min(gamma, nextGamma);
		if(gamma > 0)
			w = v - gamma * w;
		else
			w = v;
	} else {
		w = v;
	}
	
	
	double rhobar = alfa,          phibar = beta,          bnorm  = beta;
	double rnorm  = beta;
	double r1norm = rnorm;
	double r2norm = rnorm;
	const char head1[]  = "   Itn      x(1)       r1norm     r2norm ";
	const char head2[]  = " Compatible   LS      Norm A   Cond A";
	
	// Heading for iteration log.
	
	if(show > 3){
		double test1 = 1;
		double test2 = alfa/beta;
		printf("\n\n%s%s"      , head1 , head2   );
		printf("\n%6d %12.5e"  , itn  , x(1)   );
		printf(" %10.3e %10.3e", r1norm, r2norm );
		printf("  %8.1e %8.1e\n" , test1, test2  );
	}
	
	//------------------------------------------------------------------
	//     Main iteration loop.
	//------------------------------------------------------------------
	while(itn < itnlim){
		itn = itn + 1;
		
		// Perform the next step of the bidiagonalization to obtain the
		// next beta, u, alfa, v.  These satisfy the relations
		//      beta*u  =  A*v  - alfa*u,
		//      alfa*v  =  A'*u - beta*v.
		
		double lastalpha = alfa;
		step();
		
		
		Anorm = Eigen::Vector4d(Anorm, lastalpha, beta, damp).norm();
		
		//Use a plane rotation to eliminate the damping parameter.
		// This alters the diagonal (rhobar) of the lower-bidiagonal matrix.
		
		double cs1, sn1;
		double rhobar1 = rotation(cs1, sn1, rhobar, damp);
		double psi  = sn1*phibar;
		phibar      = cs1*phibar;
		
		// Use a plane rotation to eliminate the subdiagonal element (beta)
		// of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
		
		double cs, sn;
		double rho = rotation(cs, sn, rhobar1, beta);
		double theta = sn*alfa;
		rhobar     = - cs*alfa;
		double phi =   cs*phibar;
		phibar     =   sn*phibar;
		double tau =   sn*phi;
		
		// Update x and w.
		
		double t1  = phi  /rho;
		double t2  = theta/rho;
		ddnorm     = ddnorm + w.squaredNorm()/(rho*rho); // dk == (1/rho)*w;
		
		x  +=     t1*w;
		
		if(show > 5) {
			VectorType Jw = getJ() * w;
			double gamma = v.dot(getJ().transpose() * Jw) / Jw.squaredNorm();
			std::cout << "LSQR> gamma = " << gamma << ", t2 = " << t2 << "\n";
		}
		
		// update w after testing for stopping criteria
		
		
		//  if wantvar, var = var + dk.*dk; end
		
		// Use a plane rotation on the right to eliminate the
		// super-diagonal element (theta) of the upper-bidiagonal matrix.
		// Then use the result to estimate  norm(x).
		
		double delta  =   sn2*rho;
		double gambar = - cs2*rho;
		double rhs    =   phi - delta*z;
		double zbar   =   rhs/gambar;
		xnorm         =   std::sqrt(xxnorm + zbar*zbar);
		double gamma   =   rotation(cs2, sn2, gambar, theta);
		z       =   rhs   /gamma;
		xxnorm  =   xxnorm + z*z;
		
		// Test for convergence.
		// First, estimate the condition of the matrix  Abar,
		// and the norms of  rbar  and  Abar'rbar.
		
		Acond   =   Anorm*std::sqrt(ddnorm);
		double res1 =   phibar*phibar;
		res2    =   res2 + psi*psi;
		rnorm   =   std::sqrt(res1 + res2);
		Arnorm  =   alfa*std::abs(tau);
		
		// 07 Aug 2002:
		// Distinguish between
		//    r1norm = ||b - Ax|| and
		//    r2norm = rnorm in current code
		//           = sqrt(r1norm^2 + damp^2*||x||^2).
		//    Estimate r1norm from
		//    r1norm = sqrt(r2norm^2 - damp^2*||x||^2).
		// Although there is cancellation, it might be accurate enough.
		
		double r1sq    =   rnorm*rnorm - dampsq*xxnorm;
		r1norm  =   copysign(std::sqrt(std::abs(r1sq)), r1sq);
		r2norm  =   rnorm;
		
		// Now use these norms to estimate certain other quantities,
		// some of which will be small near a solution.
		
		double test1   =   rnorm /bnorm;
		double test2   =   Arnorm/(Anorm*rnorm);
		double test3   =   1/Acond;
		double test4   =   test1/(1 + Anorm*xnorm/bnorm);
		double rtol    =   btol + atol*Anorm*xnorm/bnorm;
		
		// The following tests guard against extremely small values of
		// atol, btol  or  ctol.  (The user may have set any or all of
		// the parameters  atol, btol, conlim  to 0.)
		// The effect is equivalent to the normal tests using
		// atol = eps,  btol = eps,  conlim = 1/eps.
		
		if(itn >= itnlim  ) istop = 7;
		if(1 + test3  <= 1) istop = 6;
		if(1 + test2  <= 1) istop = 5;
		if(1 + test4  <= 1) istop = 4;
		
		// Allow for tolerances set by the user.
		
		if(test3 <= ctol)  istop = 3;
		if(test2 <= atol)  istop = 2;
		if(test1 <= rtol)  istop = 1;
		
		// See if it is time to print something.
		
		if(show>2) {
			int prnt = 0;
			if(n     <= 40       ) prnt = 1;
			if(itn   <= 10       ) prnt = 1;
			if(itn   >= itnlim-10) prnt = 1;
			if(itn % 10 == 0     ) prnt = 1;
			if(test3 <=  2*ctol  ) prnt = 1;
			if(test2 <= 10*atol  ) prnt = 1;
			if(test1 <= 10*rtol  ) prnt = 1;
			if(istop !=  0       ) prnt = 1;
		
			if(prnt) {
				printf( "%6d %12.5e",        itn,   x(1) );
				printf( " %10.3e %10.3e", r1norm, r2norm );
				printf( "  %8.1e %8.1e",   test1,  test2 );
				printf( " %8.1e %8.1e %8.1e\n",    Anorm,  Acond, Arnorm );
				if(Arnorm < 1e-3){
					std::cout << Arnorm << "\t= " 
							<< phibar << " * " << alfa << " * |" << cs << "|\n\t= "
							<< alfa  << " * |" << tau << "|\n";
				}
			}
		}
		if(istop > 0) {
			if(show > 0 && recycle) {
				std::cout << "LSQR> next gamma would have been " << t2 << std::endl;
			}
			nextGamma = istop != 7 ? -1 : t2;
			break;
		}
		w   = v - t2*w;
	}
	
	// End of iteration loop.
	// Print the stopping condition.
	
	if(show){
		printf("\nlsqrSOL finished\n");
		printf("%s\n\n", msg[istop]);
		printf( "istop =%8d   r1norm =%8.1e  ",   istop, r1norm );
		printf( "Anorm =%8.1e   Arnorm =%8.1e\n", Anorm, Arnorm );
		printf( "itn   =%8d   r2norm =%8.1e  ",     itn, r2norm );
		printf( "Acond =%8.1e   xnorm  =%8.1e\n", Acond, xnorm  );
	}
	
	return true;
}


}
