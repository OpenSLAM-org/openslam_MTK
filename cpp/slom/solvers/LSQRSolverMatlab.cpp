/**
 * @file LSQRSolver.cpp
 * @brief Refactored version of Matlab's LSQR solver
 * 
 */

#include "LSQRSolverMatlab.hpp"

#include <cmath>

namespace SLOM {



bool LSQRSolverMatlab::solve(VectorType& x, const double & lambda) {
	
	// FIXME FATAL lambda is ignored!
	(void) lambda;
	
	int verbosity = 0;
	const double eps = 2.2e-16; // TODO get from limits.h or something
	//TODO make parameters class attributes
	//! Tolerance
	double tol = 1e-8;
	
	// Additional output parameters:
	int flag;
	double relres;
	int iter;
	std::vector<double> resvec, lsvec;
	
	// ** Init GKL ** //
	init();
	// ************** //
	
	int m = getM(), n=getN(); (void)m;
	
	int itnlim = std::min(std::min(n, m), max_it);
	
	
	double beta = getBeta();
	x.setZero(n);
	
	double n2b = beta; // Norm of rhs vector, b
	flag = 1;
	double tolb = tol * n2b;                  // Relative tolerance
	// Norm of residual r=b-A*x is estimated well by prod_i abs(sin_i)
	double normr = beta;
	if(beta == 0)
	{
		std::cerr << "Beta==0, " << n2b << " " << normr << std::endl;
	}
	
	double c = 1, s = 0;
	double phibar = beta;
	VectorType d_times_rho;// = Vector::Zero(n);
	
	double alpha = getAlpha();
	
	// norm((A*inv(M))'*r) = alpha_i * abs(sin_i * phi_i)
	double normar = alpha * beta;
	
	
	
	
	
	// Check for all zero solution
	if (normar == 0) {            // if alpha_1 == 0 | beta_1 == 0
		x.setZero(n);              // then  solution is all zeros
		flag = 0;                  // a valid solution has been obtained
		relres = 0;                // the relative residual is actually 0/0
		iter = 0;                  // no iterations need be performed
		resvec.push_back(beta);    // resvec(1) = norm(b-A*x) = norm(0)
		lsvec.clear();             // no estimate for norm(A*inv(M),'fro') yet
		return true;
	}
	
	
	
	// Poorly estimate norm(A*inv(M),'fro') by norm(B_{ii+1,ii},'fro')
	// which is in turn estimated very well by
	// sqrt(sum_i (alpha_i^2 + beta_{ii+1}^2))
	double norma = 0;
	// norm(inv(A*inv(M)),'fro') = norm(D,'fro')
	// which is poorly estimated by sqrt(sum_i norm(d_i)^2)
	double sumnormd2 = 0;
	resvec.reserve(itnlim+1);       // Preallocate vector for norm of residuals
	resvec.push_back(normr);       // resvec(1,1) = norm(b-A*x0)
	lsvec.reserve(itnlim);          // Preallocate vector for least squares estimates
	int stag = 0;                  // stagnation of the method
	iter = itnlim;                  // Assume lack of convergence until it happens
	int maxstagsteps = 3;
	
	double rho = 1;
	
	const VectorType& z = getZ(); // declare z before loop starts to avoid re-allocalization
	d_times_rho = z;
	
	double thet = -s * alpha;
	double rhot = c * alpha;
	
	for(int ii = 0 ; ii< itnlim; ++ii) {
		
		
		// ** GKL step ** //
		step();
		// ************** //
		
		beta = getBeta();
		
		norma = Eigen::Vector3d(norma, alpha, beta).norm();
		lsvec.push_back(normar / norma);
		rho = Eigen::Vector2d(rhot, beta).norm();
		c = rhot / rho;
		s = - beta / rho;
//		std::cerr << "rho: " << rho << "c: " << c << ", s: " << s << "\n";
		double phi = c * phibar;
		if (phi == 0)              // stagnation of the method
			stag = 1;
		phibar = s * phibar;
		sumnormd2 +=  d_times_rho.squaredNorm() / (rho*rho);
		
		// Check for stagnation of the method
		if(std::abs(phi)*d_times_rho.norm() < eps*x.norm() * rho)
			++stag;
		else
			stag = 0;
		
		if(normar/(norma*normr) <= tol){ // check for convergence in min{|b-A*x|}
			if(verbosity > 1) {
				std::cout << "normar/(norma*normr): " << normar << "/(" << norma << "*" << normr << ")= " 
						<< normar/(norma*normr) << " < " << tol << std::endl;
			}
			flag = 100;
			iter = ii;
			break;
		}
		if(normr <= tolb){         // check for convergence in A*x=b
			if(verbosity > 1)
				std::cout << "normr: " << normr << " < " << tolb << std::endl; 
			flag = 101;
			iter = ii;
			break;
		}
		
		if(stag >= maxstagsteps){
			flag = 3;
			iter = ii;
			break;
		}
		
		x += phi * d_times_rho / rho;   // CHANGED: replaced d by  d_times_rho/rho.
		normr = std::abs(s) * normr;
		resvec.push_back(normr);
		
		alpha = getAlpha();
		
		normar = alpha * std::abs( s * phi);
		
		thet = - s * alpha;
		rhot = c * alpha;
		
		d_times_rho = (z - thet * d_times_rho/rho);

		
		
	} // for ii = 1 : maxit
	
	if(flag == 1){
		if(normar/(norma*normr) <= tol){ // check for convergence in min{|b-A*x|}
			flag = 0;
			iter = itnlim;
		}
		
		if(normr <= tolb){          // check for convergence in A*x=b
			flag = 0;
			iter = itnlim;
		}
	}

	relres = normr/n2b;

	if(verbosity > 0) {
		std::cout << "LSQR: " << tol << " " << itnlim << " " << iter << " " << flag << " " << relres << std::endl;
	}
	
	
	return true;
}


}
