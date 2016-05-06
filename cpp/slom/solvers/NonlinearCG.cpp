/**
 * @file NonlinearCG.cpp
 * @brief Brief description
 * 
 */

#include "NonlinearCG.hpp"

namespace SLOM {

NonlinearCG::NonlinearCG() : normLastGradient (-1) {
	
}


bool NonlinearCG::solve(VectorType& delta, const double & lambda){
	
	
	// gradient and lastGradient are the preconditionend gradients.
	// The current cgDirection is stored in an unpreconditionend manner.
	
	// gradiant = M^-T * J^T * res \approx Q^T * res;
//	const Vec& grad = getGradiant();
	precGrad = getGradiant();
	applyPrecond(precGrad, precGrad, true);
	applyPrecond(precGrad, precGrad, false);
	const VectorType grad = precGrad;
	
	double gradNorm2 = precGrad.squaredNorm();
	
	if(normLastGradient < 0 || !true){
		// first cgDirection is just steepest descent:
		cgDirection = grad;
		std::cout << "beta = 0 -- RESTART --";
	} else {
		// calculate beta from current and last gradient:
		double dot = precGrad.dot(lastGradient);
		double beta = (gradNorm2 - dot)/normLastGradient;
//		double beta = (gradNorm2 - dot)/(dot - normLastGradient);
//		double beta = (gradNorm2)/normLastGradient;(void)dot;
		std::cout << "beta = " << beta;
		if(beta < 0.0){ // reset cgDirection
			std::cout << " -- RESTART --";
			beta = 0.0;
		}
		cgDirection = grad + beta * cgDirection;
	}
	
	// perform line search:
	Eigen::VectorXd v = getJ() * cgDirection;
	
	double alpha = v.dot(res)/(v.squaredNorm()+lambda*cgDirection.squaredNorm());
	
	std::cout << ", alpha = " << alpha << std::endl;
	delta = alpha * cgDirection;
	
	lastGradient.swap(precGrad);
	normLastGradient = gradNorm2;
	
	return true;
}


} /* namespace SLOM */
