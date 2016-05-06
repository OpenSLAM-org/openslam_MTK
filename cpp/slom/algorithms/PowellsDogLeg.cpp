/**
 * @file /mtk/slom/algorithms/PowellsDogLeg.cpp
 * @brief Brief description
 * 
 */

#include "PowellsDogLeg.hpp"

namespace SLOM {

double PowellsDogLeg::optimizeStep() {
	
	typedef Eigen::VectorXd Vect;
	
	double oldRSS = getRSS();
	double newRSS;
	double rho;
	
	
	const Vect& grad  = getGradiant();
	// solution of Gauss-Newton
	const Vect& h_GN = solve();
	double h_GN2 = h_GN.squaredNorm();
	
	double linearizedGain;
	
	double norm_update;
	
	
	if(h_GN2 <= Delta * Delta){
		// Gauss-Newton update fits into trust region:
		newRSS = applyDelta(h_GN, -1);
		linearizedGain = oldRSS;
		norm_update = std::sqrt(h_GN2);
	} else {
		double g2 = getSquaredNormOfGradient();
		// Scale for Steepest descent:
		double alpha = g2 / (getJ() * grad).squaredNorm();
		if(alpha * alpha * g2 >= Delta*Delta) {
			// Steepest descent vector is bigger than trust region so apply capped vector
			newRSS = applyDelta(grad, -Delta / std::sqrt(g2));
			linearizedGain = Delta * (2*alpha * std::sqrt(g2) - Delta) / (2*alpha);
			norm_update = Delta;
		} else {
			// Dogleg-step, find linear combination such that -alpha*grad + beta*(h_GN + alpha*grad) has norm Delta. 
			double dot_grad_GN = grad.dot(h_GN);
			// distance between Gauss-Newton and steepest descent vector:
			double dist2_GN_grad = h_GN2 - 2*alpha*dot_grad_GN + alpha*alpha * g2;
			
			double p = (alpha * dot_grad_GN - alpha*alpha*g2)/dist2_GN_grad;
			double q = (Delta*Delta - alpha*alpha * g2)/dist2_GN_grad;  // guaranteed to be positive, because steepest descent vector is smaller than Delta
			
			double beta = std::sqrt(q+p*p) - p;
			assert(beta == beta && 0 <= beta && beta <= 1 && "beta in dogleg step should lay within [0, 1]");
			
			Eigen::VectorXd delta = beta * h_GN + (1-beta)*alpha * grad;
			newRSS = applyDelta(delta, -1);
			linearizedGain = 0.5*alpha*(1-beta)*(1-beta)*g2+ beta*(2-beta)*oldRSS;
			norm_update = delta.norm();
		}
	}
	
	rho = (oldRSS - newRSS) / linearizedGain;
	std::cerr << "Dogleg: " << oldRSS << " " << newRSS << "; rho= " << rho << ", Delta = " << Delta << "\n";
	
	if(rho > 0){
		store();
	} else {
		// FIXME we could directly re-use h_GN, grad, etc and start a new attempt with reduced Delta
		restore();
	}
	if(rho > 0.75) {
		Delta = std::max(Delta, 3*norm_update);
	} else if(rho < 0.25){
		Delta /= 2;
	}
	
	return rho;
}


} /* namespace SLOM */
