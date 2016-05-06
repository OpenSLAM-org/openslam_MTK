/**
 * @file GaussNewton.cpp
 * @brief Brief description
 * 
 */

#include "GaussNewton.hpp"

namespace SLOM {

double GaussNewtonAlgorithm::optimizeStep() {
	
	double oldRSS = getRSS();
	const Eigen::VectorXd& delta = solve(0.0);
	double newRSS = applyDelta(delta, -1);
	store();
	
	return (oldRSS - newRSS) / oldRSS;
}


} /* namespace SLOM */
