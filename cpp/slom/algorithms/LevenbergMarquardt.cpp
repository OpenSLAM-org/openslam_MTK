/**
 * @file /mtk/slom/algorithms/LevenbergMarquardt.cpp
 * @brief Levenberg-Marquardt algorithm using parameter control as proposed by Hans Bruun Nielsen.
 * See http://www2.imm.dtu.dk/~hbn/Software/ for details.
 */

#include "LevenbergMarquardt.hpp"

namespace SLOM {

LevenbergMarquardt::LevenbergMarquardt(double lambda) : lambda(lambda), nu(2) {
}

double LevenbergMarquardt::optimizeStep() {
	
	typedef Eigen::VectorXd Vect;
	
	double oldRSS = getRSS();
	
	const Vect& grad  = getGradiant();
	const Vect& delta = solve(lambda);
	
	double newRSS = applyDelta(delta, -1);
	
	double rho = (oldRSS - newRSS) / delta.dot(lambda*delta + grad);
	
// FIXME make this optionally (i.e. depending on callBack object)
//	std::cerr << "LMA> " << oldRSS << " " << newRSS << " " << grad.norm() << " " << lambda << " " << rho ;//<< "\n";
	
	if(rho > 0) {
		store();
		lambda *= std::max(1.0/3, 1 - std::pow(2*rho - 1, 3));
		nu = 2;
	} else {
		restore();
		lambda *= nu; 
		nu *= 2;
	}
	if(lambda>1e100){
		lambda = 1e100;
		nu = 1;
	}
//	std::cerr << " " << lambda <<"\n";
	
	return rho;
}



} /* namespace SLOM */
