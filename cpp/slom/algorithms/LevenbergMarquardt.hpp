/**
 * @file /mtk/slom/algorithms/LevenbergMarquardt.hpp
 * @brief Levenberg-Marquardt algorithm using parameter control as proposed by Hans Bruun Nielsen.
 * See http://www2.imm.dtu.dk/~hbn/Software/ for details.
 */

#ifndef LEVENBERGMARQUARDT_HPP_
#define LEVENBERGMARQUARDT_HPP_

#include "../Algorithm.hpp"

namespace SLOM {

/*
 *
 */
class LevenbergMarquardt : public Algorithm {
	double lambda, nu;
	
public:
	LevenbergMarquardt(double lambda = 1e-3);
	
	double optimizeStep();
	double getLambda() const { return lambda; }
	void setLambda(const double &lambda) { this->lambda = lambda; }
};

} /* namespace SLOM */
#endif /* LEVENBERGMARQUARDT_HPP_ */
