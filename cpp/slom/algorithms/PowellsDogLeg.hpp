/**
 * @file /mtk/slom/algorithms/PowellsDogLeg.hpp
 * @brief Implementation of Powell's Dogleg
 * 
 */

#ifndef POWELLSDOGLEG_HPP_
#define POWELLSDOGLEG_HPP_

#include "../Algorithm.hpp"


namespace SLOM {

/*
 *
 */
class PowellsDogLeg : public Algorithm {
	//! Size of trust region
	double Delta;
	
public:
	PowellsDogLeg(const double& Delta = 1.0) : Delta(Delta) {};
	
	double optimizeStep();
	double getDelta() const { return Delta; }
	void setDelta(const double &Delta) { this->Delta = Delta; }

};

} /* namespace SLOM */
#endif /* POWELLSDOGLEG_HPP_ */
