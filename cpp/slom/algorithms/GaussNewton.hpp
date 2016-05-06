/**
 * @file GaussNewton.hpp
 * @brief Brief description
 * 
 */

#ifndef GAUSSNEWTON_HPP_
#define GAUSSNEWTON_HPP_

#include "../Algorithm.hpp"

namespace SLOM {

/**
 * Very simple Gauss-Newton solver
 */
class GaussNewtonAlgorithm : public Algorithm {
public:
	double optimizeStep();
};

} /* namespace SLOM */
#endif /* GAUSSNEWTON_HPP_ */
