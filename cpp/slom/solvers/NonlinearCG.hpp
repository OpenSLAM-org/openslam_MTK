/**
 * @file NonlinearCG.hpp
 * @brief Brief description
 * 
 */

#ifndef NONLINEARCG_HPP_
#define NONLINEARCG_HPP_

#include "PreConditionedSolver.hpp"

namespace SLOM {

/**
 * Non-linear Conjugate Gradient method
 */
class NonlinearCG : public PreConditionedSolver {
	
	VectorType precGrad, lastGradient;
	double normLastGradient;
	
	VectorType cgDirection;
	
	
public:
	NonlinearCG();
	virtual ~NonlinearCG() {};
	virtual bool solve(Eigen::VectorXd& delta, const double & lambda);
};

} /* namespace SLOM */
#endif /* NONLINEARCG_HPP_ */
