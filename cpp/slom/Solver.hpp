/**
 * @file slom/Solver.hpp
 * @brief Base class for linear solvers.
 * 
 */

#ifndef SOLVER_HPP_
#define SOLVER_HPP_

#include "src/SparseInterface.hpp"


namespace SLOM {



/**
 * Virtual base class for solvers.
 * For 
 */
class Solver : protected SparseInterface {
	//! Estimator is allowed to modify internal variables
	friend class Estimator;
	friend class SparseFunction;
	
	
	
	
	//! Solution of linear least squares optimization
	VectorType delta;
	
	
	
	
protected:
	
	static bool isFinite(VectorType& x) {
		return true;
		// FIXME quite expensive test if used exhaustively
		return ((x-x).array()==(x-x).array()).all();
	}
	
	
public:
	virtual ~Solver() {}
	
	/**
	 * Sets @c delta to @f$\argmin_h \norm{Jh - r}^2 + \lambda h^TDh @f$.
	 * If solving fails due to numerical reasons, this function shall return @c false rather
	 * than throw or fail an assertion. 
	 * If something is structurally wrong, something is wrong on the caller's side and asserting is encouraged
	 * 
	 * @param delta  result vector, has to be resized properly by this function
	 * @param lambda damping term
	 * @return       true on success, false otherwise.
	 */
	virtual bool solve(VectorType& delta, const double & lambda) = 0;
	
	
	/**
	 * Amount of memory required by the solver. Usually, it's sufficient to
	 * return an approximation, it should however include expressions taken 
	 * from SparseInterface. 
	 * Amount to store Jacobian and residuum is shall not be included (unless a copy is made)
	 * @return amount of memory required for solving.
	 */
	virtual size_t memory() const = 0; // { return 0; };
	
};

}  // namespace SLOM



#endif /* SOLVER_HPP_ */
