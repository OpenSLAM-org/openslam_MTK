/**
 * @file slom/Algorithm.hpp
 * @brief Virtual base class implementing an optimization algorithm.
 * 
 * The Algorithm will be called from Estimator and shall only access 
 * the solver object.
 */


#ifndef ALGORITHM_HPP_
#define ALGORITHM_HPP_


#include "src/SparseFunction.hpp"
#include "src/SparseInterface.hpp"
#include "Solver.hpp"

namespace SLOM {

/**
 * Base class for all algorithms.
 */
class Algorithm : protected SparseInterface {
	//! Estimator is allowed to change internal variables
	friend class Estimator;
	friend class SparseFunction;
	
	//! The solver object. Not to be accessed by the algorithm directly.
	Solver* solver;
	
	//! solution of the linear LS problem
	Eigen::VectorXd linearLS;
	
protected:
	Algorithm() : solver(0) {}
	
	/**
	 * Optimize the underlying linear system.
	 * The result is a minimum of @f$\norm{J\cdot h - r}^2 + \lambda\norm h@f$.
	 * @todo variant with different penalty term
	 */
	const Eigen::VectorXd& solve(double lambda=0.0) {
		assert(solver);
		
		bool success = solver->solve(linearLS, lambda); (void)success;
		assert(success && "Solving linear system failed!");
		assert(linearLS.rows() == getN() && "Refinement vector has wrong dimensions!");
		
		return linearLS;
	}
	
	double applyDelta(const Eigen::VectorXd& delta, double scale){
		return func->apply_delta(delta, scale);
	}
	void store() {
		func->store_or_restore(true);
	}
	void restore() {
		func->store_or_restore(false);
	}
	
	
public:
	
	/**
	 * Method which needs to be implemented by derived algorithms.
	 * While this method is running, no changes of structure will happen.
	 * Variables will only change by calls to add_delta and store or restore.
	 * 
	 * @return rho value. @f$ \rho = \frac{\norm{f(x)}^2 - \norm{f(x_{\text{new}})}^2}{\norm(J\cdot (x-x_{\text{new}}))^2
	 */
	virtual double optimizeStep() = 0;
	
	/**
	 * Print information on the internal state.
	 * @param strm   output stream to which information is printed
	 * @param header if true, print headers, otherwise print values
	 * @return       the output stream
	 * 
	 * FIXME Maybe better let it push into a vector of doubles, however this assumes all internals are numbers.
	 */
	virtual std::ostream& info(std::ostream& strm, bool header = false) const { (void)header; return strm; }; 
	virtual ~Algorithm() {}
};

}  // namespace SLOM



#endif /* ALGORITHM_HPP_ */
