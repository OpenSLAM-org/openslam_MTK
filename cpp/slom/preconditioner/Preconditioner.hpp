/**
 * @file /mtk/slom/preconditioner/Preconditioner.hpp
 * @brief Brief description
 * 
 */

#ifndef PRECONDITIONER_HPP_
#define PRECONDITIONER_HPP_

#include "../src/SparseFunction.hpp"

#include "../src/SparseInterface.hpp"


namespace SLOM {

/** 
 * A preconditioner simulates a matrix @f$ K @f$, such that @f$ J K @f$  has a better condition than @f$ J @f$.
 * Ideally  @f$ J K @f$  is (near) orthogonal, which is the case if  @f$ K\trans K\approx J\trans J @f$.
 *
 */
class Preconditioner : protected SparseInterface {
	friend class Estimator;
	friend class SparseFunction;
	
public:
	virtual ~Preconditioner() {}
	
	/**
	 * @c compute is called before the preconditioner is applied
	 */
	virtual void compute() {}
	
	/**
	 * Multiply @c vec by @f$ K\inv @f$ or @f$ K\trans[-] @f$.
	 * @param vec The vector to which @f$ K\inv @f$ shall be applied.
	 * @param transpose flag which decides whether @f$ K @f$ shall be transposed.
	 */
	virtual void apply(Eigen::VectorXd& vec, bool transpose) const = 0;
};

/**
 * @c JacobiPreconditioner divide each entry by the square root of @f$ (J\trans J)_{ii} @f$.
 */
class JacobiPreconditioner : public Preconditioner {
	virtual void apply(Eigen::VectorXd& vec, bool /*transpose*/ ) const {
		const Eigen::VectorXd& chol = getCholCov();
		vec.array() /= chol.array();
	}
};

/**
 * Similar to @c JacobiPreconditioner but uses a block for each variable multiplies 
 * by the inverse of the corresponding Cholesky factor.
 */
class BlockJacobiPreconditioner : public Preconditioner {
	virtual void apply(Eigen::VectorXd& vec, bool transpose) const {
		applyJacobiPreconditioner(vec, transpose);
	}
};


class CholIncPreconditioner : public Preconditioner {
	virtual void apply(Eigen::VectorXd& vec, bool transpose) const {
		assert(false && "currently not implemented, also transpose flags needs to be checked");
		applyCholIncPreconditioner(vec, transpose);
	}
};


} /* namespace SLOM */
#endif /* PRECONDITIONER_HPP_ */
