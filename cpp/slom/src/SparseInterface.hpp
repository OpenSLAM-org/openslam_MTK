/**
 * @file SparseInterface.hpp
 * @brief Brief description
 * 
 */


#ifndef SPARSEINTERFACE_HPP_
#define SPARSEINTERFACE_HPP_

// FIXME Avoid this include
//#include "SparseFunction.hpp"


#include "Sparse.hpp"


namespace SLOM {


class SparseFunction;
class CallBack;

// FIXME make this internal?
/**
 * SparseInterface is the common base class for classes which need read-access to the SparseFunction
 */
class SparseInterface {
	
	// Estimator is allowed to set @c func variable
	friend class Estimator;
	friend class SparseFunction;
	friend class Algorithm;
	
	//! Pointer to SparseFunction, not to be accessed directly by derived algorithms.
	SparseFunction *func;
	
	//! did the structure change since last call to hasStructureChanged();
	bool structure_changed, values_changed;
	
	void valuesChanged() {
		values_changed = true;
	}
	void structureChanged() {
		structure_changed = true;
	}
	
public:
	typedef SLOM::internal::SparseType SparseType;
	typedef Eigen::VectorXd VectorType;
	
protected:
	
	SparseInterface() : func(0) {
		structureChanged();
	};
	
	/**
	 * Virtual destructor, because all derived classes will be virtual.
	 */
	virtual ~SparseInterface() { }
	/**
	 * Determine whether the structure changed since the last call.
	 * If the structure did not change, things such as symbolic decompositions can be reused.
	 */
	bool hasStructureChanged() {
		bool ret = structure_changed;
		structure_changed = false;
		return ret;
	}
	/**
	 * Determine whether the structure changed since the last call.
	 * If the structure did not change, things such as symbolic decompositions can be reused.
	 */
	bool haveValuesChanged() {
		bool ret = values_changed;
		structure_changed = values_changed = false;
		return ret;
	}
	
	CallBack* getCallback() const;
	
	//! \name Access methods to be used by the Solver
	//@{
	//! Get the Jacobian
	const SparseType& getJ() const;
	
	
	/** @deprecated */
	//MTK_DEPRECATED
	const VectorType& getCholCov() const;
	
	void applyPreconditioner(VectorType &vec, bool transpose) const;
	
	
	/**
	 * 
	 * @deprecated Call applyPreconditioner
	 */
	void applyJacobiPreconditioner(VectorType & vec, bool transpose) const;
	
	/**
	 * 
	 * @deprecated Call applyPreconditioner
	 */
	void applyCholIncPreconditioner(VectorType & vec, bool transpose) const;
	
	
	/**
	 * res += J * v
	 */
	void add_Jv(VectorType& res, const VectorType& v) const;
	/**
	 * res += J' * u
	 */
	void add_Jtu(VectorType& res, const VectorType& u) const;
	
	/**
	 * Returns upper half of JtJ. First entry per column corresponds to diagonal entry.
	 * Although it is not good style, one can modify the diagonal entries, 
	 * e.g. to perform Levenberg-style optimizations.
	 * If the structure of the problem changed since last call structure_changed is set to true.
	 */
	const SparseType& get_JtJ() const;
	
	std::pair<int,int> getSize() const;
	
	int getM() const;
	int getN() const;
	
	double getRSS() const;
	
	const VectorType& getResiduum() const;
	
	const VectorType& getGradiant() const;
	double getSquaredNormOfGradient() const;
	
	//@}
};


}  // namespace SLOM


#endif /* SPARSEINTERFACE_HPP_ */
