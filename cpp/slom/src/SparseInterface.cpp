/**
 * @file /mtk-trunk/slom/src/SparseInterface.cpp
 * @brief Brief description
 * 
 */



#include "SparseInterface.hpp"
#include "SparseFunction.hpp"

namespace SLOM {



CallBack* SparseInterface::getCallback() const {
	return func->getCallback();
}

//! \name Access methods to be used by the Solver
//@{
//! Get the Jacobian
const SparseInterface::SparseType& SparseInterface::getJ() const {
	return func->getJ();
}


/** @deprecated */
//MTK_DEPRECATED
const SparseInterface::VectorType& SparseInterface::getCholCov() const {
	return func->getCholCov();
}

void SparseInterface::applyPreconditioner(VectorType &vec, bool transpose) const {
	func->applyPreconditioner(vec, transpose);
}



/**
 * 
 * @deprecated Call applyPreconditioner
 */
void SparseInterface::applyJacobiPreconditioner(VectorType & vec, bool transpose) const {
	func->applyJacobiPreconditioner(vec, transpose);
}

/**
 * 
 * @deprecated Call applyPreconditioner
 */
void SparseInterface::applyCholIncPreconditioner(VectorType & vec, bool transpose) const {
	(void) vec; (void) transpose;
	assert(false);
	// FIXME re-implement Incomplete Cholesky
	//		func->applyCholIncPreconditioner(vec, transpose);
}


/**
 * res += J * v
 */
void SparseInterface::add_Jv(VectorType& res, const VectorType& v) const {
	func->add_Jv(res, v);
}
/**
 * res += J' * u
 */
void SparseInterface::add_Jtu(VectorType& res, const VectorType& u) const {
	func->add_Jtu(res, u);
}

/**
 * Returns upper half of JtJ. First entry per column corresponds to diagonal entry.
 * Although it is not good style, one can modify the diagonal entries, 
 * e.g. to perform Levenberg-style optimizations.
 * If the structure of the problem changed since last call structure_changed is set to true.
 */
const SparseInterface::SparseType& SparseInterface::get_JtJ() const {
	return func->get_JtJ();
}

std::pair<int,int> SparseInterface::getSize() const {
	return func->getSize();
}

int SparseInterface::getM() const {
	return func->getM();
}
int SparseInterface::getN() const {
	return func->getN();
}

double SparseInterface::getRSS() const {
	return func->getRSS();
}

const SparseInterface::VectorType& SparseInterface::getResiduum() const {
	return func->getResiduum();
}

const SparseInterface::VectorType& SparseInterface::getGradiant() const {
	return func->getGradient();
}
double SparseInterface::getSquaredNormOfGradient() const {
	return func->getSquaredNormOfGradient();
}

}  // namespace SLOM


