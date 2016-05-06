/**
 * @file CholeskySolver.hpp
 * @brief Brief description
 * 
 */

#ifndef CHOLESKYSOLVER_HPP_
#define CHOLESKYSOLVER_HPP_

#include "../Solver.hpp"

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

/*
 *
 */
namespace SLOM {


/**
 * Solves (J'*J + lambda*I) x = J'*b using sparse Cholesky decomposition
 * @tparam CholeskyDecomposer right now can be Eigen::SimplicialLDLT, Eigen::SimplicialLLT or Eigen::CholmodDecomposition
 */
template<class CholeskyDecomposer = Eigen::SimplicialLDLT<SparseInterface::SparseType, Eigen::Upper> >
class CholeskySolver : public Solver{
	
	typedef CholeskyDecomposer Decomposer;

	Decomposer cholesky;
	bool initialized;

public:
	CholeskySolver() : initialized(false) {}
	
	template<class Mode>
	CholeskySolver(Mode m) : initialized(false) {
		cholesky.setMode(m);
	}
	
	
	/**
	 * Implements Solver::solve
	 * res is ignored
	 */
	bool solve(VectorType& delta , const double & lambda) {
//		int m = getM(), n=getN(); (void)m;
		
		bool structureChanged = hasStructureChanged();
		const SparseType& JtJ = get_JtJ();

		
		if(structureChanged){
//			std::cerr << "Analyzing Pattern\n";
			cholesky.analyzePattern(JtJ);
			//cholesky.setMode(Eigen::CholmodSupernodalLLt);
			//cholesky.setMode(Eigen::CholmodLDLt);
		}
		cholesky.setShift(lambda).factorize(JtJ);
		initialized = true;
		if(cholesky.info() != Eigen::Success) return false;
		delta = cholesky.solve(getGradiant());
		
		return true;
	}
	
	size_t memory() const {
		size_t size = 0;
		size += (get_JtJ().nonZeros() + 0* getJ().nonZeros());
		if(initialized){// FIXME does not seem to work
			SparseType L = cholesky.matrixL();
			//size +=cholesky.matrixL().nonZeros();
//			size += get_JtJ().nonZeros();
			size+=L.nonZeros();
		}
		size*= (sizeof(SparseType::Scalar) + sizeof(SparseType::Index));
		size += 2*sizeof(SparseType::Index) * getN();
		size += sizeof(SparseType::Scalar) * (getN() + getM());
		return size; 
	}
};




}

#endif /* CHOLESKYSOLVER_HPP_ */
