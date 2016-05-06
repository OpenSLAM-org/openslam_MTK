/**
 * @file /mtk/slom/preconditioner/BlockJacobiPreconditioner.hpp
 * @brief Brief description
 * 
 */

#ifndef BLOCKJACOBIPRECONDITIONER_HPP_
#define BLOCKJACOBIPRECONDITIONER_HPP_

#include "Preconditioner.hpp"

#include <Eigen/Cholesky>

namespace SLOM {

/*
 *
 */
class BlockJacobiPreconditioner : public Preconditioner {
	
	template<int DOF>
	struct Block {
		typedef Eigen::Matrix<double, DOF, DOF> Matrix;
		Eigen::LLT<Matrix> LLT;
		int start;
		
		
		template<bool transpose>
		void apply(Eigen::VectorXd& vect) const {
			// very future TODO invert matrix once and multiply in-place instead of solve
			if(transpose) // TODO check when to transpose
				LLT.matrixL().solveInPlace(vect.segment<DOF>(start));
			else
				LLT.matrixU().solveInPlace(vect.segment<DOF>(start));
		}
		
	};
	
	
	
	
	template<bool transpose>
	void apply(Eigen::VectorXd& vect) const {
		
	}
	
	
public:
	void compute();
	
	void apply(Eigen::VectorXd& vect, bool transp) const {
		if(transp){
			apply<true> (vect);
		} else {
			apply<false>(vect);
		}
	}
	
	
};

} /* namespace SLOM */
#endif /* BLOCKJACOBIPRECONDITIONER_HPP_ */
