/**
 * @file /mtk/slom/preconditioner/SchurPrecond.hpp
 * @brief Brief description
 * 
 */

#ifndef SCHURPRECOND_HPP_
#define SCHURPRECOND_HPP_

#include <Eigen/Core>

#include <map>

namespace SLOM {
// FIXME rename type1, type2 to typeA, typeB?
/*
 *
 */
template<class Scalar, int dim1, int dim2>
class SchurPrecond {
	
	typedef Eigen::Matrix<Scalar, dim1, dim1> MBlock11;
	typedef Eigen::Matrix<Scalar, dim1, dim2> MBlock12;
	typedef Eigen::Matrix<Scalar, dim2, dim2> MBlock22;
	
	
	/*
	 * Matrix is stored as
	 * 
	 * [ Block11, Block12 ]
	 * [ (sym),   Block22 ]
	 * 
	 * Where Block11 is block-tridiagonal, Block12 is random-Block-Sparse and Block22 is Block-Diagonal
	 * Block21 is not stored because it is symmetric to Block12.
	 * 
	 */
	
	
	
	typedef std::map<int, MBlock22> Block22;
	
	struct Block2_obsolete_________ {
		MBlock22 block;
		int idx;
	};
	
	
	struct Block12 {
		MBlock12 block;
		typename Block22::const_iterator diag; // iterator pointing to corresponding diagonal in Block22
//fixme obsolete		const Block2* diag; // link to corresponding diagonal (required for normalization and to obtain the index)
		
	};
	
	/**
	 * Block for a type1 variable.
	 * Contains the actual Matrix, the Matrix describing the transition and a link to the parent block.
	 * Also index of the Block and a vector of all Blocks interacting with type 2 variables.
	 */
	struct Block1 {
		MBlock11 diag;
		MBlock11 trans;
		// Block1* parent; //FIXME future work for more complex structure
		int idx; // absolute index of corresponding variable
		
		std::vector<Block12> block12;
	};
	
	
	// fixme needs indexes, maybe separate deque per type1 variable?
	std::deque<MBlock12> block12;
	
	std::deque<Block1> block11;
	
	Block22 block22;
	
	
public:
	
	//! insert block into preconditioner
	template<int height>
	void insert(const Eigen::Matrix<Scalar, height, dim1>& mat1, int idx2, const Eigen::Matrix<Scalar, height, dim2>& mat2) {
		block22[idx2] = mat2.transpose() * mat2;
		
		;
	}
	
	//! go to next type1 variable
	void transition() {
		
	}
	/**
	 * Go to next type1 variable, including a transition function.
	 * @param trans1 Jacobian w.r.t. the first variable.
	 * @param trans2 Jacobian w.r.t. the second variable.
	 */
	void transition(const MBlock11& trans1, const MBlock11& trans2);
	
	//! calculate preconditioner
	void calculate() {
		// block-Jacobi + Cholesky of Block22 
		
		// apply on Block12
		
		// calculate Schur-Complement (only diagonal and first off-diagonal)
		
		// Cholesky decomposition of tridiagonal Matrix (could happen on the fly)
		// If underflow occurs, just restart decomposition from that point
		
	}
	
	//! resize preconditioner (could happen on-the-fly?)
	void resize(int num1, int num2);
	//! reset preconditioner w/o resizing
	void reset();
	
	//! apply preconditioner to \c vec (in-place)
	void apply(Eigen::VectorXd& vec, bool trans) const;
	
	
public:
	SchurPrecond();
};

} /* namespace SLOM */
#endif /* SCHURPRECOND_HPP_ */
