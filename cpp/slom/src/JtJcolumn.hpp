/**
 * @file JtJcolumn.hpp
 * @brief Block storage of JtJ matrix
 * 
 */

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include "Sparse.hpp"

#include <vector>
#include <algorithm>


#ifndef JTJCOLUMN_HPP_
#define JTJCOLUMN_HPP_

namespace SLOM {

namespace internal {



/**
 * Storage container for one block-column of @f$ J^\top J@f$ matrix.
 * 
 * There is one block for the diagonal entries, as well as a block entry for each 
 * correspondence above the diagonal, i.e. with variables having a smaller index.
 * 
 * @todo maybe storing lower half instead of upper half?
 */
template<int _DOF, int _UpLo = Eigen::Upper>
struct JtJColumn {
	enum {DOF = _DOF, UpLo = _UpLo};
	typedef Eigen::Matrix<double, DOF, 1> VectorType;
	typedef Eigen::Matrix<double, DOF, DOF> MatrixDiag;
	typedef Eigen::Matrix<double, Eigen::Dynamic, DOF> MatrixBlock;
	typedef Eigen::LLT<MatrixDiag, UpLo> CholDiag; 
	struct JtJBlock {
		int start_idx, count;
		MatrixBlock mat;
		
		bool operator<(const JtJBlock& oth) const { return start_idx < oth.start_idx; }
		bool operator<(const int& idx) const { return start_idx < idx; }
		friend bool operator<(const int idx, const JtJBlock& oth) { return idx < oth.start_idx; }
		
		JtJBlock(int start_idx, int rows) : start_idx(start_idx), count(1), mat(rows, int(DOF)) { };
		MatrixBlock& operator*() { return mat; }
		const MatrixBlock& operator*() const { return mat; }
		int rows() const { return mat.rows(); }
		int start() const { return start_idx; }
		int end() const { return start_idx + rows(); } 
	};
	typedef std::vector<JtJBlock> Container;
	
	MatrixDiag diag;
	CholDiag cholDiag;
	Container column;
	typedef typename Container::iterator iterator;
	typedef typename Container::const_iterator const_iterator;
	
	
	//! Iterator pointing after last write position 
	iterator current_write;
	//! Iterator pointing after last read position
	const_iterator current_read;
	
	//! number of actually filled rows:
	int rows;
	bool cholCalculated;
	
	JtJColumn() : rows(0), cholCalculated(false) {}
	
	MatrixBlock& insert(int startIdx, int block_rows){
		iterator pos;
		if(column.empty() || column.back().end() <= startIdx) {
			// usually insertions will take place at the end:
			pos = column.end();
		} else {
			// otherwise do binary search:
			pos = std::lower_bound(column.begin(), column.end(), startIdx);
		}
		if(pos == column.end() || pos->start() != startIdx){
			pos = column.insert(pos, JtJBlock(startIdx, block_rows));
			assert((pos == column.begin() || pos[-1].end()<=pos->start()) 
					&& (pos == column.end()-1 || pos->end()<=pos[1].start())
					&& "Tried to insert overlapping blocks!");
			rows+=block_rows;
		} else {
			++pos->count;
		}
		assert(pos->rows() == block_rows);
		return pos->mat;
	}
	
	void setZero(){
		diag.setZero();
		for(iterator it = column.begin(); it != column.end(); ++it){
			it->mat.setZero();
		}
	}
	
	bool remove(int startIdx){
		iterator pos = std::lower_bound(column.begin(), column.end(), startIdx);
		assert(pos != column.end() && pos->start_idx == startIdx);
		if((--pos->count) == 0){ // removed last count
			rows -= pos->rows();
			column.erase(pos);
			return true;
		}
		return false;
	}
	
	MatrixBlock& find(int startIdx){
		iterator pos = std::lower_bound(column.begin(), column.end(), startIdx);
		assert(pos != column.end() && pos->start_idx == startIdx);
		return pos->mat;
	}
	MatrixBlock& operator()(int startIdx){
		return find(startIdx);
	}
	
	MatrixDiag& operator()(){
		// assume calling this method implies modification of diag.
		// FIXME not the best interface imaginable
		cholCalculated = false;
		return diag;
	}
	
	void updateChol() {
		if(cholCalculated) return;
		
		cholDiag.compute(diag);
		if(cholDiag.info() != Eigen::Success) {
			cholDiag.compute((diag.diagonal().array() + 1e-6).matrix().asDiagonal());
			assert(cholDiag.info() == Eigen::Success && "Failed to Cholesky decompose diagonal block!");
		}
		
	}
	
	int nnz() const {
		return rows * DOF + (DOF*(DOF+1))/2;
	}
	
	
	/**
	 * Assuming v is already updated, calculate 
	 * @f$ u = L^{-\top}(u - L^{-1} A^\top v) @f$
	 * note that u usually is part of v, but A just has entries above u.
	 * @param u location of sub-vector u
	 * @param v entire vector
	 */
	void backSubstitution(double * u, const double *v) {
		updateChol();
		VectorType Atv = multTranspose(v);
		cholDiag.matrixL().solveInPlace(Atv);
		Eigen::Map<VectorType> u_(u);
		u_ -= Atv;
		cholDiag.matrixU().solveInPlace(u_);
	}
	
	/**
	 * Calculates new @f$u = L^{-1} u @f$, then updates v as 
	 * @f$ v = (v - A L^{-\top} u) @f$
	 * note that u usually is part of v, but A just has entries above u.
	 * @param u location of sub-vector u
	 * @param v entire vector
	 */
	void forwardSubstitution(double * u, double *v) {
		updateChol();
		Eigen::Map<VectorType> u_(u);
		cholDiag.matrixL().solveInPlace(u_);
		VectorType minus_u = - u_;
		cholDiag.matrixU().solveInPlace(minus_u);
		addMult(v, minus_u);
	}

	
	
	/**
	 * Calculates A^T*v where A is the upper-triangular part of J^T*J.
	 * @warning no sanity-checks w.r.t. size of v
	 * @param   v  pointer to vector elements
	 * @return  a vector containing the result. Size is usually small, so passing by value should not be a problem
	 */
	VectorType multTranspose(const double *v) const {
		EIGEN_ASM_COMMENT("start multTranspose");
		VectorType result = VectorType::Zero(); // usually better vectorization with native vector?
		for(const_iterator it = column.begin(); it != column.end(); ++it) {
			result.noalias() += (**it).transpose() * Eigen::VectorXd::Map(v + it->start(), it->rows());
		}
		EIGEN_ASM_COMMENT("end multTranspose");
		return result;
	}
	
	/**
	 * calculates dest += A * vect, where A is the upper-triangular part of J^T*J.
	 * @warning no sanity-checks w.r.t. size of dest
	 * @param dest  pointer to elements of destination vector
	 * @param vect  vector with which to multiply
	 */
	void addMult(double * dest, const VectorType& vect) const {
		EIGEN_ASM_COMMENT("start addMult");
		for(const_iterator it = column.begin(); it != column.end(); ++it) {
			Eigen::VectorXd::Map(dest + it->start(), it->rows()).noalias() += (**it) * vect;
		}
		EIGEN_ASM_COMMENT("end addMult");
	}
	
	/**
	 * Export indices to SparseMatrix
	 * @param A   destination matrix. It is assumed, that enough space is available
	 * @param col start column number
	 */
	void storeIndexes(Eigen::SparseMatrix<double> & A, int col) const{
		int * Ap = A.outerIndexPtr() + col; // p[col] points is first index of col in i
		int * Ai = A.innerIndexPtr(); // first entry to be written to
		int* row = Ai + *Ap;
		
		// calculate indexes for first column:
		for(const_iterator it = column.begin(); it != column.end(); ++it){
			for(int r = it->start(); r<it->end(); ++r){
				*row++ = r;
			}
			assert(row[-1]<col && "There are entries below the diagonal in JtJ");
		}
		// entry for diagonal block
		*row++ = col;
		for(int j=1; j<DOF; ++j){
			Ap[j] = Ap[j-1] + rows + j;
			int* next = std::copy(Ai + Ap[j-1], Ai + Ap[j], Ai + Ap[j]);
			*next = col + j;
		}
		Ap[DOF] = Ap[DOF-1] + rows + DOF;
		
	}
	
	/**
	 * Export values to a Matrix previously set up using storeIndexes
	 * @param x  First address to store a value
	 * @return   pointer directly behind last write location
	 */
	double* storeValues(double *x) const {
		for(int j=0; j<DOF; ++j){
			for(const_iterator it = column.begin(); it!=column.end(); ++it){
				for(int i=0; i< it->rows(); ++i){
					*x++ = (**it)(i, j);
				}
			}
			for(int i=0; i<=j; ++i){
				*x++ = diag(i, j);
			}
		}
		return x;
	}
};


}  // namespace internal



}  // namespace SLOM

#endif /* JTJCOLUMN_HPP_ */
