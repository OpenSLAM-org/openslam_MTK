/**
 * @file Sparse.hpp
 * @brief Eigen Sparse Matrix compatibility file
 * 
 */


#ifndef SPARSE_HPP_
#define SPARSE_HPP_



#include <Eigen/Sparse>

#include <set>
#include <iostream>



namespace SLOM {

namespace internal {


// TODO implement everything generic:
// TODO put in file that is not included everywhere
typedef Eigen::SparseMatrix<double> SparseType;


typedef std::set<int> IndexSet;

/* x(0:mark-1) = x(0:mark-1) + beta * A(0:mark-1,j), where x is a dense vector and A(:,j) is sparse */
inline
void cs_scatter_upper (const SparseType &A, int j, double beta, int *w, double *x, int mark,
		SparseType *C)
{
	
	for(SparseType::InnerIterator it(A, j); it; ++it)
	{
		int i = it.index() ;                            /* A(i,j) is nonzero */
		if(i >= mark) break; // FIXME: this assumes columns of A are sorted, so we can break here
		if (w [i] < mark)
		{
			w [i] = mark ;                      /* i is new entry in column j */
			if(C){
				C->insertBackByOuterInnerUnordered( mark-1, i);
			}
			x [i] = beta * it.value() ;             /* x(i) = beta*A(i,j) */
		}
		else
		{
			x [i] += beta * it.value() ;    /* i exists in C(:,j) already */
		}
	}
}


/**
 * calculate upper half of A * B, (assuming columns of A are sorted)
 * first entry in result is the diagonal entry
 * For recalculate == true it is assumed that C already contains the correct structure.
 * 
 * @param C result
 * @param A left matrix
 * @param B right matrix
 * @param nnz estimate for nnz in result
 * @deprecated shall be replaced by corresponding Eigen function, or by block-matrix multiplication
 */
template<bool structure_changed, int UpLo>
//__attribute__((deprecated))
void cs_JtJ (SparseType & C, const SparseType &A, const SparseType &B, int nnz = -1)
{
	assert(UpLo == Eigen::Upper && "Currently, only upper half of JtJ can be calculated");
	eigen_assert(A.outerSize() == B.innerSize());
	int m = A.innerSize(), anz = A.nonZeros();
	int n = B.outerSize(), bnz = B.nonZeros() ;
	
	// get workspace:
	Eigen::VectorXi w(m);
	Eigen::VectorXd x(m);
	if(structure_changed){
		C.resize(m, n);
		C.reserve(nnz > 0 ? nnz : anz + bnz);
	} else {
		eigen_assert(C.rows() == m && C.cols() == n && "Matrix size must agree for recalculation");
	}
	for (int j = 0 ; j < n ; j++)
	{
		w[j] = j+1;    // mark entry C(j,j)
		if(structure_changed){
			C.startVec(j);
			C.insertBackByOuterInnerUnordered(j,j);
		}
		x[j] = 0;      // reset entry
		for (SparseType::InnerIterator it(B, j); it; ++it)
		{
			cs_scatter_upper(A, it.index(), it.value(), w.data(), x.data(), j+1, structure_changed ? &C : 0) ;
		}
		for (SparseType::InnerIterator it(C, j); it; ++it) {
			it.valueRef() = x[it.index()];
		}
	}
	if(structure_changed)
		C.finalize();
}

}  // namespace internal

}  // namespace SLOM




#endif /* SPARSE_HPP_ */
