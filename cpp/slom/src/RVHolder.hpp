/**
 * @file RVHolder.hpp
 * @brief Brief description
 * 
 */


#ifndef RVHOLDER_HPP_
#define RVHOLDER_HPP_

#include "internal.hpp"
#include "JtJcolumn.hpp"

namespace SLOM {

namespace internal {




template<class M>
class RVHolder : public IRVHolder
{
	template<class, int>
	friend class IMeasurement_Holder::VarRef;
	friend class VarID<M>;
	M var;
	M backup;
	
public:
	enum {DOF = M::DOF, DIM = DOF}; // DIM is needed for indexed_list
	typedef Eigen::Matrix<double, DOF, DOF> MatrixType;
	typedef Eigen::Matrix<double, DOF, 1> VectorType;
	
	RVHolder(const M& v=M(), bool optimize=true) : IRVHolder(optimize), var(v), backup(v) {}
	int getDOF() const {return DOF;}
	int getDim() const {return DOF;}
	const double* boxplus(const double* vec, double scale=1) {
		var = backup;
		var.boxplus(MTK::vectview<const typename M::scalar, DOF>(vec), scale);
		return vec + DOF;
	}
	void store() {backup = var;}
	void restore() {var = backup;}
	
	void reset(const M& m){
		var = backup = m;
	}
	
	// maybe make alignment conditional?
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
#ifdef SLOM_JACOBI_BLOCKS
	void calcBlockJacobi() {
		MatrixType temp = MatrixType::Zero();
		
		for(measurement_container::const_iterator it = measurements.begin(); it!=measurements.end(); ++it) {
			temp.template selfadjointView<Eigen::Lower>().rankUpdate(it->JacBlock<DOF>().transpose());
		}
		
		
		jacobiBlock.compute(temp);
		if(jacobiBlock.info() != Eigen::Success) {
			std::cerr << "Warning: (Near) singular Jacobi block\n";
//			std::cerr << map << std::endl;
//			assert(false && "Computation of jacobiBlock failed!");
			temp = temp.diagonal().asDiagonal(); 
			jacobiBlock.compute(temp + Eigen::Matrix<double, DOF, DOF>::Identity()*1e-9);
			assert(jacobiBlock.info() == Eigen::Success && "Computation of jacobiBlock failed!");
		}
//		return jBlock + dim*DOF;
	}
	
	double * applyBlockJacobi(double * vec, bool transpose) const {
		typename VectorType::MapType res(vec);
		
		if(transpose)
			jacobiBlock.matrixL().solveInPlace(res);
		else
			jacobiBlock.matrixU().solveInPlace(res);
		
		return vec + DOF;
	}
#endif

#ifdef SLOM_EXPERIMENTAL
	
	void applyCholInc(double* vec, bool transpose) {
		double * u = vec + idx;
		if(transpose) {
			columnJtJ.forwardSubstitution(u, vec);
		} else {
			columnJtJ.backSubstitution(u, vec);
		}
	}
	
	void resetJtJ() {
		columnJtJ.setZero();
	}
	
	void getJtJ(SparseType &JtJ, bool _JtJ_structureChanged) const {
		if(_JtJ_structureChanged) {
			columnJtJ.storeIndexes(JtJ, idx);
		}
		columnJtJ.storeValues(JtJ.valuePtr() + JtJ.outerIndexPtr()[idx]);
	}

	
protected:
	template<class Measurement>
	friend class Measurement_Holder<Measurement>::JtJAdder;
	
	template<class Measurement>
	template<bool insert>
	friend class Measurement_Holder<Measurement>::JtJInserter;
// TODO WORK IN PROGRESS
	// FIXME change to Eigen::Lower as soon as implemented
	typedef internal::JtJColumn<DOF, Eigen::Upper> JtJColumnType;
	JtJColumnType columnJtJ;
	
#endif

#ifdef SLOM_JACOBI_BLOCKS
	Eigen::LLT<MatrixType, Eigen::Lower> jacobiBlock;	// FIXME JtJColumn essentially saves/computes the same
#endif
};

}  // namespace internal

}  // namespace SLOM


#endif /* RVHOLDER_HPP_ */
