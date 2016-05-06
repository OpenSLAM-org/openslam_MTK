/**
 * @file MeasurementHolder.hpp
 * @brief Brief description
 * 
 */


#ifndef MEASUREMENTHOLDER_HPP_
#define MEASUREMENTHOLDER_HPP_


#include "internal.hpp"
#include "RVHolder.hpp"

#include "../Covariance.hpp"

namespace SLOM {

namespace internal {

/**
 * Actual implementation of an \c IMeasurement_Holder
 * @tparam Measurement determines the type of the measurement
 */
template<class Measurement>
struct Measurement_Holder : IMeasurement_Holder {
	enum {DIM = Measurement::DIM, DEPEND = Measurement::DEPEND};
	
	typedef Eigen::Matrix<double, DIM, 1> ResultType;
	typedef Eigen::Matrix<double, DEPEND, 1> VariableVecType;
	typedef Eigen::Matrix<double, DEPEND, DEPEND> JtJBlockType;
	
	
	Measurement m_;
	Measurement_Holder(const Measurement & m) : m_(m) {
#ifdef SLOM_JACOBI_BLOCKS
		 jac.resize(m_.dim, DEPEND);
#endif
	}
	
	int getDim() const { return m_.dim; }
	double* eval(double* ret, bool numerical_jacobian = false) const {
		MTK::vectview<double, DIM> retvec(ret, m_.dim);
		m_.eval(retvec, numerical_jacobian);
		if(!std::isfinite(retvec.squaredNorm())) { // FIXME mark as SLOM_DEBUG
			std::cerr << __PRETTY_FUNCTION__ << ": " << retvec << std::endl;
			assert(false);
		}
		return ret + m_.dim;
	}
	
	template<bool register_> 
	struct Registrator {
		Measurement_Holder* ptr;
		int &count;
		Registrator(Measurement_Holder* ptr, int &count) : ptr(ptr), count(count) {}
		
		template<class Type, int IDX>
		MTK_ALWAYS_INLINE void operator()(VarRef<Type, IDX> var) {
//			if(!var.getsOptimized()) return; // FIXME fixed this was fatal if variables later change from fixed to non-fixed
			IRVHolder* v = var;
			int dim = ptr->getDim();
			if(register_)
				count += v->registerMeasurement(ptr
#ifdef SLOM_JACOBI_BLOCKS
						, ptr->jac.data()+IDX * dim, dim
#else
						, 0, dim
#endif
						);
			else
				count += v->unregisterMeasurement(ptr
						, dim
						);
		}
	};

	
	int registerVariables() {
#ifdef SLOM_EXPERIMENTAL
		m_.traverse_variables(JtJInserter<true>(m_));
#endif
		int count = 0;
		m_.traverse_variables(Registrator<true>(this, count));
		return count * m_.dim;
	}
	
	int unregisterVariables() {
#ifdef SLOM_EXPERIMENTAL
		m_.traverse_variables(JtJInserter<false>(m_));
#endif
		int count = 0;
		m_.traverse_variables(Registrator<false>(this, count));
		return count * m_.dim;
	}
	
#ifdef SLOM_JACOBI_BLOCKS
	// ColMajor is strictly required here: FIXME use default and check that EIGEN_DEFAULT_ROWMAJOR is not defined
	typedef Eigen::Matrix<double, DIM, DEPEND, (DIM != 1 ? Eigen::ColMajor : Eigen::RowMajor) || Eigen::AutoAlign> JacobiBlock;
	JacobiBlock jac;
	
	template<class Type, int IDX>
	typename JacobiBlock::template NColsBlockXpr<Type::DOF>::Type jacSubblock(VarRef<Type, IDX> Measurement::*) {
		return jac.template middleCols<Type::DOF>(IDX);
	}
	template<class Type, int IDX>
	const typename JacobiBlock::template ConstNColsBlockXpr<Type::DOF>::Type jacSubblock(VarRef<Type, IDX> Measurement::*) const {
		if(0 <= IDX && Type::DOF + IDX <= DEPEND){} else {
			std::cerr << "Out of bound error " << jac.rows() << "x" << jac.cols() << ", " << Type::DOF << " + " << IDX << std::endl;
		}
			
		return jac.template middleCols<Type::DOF>(IDX);
	}
	
	struct JacobiUpdater{
		Measurement_Holder& m;
		
		JacobiUpdater(Measurement_Holder &m) : m(m) {}
		
		// TODO This method can be specialized to implement symbolic derivations
		template<class Type, int IDX>
//		MTK_ALWAYS_INLINE // inlining seems to hurt performance here
		void operator()(const VarRef<Type, IDX> var_){
			if(!var_.getsOptimized()) return;
			
			Eigen::Matrix<double, DIM, 1> temp(m.getDim());
			// FIXME replace delta by 1e-6 * Eigen::Matrix<...>::Unit(i);
			Eigen::Matrix<double, Type::DOF,1> delta;
			Type& var = var_.var();
			delta.setZero();
			Type backup = var;
			EIGEN_ASM_COMMENT("JacobiUpdate_start");
			for(int i=0; i<Type::DOF; ++i){
				delta[i] = 1e-6;
				var.boxplus(delta, 1);
				m.m_.eval(temp, true);
				var = backup;
				var.boxplus(delta, -1);
				m.m_.eval(m.jac.col(IDX + i), true);
				m.jac.col(IDX + i) =  0.5e6 *( temp - m.jac.col(IDX + i));
//				assert((m.jac.col(IDX+i).norm())<1e12);
				var = backup;
				delta[i] = 0;
			}
			EIGEN_ASM_COMMENT("JacobiUpdate_end");
		}
		
	};
	
	void updateJacobian() {
		// TODO pass numerical delta ?
		EIGEN_ASM_COMMENT("updateJacobian");
		JacobiUpdater ju(*this);
		m_.traverse_variables(ju);
		EIGEN_ASM_COMMENT("updateJacobian_end");
//		std::cout << __PRETTY_FUNCTION__ << idx << "\n" << jac << std::endl;
	};
	
	
	struct FetchData {
		VariableVecType & data;
		const double * source;
		
		FetchData(VariableVecType & data, const double * source) : data(data), source(source) {}
		
		template<class Type, int IDX>
		MTK_ALWAYS_INLINE void operator()(VarRef<Type, IDX> var) {
			if(var.getsOptimized())
				data.template segment<Type::DOF>(IDX) =
						MTK::vectview<const double, Type::DOF>(source + var.ptr->idx);
			else
				data.template segment<Type::DOF>(IDX).setZero();
		}
	};
	template<bool Add_instead_of_Store>
	struct SpreadData {
		const VariableVecType & data;
		double * dest;
		
		SpreadData(const VariableVecType & data, double * dest) : data(data), dest(dest) {}
		
		template<class Type, int IDX>
		MTK_ALWAYS_INLINE void operator()(VarRef<Type, IDX> var) {
			if(!var.getsOptimized()) return;
			MTK::vectview<double, Type::DOF> res(dest + var.ptr->idx);
			if(Add_instead_of_Store)
				res+= data.template segment<Type::DOF>(IDX);
			else
				res = data.template segment<Type::DOF>(IDX);
		}
	};
#endif
	
#ifdef SLOM_EXPERIMENTAL
	double* addJv(double *res, const double *v) const {
		VariableVecType data;
		m_.traverse_variables(FetchData(data, v));
		
		ResultType::Map(res) += jac * data;
		
		return res + DIM;
	};
	
	const double* addJtu(double *res, const double *u) const {
		VariableVecType result = jac.transpose() * ResultType::Map(u);
		m_.traverse_variables(SpreadData<true>(result, res));
		
		return u + DIM;
	};
	
	template<bool insert>
	struct JtJInserter {
		template<class Type1, int IDX1>
		struct JtJColumnInserter {
			typedef typename RVHolder<Type1>::JtJColumnType JtJColumnType;
			JtJColumnType &jtjCol;
			int idx;
			JtJColumnInserter(JtJColumnType& jtjCol, int idx) : jtjCol(jtjCol), idx(idx) {};
			template<class Type2, int IDX2>
			MTK_ALWAYS_INLINE void operator()(VarRef<Type2, IDX2> var) {
				if(!var.getsOptimized()) return;
				int idx2 = var.ptr->idx;
				if(idx2 < idx) {
					if(insert) {
						jtjCol.insert(idx2, Type2::DOF);
					} else {
						jtjCol.remove(idx2);
					}
				}
			}
		};
		const Measurement &m;
		JtJInserter(const Measurement &m) : m(m) {}
		template<class Type, int IDX>
		MTK_ALWAYS_INLINE void operator()(VarRef<Type, IDX> var) {
			if(!var.getsOptimized()) return;
			m.traverse_variables(JtJColumnInserter<Type, IDX>(var.ptr->columnJtJ, var.ptr->idx));
		}
	};
	
	struct JtJAdder {
		template<class Type1, int IDX1>
		struct JtJColumnAdder {
			typedef typename RVHolder<Type1>::JtJColumnType JtJColumnType;
			const JtJBlockType &block;
			JtJColumnType &jtjCol;
			int idx;
			JtJColumnAdder(const JtJBlockType & block, JtJColumnType& jtjCol, int idx) : block(block), jtjCol(jtjCol), idx(idx) {};
			template<class Type2, int IDX2>
			MTK_ALWAYS_INLINE void operator()(VarRef<Type2, IDX2> var) {
				if(!var.getsOptimized()) return;
				if(var.ptr->idx < idx) {
					jtjCol.find(var.ptr->idx) += block.template block<Type2::DOF, Type1::DOF>(IDX2, IDX1);
				}
			}
		};
		const Measurement &m;
		const JtJBlockType &block;
		JtJAdder(const Measurement &m, const JtJBlockType& block) : m(m), block(block) {}
		template<class Type, int IDX>
		MTK_ALWAYS_INLINE void operator()(VarRef<Type, IDX> var) {
			if(!var.getsOptimized()) return;
			var.ptr->columnJtJ() += block.template block<Type::DOF, Type::DOF>(IDX, IDX);
			m.traverse_variables(JtJColumnAdder<Type, IDX>(block, var.ptr->columnJtJ, var.ptr->idx));
		}
	};
	
	void addJtJ() const {
		JtJBlockType block;
		block.template triangularView<Eigen::Lower>()= jac.transpose() * jac;
		block.template triangularView<Eigen::StrictlyUpper>() = block.transpose();
		m_.traverse_variables(JtJAdder(m_, block));
	};
	

	
#endif /* SLOM_EXPERIMENTAL */
	
	
	
	// this is necessary, if Measurement has variables, which need to be aligned:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};


/**
 * Measurement_Holder which automatically normalizes the covariance by applying the inverse square root covariance.
 * This requires slom/Covariance.hpp to be included.
 */
template<class Measurement, class invDevType>
struct Measurement_Holder_with_Dev : Measurement_Holder<Measurement> {
	InvDeviation<invDevType> dev;
	typedef Measurement_Holder<Measurement> base;
	
	Measurement_Holder_with_Dev(const Measurement & m, InvDeviation<invDevType> dev) : base(m), dev(dev) { }
	
	
	
	double* eval(double* ret, bool numerical_jacobian = false) const {
		int dim = base::getDim();
		MTK::vectview<double, base::DIM> ret_(ret, dim);
		base::m_.eval(ret_, numerical_jacobian);
		dev.apply(ret_);
		if(!std::isfinite(ret_.squaredNorm())) { // TODO mark as SLOM_DEBUG
			std::cerr << __PRETTY_FUNCTION__ << ": " << ret_ << std::endl;
			assert(false);
		}
		return ret + dim;
	}
	
#ifdef SLOM_JACOBI_BLOCKS
	void updateJacobian() {
		// FIXME FATAL: This needs to be implemented. 
		base::updateJacobian();
		// FIXME apply InvDeviation
	}
#endif
};



}  // namespace internal

}  // namespace SLOM


#endif /* MEASUREMENTHOLDER_HPP_ */
