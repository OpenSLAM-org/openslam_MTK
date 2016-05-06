/*
 *  Copyright (c) 2008--2011, Universitaet Bremen
 *  All rights reserved.
 *
 *  Author: Christoph Hertzberg <chtz@informatik.uni-bremen.de>
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the Universitaet Bremen nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */
/**
 * @file  slom/src/internal.hpp
 * @brief Internal classes for data management.
 */

#ifndef SLOM_INTERNAL_HH_
#define SLOM_INTERNAL_HH_

#include <Eigen/Core>
#include "indexed_list.hpp"
#include <deque>
#include <boost/bind.hpp>
#include "../../mtk/src/vectview.hpp"

namespace SLOM {

template<class Manifold>
class VarID;
class SparseFunction;

namespace internal {



struct IMeasurement_Holder;

struct IRVHolder : public indexed_list_hook
{
	struct measurement_ptr {
		const IMeasurement_Holder* ptr;
		
#ifdef SLOM_JACOBI_BLOCKS
		// FIXME back-propagate constness
		const double* jac_block; // first entry of the variable's subblock in the jacobian of the measurement
#endif
		int measurement_dim;     //dimension of the measurement
		
		// FIXME for extreme optimization: function pointers to optimized matrix multiplication routines
		
		// FIXME remove default values
		measurement_ptr(const IMeasurement_Holder* ptr, double* jac_block = 0, int measurement_dim = 0)
			: ptr(ptr), 
#ifdef SLOM_JACOBI_BLOCKS
			  jac_block(jac_block), 
#endif
			  measurement_dim(measurement_dim) { (void) jac_block;}
		
		inline double * eval(double* res, bool numerical_jacobian = false) const;
		
		operator const IMeasurement_Holder*() const {return ptr;}
		
		const IMeasurement_Holder * operator->() const {
			return ptr;
		}
		const IMeasurement_Holder & operator*() const {
			return *ptr;
		}
		
		bool operator==(const IMeasurement_Holder *other) const {
			return (ptr == other);
		}
		
#ifdef SLOM_JACOBI_BLOCKS
		template<int DimVar>
		Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, DimVar> > JacBlock() const {
			return Eigen::Matrix<double, Eigen::Dynamic, DimVar>::Map(jac_block, measurement_dim, DimVar);
		}
#endif
	};
	typedef std::deque<measurement_ptr> measurement_container;
//	typedef std::deque<const IMeasurement_Holder*> measurement_container;
	measurement_container measurements; // maybe a single linked pointer list?
	int measurementDimension;
	bool optimize;
	bool isRegistered;
	
	IRVHolder(bool optimize) : optimize(optimize), isRegistered(false) {}
	
	bool has_measurements() const {
		return !measurements.empty();
	}
	
	
	/**
	 * Registers a Measurement and returns the DOFs if this variable is to be optimized.
	 */
	inline
	int registerMeasurement(const IMeasurement_Holder* m, double* jac, int dim);
	
	/**
	 * Unregisters a Measurement and returns the DOFs if this variable is to be optimized.
	 */
	inline
	int unregisterMeasurement(const IMeasurement_Holder* m, int dim);

	
	/**
	 * Return the total dimensionality of measurements depending on this variable.
	 */
	inline int measurementDimensions() const;
	
	virtual ~IRVHolder() {}
	virtual int getDOF() const = 0;
	/**
	 * Set working copy to backup [+] (scale*vec), i.e. no previous call to restore() necessary.
	 */
	virtual const double* boxplus(const double* vec, double scale=1) = 0;
	virtual void store() = 0;
	virtual void restore() = 0;
	
	
	
#ifdef SLOM_JACOBI_BLOCKS
	/**
	 * Calculate the Jacobi-Block preconditioner for this variable.
	 */
	virtual void calcBlockJacobi() {};
	/**
	 * Apply the the Jacobi-Block preconditioner to the vector starting at vec.
	 * Returns vec+DOF, to simplify iteration.
	 */
	virtual double * applyBlockJacobi(double * vec, bool transpose) const = 0;
#endif

};

template<class>
class RVHolder;

typedef indexed_list<IRVHolder, RVHolder> RVList;

/**
 * IMeasurement_Holder is (the node of) an intrusive list of all measurements.
 * Actual measurements hold (i.e. aggregate) an arbitrary measurement object, which implements
 * void eval(double*) const or a compatible function and is able to traverse the variables it depends on.
 */
struct IMeasurement_Holder : public indexed_list_hook
{
	virtual ~IMeasurement_Holder() {}
	virtual int getDim() const = 0;
	virtual double* eval(double*, bool numerical_jacobian) const = 0;
	
	template<class Type, int idx>
	struct VarRef;
protected:
	friend class SLOM::SparseFunction;
	
#ifdef SLOM_EXPERIMENTAL
	// TODO addJ** methods should be const
	//! adds JtJ blocks to corresponding blocks in referenced Variables.
	virtual void addJtJ()= 0;
	virtual double* addJv(double *res, const double *v) = 0;
	virtual const double* addJtu(double *res, const double *u) = 0;
#endif /* SLOM_EXPERIMENTAL */

#ifdef SLOM_JACOBI_BLOCKS
	virtual void updateJacobian() = 0;
#endif
	
//	void register_at_variables(int& count, IRVHolder* rv) const {
//		count += rv->registerMeasurement(this);
//	}
//	void unregister_at_variables(int& count, IRVHolder* rv) const {
//		count += rv->unregisterMeasurement(this);
//	}
};







template<class Measurement>
struct Measurement_Holder;

template<class M, int idx_in_Measurement>
struct IMeasurement_Holder::VarRef : public SLOM::VarID<M>{
	enum { IDX = idx_in_Measurement };
	typedef SLOM::VarID<M> base;
	using base::ptr;
	VarRef(const base &m) : base(m) {
		assert(m.ptr && m.ptr->isRegistered && "You must register variables before passing them to measurements");
	}
	
//	VarRef() {} // no standard constructor necessary!
	
	const M& operator*()  const { return  ptr->var;}
	const M* operator->() const { return &ptr->var;}
	
	operator IRVHolder*() const { return ptr; }
	
private:
#ifdef SLOM_JACOBI_BLOCKS
	template<class Measurement>
	friend class Measurement_Holder<Measurement>::JacobiUpdater;
#endif
	
	// access to non-const var:
	M& var() const {return ptr->var;}

	/* TODO
	 * maybe also make VarRef a linked list (for the measurements container in IRVHolder)
	 * using also a pointer to the encapsulating IMeasurement.
	 * --> in this case a std::list<IMeasurement*>::iterator might be sufficient?
	 * 
	 * Also maybe put a matrix block for the Jacobian in it (future work, would 
	 * require another template parameter for dim of measurement)
	 */
};



int IRVHolder::measurementDimensions() const {
	// This could be more efficient if this is accumulated within registerMeasurement
	int dim = 0;
	for(measurement_container::const_iterator it = measurements.begin();
			it != measurements.end(); ++it)
		dim += (**it).getDim();
	
	return dim;
}

typedef indexed_list<IMeasurement_Holder, Measurement_Holder> MeasurementList;

int IRVHolder::registerMeasurement(const IMeasurement_Holder* m, double * jac, int dim){
	assert (isRegistered && "Variable is not registered");
	assert(dim == m->getDim() && "Dimensions do not match");
	measurements.push_back(measurement_ptr(m, jac, dim));
	measurementDimension += dim;
	// FIXME This is now called with the actual type of the Variable known, so no virtual call is required.
	return optimize ? getDOF() : 0;
}

/**
 * Unregisters a Measurement and returns the DOFs if this variable is to be optimized.
 */
int IRVHolder::unregisterMeasurement(const IMeasurement_Holder* m, int dim){
	assert (isRegistered);
	assert(dim == m->getDim() && "Dimensions do not match");
	measurement_container::iterator it = std::find(measurements.begin(), measurements.end(), m);
	assert(it != measurements.end() && "Tried to remove an unregistered measurement!"); 
	measurements.erase(it);
	measurementDimension -= dim;
	return optimize ? getDOF() : 0;
}

inline double * IRVHolder::measurement_ptr::eval(double* res, bool numerical_jacobian) const {
	return ptr->eval(res, numerical_jacobian);
}




}  // namespace internal
}  // namespace SLOM


#endif /* SLOM_INTERNAL_HH_ */
