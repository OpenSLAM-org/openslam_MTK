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
	typedef std::deque<const IMeasurement_Holder*> measurement_container;
	measurement_container measurements; // maybe a single linked pointer list?
	bool optimize;
	bool isRegistered;
	
	template<class M>
	struct holder;
	
	IRVHolder(bool optimize) : optimize(optimize), isRegistered(false) {}
	
	bool has_measurements() const {
		return !measurements.empty();
	}
	
	/**
	 * Registers a Measurement and returns the DOFs if this variable is to be optimized.
	 */
	int registerMeasurement(const IMeasurement_Holder* m){
		assert (isRegistered);
		measurements.push_back(m);
		return optimize ? getDOF() : 0;
	}
	
	/**
	 * Unregisters a Measurement and returns the DOFs if this variable is to be optimized.
	 */
	int unregisterMeasurement(const IMeasurement_Holder* m){
		assert (isRegistered);
		measurement_container::iterator it = std::find(measurements.begin(), measurements.end(), m);
		assert(it != measurements.end() && "Tried to remove an unregistered measurement!"); 
		measurements.erase(it);
		return optimize ? getDOF() : 0;
	}
	
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
};
typedef indexed_list<IRVHolder> RVList;

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
	
	template<class M>
	struct holder;
	template<class Type, int idx>
	struct VarRef;
protected:
	friend class SLOM::SparseFunction;
	
	
	void register_at_variables(int& count, IRVHolder* rv) const {
		count += rv->registerMeasurement(this);
	}
	void unregister_at_variables(int& count, IRVHolder* rv) const {
		count += rv->unregisterMeasurement(this);
	}
};

template<class M>
class IRVHolder::holder : public IRVHolder
{
	template<class, int>
	friend class IMeasurement_Holder::VarRef;
	friend class VarID<M>;
	M var;
	M backup;
public:
	enum {DOF = M::DOF, DIM = DOF}; // DIM is needed for indexed_list
	holder(const M& v=M(), bool optimize=true) : IRVHolder(optimize), var(v), backup(v) {}
	int getDOF() const {return DOF;}
	const double* boxplus(const double* vec, double scale=1) {
		var = backup;
		var.boxplus(vec, scale);
		return vec + DOF;
	}
	void store() {backup = var;}
	void restore() {var = backup;}
	
	void reset(const M& m){
		var = backup = m;
	}
	
	// maybe make alignment conditional?
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};


template<class M, int idx_in_Measurement>
struct IMeasurement_Holder::VarRef : public SLOM::VarID<M>{
	enum { IDX = idx_in_Measurement };
	typedef SLOM::VarID<M> base;
	using base::ptr;
	VarRef(const base &m) : base(m) {assert(m.ptr->isRegistered);}
	
//	VarRef() {} // no standard constructor necessary!
	
	const M& operator*()  const { return  ptr->var;}
	const M* operator->() const { return &ptr->var;}
	
	operator IRVHolder*() const { return ptr; }
	
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

/**
 * Actual implementation of an \c IMeasurement_Holder
 * @tparam Measurement determines the type of the measurement
 */
template<class Measurement>
struct IMeasurement_Holder::holder : IMeasurement_Holder {
	Measurement m_;
	holder(const Measurement & m) : m_(m) { }
	
	enum {DIM = Measurement::DIM};
	
	int getDim() const { return DIM; }
	double* eval(double* ret, bool numerical_jacobian = false) const {
		m_.eval(ret, numerical_jacobian);
		return ret + DIM;
	}
	
	int registerVariables() { // FIXME this method could actually be const, but boost::bind doesn't support that
		int count = 0;
		m_.traverse_variables(boost::bind(&IMeasurement_Holder::register_at_variables, this, boost::ref(count), _1));
		return count * DIM;
	}
	
	int unregisterVariables() {
		int count = 0;
		m_.traverse_variables(boost::bind(&IMeasurement_Holder::unregister_at_variables, this, boost::ref(count), _1));
		return count * DIM;
	}
	
	// this is necessary, if Measurement has variables, which need to be aligned:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

typedef indexed_list<IMeasurement_Holder> MeasurementList;


}  // namespace internal
}  // namespace SLOM


#endif /* SLOM_INTERNAL_HH_ */
