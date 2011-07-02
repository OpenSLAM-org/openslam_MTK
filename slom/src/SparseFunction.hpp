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
 * @file slom/src/SparseFunction.hpp
 * 
 * @brief Declaration of SparseFunction, VarID and MeasID.
 */

#ifndef SPARSEFUNCTION_HH_
#define SPARSEFUNCTION_HH_

#include "internal.hpp"
#include <Eigen/Core>

#include <cs.h>

#include <deque>
#include <iostream>
#include <numeric>


namespace SLOM {


/**
 * VarID holds a pointer to a Random Variable on user level.
 * It can be used read only, to access the current value of a variable, 
 * or when creating measurements to determine data association.
 */
template<class M>
class VarID : private internal::RVList::id<M>{
	typedef internal::RVList::id<M> base;
	//	using base::ptr; // this line makes ptr public in g++ 4.3.2 (which seems to be a bug)
	friend class SparseFunction;
	//VarRef needs to access base (actually only VarRef<M, int>):
	template<class, int>
	friend class internal::IMeasurement_Holder::VarRef;
public:
	VarID(const base &m) : base(m) {}
	VarID() : base(0) {}
	
	//! checks, if this is a valid ID.
	operator bool() const
	{
		return base::ptr != 0;
	}
	
	bool getsOptimized() const {
		return base::ptr && base::ptr->optimize;
	}
	template<class M2>
	bool operator<(const VarID<M2>& v) const { return base::ptr < v.ptr; }
	
	const M& operator*()  const { return  base::ptr->backup;}
	const M* operator->() const { return &base::ptr->backup;}
};


template<class Measurement>
class MeasID : private internal::MeasurementList::id<Measurement> {
	typedef internal::MeasurementList::id<Measurement> base;
	MeasID(internal::IMeasurement_Holder::holder<Measurement> &m) : base(&m) {}
	MeasID(const internal::MeasurementList::id<Measurement> &id) : base(id) {}
	friend class SparseFunction;
public:
	MeasID() : base(0) {}
	const Measurement& operator*() const
	{
		return base::ptr->m_;
	}
	const Measurement* operator->() const
	{
		return &**this;
	}
};


class SparseFunction {
	
	friend class Estimator;
	
	typedef internal::RVList RVList;
	typedef internal::IRVHolder IRVHolder;
	
	typedef internal::MeasurementList MeasurementList;
	typedef internal::IMeasurement_Holder IMeasurement_Holder;
	
	
	//! \name Storage of variables and measurements
	//@{
	RVList variables;
	RVList fixed_variables;
	
	internal::MeasurementList measurements;
	//@}
	
	//! \name Internal helper variables
	//@{
	//! keep track whether the structure or values have changed
	bool _structureChanged, _addDiagonal, _valuesChanged;
	
	//! delta for numerical calculation of Jacobian
	double numerical_delta, half_delta_inv;
	
	//! Number of non-zeroes in the jacobian
	int nnz;
	//! sparse Jacobian of the function
	cs* jacobian;
	
	// diagonal of $J^\top J$
	Eigen::VectorXd cholCovariance;
	
	// helper variables
	cs* Jt; // transpose of jacobian
	cs* JtJ;  // Jt * J
	//@}
	
	//! \name Construction and destruction
	//@{
	/**
	 * Construct a SparseFunction object
	 * @param delta delta for numerical differentiation
	 */
	SparseFunction(double delta = 1e-6) 
	: _structureChanged(true),
	  numerical_delta(delta), half_delta_inv(0.5/delta),
	  nnz(0),
	  jacobian(0), Jt(0), JtJ(0) 
	{
		assert(delta > 0 && "delta must be positive");
	}
	
	void freeWorkspace() {
		jacobian = cs_di_spfree(jacobian);
		freeTemporaries();
	}
	
	void freeTemporaries() {
		Jt =  cs_di_spfree(Jt);
		JtJ =  cs_di_spfree(JtJ);
	}
	~SparseFunction() {
		freeWorkspace();
	}
	//@}
	
	/**
	 * \name Variable Managment
	 * Methods for inserting and removing of variables and measurements.
	 * Methods are passed encapsulated by Estimator.
	 * @{
	 */ 
	
	template<class Manifold>
	VarID<Manifold> insertRV(const Manifold& m, bool optimize = true){
		RVList &list = optimize ? variables : fixed_variables; 
		IRVHolder::holder<Manifold>* ptr = new IRVHolder::holder<Manifold>(m, optimize);
		VarID<Manifold> id = list.insert(ptr);
		ptr->isRegistered = true;
		_structureChanged |= optimize;
		return id;
	}
	
	template<class Manifold>
	VarID<Manifold> reinitRV(VarID<Manifold> id, const Manifold& m) {
		id.ptr->reset(m);
		_valuesChanged = true;
		return id;
	}
	
	template<class Manifold>
	VarID<Manifold> optimizeRV(VarID<Manifold> id, bool optimize=true) {
		if(optimize == id.ptr->optimize)
			return id;
		
		RVList &old_list = !optimize ? variables : fixed_variables;
		RVList &new_list =  optimize ? variables : fixed_variables;
		
		int dim = id.ptr->measurementDimensions();
		nnz += Manifold::DOF * (optimize ? dim : -dim);
		id.ptr->optimize = optimize;
		
		old_list.unhook(id);
		return new_list.insert(id.ptr);
		
	}
	
	template<class Manifold>
	bool removeRV(const VarID<Manifold> &id){
		VarID<Manifold> id_(id);
		if(id_.ptr->has_measurements()) {
			std::cerr << "Tried to remove a variable still in use!" << std::endl;
			return false;
		}
		bool optimize = id_.ptr->optimize;
		RVList &list = optimize ? variables : fixed_variables;
		//isInitialized &= !optimize; // if variable had to be optimized, re-initialization is necessary
		_structureChanged |= optimize;
		
		list.remove(id_);
		return true;
	}
	
	template<class Measurement>
	MeasID<Measurement> insertMeasurement(const Measurement &m) {
		MeasID<Measurement> id = measurements.insert(new IMeasurement_Holder::holder<Measurement>(m));
		nnz += id.ptr->registerVariables();
		//isInitialized = false;
		_structureChanged = true;
		
		return id;
	}
	
	template<class Measurement>
	bool removeMeasurement(const MeasID<Measurement> &id){
		nnz -= id.ptr->unregisterVariables();
		_structureChanged = true;
		
		return true;
	}
	//@}
	
	
	
	//! \name Methods needed for optimization algorithms
	//@{
	/**
	 * Creates "the big matrix", 
	 */
	void createSparse();
	void updateDiagonal(double lamda, double* cholCovariance);
	void initCovariance();
	
	/**
	 * Numerically calculates the Jacobian of the function with the 
	 * and updates cholCovariance. The result is stored in matrix.
	 */
	void calculateJacobian();
	
	double evaluate(double *result) const;
	
	
	/**
	 * get the current Jacobian. If neither values nor structure has changed since last call,
	 * no recalculation is done.
	 * If add_diagonal = true, each column contains an extra entry for diagonals
	 */ 
	const cs* getJ(){
		initialize(_addDiagonal);
		if(_valuesChanged){
			calculateJacobian();
			freeTemporaries(); // remove Jt and JtJ if values changed
			_valuesChanged = false;
		}
		return jacobian;
	}
	
	const cs* getJt(){
		if(!Jt){
			const cs* J = getJ();
			Jt = cs_transpose(J, true);
		}
		return Jt;
	}
	
	/**
	 * res += J * v
	 */
	void add_Jv(double *res, const double* v){
		cs_gaxpy(getJ(), v, res);
	}
	/**
	 * res += J' * u
	 */
	void add_Jtu(double *res, const double* u){
		cs_gaxpy(getJt(), u, res);
	}
	
	//! returns upper half of JtJ. First entry per column corresponds to diagonal entry
	const cs* get_JtJ();
	
	std::pair<int,int> getSize() const {
		int M = measurements.getDim();
		int N = variables.getDim();
//		if(_addDiagonal) M += N;
		
		return std::make_pair(M, N);
	}
	
	bool initialize(bool addDiagonal) {
		bool ret = _structureChanged;
		_addDiagonal = addDiagonal;
		_structureChanged = false;
		if(ret) {
			freeWorkspace();
			createSparse();
			initCovariance();
			_valuesChanged = true;
		}
		return ret;
	}
	
	void apply_delta(const Eigen::VectorXd &delta, double scale);
	
	void store_or_restore(bool store){
		if(store){
			for(RVList::iterator v = variables.begin(); v!= variables.end(); ++v){
				v->store();
			}
			_valuesChanged = true;
			freeTemporaries();
		}
		else
		{
			for(RVList::iterator v = variables.begin(); v!= variables.end(); ++v){
				v->restore();
			}
		}
	}
	//@}
	
public:
	//! \name Inspection functions
	//@{
	/// Returns the dimension of the measurement space
	int getM() const {
		return measurements.getDim();
	}
	
	/// Returns the dimension of the variable space
	int getN() const {
		return variables.getDim();
	}
	
	/// Returns the RSS of measurement m
	template<class Measurement>
	static double measurementRSS(const Measurement &m)
	{
		Eigen::Matrix<double, Measurement::DIM, 1> res;
		m.eval(res);
		return res.squaredNorm();
	}
	
	/// applies Functor f to all measurements which depend on variable with given id
	template<class Manifold, class Functor>
	static void traverse_measurements(const VarID<Manifold> &id, Functor &f)
	{
		
		const IRVHolder::measurement_container &measurements = id.ptr->measurements;
		for(IRVHolder::measurement_container::const_iterator it = measurements.begin(); it != measurements.end(); ++it)
		{
			f(*it);
		}
	}
	
	/// gets the RSS of all Measurements of given type and returns their total dimension
	template<class Measurement>
	int get_RSS(double &rss) const
	{
		rss = 0;
		int sum_dim = 0;
		for(MeasurementList::const_iterator it = measurements.begin(); it != measurements.end(); ++it)
		{
			const IMeasurement_Holder &mh = *it;
			typedef typename IMeasurement_Holder::holder<Measurement> m_type;
			const m_type *m = dynamic_cast<const m_type*>(&mh);
			if(m!=0)
			{
				Eigen::Matrix<double, Measurement::DIM, 1> res;
				m->eval(res.data(), false);
				rss += res.squaredNorm();
				sum_dim += Measurement::DIM;
			}
		}
		return sum_dim;
	}
	
	/// gets the RSS of all Measurements depending on variable with given ID.
	template<class Manifold>
	static int get_RSS(double &rss, const VarID<Manifold> &id) {
		const internal::IRVHolder::measurement_container& meas = id.ptr->measurements;
		
		rss = 0;
		int sum_dim = 0;
		for(internal::IRVHolder::measurement_container::const_iterator m = meas.begin(); m!= meas.end(); ++m){
			int dim = (*m)->getDim();
			sum_dim += dim;
			double res[dim];
			(*m)->eval(res, false);
			rss = std::inner_product(res, res + dim, res, rss);
		}
		return sum_dim;
	}
	
	/// returns the index of given variable
	template<class Manifold>
	static int get_index(const VarID<Manifold> &id)
	{
		return id.ptr->idx;
	}
	//@}
};



}  // namespace SLOM

#endif /* SPARSEFUNCTION_HH_ */
