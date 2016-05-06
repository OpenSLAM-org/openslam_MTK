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
#include "RVHolder.hpp"
#include "MeasurementHolder.hpp"
#include <Eigen/Core>

#include "Sparse.hpp"

#include <deque>
#include <iostream>
#include <numeric>

#include <bitset>


namespace SLOM {

class CallBack;
class Solver;
class Algorithm;
class Preconditioner;


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
	
	bool operator!() const {
		return base::ptr == 0;
	}
	
	// deprecated, because it causes trouble with automatically promoting to bool->int->double
	// use !! if required
	MTK_DEPRECATED
	operator bool() const {
		return !(!*this);
	}
	
	using base::index;
	bool operator==(const VarID<M>& y) const {
		return base::ptr == y.ptr;
	}
	bool operator!=(const VarID<M>& y) const {
		return base::ptr != y.ptr;
	}
	
	bool getsOptimized() const {
		return base::ptr && base::ptr->optimize;
	}
	template<class M2>
	bool operator<(const VarID<M2>& v) const { return base::index() < v.index(); }
	
	const M& operator*()  const { return  base::ptr->backup;}
	const M* operator->() const { return &base::ptr->backup;}
};


template<class Measurement>
class MeasID : private internal::MeasurementList::id<Measurement> {
	typedef internal::MeasurementList::id<Measurement> base;
	MeasID(internal::Measurement_Holder<Measurement> &m) : base(&m) {}
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
	
	const internal::Measurement_Holder<Measurement>& getHolder() const {
		return *base::ptr;
	}
	
	// This allows to change data members of a measurement
	template<class Data>
	Data& operator->*(Data Measurement::* data) {
		return (base::ptr->m_).*data;
	}
	
private:
	template<class Var, int idx>
	void operator->*(internal::IMeasurement_Holder::VarRef<Var, idx> Measurement::* var) const {
		assert(false && "changing variables not yet supported! Remove measurement and add it again");
		(void) var;
		// TODO return some proxy, which allows changing VarRefs
	}
	
};


class SparseFunction {
	
	friend class Estimator;
	friend class Solver;
	friend class SparseInterface;
	friend class Algorithm;
	
	CallBack *callback;
	// FIXME solver and algo and precond could be replaced by a list of SparseInterface pointers
	Solver * solver;
	Algorithm * algo;
	Preconditioner * precond;

	typedef internal::RVList RVList;
	typedef internal::IRVHolder IRVHolder;
	
	typedef internal::MeasurementList MeasurementList;
	typedef internal::IMeasurement_Holder IMeasurement_Holder;
	
	// TODO Eventually this shall get replaced by some kind of block-matrix
	typedef Eigen::SparseMatrix<double> SparseType;
	typedef Eigen::VectorXd             VectorType;
	
	//! \name Storage of variables and measurements
	//@{
	RVList variables;
	RVList fixed_variables;
	
	internal::MeasurementList measurements;
	//@}
	
	//! \name Internal helper variables
	//@{
	//! List of entities which might or might not be up to date:
	enum entities {
		ent_residuum, ent_newRes, ent_gradient,
		ent_J, ent_Jt, 
		ent_JtJ, 
		ent_BlockJ, ent_BlockJtJ, 
		
		ent_jacobiPreconditioner /* might get obsolete */, 
		
		ent_numberOfEntities
	};
	//! keep track whether certain structure and/or value entities are up to date
	std::bitset<ent_numberOfEntities> structureUpToDate, valuesUpToDate;
	
	
	bool backup_valid;
	
	void updated(entities ent) {
		structureUpToDate[ent] = valuesUpToDate[ent] = true;
	}
	
	/**
	 * This method is to be called whenever any value changed.
	 */
	void valuesChanged();
	
	/**
	 * This method is to be called whenever the structure changes.
	 * This requires re-computation of, e.g., symbolical decompositions
	 */
	void structureChanged();
	
	//! delta for numerical calculation of Jacobian
	double numerical_delta, half_delta_inv;
	
	//! Number of non-zeroes in the jacobian
	int nnz;
	//! sparse Jacobian of the function
	SparseType jacobian;
	
	// diagonal of $J^\top J$
	VectorType cholCovariance;
	
	
	//! current residuum (at backup position)
	VectorType residuum;
	double RSS;
	
	//! residuum at working position:
	VectorType newResiduum;
	double newRSS;
	
	//! current gradient;
	VectorType gradient;
	double norm2_gradient;
	
	
	// helper variables
	//! over-approximation of non-zeroes in JtJ
	int nnz_JtJ;
	// FIXME Jt might not be needed anymore
	SparseType Jt; // transpose of jacobian
	SparseType JtJ;  // Jt * J
	//@}
	
	//! \name Construction and destruction
	//@{
	/**
	 * Construct a SparseFunction object
	 * @param delta delta for numerical differentiation
	 */
	explicit
	SparseFunction(double delta = 1e-6) : 
		callback(0),
		solver(0),
		algo(0),
		precond(0),
		numerical_delta(delta), half_delta_inv(0.5/delta),
		nnz(0), nnz_JtJ(0)//,
//	  jacobian(0), Jt(0), JtJ(0) 
	{
		structureChanged();
		assert(delta > 0 && "delta must be positive");
		backup_valid = true; // true because there are no entries
	}
	
	
	~SparseFunction() {
	}
	
	
	void cleanWorkspace() {
		// mark all structure relevant things as invalid:
		structureChanged();
		
		newResiduum = residuum = gradient = cholCovariance = VectorType();
		jacobian = Jt = JtJ = SparseType();
	}
	//@}
	
	CallBack* getCallback() const {
		return callback;
	}
	
	/**
	 * \name Variable Managment
	 * Methods for inserting and removing of variables and measurements.
	 * Methods are passed encapsulated by Estimator.
	 * @{
	 */ 
	
	template<class Manifold>
	VarID<Manifold> insertRV(const Manifold& m, bool optimize = true){
		RVList &list = optimize ? variables : fixed_variables; 
		internal::RVHolder<Manifold>* ptr = new internal::RVHolder<Manifold>(m, optimize);
		VarID<Manifold> id = list.insert(ptr);
		ptr->isRegistered = true;
		if(optimize) {
			structureChanged();
			// JtJ (if calculated) gets a new diagonal block:
			nnz_JtJ += Manifold::DOF * (Manifold::DOF + 1) / 2; 
		}
		return id;
	}
	
	template<class Manifold>
	VarID<Manifold> reinitRV(VarID<Manifold> id, const Manifold& m) {
		valuesChanged();
		id.ptr->reset(m);
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
		
		structureChanged();
		return new_list.insert(id.ptr);
		
	}
	
	template<class Manifold>
	bool removeRV(const VarID<Manifold> &id){
		if(id.ptr->has_measurements()) {
			return false;
		}
		// move variable to fixed_variables (if not done already) this also updates nnz counts
		VarID<Manifold> id_ = optimizeRV(id, false);
		
		fixed_variables.remove(id_);
		return true;
	}
	
	
	template<class Measurement, class DevType>
	MeasID<Measurement> insertMeasurement(const Measurement &m, const internal::InvDeviation<DevType>& dev) {
		MeasID<Measurement> id = measurements.insert(new internal::Measurement_Holder_with_Dev<Measurement, DevType>(m, dev));
		nnz += id.ptr->registerVariables();
		
		// over-approximation for nnz in JtJ, TODO if it really matters, calculate correct number when defining Measurements
		nnz_JtJ += (Measurement::DEPEND - 1) * Measurement::DEPEND / 2;
		
		structureChanged();
		
		return id;
	}

	
	
	/**
	 * @warning NOT PROPERLY TESTED
	 * @param id
	 * @return
	 */
	template<class Measurement>
	bool removeMeasurement(MeasID<Measurement> id){
		nnz -= id.ptr->unregisterVariables();
		nnz_JtJ -= (Measurement::DEPEND - 1) * Measurement::DEPEND / 2;
		measurements.remove(id);
		
		structureChanged();
		
		return true;
	}
	//@}
	
	
	
	//! \name Methods needed for optimization algorithms
	//@{
	/**
	 * Creates "the big matrix",
	 * 
	 * @deprecated this shall be replace by block matrix entirely 
	 */
	void createSparse();
	
	//! @deprecated mark as __attribute__((deprecated))
	void initCovariance();
	
	/**
	 * Numerically calculates the Jacobian of the function with the 
	 * and updates cholCovariance. The result is stored in matrix.
	 */
	void calculateJacobian();
	
	void calculateBlockJacobian();
	
	/**
	 * Evaluate the function, store result to vector starting at result and return the RSS.
	 */
	double evaluate(VectorType & result) const;
	
	
	/**
	 * get the current Jacobian. If neither values nor structure has changed since last call,
	 * no recalculation is done.
	 * If add_diagonal = true, each column contains an extra entry for diagonals
	 * @deprecated This MIGHT get deprecated in a very future version
	 */ 
	const SparseType& getJ(){
		calculateJacobian();
		
		return jacobian;
	}
	
	//! returns the residuum at the back-up position
	const VectorType& getResiduum() {
		computeResiduum();
		return residuum;
	}
	//! gets the residual sum of squares (i.e. squared norm of residuum)
	double getRSS(){
		computeResiduum();
		return RSS;
	}
	
	//! gets the RSS at working position
	double getNewRSS() {
		computeNewRes();
		return newRSS;
	}
	
	void computeNewRes() {
		if(valuesUpToDate[ent_newRes]) return;
		newResiduum.resize(getM());
		newRSS = evaluate(newResiduum);
		updated(ent_newRes);
	}
	
	void computeResiduum() {
		if(valuesUpToDate[ent_residuum]) return;
		assert(backup_valid && "Can't compute RSS with invalid working copy");
		residuum.resize(getM());
		RSS = evaluate(residuum);
		updated(ent_residuum);
	}
	
	const VectorType& getGradient(){
		computeGradient();
		return gradient;
	}
	double getSquaredNormOfGradient(){
		computeGradient();
		return norm2_gradient;
	}
	
	void computeGradient() {
		if(valuesUpToDate[ent_gradient]) return;
		gradient.setZero(getN());
		add_Jtu(gradient, getResiduum());
		norm2_gradient = gradient.squaredNorm();
		updated(ent_gradient);
	}
	
	void updateJacobiPreconditioner();
	
	void applyJacobiPreconditioner(VectorType &vec, bool transpose);
	
	void updateBlockJtJ();
	
	void applyCholIncPreconditioner(VectorType &vec, bool transpose);
	
	void applyPreconditioner(VectorType &vec, bool transpose);
	
	//! @deprecated This feature might be removed
	const VectorType& getCholCov() {
		// make sure cholCovariance is up-to-date:
		getJ();
		return cholCovariance;
	}
	/**
	 * @todo hopefully, this function will not be necessary
	 * @deprecated This function might be removed in future versions.
	 */
	//MTK_DEPRECATED
	const SparseType& getJt() {
		if(!valuesUpToDate[ent_Jt]) 
			Jt = getJ().transpose();
		return Jt;
	}
	
	/**
	 * res += J * v
	 */
	void add_Jv(VectorType& res, const VectorType& v){
		res += getJ() * v;
	}
	/**
	 * res += J' * u
	 */
	void add_Jtu(VectorType& res, const VectorType& u_){
		assert(u_.rows() == getM() && res.rows() == getN() && "Vector dimensions must agree with Jacobian");
#ifdef SLOM_EXPERIMENTAL
		getJt();
		for(MeasurementList::iterator it= measurements.begin(); it != measurements.end(); ++it){
			u = it->addJtu(res, u);
		}
#else /* SLOM_EXPERIMENTAL */
		res += getJ().transpose() * u_;
//		cs_gaxpy(getJt(), u, res);
#endif /* SLOM_EXPERIMENTAL */
	}
	
	//! returns upper half of JtJ. First entry per column corresponds to diagonal entry
	const SparseType& get_JtJ();
	
	std::pair<int,int> getSize() const {
		int M = measurements.getDim();
		int N = variables.getDim();
//		if(_addDiagonal) M += N;
		
		return std::make_pair(M, N);
	}
	
	/**
	 * Applies the delta vector (multiplied by scale) to the backup and stores it to the working copy
	 * @param delta
	 * @param scale
	 * @return  Residual sum of squares at the new working position
	 */
	double apply_delta(const VectorType &delta, double scale);
	
	
	/**
	 * If store is true, it stores the working copy to backup, otherwise
	 * copies the backup to the working copy.
	 * @param store
	 */
	void store_or_restore(bool store){
		if(store){
			for(RVList::iterator v = variables.begin(); v!= variables.end(); ++v){
				v->store();
			}
			RSS = getNewRSS();
			assert(valuesUpToDate[ent_newRes] && "new residuum was not evaluated before storing");
			valuesChanged();
			residuum.swap(newResiduum);
			updated(ent_residuum);
		}
		else
		{
			for(RVList::iterator v = variables.begin(); v!= variables.end(); ++v){
				v->restore();
			}
		}
		backup_valid = true;
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
			typedef internal::Measurement_Holder<Measurement> m_type;
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
	
	template<class Measurement>
	MTK::vectview<const double, Measurement::DIM>
	getRes(const SLOM::MeasID<Measurement>& id) {
		return MTK::vectview<const double, Measurement::DIM>(getResiduum().data() + id.index(), id->dim);
	}
	
	/// gets the RSS of all Measurements depending on variable with given ID.
	template<class Manifold>
	int get_RSS(double &rss, const VarID<Manifold> &id) const {
		const internal::IRVHolder::measurement_container& meas = id.ptr->measurements;
		
		// allocate temporary vector for evaluation:
		VectorType res_tmp(getM());
		double* res = res_tmp.data();
		for(internal::IRVHolder::measurement_container::const_iterator m = meas.begin(); m!= meas.end(); ++m) {
			res = (*m)->eval(res, false);
		}
		int dim = res - res_tmp.data();
		rss = res_tmp.head(dim).squaredNorm();
		return dim;
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
