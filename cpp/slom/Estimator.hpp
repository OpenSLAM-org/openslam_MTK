/*
 *  Copyright (c) 2008--2015, Universitaet Bremen
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
 * @file slom/Estimator.hpp
 * @brief Defines the Estimator class
 */
#ifndef ESTIMATOR_H_
#define ESTIMATOR_H_

#include "src/SparseFunction.hpp"

#include "solvers/CholeskySolver.hpp"

#include "algorithms/GaussNewton.hpp"
#include "algorithms/LevenbergMarquardt.hpp"

#include "preconditioner/Preconditioner.hpp"

#include <boost/scoped_ptr.hpp>


#include "Covariance.hpp"

namespace SLOM {

#define SLOM_ESTIMATOR_CHANGE(what, name) \
		void change ## what (boost::scoped_ptr<what> &new ## what) { \
			name.swap(new ## what); \
			resetAlgoAndSolver(); \
		} \
		void change ## what(what *new ## what) { \
			boost::scoped_ptr<what> scoped ## what(new ## what); \
			change ## what (scoped ## what); \
		} \
		const what & get ## what() const { \
			return *name; \
		}



class CallBack;
class scoped_callback;

/**
 * Estimator class, which provides the user interface for variable and measurement handling as well as optimization.
 */
class Estimator
{
	/**
	 * Data storage of the variables and the measurement functions.
	 * Combined they form a sparse function. All modifications of variables 
	 * and function evaluation shall be done via this member.
	 */
	SparseFunction func;
	
	typedef SparseFunction::SparseType SparseType;
	
public:
	/** use the following "damping term" @f$ N @f$ in @f$ (J\trans J + N)\delta = J\trans (y - f(\beta)) @f$
	 *  @deprecated Algorithms are specified by separate classes now
	 */
	enum AlgorithmE { 
		GaussNewton        //! @f$N=0@f$
		,Levenberg         //! @f$N= \lambda*I@f$
		,LevenbergMarquardt //! @f$N= \lambda*diag(J^T J)@f$
	};
	
private:
	CallBack *callback; // TODO make this a scoped_ptr?

	// TODO maybe these can be members of SparseFunction
	boost::scoped_ptr<Solver> solver;
	boost::scoped_ptr<Algorithm> algo;
	boost::scoped_ptr<Preconditioner> precond;
	
	
public:
	
	//! \name Construction and initialization
	//@{
	/**
	 * Construct an Estimator.
	 * 
	 * @param alg   Choose Algorithm. Currently SLoM supports GaussNewton, Levenberg and LevenbergMarquardt
	 * @param lamda0 Initial lamda for Levenberg(-Marquardt)
	 * @deprecated Pass Algorithm as an object pointer instead
	 */
	MTK_DEPRECATED
	Estimator(AlgorithmE alg, double lamda0=1e-3, Solver* solver_ = new CholeskySolver<>);
	
	Estimator(Algorithm* alg = new SLOM::GaussNewtonAlgorithm(), Solver* solver = new CholeskySolver<>()) :
		callback(0),
		solver(solver),
		algo(alg)
	{
		func.callback = 0;
		resetAlgoAndSolver();
	}
	
	/**
	 * Set the CallBack object, to get status informations
	 */
	void setCallBack(CallBack &cb){
		//FIXME also set algo and solver callback?
		func.callback = callback = &cb;
	}
	
	CallBack* getCallback() const {
		return callback;
	}
	
	//! Destructor (what else to say ...)
	~Estimator();
	
	
	/**
	 * Deletes all temporary data. 
	 * Measurements, random variables and the algorithms state are kept unmodified,
	 * i.e. it is possible to continue the algorithm afterwards.
	 */
	void cleanWorkspace() {
		func.cleanWorkspace();
	}
	
	/**
	 * Call this method if you manually changed data, which is only indirectly used by the measurements.
	 */
	void dataChanged() {
		func.valuesChanged();
	}
	
	/**
	 * Initialize Jacobian and symbolical decomposition.
	 * @deprecated You do not need to call this method yourself, because the 
	 *             algorithm decides itself when initialization is necessary
	 */
	MTK_DEPRECATED
	void initialize() { }
	//@}
private:
	
	
	/**
	 * re-init pointers for interaction of algorithm and solver with each other and with function
	 */
	void resetAlgoAndSolver() {
		assert(solver && "Require valid Solver");
		assert(algo && "Require valid Algorithm");
		
		solver->func = algo->func = &func;
		
		if(precond)
			precond->func = &func;
		
		algo->solver = func.solver = solver.get();
		func.precond = precond.get();
		func.algo = algo.get();
	}
	
	
public:
	//! \name Optimization functions
	//@{
	/**
	 * Performs a single optimization step.
	 * @param linear_refine refine step length after finding optimization delta.
	 *        This parameter is highly experimental
	 */
	double optimizeStep(bool linear_refine = false);
	
	/**
	 * Function to change from between GaussNewton, Levenberg and LevenbergMarquardt
	 * @deprecated The algorithm is out-factored into the Algorithm classes
	 */
	MTK_DEPRECATED
	void changeAlgorithm(AlgorithmE algo, double lamdaNew=-1){
		switch(algo){
		case LevenbergMarquardt:
			std::cerr << "Warning only Levenberg is supported\n";
			// no break
		case Levenberg:
			if(lamdaNew < 0) {
				std::cerr << "Warning, you need to set lambda when using deprecated changeAlgorithm, setting lambda=1.0";
				lamdaNew = 1.0;
			}
			changeAlgorithm(new SLOM::LevenbergMarquardt(lamdaNew));
			break;
		case GaussNewton:
			changeAlgorithm(new SLOM::GaussNewtonAlgorithm());
			break;
		}
		resetAlgoAndSolver();
	}
	
	SLOM_ESTIMATOR_CHANGE(Algorithm, algo)
	
	/**
	 * Change the Solver, by swapping scoped pointers. The old Solver can by analyzed afterwards.
	 * @param newSolver
	 */
	SLOM_ESTIMATOR_CHANGE(Solver, solver)

	SLOM_ESTIMATOR_CHANGE(Preconditioner, precond)
	//@}
	
	//! \name Inspection functions of the optimization state
	//@{
	/**
	 * returns the last residual sum of squares.
	 * @deprecated Usage is discouraged, use getRSS() instead
	 */
	MTK_DEPRECATED
	double getLastRSS() const {
		assert(false && "called deprecated function");
		return 0;
	}
	/**
	 * Gets the current residual sum of squares.
	 * Internally it updates the residual vector if necessary.
	 */
	double getRSS()
	{
		return func.getRSS();
	}
	
	double getRMS() {
		return std::sqrt(getRSS()/getM());
	}
	
	/**
	 * Gets the squared norm of the gradient
	 * @return
	 */
	double getSquaredNormOfGradient() {
		return func.getSquaredNormOfGradient();
	}
	
	const Eigen::VectorXd& getGradient() {
		return func.getGradient();
	}
	
	/**
	 * Return the current residuum vector.
	 * This is a reference to the internal residuum, so it's value is likely to 
	 * change while the algorithm is running.
	 */
	const Eigen::VectorXd& getRes(){
		return func.getResiduum();
	}
	
	/**
	 * Return @f$diag(J\trans J) @f$.
	 * This is a reference to an internal variable, so it's value is likely to 
	 * change while the algorithm is running.
	 */
	const Eigen::VectorXd& getCholCovariance() const {
		return func.cholCovariance;
	}
	
	/**
	 * Current lamda of Levenberg(-Marquardt) algorithm.
	 */
	MTK_DEPRECATED
	double getLamda () const {
		assert(false && "You have to obtain lambda from your algorithm directly now");
		return 0.0/0.0;
	}
	
	/**
	 * Dump Jacobian as triplet form.
	 */
	std::ostream& dumpJacobian(std::ostream& out);
	
	const SparseType& get_JtJ() {
		return func.get_JtJ();
	}

	const SparseType& getJ() {
		return func.getJ();
	}
	
	/**
	 * prints the current Jacobian
	 * @deprecated Use dumpJacobian(std::ostream&) instead
	 */
	void printJacobian(bool /*brief*/=false) const __attribute__((deprecated)) {
		//cs_print(func.jacobian, brief);
	}

	//@}
	
	
	/**
	 * \name Register and Modify Variables and Measurements
	 * These functions are passed through to @ref SparseFunction
	 */ 
	//@{
	
	/**
	 * Adds a new RandomVariable based on Manifold m into the Estimator. 
	 * Returns an id. The id needs to be passed to the corresponding measurements.
	 * Via the id it is also possible, to obtain the current value of the variable (using operators * and -> )
	 */ 
	template<class Manifold>
	VarID<Manifold> insertRV(const Manifold& m, bool optimize = true) {
		return func.insertRV(m, optimize);
	}

	/**
	 * Reinitialize value of RV with given id.
	 * Manifold2 must be the same as or convertible to Manifold
	 */
	template<class Manifold, class Manifold2>
	VarID<Manifold> reinitRV(VarID<Manifold> id, const Manifold2& m) {
		return func.reinitRV(id, Manifold(m));
	}
	
	/**
	 * Change the optimization status of a variable later on.
	 * It is possible to optimize a variable which was previously fixed (with optimize==true)
	 * or to fix a variable which was previously optimized (with optimize == false).
	 */
	template<class Manifold>
	VarID<Manifold> optimizeRV(VarID<Manifold> id, bool optimize=true) {
		return func.optimizeRV(id, optimize);
	}

	/**
	 * Removes the RandomVar with @a id from the estimator.
	 * Returns true on success. If the RV is still used by measurements, 
	 * it returns false.
	 * If removing is not possible one can always disable optimization using @ref optimizeRV
	 * 
	 * @todo maybe add a version which throws instead of returning @c false
	 */
	template<class Manifold>
	bool removeRV(const VarID<Manifold> &id){
		bool success = func.removeRV(id);
		if(!success)
			std::cerr << "Tried to remove a variable still in use!" << std::endl;
		return success;
	}
	
	/**
	 * Inserts a new measurement. 
	 * The Measurement itself registers the variables it depends on. 
	 */
	template<class Measurement>
	MeasID<Measurement> insertMeasurement(const Measurement &m) {
		return func.insertMeasurement(m, internal::InvDeviation<void>());
	}
	
	/**
	 * Inserts a new measurement with additional covariance information. 
	 */
	template<class Measurement, class DevType>
	MeasID<Measurement> insertMeasurement(const Measurement &m, const internal::InvDeviation<DevType>& cov) {
		return func.insertMeasurement(m, cov);
	}
	
	/**
	 * Removes the Measurement with id from the estimator.
	 * Returns true on success (should always be the case).
	 */
	template<class Measurement>
	bool removeMeasurement(const MeasID<Measurement> &id){
		return func.removeMeasurement(id);
	}
	//@}
	
	
	/**
	 * \name Get status of the sparse function.
	 */
	//@{
	
	/**
	 * Get Number of Non-Zeroes in the Jacobian.
	 */
	int getNNZ() const {
		return func.nnz;
	}
	
	/**
	 * Get Dimension of Measurement Space
	 */
	int getM() const {
		return func.getM();
	}
	
	/**
	 * Get Dimension of Variable Space
	 */
	int getN() const {
		return func.getN();
	}
	/**
	 * get the RSS for all measurements of type Measurement (template parameter)
	 */
	template<class Measurement>
	int getRSS(double &rss) const {
		return func.get_RSS<Measurement>(rss);
	}
	template<class Measurement>
	MTK::vectview<const double, Measurement::DIM>
	getRes(const SLOM::MeasID<Measurement>& id) {
		return func.getRes(id);
	}
	template<class Measurement>
	double
	getRSS(const SLOM::MeasID<Measurement>& id) {
		return func.getRes(id).squaredNorm();
	}
	/**
	 * get the RSS of all measurements depending on Variable id
	 */
	template<class Manifold>
	int getRSS(double &rss, const VarID<Manifold> &id) const {
		return func.get_RSS(rss, id);
	}
	
	/// call Functor on every measurement depending on id
	template<class Functor, class Manifold>
	void traverse_measurements(const SLOM::VarID<Manifold>& id, Functor f) const {
		return SparseFunction::traverse_measurements(id, f);
	}

	//@}
};

}  // namespace SLOM

#endif /*ESTIMATOR_H_*/
