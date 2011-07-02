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
 * @file slom/Estimator.hpp
 * @brief Defines the Estimator class
 */
#ifndef ESTIMATOR_H_
#define ESTIMATOR_H_

#include "src/SparseFunction.hpp"



namespace SLOM {

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
	
public:
	//! use the following "damping term" @f$ N @f$ in @f$ (J\trans J + N)\delta = J\trans (y - f(\beta)) @f$
	enum Algorithm{ 
		GaussNewton        //! @f$N=0@f$
		,Levenberg         //! @f$N= \lambda*I@f$
		,LevenbergMarquardt //! @f$N= \lambda*diag(J^T J)@f$
	};
	
	/**
	 * Linear solver
	 * @deprecated currently only Cholesky is working
	 */
	enum Solver{
		QR, // QR is not supported at the moment (and maybe will never be again)
		Cholesky
	};
private:
	friend class CallBack;
	CallBack *callback;
	
	
	enum Algorithm usedAlgorithm;
	enum Solver usedSolver;
	
	css* symbolic; // symbolic decomposition of jacobian or JtJ 
	csn* numeric;  // numeric decomposition of jacobian or JtJ
	
	// the current residuum:
	Eigen::VectorXd res, newRes;
	
	// workspace for solving:
	Eigen::VectorXd workspace;
	
	// current delta vector
	Eigen::VectorXd delta;
	
	
	// lamda parameter for LMA:
	double lamda; 
	
	
	/** the last Residual Sum of Squares
	 */
	double lastRSS;
	
public:
	
	//! \name Construction and initialization
	//@{
	/**
	 * Construct an Estimator.
	 * 
	 * @param alg   Choose Algorithm. Currently SLoM supports GaussNewton, Levenberg and LevenbergMarquardt
	 * @param lamda0 Initial lamda for Levenberg(-Marquardt)
	 */
	Estimator(Algorithm alg=GaussNewton, double lamda0=1e-3) : 
		callback(0),
		usedAlgorithm(alg), usedSolver(Cholesky),
		symbolic(0), numeric(0), lamda(lamda0)
	{};
	
	/**
	 * Construct Estimator, also choosing linear solver.
	 * @deprecated QR Solving is not supported anymore
	 */
	Estimator(Solver solver, Algorithm alg=GaussNewton, double lamda0=1e-3) :
		callback(0),
		usedAlgorithm(alg), usedSolver(solver),
		symbolic(0), numeric(0), lamda(lamda0)
	{
		assert(usedSolver != QR && "usage of QR is deprecated now, please use Cholesky instead");
	}
	/**
	 * Set the CallBack object, to get status informations
	 */
	void setCallBack(CallBack &cb){
		callback = &cb;
	}
	//! Destructor (what else to say ...)
	virtual ~Estimator();
	
	/**
	 * Initialize Jacobian and symbolical decomposition.
	 * @deprecated You do not need to call this method yourself, because the 
	 *             algorithm decides itself when initialization is necessary
	 */
	void initialize()  __attribute__((deprecated)) {
		init();
	}
	//@}
private:
	
	//! \name internal methods of the algorithm
	//@{
	/**
	 * Evaluate the function, store result to vector starting at result and return the RSS.
	 */
	double evaluate(double *result) const;
	
	/**
	 * Internal method for initialization.
	 */
	void init();
	
	
	
	/**
	 * calculates delta = (J^T*J + N)^-1 * J^T * res, 
	 * with damping term N depending on current algorithm.
	 */
	void calculate_delta();
	
	
	/**
	 * add scale * delta to variables and evaluate new residuum to workspace
	 */
	double apply_delta(double scale = -1);
	
	/**
	 * Stores the modified variables or restores the old ones depending on newRSS
	 * and algorithm. Returns gain = (oldRSS - newRSS)/(newRSS)
	 */
	double store_or_restore(double newRSS, scoped_callback &cb);
	
	void freeWorkspace();
	bool qrSolve(double* delta, const double *res);
	bool choleskySolve(double* delta, const double *res);
	
	/**
	 * Update the sparse Matrix depending on the current algorithm
	 */
	void updateSparse();
	//@}
	
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
	 * @note The algorithm will be out-factored in a future versions
	 */
	void changeAlgorithm(Algorithm algo, double lamdaNew=-1){
		usedAlgorithm = algo;
		if(lamdaNew > 0) lamda = lamdaNew;
		freeWorkspace();
		init();
	}
	//@}
	
	//! \name Inspection functions of the optimization state
	//@{
	/**
	 * returns the last residual sum of squares.
	 * @deprecated Usage is discouraged, use getRSS() instead
	 */
	double getLastRSS() const {
		return lastRSS;
	}
	/**
	 * Gets the current residual sum of squares.
	 * Internally it updates the residual vector if necessary.
	 */
	double getRSS()
	{
		init();
		if(lastRSS < 0){
			lastRSS = evaluate(res.data());
		}
		return lastRSS;
	}
	
	/**
	 * Return the current residuum vector.
	 * This is a reference to the internal residuum, so it's value is likely to 
	 * change while the algorithm is running.
	 */
	const Eigen::VectorXd& getRes(){
		getRSS();
		return res;
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
	double getLamda () const {
		return lamda;
	}
	
	/**
	 * Dump Jacobian as triplet form.
	 */
	std::ostream& dumpJacobian(std::ostream& out);
	
	/**
	 * prints the current Jacobian
	 * @deprecated Use dumpJacobian(std::ostream&) instead
	 */
	void printJacobian(bool brief=false) const __attribute__((deprecated)) {
		cs_print(func.jacobian, brief);
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
	 */
	template<class Manifold>
	VarID<Manifold> reinitRV(VarID<Manifold> id, const Manifold& m) {
		return func.reinitRV(id, m);
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
		return func.insertMeasurement(m);
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
	/**
	 * get the RSS of all measurements depending on Variable id
	 */
	template<class Manifold>
	static int getRSS(double &rss, const VarID<Manifold> &id) {
		return SparseFunction::get_RSS(rss, id);
	}
	//@}
};

}  // namespace SLOM

#endif /*ESTIMATOR_H_*/
