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
 * @file slom/src/SparseFunction.cpp
 * 
 * @brief Implementation of SparseFunction.
 */

#include "SparseFunction.hpp"

#include "../TicToc.hpp"
#include "../CallBack.hpp"
#include "../Solver.hpp"
#include "../Algorithm.hpp"
#include "../preconditioner/Preconditioner.hpp"
#include <typeinfo>

namespace SLOM {



double SparseFunction::evaluate(VectorType& result) const{
	CALLBACK_START(evaluate);
	MeasurementList::const_iterator it = measurements.begin();
	assert(result.rows()>=getM() && "Vector not big enough to hold result");
	double sum = 0;
	double *ptr = result.data();
	for(; it != measurements.end(); it++){
		assert(result.data() + it->idx == ptr);
		double * end = it->eval(ptr, false);
		for( ; ptr < end; ptr++){
			sum += std::pow(*ptr, 2);
		}
	}
	return sum;
}




void SparseFunction::createSparse(){
	if(structureUpToDate[ent_J]) return;
	CALLBACK_START(createSparse);
	initCovariance();
	// FIXME this method is not required anymore for new Block-Jacobian
	int M = measurements.getDim();
	int N = variables.getDim();
	int size = nnz;
	//allocate a compressed column matrix:
	jacobian.resize(M, N);
	jacobian.resizeNonZeros(size);
	
	//create the structure of the matrix:
	int *cIdx=jacobian.outerIndexPtr(); // column pointer
	int *rIdx=jacobian.innerIndexPtr(); // row indices
	*cIdx = 0; //first column starts at 0.
	
	for(RVList::const_iterator var = variables.begin(); var!= variables.end(); ++var){
		int vDOF = var->getDOF();
		assert(vDOF>0);
		if(!var->has_measurements()){
			std::cerr << "Warning: No measurement for Variable " 
					<< var->idx << " " << typeid(*var).name()<< std::endl;
		}
		for(IRVHolder::measurement_container::const_iterator meas= var->measurements.begin(); meas!= var->measurements.end(); ++meas){
			int row = (*meas)->idx;
			for(int mDim = (*meas)->getDim(); mDim>0; mDim--){
				*rIdx++ = row++;
			}
		}
		*++cIdx = rIdx - jacobian.innerIndexPtr(); //start of next column
		while(--vDOF > 0){ // copy the current column for each DOF of the variable
			rIdx = std::copy(jacobian.innerIndexPtr() + cIdx[-1], // start of previous column
					rIdx, rIdx);
			*++cIdx = rIdx - jacobian.innerIndexPtr();
		}
	}
	assert(*cIdx == size);
	structureUpToDate[ent_J] = true;
}


void SparseFunction::initCovariance(){
	int n=variables.getDim();
	cholCovariance.resize(n);
	cholCovariance.fill(0);
}


void SparseFunction::calculateBlockJacobian(){
#ifdef SLOM_JACOBI_BLOCKS
	if(!valuesUpToDate[ent_BlockJ]) {
		CALLBACK_START(calculateBlockJacobian);
		MeasurementList::iterator it = measurements.begin();
		for( ; it != measurements.end(); ++it){
			it->updateJacobian();
		}
		updated(ent_BlockJ);
	}
#endif
}


void SparseFunction::calculateJacobian(){
	calculateBlockJacobian(); // FIXME make it possible to update block Jacobian from outside

	// FIXME this shall be replaced by assembling the entries of the block-Jacobian
	if(valuesUpToDate[ent_J]) return;
//	std::cerr << "Warning: called " << __PRETTY_FUNCTION__ << "\n";
	CALLBACK_START(calculateJacobian)
	createSparse();  // update structure if necessary
	int m = jacobian.rows(), n = jacobian.cols();
	double *x = jacobian.valuePtr();
	
	assert(m == measurements.getDim());
	assert(n == variables.getDim());
	
	double *chol = cholCovariance.data();
	// FIXME a little bit conservative workspaces; will be redundant with symbolic derivations/block-Jacobians
	VectorType delta(n);
	VectorType workspace(m);
	delta.fill(0);
	double *add = delta.data(); // temp-array for adding
	
	for(RVList::iterator var = variables.begin(); var!= variables.end(); ++var){
		assert(var->optimize);
		int vDOF = var->getDOF();
		assert(vDOF>0);
		for(int k=0; k<vDOF; k++){
			add[k] = numerical_delta;
			
			// store $f(\mu \mplus 1/d)$ in res:
			var->boxplus(add);
			double *temp=workspace.data();
			for(IRVHolder::measurement_container::const_iterator meas= var->measurements.begin(); meas!= var->measurements.end(); meas++){
				double *start = temp; (void) start;
				temp=(*meas)->eval(temp, true);
				assert(temp - start == (*meas)->getDim() && "Dimension mismatch");
			}
			
			// store $f(\mu \mplus -1/d)$ directly in the matrix:
			var->boxplus(add, -1);
			temp = x;
			for(IRVHolder::measurement_container::const_iterator meas= var->measurements.begin(); meas!= var->measurements.end(); meas++){
				double *start = temp; (void) start;
				temp=(*meas)->eval(temp, true);
				assert(temp - start == (*meas)->getDim() && "Dimension mismatch");
			}
			
			// calculate difference and multiply by $0.5d$
			// also accumulate results for new inverse covariance
			double c = 0;
			for(double *xP=workspace.data() ;x<temp; x++){
				if(!std::isfinite(*x) || !std::isfinite(*xP)) {
					std::cerr << "\t" << xP - workspace.data() << " : " << temp-x << " " << *xP << " " << *x << " Variable: " << typeid(*var).name() << std::endl;
				}
				*x = half_delta_inv*(*xP++ - *x);
//				assert(std::isfinite(*x));
				c += std::pow(*x,2);
			}
			*chol++ = std::sqrt(c);
			add[k] = 0; // reset delta-vector
		}
		var->restore();
	}
	valuesUpToDate[ent_J] = true;
}


double SparseFunction::apply_delta(const VectorType &delta, double scale){
	CALLBACK_START(apply_delta);
	assert(valuesUpToDate[ent_residuum] && "Last residuum is not up to date");
	backup_valid = false;
	assert(delta.data() && delta.rows() == variables.getDim() && "Invalid delta vector!");
	// FIXME this function must block read-access to the variables
	const double *temp = delta.data();
	// Add delta-vector to the variables:
	for(RVList::iterator var = variables.begin(); var!= variables.end(); ++var){
		assert(var->optimize && "Fixed Variable in variables list");
		temp = var->boxplus(temp, scale);
	}
	valuesUpToDate[ent_newRes] = false;
	
	return getNewRSS();
}


const SparseFunction::SparseType& SparseFunction::get_JtJ(){
	if(valuesUpToDate[ent_JtJ]) return JtJ;  // if nothing changed, just return
	CALLBACK_START("get_JtJ");
	// TODO exploit block-structure of jacobian [X] mostly done
	// TODO make UpLo customizable
	enum { UpLo = Eigen::Upper };
		const SparseType &J = getJ(); // makes sure Jacobian is up to date
		if(!valuesUpToDate[ent_JtJ]) {
			if(structureUpToDate[ent_JtJ]){
				CALLBACK_START("re-calculateJtJ")
				internal::cs_JtJ<false, UpLo>(JtJ, getJt(), J);
			} else {
				CALLBACK_START("calculateJtJ")
				internal::cs_JtJ<true, UpLo>(JtJ, getJt(), J);
			}
			updated(ent_JtJ);
		}
	
	return JtJ;
}

void SparseFunction::valuesChanged() {
	valuesUpToDate.reset();
	if(solver){
		solver->valuesChanged();
	}
	if(algo) {
		algo->valuesChanged();
	}
	if(precond) {
		precond->valuesChanged();
	}
}

void SparseFunction::structureChanged() {
	{
		// a changed structure also implies that somehow the values changed
		valuesChanged();
		structureUpToDate.reset();
		if(solver)
			solver->structureChanged();
		if(algo) {
			algo->structureChanged();
		}
		if(precond) {
			precond->structureChanged();
		}
	}
}


void SparseFunction::updateJacobiPreconditioner(){
	if(valuesUpToDate[ent_jacobiPreconditioner]) return;
	CALLBACK_START(updateJacobiPreconditioner);
#ifdef SLOM_JACOBI_BLOCKS
	calculateBlockJacobian();
//	const double * jPointer = getJ().valuePtr();
	for(RVList::iterator var = variables.begin(); var!= variables.end(); ++var){
		assert(var->optimize);
		var->calcBlockJacobi();
	}
#else
	assert(false && "You must enable SLOM_JACOBI_BLOCKS for this to work");
#endif
	updated(ent_jacobiPreconditioner);
}


void SparseFunction::applyJacobiPreconditioner(VectorType &vec, bool transpose) {
#ifdef SLOM_JACOBI_BLOCKS
	updateJacobiPreconditioner();
	assert(vec.rows() == variables.getDim() && "Vector dimension does not equal variable dimension");
	double * vecPtr = vec.data();
	for(RVList::const_iterator var = variables.begin(); var!= variables.end(); ++var){
		assert(var->optimize);
		vecPtr = var->applyBlockJacobi(vecPtr, transpose);
	}
#else
	(void)vec; (void) transpose;
	assert(false && "You must enable SLOM_JACOBI_BLOCKS for this to work");
#endif
	
}


void SparseFunction::applyPreconditioner(VectorType & vec, bool transpose) {
	if(precond) {
		precond->compute(); // FIXME reduce number of calls if virtual overhead is significant
		precond->apply(vec, transpose);
	}
}

}  // namespace SLOM
