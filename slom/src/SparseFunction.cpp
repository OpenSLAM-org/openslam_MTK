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

#include "cs_helper.hpp"

namespace SLOM {


double SparseFunction::evaluate(double * result) const{
	MeasurementList::const_iterator it = measurements.begin();
	double sum = 0;
	double *ptr = result;
	for(; it != measurements.end(); it++){
		assert(result + it->idx == ptr);
		it->eval(ptr, false);
		double *end = ptr+it->getDim();
		for( ; ptr < end; ptr++){
			sum += std::pow(*ptr, 2);
		}
	}
	return sum;
}




void SparseFunction::createSparse(){
	freeWorkspace();
	int M = measurements.getDim();
	int N = variables.getDim();
	int size = nnz;
	int skip = 0;
//	if(_addDiagonal)
//	{
//		size += N;
//		M += N;
//		skip = 1;
//	}
	//allocate a compressed column matrix:
	jacobian = cs_spalloc(M, N,   // size
			size,    // number of non-zeros
			true,   // allocate space for values
			false); // compressed-column
	
	//create the structure of the matrix:
	int *cIdx=jacobian->p; // column pointer
	int *rIdx=jacobian->i; // row indices
	*cIdx = 0; //first column starts at 0.
//	int n=measurements.getDim();
	for(RVList::const_iterator var = variables.begin(); var!= variables.end(); ++var){
		int vDOF = var->getDOF();
		assert(vDOF>0);
		if(!var->has_measurements()){
			std::cerr << "No measurement for Variable " 
					<< var->idx << std::endl;
			assert(_addDiagonal);
		}
		for(IRVHolder::measurement_container::const_iterator meas= var->measurements.begin(); meas!= var->measurements.end(); ++meas){
			int row = (*meas)->idx;
			for(int mDim = (*meas)->getDim(); mDim>0; mDim--){
				*rIdx++ = row++;
			}
		}
//		if(_hasDiagonal){
//			*rIdx++ = n++;
//		}
		*++cIdx = rIdx - jacobian->i; //start of next column
		while(--vDOF > 0){ // copy the current column for each DOF of the variable
			rIdx = std::copy(jacobian->i + cIdx[-1], // start of previous column
					rIdx-skip, rIdx);
//			if(_hasDiagonal){
//				*rIdx++ = n++;
//			}
			*++cIdx = rIdx - jacobian->i;
		}
	}
	assert(*cIdx == size);
}


void SparseFunction::initCovariance(){
	int n=variables.getDim();
	cholCovariance.resize(n);
	cholCovariance.fill(0);
}

void SparseFunction::calculateJacobian(){
	assert(jacobian);   // matrix is allocated
	int m = jacobian->m, n = jacobian->n;
	double *x = jacobian->x;
	
	assert(m == measurements.getDim());
	assert(n == variables.getDim());
	
	double *chol = cholCovariance.data();
	// FIXME a little bit conservative workspaces; will be redundant with symbolic derivations/block-Jacobians
	Eigen::VectorXd delta(n);
	Eigen::VectorXd workspace(m);
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
				temp=(*meas)->eval(temp, true);
			}
			
			// store $f(\mu \mplus -1/d)$ directly in the matrix:
			var->boxplus(add, -1);
			temp = x;
			for(IRVHolder::measurement_container::const_iterator meas= var->measurements.begin(); meas!= var->measurements.end(); meas++){
				temp=(*meas)->eval(temp, true);
			}
			
			// calculate difference and multiply by $0.5d$
			// also accumulate results for new inverse covariance
			double c = 0;
			for(double *xP=workspace.data() ;x<temp; x++){
				*x = half_delta_inv*(*xP++ - *x);
				assert(std::isfinite(*x));
				c += std::pow(*x,2);
			}
			*chol++ = std::sqrt(c);
			add[k] = 0; // reset delta-vector
		}
		var->restore();
	}
}

void SparseFunction::updateDiagonal(double lamda, double* cholCovariance) {
	if(!_addDiagonal) return;
	// it is kind of dirty to change values in a const object ...
	const cs* JtJ = get_JtJ();
	// for LMA set the last entry of each column to lamda or lamda*cholCovariance;
	int n=JtJ->n;
	int *p = JtJ->p;
	double *x = JtJ->x;
	
	for(int k=0; k<n; k++){
		// Diagonal entry is always first entry
		x[p[k]] += (!cholCovariance ?
				lamda : lamda*cholCovariance[k]);
	}
	
}

void SparseFunction::apply_delta(const Eigen::VectorXd &delta, double scale){
	const double *temp = delta.data();
	// Add delta-vector to the variables:
	for(RVList::iterator var = variables.begin(); var!= variables.end(); ++var){
		if(var->optimize){
			temp = var->boxplus(temp, scale);
		} else {
			temp += var->getDOF();
		}
	}
}


CS_INT filterUpper (CS_INT i, CS_INT j, CS_ENTRY, void *){
	return i <= j;
}

const cs* SparseFunction::get_JtJ(){
	if(!JtJ){
		getJt();
//		JtJ = cs_di_multiply(Jt, jacobian);
		JtJ = internal::cs_JtJ(Jt, jacobian);
//		std::cerr << JtJ->nzmax << " ";
//		cs_di_fkeep(JtJ, filterUpper, 0);
//		std::cerr << JtJ->nzmax << "\n";
	}
	
	return JtJ;
}



}  // namespace SLOM
