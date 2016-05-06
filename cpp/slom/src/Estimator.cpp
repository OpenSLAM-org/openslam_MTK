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
 * @file slom/src/Estimator.cpp
 * @brief implements the Estimator class
 */

#include "../Estimator.hpp"

#include <cmath>
#include <cassert>

#include <iostream>
#include <fstream>



#include "../CallBack.hpp"


namespace SLOM {



Estimator::~Estimator(){
	// every temporary destructs itself
}






double Estimator::optimizeStep(bool linear_refine){
	CALLBACK_START(optimizeStep);
	(void) linear_refine;
	return algo->optimizeStep();
#if false
	// TODO better parameter control for LMA
	// TODO factor out into Algorithm class
	
	
	// calculate delta = (J^T * J + lambda*I) \ (J^T * res):
	calculate_delta();
	
	double newRSS;
	if(!linear_refine)
	{
		newRSS = apply_delta(-1);
	}
	else
	{
		// try linear refinement (highly experimental and not optimized)
		cb.output(5) << "RSS(scale): ";
		double bestRSS = lastRSS;
		double bestScale = 0;
		for(double scale = -2; scale<0; scale+=1.0/16)
		{
			newRSS = apply_delta(scale);
			cb.output(5) << newRSS << " ";
			if(newRSS < bestRSS)
			{
				bestRSS = newRSS;
				bestScale = scale;
			}
		}
		cb.output(5) << lastRSS << "\n";
		
		newRSS = apply_delta(bestScale);
	}
	
	cb.output(2) << "lastRSS: " << lastRSS; 
	double gain = store_or_restore(newRSS, cb);
	
	cb.output(2) << std::endl;

	return gain;
#endif
}

std::ostream& Estimator::dumpJacobian(std::ostream& out)
{
	const SparseFunction::SparseType& jacobian = func.getJ();
	
	for(int k=0; k<jacobian.outerSize(); ++k) {
		for(SparseFunction::SparseType::InnerIterator it(jacobian, k); it; ++it) {
			out << it.row() << " " << it.col() << " " << it.value() << "\n";
		}
	}
	
	return out;
}


}  // namespace SLOM
