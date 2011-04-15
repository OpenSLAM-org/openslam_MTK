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

// MTK's pose and orientation definition:
#include <mtk/types/pose.hpp>
#include <mtk/types/SOn.hpp>

// Main SLoM class:
#include <slom/Estimator.hpp>
// CallBack class for progress/debug output
#include <slom/CallBack.hpp>
// AutoConstruct macros to define measurements
#include <slom/BuildMeasurement.hpp>
// CholeskyCovariance to apply Covariance/Information-matrix
#include <slom/CholeskyCovariance.hpp>
// Simple timing class
#include <slom/TicToc.hpp>

#include <vector>

// We can't use types having a comma inside AutoConstruct macros :(
typedef MTK::vect<3, double> vec3;
typedef MTK::SO3<double> SO3;



// Define a Compund manifold
MTK_BUILD_MANIFOLD(PoseCompund, /* Name of the manifold */
	((vec3, pos))      /* each sub-variable in form ((type, name)) */
	((SO3, orient))    /* no comma between variables! */
);	                   /* semicolon is optional */

// Alternatively use predefined trafo-Type (already implements transformations):
typedef MTK::trafo<MTK::SO2<double> > Pose;



// Define a measurement
SLOM_BUILD_MEASUREMENT(Odo, Pose::DOF, /* Name and dimension of measurement */
		((Pose, t0))    /* ((type, name)) for each variable */ 
		((Pose, t1))    /* type must be an MTK manifold */
		,               /* Comma between variables and data */
		((Pose, odo))   /* ((type, name)) for each data element */
		/* Type can be arbitrary, i.e. also pointers or references */
		((SLOM::CholeskyCovariance<Pose::DOF>, cov)) 
);		                /* semicolon is optional */

// Implementation of measurement:
SLOM_IMPLEMENT_MEASUREMENT(Odo, ret) // Odo is the name of the measurement, 
{                                    // ret is the return variable
	// variables need to be accessed with * or -> operator
	// data is stored by value (can be accessed directly)
	
	// calculate relative transformation, diff = t0^-1 * t1
	Pose diff = t0->world2Local(*t1);
	// calculate ret = diff [-] odo
	diff.boxminus(ret, odo);
	// multiply ret with Cholesky factor of information matrix
	// this normalizes the measurement to have unit-covariance
	cov.apply(ret);
}


// Method to read odometry (needs to be implemented somewhere)
bool getOdometry(Pose::vect_type& trans, Pose::rotation_type& rot, 
                 SLOM::CholeskyCovariance<Pose::DOF>& cov){
	return false;
}


int main()
{
	// create Estimator, usually GaussNewton is a good choice;
	// use Parameters (SLOM::Levenberg, initalLambda) for rank-deficient problems
	SLOM::Estimator est(SLOM::Estimator::GaussNewton);
	typedef SLOM::VarID<Pose> PoseID;
	
	// Variable-storage can now be done by Estimator:
	PoseID lastPose = est.insertRV(Pose(), false); // don't optimize first pose
	
	// just a list of references (i.e. pointers) to observe result afterwards:
	// could also be used to remove variables or change correspondences (future work ...)
	std::vector<PoseID> poses; 
	poses.push_back(lastPose);
	
	Pose::rotation_type rot;
	Pose::vect_type trans;
	SLOM::CholeskyCovariance<Pose::DOF> cov;
	while( getOdometry(trans, rot, cov) ) // get Odometry from somewhere ...
	{
		Pose odo(trans, rot); // Construct odometry from rotation and translation
		// Calculate next pose from last pose and odometry and insert into Estimator
		PoseID currentPose = est.insertRV(lastPose->local2World(odo));
		// MeasID can be used to remove measurements (otherwise it can be ignored)
		//Odo::id x;
		SLOM::MeasID<Odo> odoMeas = 
			est.insertMeasurement(
				// Parameters for Odometry constructor in order of definition:
				Odo(lastPose, currentPose, odo, cov)
			);
		// store ID of current pose (only required for data output)
		poses.push_back(currentPose);
		lastPose = currentPose;
	}
	
	// Run optimize-step till convergence:
	int kMax = 50;
	for(int k=1; k<=kMax; k++){
		std::cout << "Step " << k << ": ";
		double gain = est.optimizeStep();
		
		if( 0 <= gain && gain < 1e-9) break;
	}
	
	// output result
	for(std::vector<PoseID >::const_iterator it = poses.begin(); it!= poses.end(); ++it)
	{
		const Pose& p  = **it; // derefence from iterator to PoseID& and to Pose&
		// output pose (preferably into a file ...) using automatically generated streaming operator
		std::cout << p << std::endl;
	}
	
	return 0;
}
