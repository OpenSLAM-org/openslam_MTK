/*
 *  Copyright (c) 2010--2011, Universitaet Bremen
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
 * @file relation_parser.hpp
 */

#ifndef RELATION_PARSER_HPP_
#define RELATION_PARSER_HPP_

#include <boost/array.hpp>
#include <map>

#include <slom/Estimator.hpp>
#include <slom/CholeskyCovariance.hpp>
#include <slom/TicToc.hpp>
#include "../tools.h"

#include "../../unit_test/random_vector.hh"

#include <slom/BuildMeasurement.hpp>
//#include <mtk/types/pose.hpp>
//#include <mtk/types/SOn.h>



// AutoConstruct does not work properly with templated types
template<class PoseType, class PoseDiffType=PoseType>
struct OdoT {
	SLOM_VAR_REF_TYPE(PoseType, 0) t0;
	SLOM_VAR_REF_TYPE(PoseType, PoseType::DOF) t1;
	enum {DIM = PoseDiffType::DOF, DEPEND=2*PoseType::DOF}; 
	
	PoseDiffType odo;
	OdoT( SLOM::VarID<PoseType> t0, SLOM::VarID<PoseType> t1, const PoseDiffType &odo ) : 
		t0(t0), t1(t1), odo(odo) {}
	template<class Func> 
	void traverse_variables(Func __function) { 
		__function(t0); __function(t1);
	}
	void eval(MTK::vectview<double, DIM> __ret, bool numerical_jacobian) const
	{ eval(__ret);}
	void eval(MTK::vectview<double, DIM> __ret) const {
		PoseDiffType diff = t0->world2Local(*t1);
		diff.boxminus(__ret, odo);
	}
};

template<class PoseType, int cov_inverse=1, class PoseDiffType=PoseType>
struct OdoTCov : public OdoT<PoseType, PoseDiffType> {
	typedef OdoT<PoseType, PoseDiffType> base;
	SLOM::CholeskyCovariance<PoseDiffType::DOF> cov;
	using base::DIM;
	OdoTCov( SLOM::VarID<PoseType> t0, SLOM::VarID<PoseType> t1, 
			const PoseDiffType &odo, const SLOM::CholeskyCovariance<DIM> &cov) : 
				base(t0, t1, odo), cov(cov) {}
	void eval(MTK::vectview<double, DIM> __ret, bool numerical_jacobian) const
	{ eval(__ret);}
	void eval(MTK::vectview<double, DIM> __ret) const {
		base::eval(__ret);
		if(cov_inverse > 0)
			cov.apply(__ret);
		else
			cov.invApply(__ret);
	}
};



// TODO also support Landmark measurements?
template<class PoseType, class PoseDiffType=PoseType>
class RelationParser 
{
	SLOM::Estimator &e;
	
	typedef PoseType Pose;
	typedef SLOM::VarID<Pose> PoseID;
	typedef std::map<int, PoseID> Poses;
	
	Poses poses;
	
	random_vector<double> rng;
public:
	typedef void (*ParsePose)(std::istream&, PoseDiffType&);
	typedef void (*ParseCov)(std::istream&, SLOM::CholeskyCovariance<PoseDiffType::DOF>&);
	typedef void (*OutputPose)(std::ostream &, const PoseDiffType &);
	
private:
	template<class Container, class ValueType>
	bool contains(const Container &c, const ValueType &value)
	{
		return std::find(c.begin(), c.end(), value) != c.end();
	}
	
	// no Covariance parsing
	void insertPoseRel(PoseID start, PoseID end, const PoseDiffType& delta, std::istream& inp, int, int){
		e.insertMeasurement(OdoT<Pose, PoseDiffType>(start, end, delta));
	}
	// with Covariance parsing
	template<class CovParser>
	void insertPoseRel(PoseID start, PoseID end, const PoseDiffType& delta, std::istream& inp, CovParser covParser, int inverseCov){
		if(!covParser){
			e.insertMeasurement(OdoT<Pose, PoseDiffType>(start, end, delta));
			return;
		}
		SLOM::CholeskyCovariance<PoseDiffType::DOF> cov;
		covParser(inp, cov);
		if(inverseCov > 0)
			e.insertMeasurement(OdoTCov<Pose, +1, PoseDiffType>(start, end, delta, cov));
		else
			e.insertMeasurement(OdoTCov<Pose, -1, PoseDiffType>(start, end, delta, cov));
	}
	
	
public:
	RelationParser(SLOM::Estimator &e, bool fix1stPose = true, unsigned int seed=5489) : e(e), rng(seed) {
		fixFirstPose(fix1stPose);
	}
	
	void fixFirstPose(bool fix = true){
		if(!poses[0]){
			if(fix) {
				poses[0] = e.insertRV(Pose(), !fix);
			}
		} else {
			e.optimizeRV(poses[0], !fix);
		}
	}
	
	void insertPose(int id, const Pose& p) {
		poses[id] = e.insertRV(p);
	}
	
	template<size_t vs, size_t os>
	void parse(std::istream &input, 
			const boost::array<std::string, vs> &vertexStrings, const boost::array<std::string, os> &odoStrings,
			ParsePose pp, ParseCov cp, int inverseCov = 1, double addNoise=-1.0)
	{
		std::string line;
		
		while(getline(input, line)){
			std::istringstream inp(line);
			std::string tag;
			inp >> tag;
			if(contains(vertexStrings, tag)) {
				int id;
				PoseDiffType pose;
				inp >> id;
				pp(inp, pose);
				Pose pose_ = pose;
				if(poses.find(id)==poses.end()){
					insertPose(id, pose_);
				} else {
					e.reinitRV(poses[id], pose_);
				}
			} else if (contains(odoStrings, tag)) {
				int frameA, frameB;
				PoseDiffType delta;
				inp >> frameA >> frameB;
				pp(inp, delta);
				if(addNoise > 0.0){
					delta.boxplus(rng.get_normal<PoseDiffType::DOF>(), addNoise);
				}
				
				// if either pose is not available yet, compute its position from the other one.
				if(poses.find(frameB)==poses.end()){
					// if neither pose is available, initialize first with standard pose
					if(poses.find(frameA)==poses.end()){
						insertPose(frameA, Pose());
//						poses[frameA] = e.insertRV(Pose());
						std::cerr << "Warning: unconnected " << frameA << "->" << frameB << "   " << std::endl;
					}
					insertPose(frameB, poses[frameA]->local2World(delta));
//					poses[frameB]= e.insertRV(poses[frameA]->local2World(delta));
				}
				if(poses.find(frameA)==poses.end()){
					insertPose(frameA, *poses[frameB] / (delta));
//					poses[frameA]= e.insertRV(*poses[frameB] / (delta));
				}
				
				insertPoseRel(poses[frameA], poses[frameB], delta, inp, cp, inverseCov);
			}
		}
	}
	
	void outputPoses(const std::string& filename) const {
		std::ofstream out(filename.c_str());
		outputPoses(out);
	}
	void outputPoses(std::ostream& out) const {
		for (typename Poses::const_iterator it = poses.begin(); it != poses.end(); ++it) {
			out << *it->second << std::endl;
			out << std::endl;
		}
	}
	
	void iterate(int max_steps, const std::string& outname, OutputPose outputPose, std::ostream& statstream = std::cout) {
		TicToc timer;
		
		statstream << "# Step\ttime\tRSS\n";
		double gain = 1e9;
		for(int k=0; k<std::abs(max_steps); ++k){
			if(!outname.empty())
				outputPoses(make_filename(outname, k, ".pos"));
			double rss = e.getRSS(); // evaluate RSS before getting the time
			statstream << k << "\t" << timer() << "\t" << rss << std::endl;
			
			if( max_steps > 0 && 0 <= gain && gain < 1e-9){
				max_steps = k+1; // output stats, then quit
			} else {
				try{
					gain = e.optimizeStep();
				} catch(const char* err) {
					for(++k ; k < -max_steps; ++k){ // for negative max_steps fill statstream with inf entries
						statstream << k << "\t" << timer() << "\tinf" << std::endl;
					}
					throw err;
				}
			}
		}
		
	}
};


template<SLOM::CholeskyMode::CM cm, int DOF>
void parseCov(std::istream &inp, SLOM::CholeskyCovariance<DOF> &cov) {
	Eigen::Matrix<double, DOF, DOF, Eigen::RowMajor> m;
	for(int i=0; i<DOF; ++i){
		for(int j=i; j<DOF; ++j){
			inp >> m(i,j);
		}
	}
	cov = SLOM::CholeskyCovariance<DOF>(m.data(), cm);
}



#endif /* RELATION_PARSER_HPP_ */
