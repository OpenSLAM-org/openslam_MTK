/*
 *  Copyright (c) 2011, Universitaet Bremen
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
 * @file relation_slameval.cpp
 * @brief Parser hand-crafted for SLAM-Eval workshop at RSS 2011
 * 
 */



#include <iostream>
#include <sstream>

#include <boost/algorithm/string.hpp>
#include <boost/scoped_ptr.hpp>

#include <vector>

// Just using OdoT from here:
#include "relation_parser.hpp"

#include <mtk/types/SOn.hpp>
#include <mtk/types/pose.hpp>

#include "euler_stable.hpp"
#include "../tools.h"


typedef MTK::trafo<MTK::SO2<double> > Pose2D;
typedef MTK::trafo<MTK::SO3<double> > Pose3D;



class IParser {
protected:
	SLOM::Estimator &e;
public:	
	IParser(SLOM::Estimator& e) : e(e) {}
	
	virtual ~IParser() {}
	
	virtual void parseEdge(std::istream& inp) = 0;
	virtual void outputPoses(std::ostream& out, 
			const std::vector<int>& indexes = std::vector<int>()) const = 0;
};


template<class PoseType>
class Parser : public IParser {
	
	
	char buffer_mem[1024];
	char *buffer_start;
	
	typedef PoseType Pose;
	typedef SLOM::VarID<Pose> PoseID;
	typedef std::map<int, PoseID> Poses;
	
	Poses poses;
	
	// Can be specialized to allow Euler input
	static std::istream& parsePose(std::istream &inp, Pose& pose) {
		return inp >> pose;
	}
	
	
	static char* outputPose(char*, const Pose& pose); 
	
	char* outputPose(const int& id, const Pose& pose) const {
		char* buf = itoa(id, buffer_start);
		*buf++ = ' ';
		buf = outputPose(buf, pose);
		*buf++ = '\n';
		return buf;
		
	}
	
	void insertPoseRel(PoseID start, PoseID end, const Pose& delta, std::istream& inp){
		SLOM::CholeskyCovariance<Pose::DOF> cov;
		parseCov<SLOM::CholeskyMode::CHOLESKY_FULL, Pose::DOF>(inp, cov);
		e.insertMeasurement(OdoTCov<Pose, +1, Pose>(start, end, delta, cov));
	}
	void insertPose(int id, const Pose& p) {
		poses[id] = e.insertRV(p);
	}

public:
	Parser(SLOM::Estimator &e, const::std::string &vertex_tag) : IParser(e) {
		// fix first pose
		poses[0] = e.insertRV(Pose(), false);
		buffer_start = std::copy(vertex_tag.begin(), vertex_tag.end(), buffer_mem);
		*buffer_start++ = ' ';
	}
	
	void parseEdge(std::istream &inp){
		int edgeID, frameA, frameB;
		Pose delta;
		inp >> edgeID >> frameA >> frameB;
		parsePose(inp, delta);
		
		// if either pose is not available yet, compute its position from the other one.
		if(poses.find(frameB)==poses.end()){
			// if neither pose is available, initialize first with standard pose
			if(poses.find(frameA)==poses.end()){
				insertPose(frameA, Pose());
				std::cerr << "Warning: unconnected " << frameA << "->" << frameB << "   " << std::endl;
			}
			insertPose(frameB, poses[frameA]->local2World(delta));
		}
		if(poses.find(frameA)==poses.end()){
			insertPose(frameA, *poses[frameB] / (delta));
		}
		
		insertPoseRel(poses[frameA], poses[frameB], delta, inp);
	}
	
	void outputPoses(std::ostream& out, const std::vector<int>& indexes) const {
		if(indexes.empty()){
			
			for(typename Poses::const_iterator it = poses.begin(); it != poses.end(); ++it){
				char *end = outputPose(it->first, *it->second);
				out.write(buffer_mem, end - buffer_mem);
			}
		} else {
			for(std::vector<int>::const_iterator it = indexes.begin(); it != indexes.end(); ++it){
				char *end = outputPose(*it, *poses.at(*it));
				out.write(buffer_mem, end - buffer_mem);
//				out << tag << " " << *it << " ";
//				Parser::outputPose(out, *poses.at(*it)) << "\n";
			}
		}
	}
};


void old_euler_to_quat(Eigen::Quaterniond& res, double roll, double pitch, double yaw) {
	double sr = std::sin(roll*0.5);
	double cr = std::cos(roll*0.5);
	double sp = std::sin(pitch*0.5);
	double cp = std::cos(pitch*0.5);
	double sy = std::sin(yaw*0.5);
	double cy = std::cos(yaw*0.5);
	res.w() = cr*cp*cy + sr*sp*sy;
	res.x() = sr*cp*cy - cr*sp*sy;
	res.y() = cr*sp*cy + sr*cp*sy;
	res.z() = cr*cp*sy - sr*sp*cy;
	
}

void old_quat_to_euler(double& roll, double& pitch, double& yaw, const Eigen::Quaterniond &q) {
  const double q0 = q.w();
  const double q1 = q.x();
  const double q2 = q.y();
  const double q3 = q.z();
  roll = atan2(2*(q0*q1+q2*q3), 1-2*(q1*q1+q2*q2));
  pitch = asin(2*(q0*q2-q3*q1));
  yaw = atan2(2*(q0*q3+q1*q2), 1-2*(q2*q2+q3*q3));
}

void wRo_to_euler(double& roll, double& pitch, double& yaw, const Eigen::Matrix3d& wRo) {
	yaw = std::atan2(wRo(1,0), wRo(0,0));
  double c = std::cos(yaw);
  double s = std::sin(yaw);
  pitch = std::atan2(-wRo(2,0), wRo(0,0)*c + wRo(1,0)*s);
  roll  = std::atan2(wRo(0,2)*s - wRo(1,2)*c, -wRo(0,1)*s + wRo(1,1)*c);
}

template<>
std::istream& Parser<Pose3D>::parsePose(std::istream & inp, Pose3D& pose){
	inp >> pose.pos;
	double roll, pitch, yaw;
	inp >> roll >> pitch >> yaw;
	euler_to_quat(pose.orient, roll, pitch, yaw);
	return inp;
}

// Output using std::strings turns out to be very slow
template<>
char* Parser<Pose3D>::outputPose(char* out, const Pose3D& pose){
	for(int i=0; i<3; ++i){
		out = ftoa(pose.pos[i], out);
		*out++ = ' ';
	}
	double roll, pitch, yaw;
	quat_to_euler(roll, pitch, yaw, pose.orient);
	out = ftoa(roll , out); *out++=' ';
	out = ftoa(pitch, out); *out++=' ';
	out = ftoa(yaw  , out);
	
	return out;
}
template<>
char* Parser<Pose2D>::outputPose(char* out, const Pose2D& pose){
	for(int i=0; i<2; ++i){
		out = ftoa(pose.pos[i], out);
		*out++ = ' ';
	}
	out = ftoa(pose.orient, out);
	
	return out;
}



static std::istream & getline_stripped(std::istream& input, std::string &line){
	std::getline(input, line);
	
	// remove trailing ';'
	std::string::iterator it = line.end() - 1;
	if(*it == ';'){
		line.erase(it);
	}
	return input;
}



int main() {
	SLOM::Estimator est;
	
	// virtual optimizer:
	boost::scoped_ptr<IParser> opt;
	
	
	std::istream &input = std::cin;
	std::ostream &output = std::cout;
	
	output.precision(6);
	
	std::string line, edgetag, vrtxtag;
	
	// First loop to find out if 2d or 3d problem:
	do{
		std::vector<std::string> tags;
		getline_stripped(input, line);
		boost::split(tags, line, boost::is_any_of(" _;"));
		
		// First tag should be ADD:
		if(tags[0] != "ADD"){
			std::cerr << "ERROR: Expected ADD ... as first command. Ignoring:\n"
					<< line << std::endl;
			continue;
		}
		
		if(tags[1] != "VERTEX" && tags[1] != "EDGE") {
			std::cerr << "ERROR: Expected ADD VERTEX_... or ADD EDGE_... Ignoring:\n"
					<< line << std::endl;
			continue;
		}
		
		if(tags[2] == "XYT"){
			// 2D problem:
			vrtxtag = "VERTEX_XYT";
			opt.reset(new Parser<Pose2D>(est, vrtxtag));
			edgetag =   "EDGE_XYT";
		} else if(tags[2] == "XYZRPY") {
			// 3D problem:
			vrtxtag = "VERTEX_XYZRPY";
			opt.reset(new Parser<Pose3D>(est, vrtxtag));
			edgetag =   "EDGE_XYZRPY";
		} else {
			std::cerr << "Unknown Pose format. Currently only XYT and XYZRPY are supported\n";
			exit(1);
		}
	} while(opt == 0);
	
	
	double steps = 5;
	//infinite loop:
	for(;;){
		// process last line:
		std::istringstream inp(line);
		std::string tag;
		inp >> tag;
		if(tag=="ADD"){
			inp >> tag;
			steps += 0.02;
			if(tag==edgetag){
				opt->parseEdge(inp);
				steps-=.019;
			}
		} else if(tag=="QUERY_STATE") {
			// Do optimization:
			for( ; steps > 0; steps-=1.0){
				double gain = est.optimizeStep();
				//std::cerr << "\noptimizeStep " << steps << " gain: " << gain << " RSS: " << est.getRSS() << std::endl; 
				if((0 < gain && gain < 1e-5) || est.getRSS() < 1e-9) break;
			}
			std::vector<int> indexes;
			while(inp) {
				int idx;
				inp >> idx;
				if(!inp) break;
				indexes.push_back(idx);
			}
			output << "BEGIN\n";
			opt->outputPoses(output, indexes);
			output << "END" << std::endl;
		}
		
		getline_stripped(input, line);
	}
	
	
	return 0;
}
