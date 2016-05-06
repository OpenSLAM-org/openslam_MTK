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
 * @file relation.cpp
 * @brief Brief description
 * 
 */


#include <mtk/types/pose.hpp>
#include <mtk/types/SOn.hpp>
#include <slom/CallBack.hpp>
#include <slom/Optimizer.hpp>

#include "SO3singular.hpp"
#include "relation_parser.hpp"


#include <fstream>

#include <boost/program_options.hpp>
namespace po=boost::program_options;

#include <boost/scoped_ptr.hpp>



typedef MTK::trafo<MTK::SO3<double> > Pose3Dexp;

// Just experimental representations, proven to be inferior:
typedef MTK::trafo<MTK::SO3singular<MTK::scaledAxis<double> > > Pose3Dscaled;
typedef MTK::trafo<MTK::SO3singular<MTK::euler<double> > > Pose3Deuler;
typedef MTK::trafo<MTK::SO3singular<MTK::vec4SO3<double> > > Pose3Dquat4;

// When using 4d-quaternions as representation, an additional measurement might be required:
SLOM_BUILD_MEASUREMENT(UnitQuat, 1, ((Pose3Dquat4, pose)), )
SLOM_IMPLEMENT_MEASUREMENT(UnitQuat, ret){
	ret(0) = pose->orient.squaredNorm() - 1;
}

// template specialization for 4d-quaternions:
template<>
void RelationParser<Pose3Dquat4, Pose3Dexp>::insertPose(int id, const Pose3Dquat4& p) {
	poses[id] = e.insertRV(p);
	if(optim.getLambda() <=0 ) // add unit-length measurement for Gauss-Newton mode
		e.insertMeasurement(UnitQuat(poses[id]));
}


typedef MTK::trafo<MTK::SO2<double> > Pose2D;



/**
 * Convert Euler angles in XYZ-roll-pitch-yaw convention to a Quaternion.
 * Conversion is not done in the most efficient way, but the most obvious.
 */
template<class scalar>
void euler2quaternion(Eigen::Quaternion<scalar> &res,
		const scalar& phi, const scalar& theta, const scalar& psi) {
	Eigen::Quaternion<scalar> 
			qx(std::cos(phi/2), std::sin(phi/2), 0, 0),
			qy(std::cos(theta/2), 0, std::sin(theta/2), 0),
			qz(std::cos(psi/2), 0, 0, std::sin(psi/2));
	res = qz * qy * qx;
}

void parsePoseEuler(std::istream & inp, Pose3Dexp & pose){
	inp >> pose.pos;
	double phi, theta, psi;
	inp >> phi >> theta >> psi;
	
	euler2quaternion(pose.orient, phi, theta, psi);
}

void parsePoseQuat(std::istream &inp, Pose3Dexp & pose)
{
	inp >> pose;
}

enum Indexes {RGX, RGY, RGPHI};

void parsePose2D(std::istream& inp, Pose2D &p){
	inp >> p;
}

// legacy code to support strange toro convention to store the covariance:
template<SLOM::CholeskyMode::CM cm>
void parseCovToro(std::istream &inp, SLOM::CholeskyCovariance<3>& cov){
	double xCov[3][3];
   inp >> xCov[RGX][RGX] >> xCov[RGX][RGY] >> xCov[RGY][RGY]
       >> xCov[RGPHI][RGPHI] >> xCov[RGX][RGPHI] >> xCov[RGY][RGPHI];
   cov = SLOM::CholeskyCovariance<3>(xCov[0], cm);
}


void outputPose3D(std::ostream &out, const Pose3Dexp &p){
	out << p;
}



const boost::array<std::string, 0> empty = {};
const boost::array<std::string, 4> vertexes = {{"VERTEX", "VERTEX2", "VERTEX3", "VERTEX_SE3:QUAT"}};
const boost::array<std::string, 5> edges = {{"EDGE", "EDGE2", "EDGE3", "EDGE_SE3:QUAT", "ODOMETRY"}};

template<class MeasurementType>
struct IRPHolder{
	typedef RelationParser<MeasurementType, MeasurementType> RP;
	
	typedef  typename RP::ParsePose ParsePose;
	typedef  typename RP::ParseCov  ParseCov;
	
	virtual ~IRPHolder() {};
	virtual void parse(std::istream &logfile, bool init,
			ParsePose pp, ParseCov cp, int inverseCov = 1, double addNoise=-1.0) = 0;
	virtual void iterate(int max_steps, const std::string& outname, std::ostream& statstream = std::cout) = 0;
};

template<class PoseType, class MeasurementType=PoseType>
struct RPHolder : public IRPHolder<MeasurementType> {
	typedef RelationParser<PoseType, MeasurementType> RP;
	
	typedef typename RP::ParsePose ParsePose;
	typedef typename RP::ParseCov  ParseCov;
	
	RP rp;
	
	RPHolder(SLOM::Optimizer &optim, int seed) : rp(optim, true, seed) {}
	
	void parse(std::istream &logfile, bool init,
			ParsePose pPose, ParseCov pCov, int covariance, double noise=-1.0) {
		if(init)
			rp.parse(logfile, vertexes, edges, pPose, pCov, covariance, noise);
		else
			rp.parse(logfile, empty, edges, pPose, pCov, covariance, noise);
		
	}
	
	static
	void outputPose(std::ostream &out, const MeasurementType &p){
		out << p;
	}

	
	void iterate(int max_steps, const std::string& outname, std::ostream& statstream = std::cout) {
		try{
			rp.iterate(max_steps, outname, outputPose, statstream);
		} catch (const char* err) {
			std::cerr << err << std::endl;
		}
	}

};


int main(int argc, char* argv[])
{
	
	SLOM::OptimizerOptions optiOpts;
	
	int dim;
	char input_angles;
	bool toro;
	
//	int max_steps;
	int covariance;
	bool init;
	std::string output;
	double noise;
	unsigned int seed;
	
	char orientation_type;
	
	std::string measurements, dump_jacobian, statsfile;
//	int verbosity;
//	double lambda;
	po::options_description desc("Options");
	desc.add(optiOpts);
	desc.add_options()
		("help,h", "print this help message")
//		("lambda,l", po::value(&lambda)->default_value(0), "initial lambda for Levenberg-Marquard optimization, if <= 0 Gauss-Newton is performed")
//		("steps,s", po::value(&max_steps)->default_value(10), "maximal number of optimization steps, for negative values run exactly abs(steps)")
		("dim,d", po::value(&dim)->default_value(3), "2D or 3D problem")
		("cov,C", po::value(&covariance)->default_value(0)->implicit_value(-1), 
				"Covariance information is given as Covariance (>0) or Information (<0) and is already (+/-2) or needs to be (+/-1) decomposed")
		("toro", po::value(&toro)->default_value(false)->implicit_value(true), "Use strange toro convention for 2D covariance")
		("input-angles,a", po::value(&input_angles)->default_value('Q'), "input angles are represented as 'Q'uaternion, or 'E'uler angles")
		("init,I",  po::value(&init)->default_value(false)->implicit_value(true), "Use init values from logfile if present")
		("measurements,m",po::value(&measurements))
		("output,o",      po::value(&output)->implicit_value("output"), "output filename, will be extended by 3-digits and type extensions")
		("stats,O", po::value(&statsfile)->implicit_value("stats.dat"), "output statistics to this file")
		("dump-jacobian,J", po::value(&dump_jacobian)->implicit_value("jacobian.dat"), "file name to dump last jacobian")
//		("verbosity,v", po::value(&verbosity)->default_value(2), "Verbosity of optimizer, 0=off")
		("noise,N", po::value(&noise)->default_value(-1.0), 
				"add artificial noise to measurements (non-positive values are ignored, just works with -C 0, i.e. only N parameter is used")
		("seed", po::value(&seed)->default_value(5489), "seed for random generator")
		("orientation-type,T", po::value(&orientation_type)->default_value('M'), "internal orientation type M=Manifold, E=Euler, S=Scaled-axis, Q=Quaternion as R^4")
		;
	po::variables_map vm;
	try{
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);    
		
	} catch (std::exception &e) {
		std::cerr << e.what() << std::endl;
		std::cout << desc << "\n";
		return -1;
	}
	if (vm.count("help")) {
		std::cout << desc << "\n";
		return 1;
	}
	
	std::ifstream logfile(measurements.c_str());
	if( !logfile.good()) {
		std::cerr << "Bad input file\n";
		return -1;
	}
	
	SLOM::Optimizer optim(optiOpts, std::cerr);
	
	int max_steps = optiOpts.optimization_steps;
	
	SLOM::Estimator &e = optim.getEstimator();
	
	if(noise > 0) {
		if(covariance != 0) {
			std::cerr << "For artificial noise, original noise is ignored";
			covariance = 0;
		}
	}
	
	if(dim == 3) {
		typedef IRPHolder<Pose3Dexp> IRPHolder3D;
		
		if(toro)
			std::cerr << "--toro option is ignored in 3D mode";
		input_angles = std::toupper(input_angles);
		IRPHolder3D::ParsePose pPose;
		if(input_angles == 'Q')
			pPose = parsePoseQuat;
		else if(input_angles == 'E')
			pPose = parsePoseEuler;
		else {
			std::cerr << "--input-angles must be Q or E (quaternion or euler)\n";
			return -1;
		}
		IRPHolder3D::ParseCov  pCov = 0;
		if(std::abs(covariance) == 1){
			pCov = parseCov<SLOM::CholeskyMode::CHOLESKY_FULL,6>;
		} else if(std::abs(covariance)==2){
			pCov = parseCov<SLOM::CholeskyMode::COPY_UPPER_FULL,6>;
		}
		boost::scoped_ptr<IRPHolder3D> rp;
		
		switch(std::toupper(orientation_type)){
		case 'M': rp.reset(new RPHolder<Pose3Dexp   , Pose3Dexp>(optim, seed)); break;
		case 'E': rp.reset(new RPHolder<Pose3Deuler , Pose3Dexp>(optim, seed)); break;
		case 'S': rp.reset(new RPHolder<Pose3Dscaled, Pose3Dexp>(optim, seed)); break;
		case 'Q': rp.reset(new RPHolder<Pose3Dquat4 , Pose3Dexp>(optim, seed)); break;
		default:
			std::cerr << "Invalid internal oriention chosen. Valid are: M,E,S,Q\n";
			return -1;
		}
		
		rp->parse(logfile, init, pPose, pCov, -covariance, noise);
		try{
			std::ofstream out(statsfile.c_str());
			rp->iterate(max_steps, output, statsfile.empty() ? std::cout : out);
		} catch (const char* err) {
			std::cerr << err << std::endl;
		}
	} else if(dim == 2) {
		typedef RPHolder<Pose2D> RP;
		RP rp(optim, seed);
		RP::ParseCov  pCov = 0;
		if(std::abs(covariance) == 1){
			pCov = parseCov<SLOM::CholeskyMode::CHOLESKY_FULL,3>;
			if(toro) {
				pCov = parseCovToro<SLOM::CholeskyMode::CHOLESKY_FULL>;
			}
		} else if(std::abs(covariance)==2){
			pCov = parseCov<SLOM::CholeskyMode::COPY_UPPER_FULL,3>;
			if(toro) {
				pCov = parseCovToro<SLOM::CholeskyMode::COPY_UPPER_FULL>;
			}
		}
		
		rp.parse(logfile, init, parsePose2D, pCov, -covariance, noise);
		
		try{
			std::ofstream out(statsfile.c_str());
			rp.iterate(max_steps, output, statsfile.empty() ? std::cout : out);
		} catch (const char* err) {
			std::cerr << err << std::endl;
		}
		
	} else { // dim !=2 && dim !=3
		std::cerr << "Dimension must be 2 or 3";
		return -1;
	}
	if(!dump_jacobian.empty())
	{
		std::ofstream out(dump_jacobian.c_str());
		e.dumpJacobian(out);
	}
	
	std::cout << "\nFinished\n";

	return 0;
}
