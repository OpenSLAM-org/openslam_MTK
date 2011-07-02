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
 * @file slom-sba.cc
 * 
 * This programs reads and optimizes Sparse Bundle Adjustment problems as provided by the SBA package
 */


//#include <mtk/tools/MakePose.hpp>
#include <mtk/types/pose.hpp>
#include <mtk/types/SOn.hpp>

#include <slom/BuildMeasurement.hpp>

#include <slom/Estimator.hpp>
#include <slom/CallBack.hpp>


#include <vector>
#include <iostream>
#include <fstream>

///////// Required random variables //////////

//MTK_MAKE_POSE3D(pose, double, pos, orient, )
typedef MTK::trafo<MTK::SO3<double> > pose;
SLOM_RANDOMVAR_ID(pose)

typedef MTK::vect<3, double> vec3;
SLOM_RANDOMVAR_ID(vec3)

typedef MTK::vect<2, double> vec2;
typedef MTK::Scalar<double> real;

//////////// Simple (normalized) camera model ///////////

SLOM_BUILD_MEASUREMENT(point_meas_normalized, 2, 
                       ((pose, pose_)) ((vec3, point)), 
                       ((vec2, coordinates)))
SLOM_IMPLEMENT_MEASUREMENT(point_meas_normalized, ret)
{
	vec3 local = pose_->world2Local(*point);
	ret = local.head<2>() / local(2) - coordinates;
}




/// camera variable and measurement with arbitrary linear camera ///

SLOM_BUILD_RANDOMVAR(camera, 
                     ((real, focal_length)) 
                     ((vec2, principal_point)) 
                     ((real, aspect_ratio)) 
                     ((real, skew))
)

SLOM_BUILD_MEASUREMENT(point_meas_cam, 2, 
                       ((pose, pose_)) ((vec3, point)) ((camera, cam)), 
                       ((vec2, coordinates)))
SLOM_IMPLEMENT_MEASUREMENT(point_meas_cam, ret)
{
	vec3 p_local = pose_->local2World(*point);
//			pose_->orient * *point + pose_->pos; 
			
			//pose_->world2Local(*point);
	vec2 p_pinhole = p_local.head<2>() / p_local(2);
	
	
	double f = cam->focal_length, a_r = cam->aspect_ratio, s = cam->skew;
	const vec2& c = cam->principal_point;
	
	Eigen::Vector2d focal_xy(f, f * a_r);
#if MTK_EIGEN < 300
	vec2 p_cam = p_pinhole.cwise() * focal_xy + c;
#else
	vec2 p_cam = p_pinhole.cwiseProduct(focal_xy) + c;
#endif
	p_cam.x() += p_pinhole.y() * s;
	
	ret = p_cam - coordinates;
}




/**
 * Simulates behavior of sba_driver using SLoM
 */
void slom_sba(std::istream &cams, std::istream &pts, std::istream *calib,
              std::ostream *pose_out, std::ostream *points_out)
{
	TicToc t;
	SLOM::Estimator e;
	SLOM::CallBack cb(e, 10);
	if(0)
		e.setCallBack(cb);
	
	std::vector<std::pair<pose_id, camera_id> > poses;
	std::vector<vec3_id> landmarks;
	
	std::cout.precision(16);
	
	camera_id global_cam;
	
	if(calib)
	{
		double f_x, f_y, c_x, c_y, s, dummy;
		*calib >> f_x >> s >> c_x >> dummy >> f_y >> c_y;
		
		// insert camera, but don't optimize:
		global_cam = e.insertRV(camera(f_x, Eigen::Vector2d(c_x, c_y), f_y/f_x, s), false);
	}
	
	bool optimize = false;
	
	while(cams.good())
	{
		camera_id cam = global_cam;
		if(!cam)
		{
			double f_x, a_r, c_x, c_y, s;
			cams >> f_x >> c_x >> c_y >> a_r >> s;
			cam = e.insertRV(camera(f_x, Eigen::Vector2d(c_x, c_y), a_r, s), true);
		}
		Eigen::Quaterniond q;
		Eigen::Vector3d vec;
		cams >> q.w();
		for(int i=0; i<3; ++i) cams >> q.vec()[i];
		for(int i=0; i<3; ++i) cams >> vec[i];
		
		if(!cams.good())
			break;
		
		pose_id cur_pose = e.insertRV(pose(vec, q), optimize);
		
		poses.push_back(std::make_pair(cur_pose, cam));
		optimize = true;
	}
	
	while(pts.good())
	{
		vec3 lm;
		for(int i=0; i<3; ++i)
		{
			pts >> lm[i];
		}
		int nr_of_frames;
		pts >> nr_of_frames;
		
		if(!pts.good())
			break;
		
		vec3_id lm_id = e.insertRV(lm);
		landmarks.push_back(lm_id);
		
		for(int i=0; i<nr_of_frames; ++i)
		{
			int frame_nr;
			vec2 co;
			pts >> frame_nr >> co[0] >> co[1];
			assert(pts.good());
			std::pair<pose_id, camera_id> &frame = poses[frame_nr];
			e.insertMeasurement(point_meas_cam(frame.first, lm_id, frame.second, co));
		}
	}
	
	
	/// optimize
	e.changeAlgorithm(SLOM::Estimator::Levenberg, 1e3);
	
	std::cerr << poses.size() << " " << landmarks.size() << "\nInitialization took " << t() << std::endl;
	
	for(int i=0; i<5000; ++i)
	{
		std::cout << i << " " << t() << " " << e.getRSS() << "\n";
		double gain = e.optimizeStep(false);
		if(0<= gain && gain < 1e-8) break;
	}
	
	std::cerr << "finished optimization in " << ~t; 
	
	
	if(pose_out)
	{
		for(size_t i=0; i<poses.size(); ++i)
		{
			if(!calib)
			{
				const camera& cm = *(poses[i].second);
				*pose_out << cm.focal_length << " "; // TODO
			}
			const pose &p = *(poses[i].first);
			*pose_out << p.orient.w() << " " << p.orient.vec().transpose() << " " << p.pos.transpose() << std::endl; 
		}
	}
	
	if(points_out)
	{
		for(size_t i=0; i<landmarks.size(); ++i)
		{
			*points_out << landmarks[i]->transpose() << std::endl;
		}
		std::ofstream jacobian("jacobian.dat");
		e.dumpJacobian(jacobian);
		std::ofstream residues("residues.dat");
		residues << e.getRes().segment(0,e.getM());
	}
	
}




int main(int argc, char** argv)
{
	std::ofstream poses("poses.dat"), landmarks("points.dat");
	
	
	if(argc >= 3)
	{
		std::ifstream cam_file(argv[1]);
		std::ifstream point_file(argv[2]);
		std::string dummy;
		std::getline(point_file, dummy);
		if(argc >=4)
		{
			std::ifstream calib_file(argv[3]);
			slom_sba(cam_file, point_file, &calib_file, &poses, &landmarks);
		}
		else
		{
			slom_sba(cam_file, point_file, 0, &poses, &landmarks);
		}
	}
	else
	{
		std::cerr << "Usage is " << argv[0] << " <camera params> <point params> [<intrinsic calibration params>]\n";
		return -1;
	}
	
	return 0;
}


