% NAOCALIB

%  Copyright (c) 2011 DFKI GmbH
%  All rights reserved
%
%  Author: Oliver Birbach <oliver.birbach@dfki.de>
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions
%  are met:
%
%   * Redistributions of source code must retain the above copyright
%     notice, this list of conditions and the following disclaimer.
%   * Redistributions in binary form must reproduce the above
%     copyright notice, this list of conditions and the following
%     disclaimer in the documentation and/or other materials provided
%     with the distribution.
%   * Neither the name of DFKI GmbH nor the names of any
%     contributors may be used to endorse or promote products derived
%     from this software without specific prior written permission.
%
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
%  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
%  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
%  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
%  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
%  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
%  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
%  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
%  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
%  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%  POSSIBILITY OF SUCH DAMAGE.
%

clear classes

%read_images;

% declare parameters (and load data, since num_images is needed)
load matlab.mat

% define manifold types
intrinsics = mtk.make_compound('focal_length', @mtk.Rn, 1, ...
                               'offset', @mtk.Rn, 2, ...
                               'distortion', @mtk.Rn, 1);                                

% intrinsic parameters for each camera
i_top = intrinsics();
i_bottom = intrinsics();

% extrinsic parameters for each image
extrinsic_for_image = cell(both_cam_images);
for i=1:both_cam_images
    extrinsic_for_image{i} = mtk.SE3();
end;

extrinsic_for_top_only = cell(top_only_images);
for i=1:top_only_images
    extrinsic_for_top_only{i} = mtk.SE3();
end;

extrinsic_for_bottom_only = cell(bottom_only_images);
for i=1:bottom_only_images
    extrinsic_for_bottom_only{i} = mtk.SE3();
end;

% transformation between both cameras
top2bottom = mtk.SE3(); 

% compute initial intrinsic, extrinsic, top2bottom
compute_initial_parameters;

% create optimization probelm
o = mtk.OptimizationProblem();

% add intrinsics
h_i_top = o.add_random_var(i_top);
h_i_bottom = o.add_random_var(i_bottom);

% add stereo
h_stereo_calib = o.add_random_var(top2bottom);

% add extrinsic parameters and measurements
h_extrinsic_for_image = cell(both_cam_images);
h_extrinsic_for_top_only = cell(top_only_images);
h_extrinsic_for_bottom_only = cell(bottom_only_images);

for i=1:both_cam_images
    % get measurements
    eval(['z_top = transpose(p_image_top' num2str(i) ');']);
    eval(['p_top = transpose(p_grid_top' num2str(i) ');']);
    eval(['z_bottom = transpose(p_image_bottom' num2str(i) ');']);
    eval(['p_bottom = transpose(p_grid_bottom' num2str(i) ');']);

    % add extrinsic parameter for every image to problem
    h_extrinsic_for_image{i} = o.add_random_var(extrinsic_for_image{i});       
    
    % add measurements
    o.add_measurement(z_top(:), @project_points, {h_i_top, h_extrinsic_for_image{i}}, {mtk.SE3(eye(4)), p_top(:)}, 0.25);    
    o.add_measurement(z_bottom(:), @project_points, {h_i_bottom, h_extrinsic_for_image{i}, h_stereo_calib}, {p_bottom(:)}, 0.25)
end

for i=1:top_only_images
    % get measurements
    eval(['z = transpose(p_image_top_only' num2str(i) ');']);
    eval(['p = transpose(p_grid_top_only' num2str(i) ');']);
    
    % add extrinsic parameter for every image to problem
    h_extrinsic_for_top_only{i} = o.add_random_var(extrinsic_for_top_only{i});       
    
    % add measurements
    o.add_measurement(z(:), @project_points, {h_i_top, h_extrinsic_for_top_only{i}}, {mtk.SE3(eye(4)), p(:)}, 0.25);
end

for i=1:bottom_only_images
    % get measurements
    eval(['z = transpose(p_image_bottom_only' num2str(i) ');']);
    eval(['p = transpose(p_grid_bottom_only' num2str(i) ');']);
    
    % add extrinsic parameter for every image to problem
    h_extrinsic_for_bottom_only{i} = o.add_random_var(extrinsic_for_bottom_only{i});       
    
    % add measurements
    o.add_measurement(z(:), @project_points, {h_i_bottom, h_extrinsic_for_bottom_only{i}}, {mtk.SE3(eye(4)), p(:)}, 0.25);    
end

% solve problem
[X] = nrlm(@o.fun, o.X0);

% print results
i = X.get_random_var(h_i_top);
disp(['TOP camera -> focal length: ' num2str(i.focal_length.vec) ' offset:  ' num2str(i.offset.vec') ] );
i = X.get_random_var(h_i_bottom);
disp(['BOTTOM camera -> focal length: ' num2str(i.focal_length.vec) ' offset:  ' num2str(i.offset.vec') ] );

disp('top2bottom transformation:');
e = X.get_random_var(h_stereo_calib);
disp(num2str([e.Q e.pos; 0 0 0 1]));

%remove_outliers;

% show plots
%plot_results;







