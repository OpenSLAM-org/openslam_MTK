% COMPUTE_INITIAL_PARAMETERS

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

% compute initial intrinsic parameters
stacked_homographies_rgb = [];
stacked_homographies_ir = [];

stacked_measured_gs = zeros(3, num_images);
for k=1:num_images
    eval(['Hp_grid_left' num2str(k) ' = p_grid_left' num2str(k) ';']);
    eval(['Hp_grid_left' num2str(k) '(:,3) = 1;']);
    
    eval(['Hp_image_left' num2str(k) ' = p_image_left' num2str(k) ';']);
    eval(['Hp_image_left' num2str(k) '(:,3) = 1;']);    
    
    eval(['Hp_grid_right' num2str(k) ' = p_grid_right' num2str(k) ';']);
    eval(['Hp_grid_right' num2str(k) '(:,3) = 1;']);
    
    eval(['Hp_image_right' num2str(k) ' = p_image_right' num2str(k) ';']);
    eval(['Hp_image_right' num2str(k) '(:,3) = 1;']);    
    
    eval(['H_left = estimate_homography(transpose(Hp_image_left' num2str(k) '), transpose(Hp_grid_left' num2str(k) '));']);
    eval(['H_right = estimate_homography(transpose(Hp_image_right' num2str(k) '), transpose(Hp_grid_right' num2str(k) '));']);
    
    stacked_homographies_rgb((k-1)*3+1:(k-1)*3+3, 1:3) = H_left;
    stacked_homographies_ir((k-1)*3+1:(k-1)*3+3, 1:3) = H_right;
    eval(['measured_g = v_inertial' num2str(k) ';']);
    stacked_measured_gs(1:3, k) = measured_g;
end;

% compute and set initial intrinsic parameters
A_rgb = estimate_initial_intrinsic(num_images, stacked_homographies_rgb);
A_ir = estimate_initial_intrinsic(num_images, stacked_homographies_ir);

i_rgb.focal_length.vec = A_rgb(1,1);
i_rgb.offset.vec = [A_rgb(1,3); A_rgb(2,3)];
i_rgb.distortion.vec = 0.0;

i_ir.focal_length.vec = A_ir(1,1);
i_ir.offset.vec = [A_ir(1,3); A_ir(2,3)];
i_ir.distortion.vec = 0.0;

% set initial extrinsic parameters using H and A
for k=1:num_images
    H = stacked_homographies_rgb((k-1)*3+1:(k-1)*3+3, 1:3);
    [extrinsic_for_image{k}.orientation.Q extrinsic_for_image{k}.position.vec] = estimate_initial_extrinsic(H, A_rgb);
    stacked_world2camera_Qs((k-1)*3+1:(k-1)*3+3, 1:3) = extrinsic_for_image{k}.orientation.Q;
end;

%
% compute rgb2ir transformation (stereo)
% use only image, should be accurate enough as an initial guess
[Q_rgb t_rgb] = estimate_initial_extrinsic(stacked_homographies_rgb(1:3, 1:3), A_rgb);
[Q_ir t_ir] = estimate_initial_extrinsic(stacked_homographies_ir(1:3, 1:3), A_ir);

% rgb2acc = world2ir * world2rgb^-1
rgb2ir_hom = [Q_ir t_ir; 0. 0. 0. 1.] * inv([Q_rgb t_rgb; 0. 0. 0. 1.]);

rgb2ir.orientation.Q = rgb2ir_hom(1:3, 1:3);
rgb2ir.position.vec = rgb2ir_hom(1:3, 4);

%
% compute initial accelerometer scale
avg_meas_gs = 0;
for i=1:size(stacked_measured_gs,2)
    avg_meas_gs = avg_meas_gs + norm(stacked_measured_gs(1:3, i)) / size(stacked_measured_gs,2);
end;
acc_calibration.scale.vec = avg_meas_gs / 9.81;

%
% compute left2acc orientation
acc_calibration.rgb2acc.Q = estimate_initial_left2imu(num_images, stacked_world2camera_Qs, stacked_measured_gs ./ acc_calibration.scale.vec);

%
% depth
depth_calib.base_line.vec = 0.075;
depth_calib.disparity_offset.vec = 1090.;