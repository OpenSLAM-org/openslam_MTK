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

% intrinsic parameters for each camera
i_left = intrinsic_t();
i_right = intrinsic_t();

% extrinsic parameters for each image
left2world_for_img = cell(num_images);
for i=1:num_images
    left2world_for_img{i} = mtk.SE3();
end;

% transformation between both cameras
left2right = mtk.SE3(); 

% rotational between left camera and IMU
left2imu = mtk.SO3();

% camera2head calibration
cam2head = mtk.SE3();

% position of feature pose from last link of arm
p_hand_l = [0;0;0];
p_hand_r = [0;0;0];

% compute initial intrinsic parameters
stacked_homographies_left = [];
stacked_homographies_right = [];

stacked_measured_gs = zeros(3, num_images);
for k=1:num_images
    eval(['Hp_grid_left' num2str(k) ' = p_world_l{' num2str(k) '};']);
    eval(['Hp_grid_left' num2str(k) '(3,:) = 1;']);
    
    eval(['Hp_image_left' num2str(k) ' = p_img_l{' num2str(k) '};']);
    eval(['Hp_image_left' num2str(k) '(3,:) = 1;']);    
    
    eval(['Hp_grid_right' num2str(k) ' = p_world_r{' num2str(k) '};']);
    eval(['Hp_grid_right' num2str(k) '(3,:) = 1;']);
    
    eval(['Hp_image_right' num2str(k) ' = p_img_r{' num2str(k) '};']);
    eval(['Hp_image_right' num2str(k) '(3,:) = 1;']);    
   
    eval(['H_left = estimate_homography(Hp_image_left' num2str(k) ', Hp_grid_left' num2str(k) ');']);
    eval(['H_right = estimate_homography(Hp_image_right' num2str(k) ', Hp_grid_right' num2str(k) ');']);
    
    stacked_homographies_left((k-1)*3+1:(k-1)*3+3, 1:3) = H_left;
    stacked_homographies_right((k-1)*3+1:(k-1)*3+3, 1:3) = H_right;
    eval(['measured_g = g_imu{' num2str(k) '};']);
    stacked_measured_gs(1:3, k) = measured_g;
end;

% compute and set initial intrinsic parameters
A_left = estimate_initial_intrinsic(num_images, stacked_homographies_left);
A_right = estimate_initial_intrinsic(num_images, stacked_homographies_right);

i_left.focal_length.vec = A_left(1,1);
i_left.offset.vec = [A_left(1,3); A_left(2,3)];
i_left.distortion.vec = 0.0;

i_right.focal_length.vec = A_right(1,1);
i_right.offset.vec = [A_right(1,3); A_right(2,3)];
i_right.distortion.vec = 0.0;

% set initial extrinsic parameters using H and A
for k=1:num_images
    H = stacked_homographies_left((k-1)*3+1:(k-1)*3+3, 1:3);
    [R t] = estimate_initial_extrinsic(H, A_left);
    T = inv([R t; 0 0 0 1]);
    left2world_for_img{k}.Q = T(1:3, 1:3);
    left2world_for_img{k}.pos = T(1:3, 4);
    stacked_world2camera_Qs((k-1)*3+1:(k-1)*3+3, 1:3) = R;%left2world_for_img{k}.Q;
end;

%
% compute left2right transformation (stereo)
% use only image, should be accurate enough as an initial guess
[Q_left t_left] = estimate_initial_extrinsic(stacked_homographies_left(1:3, 1:3), A_left);
[Q_right t_right] = estimate_initial_extrinsic(stacked_homographies_right(1:3, 1:3), A_right);

% left2right = world2right * world2left^-1
left2right_hom = [Q_right t_right; 0. 0. 0. 1.] * inv([Q_left t_left; 0. 0. 0. 1.]);

left2right.Q = left2right_hom(1:3, 1:3);
left2right.pos = left2right_hom(1:3, 4);

%
% compute left2imu orientation
left2imu.Q = estimate_initial_left2imu(num_images, stacked_world2camera_Qs, stacked_measured_gs);

%
% compute initial camera2head pose
head_points = zeros(size(p_lhand_img_l,1) + size(p_rhand_img_l,1), 3);
cam_points = zeros(size(p_lhand_img_l,1) + size(p_rhand_img_l,1), 3);

for i=1:size(head_points,1)
        
    % Transform [0;0;0] in last arm coordinates into head coordinates
    if i <= size(p_lhand_img_l,1)        
        left_pt = p_lhand_img_l(i,:);
        right_pt = p_lhand_img_r(i,:);
        hand2head = lhand2head{i,1}.transform();
    else        
        j = i - size(p_lhand_img_l,1);
        hand2head = rhand2head{j,1}.transform();
        left_pt = p_rhand_img_l(j,:);
        right_pt = p_rhand_img_r(j,:);
    end;
    
    p_in_head = hand2head * [0; 0; 0; 1];    
    head_points(i, 1:3) = p_in_head(1:3,1)';            
    
    % Compute corresponding 3D point from stereo        
    [pl vl] = image2SensorRay(left_pt(1), left_pt(2), 0.0, 0.0, 0.0, 0.0, 0.0, [A_left(1,1) A_left(2,2)], [A_left(1,3) A_left(2,3)]);
    [pr vr] = image2SensorRay(right_pt(1), right_pt(2), 0.0, 0.0, 0.0, 0.0, 0.0, [A_right(1,1) A_right(2,2)], [A_right(1,3) A_right(2,3)]);       
    
    pr = inv(left2right_hom) * [pr';1];
    pr = pr(1:3, 1);
    vr = inv(left2right_hom(1:3,1:3)) * vr';    
    
    cam_point = intersectRays(pl', vl', pr, vr);
    cam_points(i, 1:3) = cam_point';      
end;

% Use corresponding point pairs to estimate 4x4 homegenous matrix estimate
for i = 1:size(head_points,1)
    K = head_points(i,1:3)' * cam_points(i,1:3);
end;

% Estimate 3x3 matrix
[V A U] = svd(K);

c = zeros(3,3);
c(1,1) = 1;
c(2,2) = 1;
c(3,3) = det(V*U');

cam2head.Q = V * c * U';

% Set 3x1 vector to zero, since camera and headlink frames are quite near
% to each other
cam2head.pos = [0;0;0];