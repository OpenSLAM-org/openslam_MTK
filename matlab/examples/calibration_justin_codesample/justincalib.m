% JUSTINCALIB

%  Copyright (c) 2011 DFKI GmbH
%  All rights reserved
%
%  Author: Oliver Birbach <oliver.birbach@dfki.de>
%          Rene Wagner <rene.wagner@dfki.de>
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

% load p_img_l, p_world_l. p_img_r, p_world_r, g_imu,
% p_lhand_img_[l,r], p_rhand_img_[l,r]
load_data;

intrinsic_t = mtk.make_compound('focal_length', @mtk.Rn, 1, ...
                               'offset', @mtk.Rn, 2, ...
                               'distortion', @mtk.Rn, 1);

% create and init i_left, i_right, left2right, left2imu,  
% cam2head, p_hand_l, p_hand_r, left2world_for_img
compute_initial_parameters;

o = mtk.OptimizationProblem();

% add intrinsic parameters
i_left_id = o.add_random_var(i_left);
i_right_id = o.add_random_var(i_right); 

% add transformations
left2right_id = o.add_random_var(left2right); 
left2imu_id = o.add_random_var(left2imu); 
cam2head_id = o.add_random_var(cam2head);

% add marker positions on hand 
p_hand_l_id = o.add_random_var(p_hand_l); 
p_hand_r_id = o.add_random_var(p_hand_r);

left2world_id_for_img = cell(num_images);
for i=1:num_images
    % add extrinsic parameter for every image to problem
    left2world_id_for_img{i} = o.add_random_var(left2world_for_img{i});       
        
    o.add_measurement(p_img_l{i}(:), @project_points, {i_left_id, left2world_id_for_img{i}}, {mtk.SE3(eye(4)), p_world_l{i}(:)}, Sigma_chk);   
    o.add_measurement(p_img_r{i}(:), @project_points, {i_right_id, left2world_id_for_img{i}, left2right_id}, {p_world_r{i}(:)}, Sigma_chk);
    
    if (g_imu{i}(1) ~= -1) % test for valid IMU measurement
        o.add_measurement(g_imu{i}(:), @g_world2g_imu, {left2world_id_for_img{i}, left2imu_id}, {}, Sigma_g);
    end;
end

for i=1:size(p_lhand_img_l, 1) % add left arm measurements
    o.add_measurement(p_lhand_img_l(i,:)', @project_marker, {i_left_id, cam2head_id, p_hand_l_id}, {mtk.SE3(eye(4)), lhand2head{i,1}}, Sigma_hand);
    o.add_measurement(p_lhand_img_r(i,:)', @project_marker, {i_right_id, cam2head_id, p_hand_l_id, left2right_id}, {lhand2head{i,1}}, Sigma_hand);
end;

for i=1:size(p_rhand_img_l, 1) % add right arm measurements    
    o.add_measurement(p_rhand_img_l(i,:)', @project_marker, {i_left_id, cam2head_id, p_hand_r_id}, {mtk.SE3(eye(4)), rhand2head{i,1}}, Sigma_hand);
    o.add_measurement(p_rhand_img_r(i,:)', @project_marker, {i_right_id, cam2head_id, p_hand_r_id, left2right_id}, {rhand2head{i,1}}, Sigma_hand);
end;

[X] = nrlm(@o.fun, o.X0); % solve problem

% print results
i = X.get_random_var(i_left_id);
disp(['left camera -> focal length: ' num2str(i.focal_length.vec) ' offset:  ' num2str(i.offset.vec') ' dist: ' num2str(i.distortion.vec) ] );
i = X.get_random_var(i_right_id);
disp(['right camera -> focal length: ' num2str(i.focal_length.vec) ' offset:  ' num2str(i.offset.vec') ' dist: ' num2str(i.distortion.vec) ] );

disp('left2right:');
e = X.get_random_var(left2right_id);
disp(num2str(e.transform()));

disp('left2imu');
e = X.get_random_var(left2imu_id);
disp(num2str(e.Q));

disp('left2head');
e = X.get_random_var(cam2head_id);
disp(num2str(e.transform()));

for i=1:num_images
    e = X.get_random_var(left2world_id_for_img{i});
    disp(['left2world for image ' num2str(i) ':' ]);
    disp(num2str(e.transform()));        
end;

% Calculate and print RMS
[r J r_orig] = o.fun(X);

g_xyz = [];
chk = [];
cnt = 0;
for i=1:num_images    
    num_chk_pts = length(p_img_l{i}(:)) + length(p_img_r{i}(:));
    chk = [chk; r_orig(cnt+1:cnt + num_chk_pts)];
    cnt = cnt + num_chk_pts;
    
    if (g_imu{i,1}(1) ~= -1)        
        g_xyz = [g_xyz; r_orig(cnt+1: cnt+3)];
        cnt = cnt + 3;
    end;
end;

hand = r_orig(cnt+1:end);

% RMS Checkerboard measurements
chk = reshape(chk, 2, [])';
rms_chk = sqrt(sum(chk.^2,1) / size(chk,1))

% RMS IMU measurements
g_xyz = reshape(g_xyz, 3, [])';
rms_g_xyz = sqrt(sum(g_xyz.^2, 1) / size(g_xyz,1))

% RMS Hand marker measurements
hand = reshape(hand, 2, [])';
rms_hand = sqrt(sum(hand.^2,1) / size(hand,1))

% Print single checkerboard image with projected world points
img_num = 4;
f = figure(100);
hold on;
set(gca,'ydir','reverse')
i_right = X.get_random_var(i_right_id);
left2world = X.get_random_var(left2world_id_for_img{img_num});
left2right = X.get_random_var(left2right_id);
z = project_points(i_right, left2world, left2right, p_world_r{img_num}(:));
A = double(imread(['r-' num2str(img_num) '.jpg']));
image(A);
colormap(gray(256));
axis off ;
axis image;
plot(z(1:2:end), z(2:2:end), 'x', 'color', [ 1. .0 .0 ], 'linewidth', 2, 'Markersize', 10);
set(gca, 'Units', 'normalized', 'Position', [0,0,1,1], 'visible', 'off');
f = getframe(gcf);
imwrite(f.cdata,'checkerboard.png');

% TODO - remove code duplication
% Print marker on hand image with projected marker point
img_num = 9;
ja_num = 3;
f = figure(101);
hold on;
set(gca,'ydir','reverse')
i_left = X.get_random_var(i_left_id);
cam2head = X.get_random_var(cam2head_id);
left2right = X.get_random_var(left2right_id);
p_hand_l = X.get_random_var(p_hand_l_id);
z = project_marker(i_left, cam2head, p_hand_l, mtk.SE3(eye(4)), lhand2head{ja_num,1});
A = double(imread(['l-' num2str(img_num) 'JA.jpg']));
image(A);
colormap(gray(256));
axis off ;
axis image;
plot(z(1:2:end), z(2:2:end), 'o', 'color', [ 1. .0 .0 ], 'linewidth', 3, 'Markersize', 20);
plot(z(1:2:end), z(2:2:end), 'x', 'color', [ 1. .0 .0 ], 'linewidth', 2, 'Markersize', 10);
set(gca, 'Units', 'normalized', 'Position', [0,0,1,1], 'visible', 'off');
f = getframe(gcf);
imwrite(f.cdata,'l-9JA_augmented.png');

% Print marker on hand image with projected marker point
img_num = 18;
ja_num = 2;
f = figure(101);
hold on;
set(gca,'ydir','reverse')
i_right = X.get_random_var(i_right_id);
cam2head = X.get_random_var(cam2head_id);
left2right = X.get_random_var(left2right_id);
p_hand_r = X.get_random_var(p_hand_r_id);
z = project_marker(i_right, cam2head, p_hand_r, left2right, rhand2head{ja_num,1});
A = double(imread(['r-' num2str(img_num) 'JA.jpg']));
image(A);
colormap(gray(256));
axis off ;
axis image;
plot(z(1:2:end), z(2:2:end), 'o', 'color', [ 1. .0 .0 ], 'linewidth', 3, 'Markersize', 20);
plot(z(1:2:end), z(2:2:end), 'x', 'color', [ 1. .0 .0 ], 'linewidth', 2, 'Markersize', 10);
set(gca, 'Units', 'normalized', 'Position', [0,0,1,1], 'visible', 'off');
f = getframe(gcf);
imwrite(f.cdata,'r-18JA_augmented.png');

% Please note that no frame plot including arm configuration is included as
% its used model is not publicly available 