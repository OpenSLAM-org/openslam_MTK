% KINECTCALIB

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

% define manifold types
intrinsics = mtk.make_compound('focal_length', @mtk.Rn, 1, ...
                               'offset', @mtk.Rn, 2, ...
                               'distortion', @mtk.Rn, 1);
pose = mtk.make_compound('orientation', @mtk.SO3, ...   
                         'position', @mtk.Rn, 3);
                                         
accelerometer_compound = mtk.make_compound('rgb2acc', @mtk.SO3, ...
                                           'scale', @mtk.Rn, 1);
                                       
kinect_depth_compound = mtk.make_compound('base_line', @mtk.Rn, 1, ...
                                         'disparity_offset', @mtk.Rn, 1);

% declare parameters (and load data, since num_images is needed)
load matlab.mat

% intrinsic parameters for each camera
i_rgb = intrinsics();
i_ir = intrinsics();

% extrinsic parameters for each image
extrinsic_for_image = cell(num_images);
for i=1:num_images
    extrinsic_for_image{i} = pose();
end;

% transformation between both cameras
rgb2ir = pose(); 

% rotational between left camera and accelerometer, scale of g
acc_calibration = accelerometer_compound();

% depth calibration
depth_calib = kinect_depth_compound();

% compute initial intrinsic, extrinsic, left2right and acc_calibration parameters
compute_initial_parameters;

% create optimization probelm
o = mtk.OptimizationProblem();

% add intrinsics
h_i_rgb = o.add_random_var(i_rgb);
h_i_ir = o.add_random_var(i_ir);

% add stereo, accelerometer, depth
h_stereo_calib = o.add_random_var(rgb2ir);
h_acc_calibration = o.add_random_var(acc_calibration);
h_depth_calib = o.add_random_var(depth_calib);

% add extrinsic parameters and measurements
h_extrinsic_for_image = cell(num_images);

for i=1:num_images
    % get measurements
    eval(['z_left = transpose(p_image_left' num2str(i) ');']);
    eval(['p_left = transpose(p_grid_left' num2str(i) ');']);
    eval(['z_right = transpose(p_image_right' num2str(i) ');']);
    eval(['p_right = transpose(p_grid_right' num2str(i) ');']);
    eval(['z_depth = transpose(depth' num2str(i) ');']);
    eval(['p_depth = transpose(p_image_depth' num2str(i) ');']);
    eval(['g_acc = v_inertial' num2str(i) ';']);    
    
    % add extrinsic parameter for every image to problem
    h_extrinsic_for_image{i} = o.add_random_var(extrinsic_for_image{i});       
    
    % add measurements
    o.add_measurement(z_left(:), @rgb_pinhole_projection, {h_i_rgb, h_extrinsic_for_image{i}}, {p_left(:)});    
    o.add_measurement(z_right(:), @ir_pinhole_projection, {h_i_ir, h_extrinsic_for_image{i}, h_stereo_calib}, {p_right(:)})

    if (g_acc(1) ~= -1) % test for valid inertial measurement
        o.add_measurement(g_acc, @g_world2g_acc, {h_extrinsic_for_image{i}, h_acc_calibration}, {});
    end;        
    
    % Add depth measurements
    o.add_measurement(double(z_depth(:)), @depth, {h_i_ir, h_extrinsic_for_image{i}, h_stereo_calib, h_depth_calib}, {double(p_depth(:,1:2))'});    
end

% solve problem
[X] = nrlm(@o.fun, o.X0);

% print results
i = X.get_random_var(h_i_rgb);
disp(['RGB camera -> focal length: ' num2str(i.focal_length.vec) ' offset:  ' num2str(i.offset.vec') ' dist: ' num2str(i.distortion.vec) ] );
i = X.get_random_var(h_i_ir);
disp(['IR camera -> focal length: ' num2str(i.focal_length.vec) ' offset:  ' num2str(i.offset.vec') ' dist: ' num2str(i.distortion.vec) ] );

disp('rgb2ir:');
e = X.get_random_var(h_stereo_calib);
disp(num2str([e.orientation.Q e.position.vec; 0 0 0 1]));

disp('rgb2acc:');
e = X.get_random_var(h_acc_calibration);
disp(num2str(e.rgb2acc.Q));

disp('accelerometer scale:');
disp(num2str(e.scale.vec));

disp('depth baseline:');
e = X.get_random_var(h_depth_calib);
disp(num2str(e.base_line.vec));

disp('depth disparity offset:');
disp(num2str(e.disparity_offset.vec));

for i=1:num_images
    e = X.get_random_var(h_extrinsic_for_image{i});
    disp(['extrinsic for image ' num2str(i) ':' ]);
    disp(num2str([e.orientation.Q e.position.vec; 0 0 0 1]));    
end;

% show plots
plot_results;

img_num = 5;
f = figure(100);
hold on;
set(gca,'ydir','reverse')
A = double(imread(['l-' num2str(img_num) '.jpg']));
image(A);
colormap(gray(256));
axis off ;
axis image;
plot(p_image_left5(:,1), p_image_left5(:,2), 'x', 'color', [ 1. .0 .0 ], 'linewidth', 2, 'Markersize', 20);
set(gca, 'Units', 'normalized', 'Position', [0,0,1,1], 'visible', 'off');
f = getframe(gcf);
%imwrite(f.cdata,'kinect_left.png');
%export_fig -native -q75 kinect_left.jpg

img_num = 5;
f = figure(101);
hold on;
set(gca,'ydir','reverse')
A = double(imread(['r-' num2str(img_num) '.jpg']));
image(A);
colormap(gray(256));
axis off ;
axis image;
plot(p_image_right5(:,1), p_image_right5(:,2), 'x', 'color', [ 1. .0 .0 ], 'linewidth', 2, 'Markersize', 20);
set(gca, 'Units', 'normalized', 'Position', [0,0,1,1], 'visible', 'off');
f = getframe(gcf);
%imwrite(f.cdata,'kinect_right.png');
%export_fig -native -q75 kinect_right.jpg

img_num = 5;
f = figure(102);
hold on;
set(gca,'ydir','reverse')
A = double(load_kinect_depth('d-5.bin')) - 505.;
A = A ./ max(max(A)) * 965.; % Augment by factor 3 so the nearby area (that is only visible in the image) is better visible
image(A');
hold on;
colormap(jet(255));
axis off ;
axis image;
plot(p_image_depth5(1,:) ,p_image_depth5(2,:), 'x', 'color', [ 1. .0 .0 ], 'linewidth', 2, 'Markersize', 20);
set(gca, 'Units', 'normalized', 'Position', [0,0,1,1], 'visible', 'off');
f = getframe(gcf);
%imwrite(f.cdata,'kinect_depth.png');
%export_fig -native -q75 kinect_depth.jpg