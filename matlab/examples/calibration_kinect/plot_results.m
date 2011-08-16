% PLOT_RESULTS

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

% Define camera cone -- RGB
i = X.get_random_var(h_i_rgb);
camera_rgb = 0.15 * [1/i.focal_length.vec 0 0;0 1/i.focal_length.vec 0;0 0 1]* [1 0 -i.offset.vec(1);0 1 -i.offset.vec(2); 0 0 1] * [0 image_size_x image_size_x 0 0 ; 0 0 image_size_y image_size_y 0;1 1 1 1 1];
for i=1:4
    camera_rgb = horzcat(camera_rgb, [0;0;0]);
    camera_rgb = horzcat(camera_rgb, camera_rgb(:,i));
end;
camera_rgb(4,:) = 1.0;

% Define camera cone -- IR
i = X.get_random_var(h_i_ir);
camera_ir = 0.15 * [1/i.focal_length.vec 0 0;0 1/i.focal_length.vec 0;0 0 1]* [1 0 -i.offset.vec(1);0 1 -i.offset.vec(2); 0 0 1] * [0 image_size_x image_size_x 0 0 ; 0 0 image_size_y image_size_y 0;1 1 1 1 1];
for i=1:4
    camera_ir = horzcat(camera_ir, [0;0;0]);
    camera_ir = horzcat(camera_ir, camera_ir(:,i));
end;
camera_ir(4,:) = 1.0;

% Plot checkerboard pattern + frame
figure(1);
hold on;
plot_checkerboard_pattern(pattern);
draw_axes(eye(4), 0.15, 3);

% Plot all camera frames and cones as well as accelerometer readings
rgb2ir_rv = X.get_random_var(h_stereo_calib);
ir2rgb = inv([rgb2ir_rv.orientation.Q rgb2ir_rv.position.vec; 0 0 0 1]);

acc_calibration = X.get_random_var(h_acc_calibration);

g_scale = 0.05; % Scale for plotting g

for i=1:num_images
    e = X.get_random_var(h_extrinsic_for_image{i});
    
    % RGB2pattern transformation
    rgb2pattern_hom = inv([e.orientation.Q e.position.vec; 0 0 0 1]);
    
    % Draw RGB camera axes and cone
    draw_axes(rgb2pattern_hom, 0.1, 2);
    cam_cone = rgb2pattern_hom * camera_rgb;
    plot3(cam_cone(1,:),cam_cone(2,:),cam_cone(3,:), 'black-', 'linewidth',1);          
    
    % Draw IR camera axes and cone
    draw_axes(rgb2pattern_hom * ir2rgb, 0.1, 2);
    cam_cone =  rgb2pattern_hom * ir2rgb * camera_ir;
    plot3(cam_cone(1,:),cam_cone(2,:),cam_cone(3,:), 'black-', 'linewidth',1);
    
    % Get accelerometer measurement and transform into pattern coordinates
    eval(['g_acc = v_inertial' num2str(i) ';']);
    g_pattern = rgb2pattern_hom(1:3,1:3) * inv(acc_calibration.rgb2acc.Q) * (g_acc ./ acc_calibration.scale.vec);    
    
    % Compute lines for world gravity and measured gravity
    gravity_line = [rgb2pattern_hom(1:3,4) rgb2pattern_hom(1:3,4) + g_scale * [0; 0; -9.81]];
    meas_gravity_line = [rgb2pattern_hom(1:3,4) rgb2pattern_hom(1:3,4) + g_scale * g_pattern];
    
    % Draw gravity lines
    line(meas_gravity_line(1,:), meas_gravity_line(2,:), meas_gravity_line(3,:), 'Color', 'magenta', 'linewidth',2);
    line(gravity_line(1,:), gravity_line(2,:), gravity_line(3,:), 'Color', 'black', 'linewidth',2);
end;

axis equal;
grid on;