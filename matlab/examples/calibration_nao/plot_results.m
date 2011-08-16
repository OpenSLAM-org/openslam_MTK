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

clf;

% Define camera cone -- top
i = X.get_random_var(h_i_top);
camera_top = 0.150 * [1/i.focal_length.vec 0 0;0 1/i.focal_length.vec 0;0 0 1]* [1 0 -i.offset.vec(1);0 1 -i.offset.vec(2); 0 0 1] * [0 image_size_x image_size_x 0 0 ; 0 0 image_size_y image_size_y 0;1 1 1 1 1];
for i=1:4
    camera_top = horzcat(camera_top, [0;0;0]);
    camera_top = horzcat(camera_top, camera_top(:,i));
end;
camera_top(4,:) = 1.0;

% Define camera cone -- bottom
i = X.get_random_var(h_i_bottom);
camera_bottom = 0.150 * [1/i.focal_length.vec 0 0;0 1/i.focal_length.vec 0;0 0 1]* [1 0 -i.offset.vec(1);0 1 -i.offset.vec(2); 0 0 1] * [0 image_size_x image_size_x 0 0 ; 0 0 image_size_y image_size_y 0;1 1 1 1 1];
for i=1:4
    camera_bottom = horzcat(camera_bottom, [0;0;0]);
    camera_bottom = horzcat(camera_bottom, camera_bottom(:,i));
end;
camera_bottom(4,:) = 1.0;

% Plot checkerboard pattern + frame
figure(1);
hold on;
pattern = [];
pattern.p_along_y = 9;
pattern.p_along_x = 7;
pattern.space_betw_p_x = 0.11125;
pattern.space_betw_p_y = 0.11125;
pattern.p_center_x = 2;
pattern.p_center_y = 3;

plot_checkerboard_pattern(pattern);
draw_axes(eye(4), 0.15, 3);

% Plot all camera frames and cones as well as accelerometer readings
top2bottom_rv = X.get_random_var(h_stereo_calib);
bottom2top = inv(top2bottom_rv.transform());

for i=1:num_images
    e = X.get_random_var(h_extrinsic_for_image{i});
    
    % top2pattern transformation
    top2pattern_hom = e.transform();
    
    % Draw top camera axes and cone
    draw_axes(top2pattern_hom, 0.1, 2);
    cam_cone = top2pattern_hom * camera_top;
    plot3(cam_cone(1,:),cam_cone(2,:),cam_cone(3,:), 'black-', 'linewidth',1);          
    
    % Draw bottom camera axes and cone
    draw_axes(top2pattern_hom * bottom2top, 0.1, 2);
    cam_cone =  top2pattern_hom * bottom2top * camera_bottom;
    plot3(cam_cone(1,:),cam_cone(2,:),cam_cone(3,:), 'black-', 'linewidth',1);       
end;

axis equal;
grid on;