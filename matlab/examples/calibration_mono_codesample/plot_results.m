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

% Define our pattern
pattern = [];
pattern.p_along_y = 12;
pattern.p_along_x = 8;
pattern.space_betw_p_x = 0.070135;
pattern.space_betw_p_y = 0.070135;
pattern.p_center_x = 2;
pattern.p_center_y = 4;

% Define camera cone
i = X.get_random_var(intrinsic_id);
camera = 0.15 * [1/i.focal_length.vec 0 0;0 1/i.focal_length.vec 0;0 0 1]* [1 0 -i.offset.vec(1);0 1 -i.offset.vec(2); 0 0 1] * [0 image_size_x image_size_x 0 0 ; 0 0 image_size_y image_size_y 0;1 1 1 1 1];
for i=1:4
    camera = horzcat(camera, [0;0;0]);
    camera = horzcat(camera, camera(:,i));
end;
camera(4,:) = 1.0;

% Plot checkerboard pattern + frame
figure(1);
hold on;
plot_checkerboard_pattern(pattern);
draw_axes(eye(4), 0.15, 3);

% Plot camera-frames and cones in the grid relative to grid
for i=[1,3,5,7,9,14,15,16,17,18]%num_images
    e = X.get_random_var(cam2world_id_for_img{i});
    
    grid2cam_hom = e.transform();
    draw_axes(grid2cam_hom, 0.1, 2);

    cam_cone = grid2cam_hom * camera;
    plot3(cam_cone(1,:),cam_cone(2,:),cam_cone(3,:), 'black-', 'linewidth',1);          
end;

axis equal;
grid on;
set(gca,'LineWidth',1); 
set(gca,'LineWidth',1); 
view(20,30);
set(gca,'ztick',[0:0.5:1.5]);
%set(gcf, 'Color', 'none'); % Sets figure background
%set(gca, 'Color', 'none'); % Sets axes background
%export_fig('mono_result.pdf');