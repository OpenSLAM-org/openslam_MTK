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
stacked_homographies_top = [];
stacked_homographies_bottom = [];

for k=1:both_cam_images
    eval(['Hp_grid_top' num2str(k) ' = p_grid_top' num2str(k) ';']);
    eval(['Hp_grid_top' num2str(k) '(:,3) = 1;']);
    
    eval(['Hp_image_top' num2str(k) ' = p_image_top' num2str(k) ';']);
    eval(['Hp_image_top' num2str(k) '(:,3) = 1;']);    
    
    eval(['Hp_grid_bottom' num2str(k) ' = p_grid_bottom' num2str(k) ';']);
    eval(['Hp_grid_bottom' num2str(k) '(:,3) = 1;']);
    
    eval(['Hp_image_bottom' num2str(k) ' = p_image_bottom' num2str(k) ';']);
    eval(['Hp_image_bottom' num2str(k) '(:,3) = 1;']);    
    
    eval(['H_top = estimate_homography(transpose(Hp_image_top' num2str(k) '), transpose(Hp_grid_top' num2str(k) '));']);
    eval(['H_bottom = estimate_homography(transpose(Hp_image_bottom' num2str(k) '), transpose(Hp_grid_bottom' num2str(k) '));']);
    
    stacked_homographies_top((k-1)*3+1:(k-1)*3+3, 1:3) = H_top;
    stacked_homographies_bottom((k-1)*3+1:(k-1)*3+3, 1:3) = H_bottom;    
end;

% compute and set initial intrinsic parameters
A_top = estimate_initial_intrinsic(num_images, stacked_homographies_top);
A_bottom = estimate_initial_intrinsic(num_images, stacked_homographies_bottom);

i_top.focal_length.vec = A_top(1,1);
i_top.offset.vec = [A_top(1,3); A_top(2,3)];

i_bottom.focal_length.vec = A_bottom(1,1);
i_bottom.offset.vec = [A_bottom(1,3); A_bottom(2,3)];

% set initial extrinsic parameters using H and A
for k=1:both_cam_images
    H = stacked_homographies_top((k-1)*3+1:(k-1)*3+3, 1:3);
    [R t] = estimate_initial_extrinsic(H, A_top);
    T = inv([R t; 0 0 0 1]);
    extrinsic_for_image{k}.Q = T(1:3, 1:3);
    extrinsic_for_image{k}.pos = T(1:3, 4);
    stacked_world2camera_Qs((k-1)*3+1:(k-1)*3+3, 1:3) = R;
end;

%
% compute rgb2ir transformation (stereo)
% use only image, should be accurate enough as an initial guess
[Q_top t_top] = estimate_initial_extrinsic(stacked_homographies_top(1:3, 1:3), A_top);
[Q_bottom t_bottom] = estimate_initial_extrinsic(stacked_homographies_bottom(1:3, 1:3), A_bottom);

% top2bottom = world2bottom * world2top^-1
top2bottom_hom = [Q_bottom t_bottom; 0. 0. 0. 1.] * inv([Q_top t_top; 0. 0. 0. 1.]);

top2bottom.Q = top2bottom_hom(1:3, 1:3);
top2bottom.pos = top2bottom_hom(1:3, 4);


% Now compute the initial extrinsic parameters for the top_only* and
% bottom_only* images
for k=1:top_only_images
    eval(['Hp_grid_top_only' num2str(k) ' = p_grid_top_only' num2str(k) ';']);
    eval(['Hp_grid_top_only' num2str(k) '(:,3) = 1;']);
    
    eval(['Hp_image_top_only' num2str(k) ' = p_image_top_only' num2str(k) ';']);
    eval(['Hp_image_top_only' num2str(k) '(:,3) = 1;']);           
    
    eval(['H = estimate_homography(transpose(Hp_image_top_only' num2str(k) '), transpose(Hp_grid_top_only' num2str(k) '));']);    
    
    [R t] = estimate_initial_extrinsic(H, A_top);
    T = inv([R t; 0 0 0 1]);
    extrinsic_for_top_only{k}.Q = T(1:3, 1:3);
    extrinsic_for_top_only{k}.pos = T(1:3,4);
end;

for k=1:bottom_only_images
    eval(['Hp_grid_bottom_only' num2str(k) ' = p_grid_bottom_only' num2str(k) ';']);
    eval(['Hp_grid_bottom_only' num2str(k) '(:,3) = 1;']);
    
    eval(['Hp_image_bottom_only' num2str(k) ' = p_image_bottom_only' num2str(k) ';']);
    eval(['Hp_image_bottom_only' num2str(k) '(:,3) = 1;']);           
    
    eval(['H = estimate_homography(transpose(Hp_image_bottom_only' num2str(k) '), transpose(Hp_grid_bottom_only' num2str(k) '));']);    
    
    [R t] = estimate_initial_extrinsic(H, A_bottom);
    T = inv([R t; 0 0 0 1]);
    extrinsic_for_bottom_only{k}.Q = T(1:3, 1:3);
    extrinsic_for_bottom_only{k}.pos = T(1:3, 4);
end;