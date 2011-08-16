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

intrinsic = intrinsic_t();
cam2world_for_img = cell(num_images,1);
for i=1:num_images
    cam2world_for_img{i,1} = mtk.SE3();
end;

% compute initial intrinsic parameters
stacked_homographies = [];
for k=1:num_images
    Hp_grid_left = M{k,1};
    Hp_grid_left(3,:) = 1;
    
    Hp_image_left = Z{k,1};
    Hp_image_left(3,:) = 1;
    
    H = estimate_homography(Hp_image_left, Hp_grid_left);
    
    stacked_homographies((k-1)*3+1:(k-1)*3+3, 1:3) = H;
end;

A = estimate_initial_intrinsic(num_images, stacked_homographies);

% set initial intrinsic parameters
intrinsic.focal_length.vec = A(1,1);
intrinsic.offset.vec = [A(1,3); A(2,3)];
intrinsic.distortion.vec = 0.0;

% set initial extrinsic parameters using H and A
for k=1:num_images
    H = stacked_homographies((k-1)*3+1:(k-1)*3+3, 1:3);
    [R t] = estimate_initial_extrinsic(H, A);
    T = inv([R t; 0 0 0 1]);
    cam2world_for_img{k}.Q = T(1:3, 1:3);
    cam2world_for_img{k}.pos = T(1:3, 4);
end;