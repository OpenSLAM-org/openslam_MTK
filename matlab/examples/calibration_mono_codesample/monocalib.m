% MONOCALIB

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

load_data;

% define custom manifold compound type
intrinsic_t = mtk.make_compound('focal_length', @mtk.Rn,1, ...
                                'offset', @mtk.Rn,2, ...
                                'distortion', @mtk.Rn,1);

% compute initial intrinsic and extrinsic parameters
compute_initial_parameters;

% create optimization problem
o = mtk.OptimizationProblem();

% add intrinsics
intrinsic_id = o.add_random_var(intrinsic);

cam2world_id_for_img = cell(num_images);
for i=1:num_images  % for each calibration image
    p_img = Z{i}; p_world = M{i};  % get measurements
    
    % add extrinsic parameter (camera pose)
    cam2world_id_for_img{i} = o.add_random_var(cam2world_for_img{i,1});
    
    % add measurements (checkerboard corner points)
    for j=1:size(p_img,2)                           
        o.add_measurement(p_img(:, j), @project_point, {intrinsic_id, cam2world_id_for_img{i}}, {mtk.SE3(eye(4)), p_world(:, j)}, Sigma);
    end;
end

[X] = nrlm(@o.fun, o.X0);  % solve problem

% print results
i = X.get_random_var(intrinsic_id);
disp(['focal length: ' num2str(i.focal_length.vec) ' offset:  ' num2str(i.offset.vec') ' dist: ' num2str(i.distortion.vec) ] )
for i=1:num_images
    e = X.get_random_var(cam2world_id_for_img{i});
    disp(['cam2world for image ' num2str(i) ':' ]);
    disp(num2str(e.transform()));    
end;

% Calculate and print RMS
[r J r_orig] = o.fun(X);

x = r_orig(1:2:end);
y = r_orig(2:2:end);

rms_x = sqrt(sum(x.^2) / length(x))
rms_y = sqrt(sum(y.^2) / length(y))

% plot results
plot_results;