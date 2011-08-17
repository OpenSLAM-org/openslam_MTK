% MONOCALIB

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

% declare parameters (and load data, since num_images is needed)
load matlab.mat

intrinsic = intrinsics();
extrinsic_for_image = cell(num_images);
for i=1:num_images
    extrinsic_for_image{i} = pose();
end;

% compute initial intrinsic and extrinsic parameters
compute_initial_parameters;

% create optimization probelm
o = mtk.OptimizationProblem();

% add intrinsics
h_intrinsic = o.add_random_var(intrinsic);

% add extrinsic parameters and measurements
h_extrinsic_for_image = cell(num_images);

for i=1:num_images
    % get measurements
    eval(['z = transpose(p_image_left' num2str(i) ');']);
    eval(['p = transpose(p_grid_left' num2str(i) ');']);
    
    % add extrinsic parameter for every image to problem
    h_extrinsic_for_image{i} = o.add_random_var(extrinsic_for_image{i});

    % add measurements     
    o.add_measurement(z(:), @pinhole_projection, {h_intrinsic, h_extrinsic_for_image{i}}, {p(:)});        
end

% solve problem
[X] = nrlm(@o.fun, o.X0);

% print results
i = X.get_random_var(h_intrinsic);
disp(['focal length: ' num2str(i.focal_length.vec) ' offset:  ' num2str(i.offset.vec') ' dist: ' num2str(i.distortion.vec) ] )
for i=1:num_images
    e = X.get_random_var(h_extrinsic_for_image{i});
    disp(['extrinsic for image ' num2str(i) ':' ]);
    disp(num2str([e.orientation.Q e.position.vec; 0 0 0 1]));    
end;

% plot results
plot_results;