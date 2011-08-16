% STEREO_PINHOLE_PROJECTION

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

function [z, sigma] = stereo_pinhole_projection(intrinsic, extrinsic, stereo_calib, x)
% Measurement function that projects a set of world points onto the image
% plane of the left camera -- computes directly the residual

% project world points into camera frame
reshaped_x = reshape(x, 3, []);

world2top = [extrinsic.orientation.Q; 0 0 0];
world2top(1:4, 4) = [extrinsic.position.vec; 1];

top2bottom = [stereo_calib.orientation.Q; 0 0 0];
top2bottom(1:4, 4) = [stereo_calib.position.vec; 1];

world2bottom = top2bottom * world2top;

x_cam = world2bottom(1:3,1:3) * reshaped_x + repmat(world2bottom(1:3,4), 1, size(reshaped_x, 2));
x_cam = x_cam';

% Project onto the image plane
f = intrinsic.focal_length.vec;
proj = camera2Image(x_cam, 0.0, 0.0, 0.0, 0.0, 0.0, [f f], intrinsic.offset.vec)';

z = proj(:);
sigma = .25; % px