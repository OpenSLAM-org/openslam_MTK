% DEPTH

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

function [z sig] = depth(intrinsic, extrinsic, stereo_calib, depth_calib, p_image_depth)

world2rgb = [extrinsic.orientation.Q extrinsic.position.vec; 0 0 0 1];
rgb2ir = [stereo_calib.orientation.Q stereo_calib.position.vec; 0 0 0 1];

world2ir = rgb2ir * world2rgb;

dist = intrinsic.distortion.vec;
f = intrinsic.focal_length.vec;

z = zeros(size(p_image_depth,2), 1);
for i=1:size(p_image_depth,2)
    [l0 l] = image2SensorRay(double(p_image_depth(1, i)) - 4.0, double(p_image_depth(2, i)) - 4.0, 0.0, -dist, 0.0, dist, 0.0, [f f], intrinsic.offset.vec);

    p0 = world2ir * [0; 0; 0; 1.0];
    n = world2ir(1:3,1:3) * [0; 0; 1.0];

    p0 = p0(1:3,1);
    n = n(1:3, 1);

    % line-plane intersection - compute scale where ray hits plane
    d = (p0 - l0')' * n / (l * n);

    % compute disparity -- equation token from ROS wiki, solved for disparity
    % using Matlab symbolic toolbox
    y = -1. * ((-1. * (d .* l(3)) * depth_calib.disparity_offset.vec + (1. / 0.125) * depth_calib.base_line.vec * f) / (d .* l(3)));
    z(i,1) = y;
end;

sig = 4.0;