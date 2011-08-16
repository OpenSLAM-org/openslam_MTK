% ESTIMATE_HOMOGRAPHY

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

function H = estimate_homography(m, M)
% m -- homogenous coordinates in the image image plane
% M -- homogenous coordinates in the world plane

% Get the number of points we are handling here
points = size(m,2);

m = m ./ (ones(3,1)*m(3,:));
M = M ./ (ones(3,1)*M(3,:));

% Normalize according to Hartley's Isotropic Scaling
% Translates the image points so that their centroid is at the origin
% Scale the points, so that their average distance is the mean of the
% distances from all points and the centroid -> NOTE: different from
% Hartley

mean_x = mean(m(1,:)); 
mean_y = mean(m(2,:)); 
m_minus_mean_x = m(1,:) - mean_x; 
m_minus_mean_y = m(2,:) - mean_y; 

mean_m_minus_mean_x = mean(abs(m_minus_mean_x)); 
mean_m_minus_mean_y = mean(abs(m_minus_mean_y));

NORM = [1/mean_m_minus_mean_x  0                      -mean_x/mean_m_minus_mean_x; ...
        0                      1/mean_m_minus_mean_y  -mean_y/mean_m_minus_mean_y; ...
        0                      0                      1];

% Normalize the image coordinates
norm_m = NORM * m;

% Set up the 2*points x 9 matrix according to Zhang
L = zeros(2*points,9);

L(1:2:2*points,1:3) = M';
L(2:2:2*points,4:6) = M';
L(1:2:2*points,7:9) = -((ones(3,1)*norm_m(1,:)).* M)';
L(2:2:2*points,7:9) = -((ones(3,1)*norm_m(2,:)).* M)';

% Solution to Lx = 0 equations, where x is defined up to a scale factor, 
% are well known to be the right singular vector of L associated with 
% the smallest singular value

[U,S,V] = svd(L);

h = V(:,9); % Get the 9th column of the right singular vector
h= h/h(9); % Divide every element by the 9th element

% Bring h back to a 3x3 shape and multiply by the inverse of our
% normalization matrix, since we estimated on normalized values
H = inv(NORM) * reshape(h,3,3)';
