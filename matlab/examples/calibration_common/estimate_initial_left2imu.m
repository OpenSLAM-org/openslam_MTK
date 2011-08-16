% ESTIMATE_INITIAL_LEFT2IMU

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

function R = estimate_initial_left2imu(num_images, stacked_world2camera_Qs, stacked_measured_gs)
% Compute the initial quaternion representing the rotational difference of
% cam and inertial sensor

Sxx=0; Sxy=0; Sxz=0;
Syx=0; Syy=0; Syz=0;
Szx=0; Szy=0; Szz=0;
for i = 1:num_images
    % Gravity sensed by the inertial sensor in inertial coordinates          
    g_acc = stacked_measured_gs(:, i);
        if (g_acc(1) ~= -1)  
        g_camera = stacked_world2camera_Qs((i-1)*3+1:(i-1)*3+3, 1:3) * [0.; 0.; -9.81];
        
        Sxx=Sxx + g_camera(1) * g_acc(1);
        Sxy=Sxy + g_camera(1) * g_acc(2);
        Sxz=Sxz + g_camera(1) * g_acc(3);
        Syx=Syx + g_camera(2) * g_acc(1);
        Syy=Syy + g_camera(2) * g_acc(2);
        Syz=Syz + g_camera(2) * g_acc(3);
        Szx=Szx + g_camera(3) * g_acc(1);
        Szy=Szy + g_camera(3) * g_acc(2);
        Szz=Szz + g_camera(3) * g_acc(3);
    end;
end;

N = [ (Sxx+Syy+Szz)     Syz-Szy         Szx-Sxz         Sxy-Syx
    Syz-Szy             (Sxx-Syy-Szz)   Sxy+Syx         Szx+Sxz
    Szx-Sxz             Sxy+Syx         (-Sxx+Syy-Szz)  Syz+Szy
    Sxy-Syx             Szx+Sxz         Syz+Szy         (-Sxx-Syy+Szz) ];

% The Eigenvctor associated with the greatest Eigenvalue is our wanted
% quaternion
[V,D] = eig(N);
[B,index] = max(diag(D));
q = V(:,index);

% Chances are high that this quaternion is not a valid unit-quaternion, so
% divide by norm
q = q/norm(q);
R = quat2mat33(q);