% TORO3D

%  Copyright (c) 2011 DFKI GmbH
%  All rights reserved
%
%  Author: Rene Wagner <rene.wagner@dfki.de>
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

function toro3d(filename)

    % import data
    [v e] = parse_toro3d(filename);
    
    % create optimization probelm
    global o;
    o = mtk.OptimizationProblem();

    % initialize and add poses
    poses = init_poses(v, e);
    o.add_random_vars(poses); % not storing any handles since they are
                              % indices that can be computed from the
                              % vertex id
    
    % add measurements
    nmeasurements = size(e, 1);
    measurements = cell(nmeasurements, 1);
    for i=1:nmeasurements
        
        observing_vertex_id = e(i, 1);
        observed_vertex_id = e(i, 2);
        
        x = e(i, 3);
        y = e(i, 4);
        z = e(i, 5);
        roll = e(i, 6);
        pitch = e(i, 7);
        yaw = e(i, 8);
        
        if size(e,2) > 8
            U = e(i,9:end);
            Omega = sym_from_U(U, 6);
        
            Sigma = inv(Omega);
        else
            Sigma = eye(6);
        end
        
        meas = mtk.SE3();
        meas.pos = [x; y; z];
        meas.Q = toro_euler_to_dcm(roll, pitch, yaw);
        
        measurements{i} = {meas, @pose_relation_3d, ...
                           {observing_vertex_id + 1, ...
                           observed_vertex_id + 1}, ...
                           {}, ...
                           Sigma};
    end
    o.add_measurements(measurements);
    
    disp('problem construction complete')
    
    t = tic;
    global X;
    global chisq;
    [X chisq] = nrlm(@o.fun, o.X0, [50, 4, 1e-3]); %, @plot_3d_poses);
    toc(t);
end