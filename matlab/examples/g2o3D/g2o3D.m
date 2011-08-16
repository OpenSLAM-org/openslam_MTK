% G2O3D

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

function g2o3D(filename)

    % import data
    [v e] = parse_g2o(filename);
    
    % create optimization problem
    global o;
    o = mtk.OptimizationProblem();

    % initialize and add poses
    poses = init_poses(v, e);
    poses_ = {poses{2:end,:}}';
    o.add_random_vars(poses_); % not storing any handles since they are
                              % indices that can be computed from the
                              % vertex id
                              
    % add measurements
    nmeasurements = size(e, 1);
    measurements = cell(nmeasurements, 1);
    for i=1:nmeasurements
        
        observing_vertex_id = e(i, 1);
        observed_vertex_id = e(i, 2);
        
        pos = e(i, 3:5)';
        q = e(i, 6:9)';
        
        U = e(i,10:end);
        Sigma = inv(sym_from_U(U, 6));
        
        meas = mtk.SE3();
        meas.pos = pos;
        meas.Q = q_to_dcm(q);

        % sanity checks
        if abs(norm(q)-1) > 1e-5
            error('not a unit quaternion');
        end
        if abs(det(meas.Q)-1) > 1e-5
            error('det not 1');
        end

        % initial pose is fixed
        if observing_vertex_id == 0
            measurements{i} = {meas, @pose_relation_3d_fixed, ...
                               {observed_vertex_id}, ...
                               {poses{1}}, ...
                               Sigma};                           
        elseif observed_vertex_id == 0
            measurements{i} = {meas, @pose_relation_3d, ...
                               {observing_vertex_id}, ...
                               {poses{1}}, ...
                               Sigma};                          
        else
            measurements{i} = {meas, @pose_relation_3d, ...
                               {observing_vertex_id, ...
                                observed_vertex_id}, ...
                               {}, ...
                               Sigma};
        end
    end
    o.add_measurements(measurements);
    
    disp('problem construction complete')
    
    t = tic;
    global X;
    global chisq;
    %[X chisq] = nrlm(@o.fun, o.X0, [40, 4, 1e-3]);
    [X chisq] = gn(@o.fun, o.X0);
    toc(t);
end