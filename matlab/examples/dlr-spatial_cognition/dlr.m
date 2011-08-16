% DLR

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

function dlr(filename, problem_only, n)

    % import data
    [odo lm] = parse_dlr_spatial_cognition(filename, n);
    
    % create optimization problem
    global o;
    o = mtk.OptimizationProblem();

    % initialize and add poses
    poses = init_poses(odo);
    poses_ = {poses{2:end,:}}';
    o.add_random_vars(poses_); % not storing any handles since they are
                              % indices that can be computed from the
                              % vertex id
    
    disp(['poses: ' num2str(length(poses))]);
                              
    % add odometry measurements
    nmeasurements = size(odo, 1);
    disp(['odometry measurements: ' num2str(nmeasurements)]);
    
    measurements = cell(nmeasurements, 1);
    for i=1:nmeasurements
        
        observing_vertex_id = odo(i, 1);
        observed_vertex_id = odo(i, 2);
        
        dx = odo(i, 3);
        dy = odo(i, 4);
        dpsi = odo(i, 5);
        
        C_xx = odo(i, 6);
        C_xy = odo(i, 7);
        C_yy = odo(i, 8);
        C_xpsi = odo(i, 9);
        C_ypsi = odo(i, 10);
        C_psipsi = odo(i, 11);

        Sigma = [C_xx   C_xy   C_xpsi;
                 C_xy   C_yy   C_ypsi;
                 C_xpsi C_ypsi C_psipsi];
        
        z = mtk.SE2();
        z.pos = [dx; dy];
        z.phi = dpsi;
        
        if observing_vertex_id == 1
        measurements{i} = {z, @pose_relation_2d_fixed, ...
                           {observed_vertex_id - 1}, ...
                           {poses{observing_vertex_id}}, ...
                           Sigma};
        else            
        measurements{i} = {z, @pose_relation_2d, ...
                           {observing_vertex_id - 1, ...
                           observed_vertex_id - 1}, ...
                           {}, ...
                           Sigma};
        end
    end
    o.add_measurements(measurements);
    
    % add circle landmarks and corresponding measurements
    k = 0;
    q = 0;
    idmap = java.util.Hashtable;

    nmeasurements = size(lm, 1);
    disp(['landmark measurements: ' num2str(nmeasurements)]);

    landmarks = cell(nmeasurements, 1);
    measurements = cell(nmeasurements, 1);
    for i=1:nmeasurements
        
        observing_vertex_id = lm(i, 1);
        observed_landmark_id = lm(i, 2);
        
        if observed_landmark_id < 1
            continue
        end
        
        px = lm(i, 3);
        py = lm(i, 4);
        
        id = idmap.get(observed_landmark_id);
        
        if isempty(id)
            
            % previously unseen landmark
            k = k + 1;
            id = k;
            idmap.put(observed_landmark_id, k);
            
            % initialize and insert
            observing = poses{observing_vertex_id};

            landmarks{id} = observing.fromlocal([px;py]);
            %p = [observing.pos; observing.phi]
            %l = landmarks{id}
        end
        
        C_xx = lm(i, 5);
        C_xy = lm(i, 6);
        C_yy = lm(i, 7);

        Sigma = [C_xx   C_xy;
                 C_xy   C_yy];
        
        z = [px; py];
        
        q = q + 1;
        if observing_vertex_id == 1
        measurements{q} = {z, @landmark_relation_fixed, ...
                           {id + length(poses_)}, ...
                           {poses{observing_vertex_id}}, ...
                           Sigma};
        else
        measurements{q} = {z, @landmark_relation, ...
                           {observing_vertex_id - 1, ...
                           id + length(poses_)}, ...
                           {}, ...
                           Sigma};
        end
    end
    o.add_random_vars(landmarks(1:k,1));
    o.add_measurements(measurements(1:q,1));
    
    disp(['valid landmark measurements: ' num2str(q)]);
    disp(['valid landmarks: ' num2str(idmap.size())]);
    
    disp('problem construction complete')
    
    if problem_only
        return;
    end
    
    t = tic;
    global X;
    global chisq;
    [X chisq] = gn(@o.fun, o.X0);
    toc(t);
end