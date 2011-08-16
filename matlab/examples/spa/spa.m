% SPA

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

function spa(filename)

    % import data
    [v e] = parse_toro2d(filename);
    
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
        
        % NOTE: order differs from TORO
        observing_vertex_id = e(i, 1);
        observed_vertex_id = e(i, 2);
        
        forward = e(i, 3);
        sideward = e(i, 4);
        rotate = e(i, 5);
        
        % NOTE: order differs from TORO
        inf_ff = e(i, 6);
        inf_fs = e(i, 7);
        inf_fr = e(i, 8);
        inf_ss = e(i, 9);
        inf_sr = e(i, 10);
        inf_rr = e(i, 11);
        
        omega = [inf_ff inf_fs inf_fr;
                 inf_fs inf_ss inf_sr;
                 inf_fr inf_sr inf_rr];
        
        sigma = inv(omega);
        
        z = mtk.SE2();
        z.pos = [forward; sideward];
        z.phi = rotate;
        
        measurements{i} = {z, @pose_relation_2d, ...
                           {observing_vertex_id + 1, ...
                           observed_vertex_id + 1}, ...
                           {}, ...
                           sigma};
    end
    o.add_measurements(measurements);
    
    disp('problem construction complete')
  
    t = tic;
    global X;
    global chisq;
    [X chisq] = nrlm(@o.fun, o.X0);
    toc(t);
end