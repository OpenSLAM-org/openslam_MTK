% INIT_POSES

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

function poses = init_poses(v, e)
    nposes = size(v,1);
    
    poses = cell(nposes, 1);
    
    start = mtk.SE3();
    start.pos = zeros(3,1);
    start.Q = eye(3);
    poses{1} = start;

    nposes_initialized = 1;
    
    nmeasurements = size(e, 1);
    for i=1:nmeasurements
        observing_vertex_id = e(i, 1);
        observed_vertex_id = e(i, 2);

        % look for odometry-like measurements
        %if (abs(observed_vertex_id - observing_vertex_id) == 1)

            %fprintf('%d -> %d\n', observing_vertex_id, observed_vertex_id);
                
            x = e(i, 3);
            y = e(i, 4);
            z = e(i, 5);
            roll = e(i, 6);
            pitch = e(i, 7);
            yaw = e(i, 8);
        
            delta_pos_local = [x; y; z];
            Q = toro_euler_to_dcm(roll, pitch, yaw);

            if (~isempty(poses{observing_vertex_id + 1}) && ...
                    isempty(poses{observed_vertex_id + 1}))
            
                observing = poses{observing_vertex_id + 1};
        
                observed = mtk.SE3();
                observed.pos = observing.fromlocal(delta_pos_local);
                observed.Q = observing.Q * Q;
            
                poses{observed_vertex_id + 1} = observed;
            
                nposes_initialized = nposes_initialized + 1;
            
            elseif (isempty(poses{observing_vertex_id + 1}) && ...
                    ~isempty(poses{observed_vertex_id + 1}))
            
                observed = poses{observed_vertex_id + 1};
        
                observing = mtk.SE3();
                observing.Q = observed.Q * ...
                    inv(Q);
                observing.pos = observed.pos - ...
                    observing.Q * delta_pos_local;
            
                poses{observing_vertex_id + 1} = observing;
                
                nposes_initialized = nposes_initialized + 1;
            
            end
        %end
    end
    
    if nposes_initialized < nposes
        init = nposes_initialized
        n = nposes
        error('initialization incomplete');
    end
end