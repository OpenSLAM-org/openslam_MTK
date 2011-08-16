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
    
    start = mtk.SE2();
    start.pos = [0;0];%[v(1, 2); v(1, 3)];
    start.phi = 0;%v(1, 4);
    poses{1} = start;

    nposes_initialized = 1;
    
    nmeasurements = size(e, 1);
    for i=1:nmeasurements
        % NOTE: order differs from TORO
        observing_vertex_id = e(i, 1);
        observed_vertex_id = e(i, 2);

        % look for odomety-like measurements
        if (observed_vertex_id - observing_vertex_id == 1) ...
                && (~isempty(poses{observing_vertex_id + 1}) && ...
                    isempty(poses{observed_vertex_id + 1}))
            
            observing = poses{observing_vertex_id + 1};
            
            forward = e(i, 3);
            sideward = e(i, 4);
            rotate = e(i, 5);
            
            delta_pos_local = [forward; sideward];
            
            observed = mtk.SE2();
            observed.pos = observing.fromlocal(delta_pos_local);
            observed.phi = observing.phi + rotate;
            
            poses{observed_vertex_id + 1} = observed;
            
            nposes_initialized = nposes_initialized + 1;
        end
    end
    
    if nposes_initialized < nposes
        error('initialization incomplete');
    end
end