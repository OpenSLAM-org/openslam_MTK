% PLOT_3D_POSE_RELATIONS

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

function plot_3d_pose_relations(o, X, nodesize)
    poses = X.rvs;
    nposes = length(poses);
    x = zeros(nposes, 1);
    y = zeros(nposes, 1);
    z = zeros(nposes, 1);
    for i=1:nposes
        p = poses{i};
        x(i) = p.pos(1);
        y(i) = p.pos(2);
        z(i) = p.pos(3);
    end
    
    %hfig = figure;
    %set(hfig,'units','centimeters','NumberTitle','off','Name','figspheremanifold');
    %pos = get(hfig,'position');
    %set(hfig,'position',[pos(1:2),9,9]);

    scatter3(x, y, z, nodesize*ones(length(x), 1), 'filled')
    hold on;
    
    ms = o.measurements;
    nms = length(ms);
    x = zeros(2, nms);
    y = zeros(2, nms);
    z = zeros(2, nms);
    for i=1:length(ms)
        m = ms{i};

        rvhandles = m{3};
        if length(rvhandles) ~= 2
            continue
        end
        rvs = cellfun(@(h) X.rvs{h}, rvhandles, 'UniformOutput', false);
        
        for k=1:2
            x(k, i) = rvs{k}.pos(1);
            y(k, i) = rvs{k}.pos(2);
            z(k, i) = rvs{k}.pos(3);
        end
    end
    plot3(x,y,z,'b');

    grid on
    axis square
    daspect([1 1 1])
    %axis([-100 100 -100 100 -100 100])
    %view([-100 -100 70])
    %set(gca,'XTick',-100:50:100);
    %set(gca,'YTick',-100:50:100);
    %set(gca,'ZTick',-100:50:100);

    drawnow;
end