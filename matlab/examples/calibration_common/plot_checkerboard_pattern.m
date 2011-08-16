% PLOT_CHECKERBOARD_PATTERN

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

function plot_checkerboard_pattern( grid )
p_grid_x = zeros(grid.p_along_y,grid.p_along_x);
p_grid_y = zeros(grid.p_along_y,grid.p_along_x);
p_grid_z = zeros(grid.p_along_y,grid.p_along_x);

p_grid = zeros(grid.p_along_x * grid.p_along_y, 2);
for i=0:grid.p_along_x-1
    x_pos = (i-grid.p_center_x) * grid.space_betw_p_x;
    for j=0:grid.p_along_y-1
        y_pos = (j-grid.p_center_y) * grid.space_betw_p_y;
        p_grid(i*grid.p_along_y+j+1,:) = [x_pos y_pos];
    end;
end;
p_grid(:,3) = 0.0;

p_grid_x(:) = p_grid(:,1);
p_grid_y(:) = p_grid(:,2);
p_grid_z(:) = p_grid(:,3);

px_b = [];
py_b = [];
pz_b = [];
px_w = [];
py_w = [];
pz_w = [];
start_black = false;
for b=1:(grid.p_along_y-1)
    inner_black = start_black;
    for c=1:(grid.p_along_x-1)
        if inner_black==true
            px_b = [px_b ;p_grid_x(b,c) p_grid_x(b,c+1) p_grid_x(b+1,c+1) p_grid_x(b+1,c)];
            py_b = [py_b ;p_grid_y(b,c) p_grid_y(b,c+1) p_grid_y(b+1,c+1) p_grid_y(b+1,c)];
            pz_b = [pz_b ;p_grid_z(b,c) p_grid_z(b,c+1) p_grid_z(b+1,c+1) p_grid_z(b+1,c)];
        else
            
            px_w = [px_w ;p_grid_x(b,c) p_grid_x(b,c+1) p_grid_x(b+1,c+1) p_grid_x(b+1,c)];
            py_w = [py_w ;p_grid_y(b,c) p_grid_y(b,c+1) p_grid_y(b+1,c+1) p_grid_y(b+1,c)];
            pz_w = [pz_w ;p_grid_z(b,c) p_grid_z(b,c+1) p_grid_z(b+1,c+1) p_grid_z(b+1,c)];
        end;
        inner_black = ~inner_black;
    end;
    start_black = ~start_black;
end;

% Plot the checkerboard
patch(px_b', py_b', pz_b', 'black');
patch(px_w', py_w', pz_w', 'white');