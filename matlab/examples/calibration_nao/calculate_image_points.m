% CALCULATE_IMAGE_POINTS

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

function [p_grid p_image] = calculate_image_points(A, grid, cal, use_cal, is_top)


figure(1);
clf;
hold on;

%image(A);
%ginput(1);
% Draw the image
%A = rgb2gray(A);
image(A);
colormap(gray(256));
axis off ;
axis image;

%
%
%
% Compute the grid points from the grid properties struct
p_grid = zeros(grid.p_along_x * grid.p_along_y, 2);

for i=0:grid.p_along_x-1
    x_pos = (i-grid.p_center_x) * grid.space_betw_p_x;
    for j=0:grid.p_along_y-1
        y_pos = (j-grid.p_center_y) * grid.space_betw_p_y;
        p_grid(i*grid.p_along_y+j+1,:) = [x_pos y_pos];
    end;
end;
% Create also a set of homgenous coordinates
p_grid_hom = p_grid;
p_grid_hom(:,3) = 1;


%figure(12);
%plot(p_grid(:,1), p_grid(:,2),'rx');
%return;
%
%
%
% Ask the user to mark four points on the grid ...
disp('Bitte die 4. Randpunkte des Musters markieren. 1. - 2. Markierung sind dann die Y- Achse, 4. - 3. Markierung die X-Achse.');

cornerpoints_I = zeros(4,3);
for count=1:4
    p = ginput(1);
    [x y] = localmarkerdetector(p(1), p(2), A, 10);
    plot(x, y, 'x', 'color', [ 1.000 .0 .0 ], 'linewidth', 2, 'Markersize',15);
    cornerpoints_I(count,:) = [x y 1];
    drawnow;
end;

% ... compute the corresponding points in the world ...
if (is_top == 0)
    cornerpoints_W = [ 0.445 -0.1112 1;
        -0.2225 -0.1112 1;
        -0.2225 -0.3338 1;
        0.445 -0.3338 1
        ];
    
    
elseif (is_top == 1)
    cornerpoints_W = [0.445 0.5562 1;
        -0.2225 0.5562 1;
        -0.2225 0.3338 1;
        0.445 0.3338 1
        ];
    
else
    cornerpoints_W = [0 0 1;
        %0 grid.space_betw_p_y 1
        0 grid.space_betw_p_y*2 1;
        %grid.space_betw_p_x grid.space_betw_p_y*2 1;
        grid.space_betw_p_x*2 grid.space_betw_p_y*2 1;
        %grid.space_betw_p_x*2 grid.space_betw_p_y 1;
        grid.space_betw_p_x*2 0 1;
        %grid.space_betw_p_x 0 1
        ];
end;



% and compute the homogrpahy between these two point sets
H_cornerpoints = estimate_homography(cornerpoints_I', cornerpoints_W');

%
%
%
% If the user provides an initial guess about the calibration, use this
% information to get much better initial points for the
% localmarkerdetector
if (use_cal==true)
    % First, undistort the four image points provided by the user...
    undist_cornerpoints = zeros(4,3);
    for i=1:4
        [unx uny] = undistort((cornerpoints_I(i,1)-cal.offset(1)) / cal.f(1), (cornerpoints_I(i,2)-cal.offset(2)) / cal.f(2), 0.0, -cal.kappa, 0.0, cal.kappa, 0.0);
        undist_cornerpoints(i,:) = [unx uny 1];
    end;
    
    % ... and compute a homography between these and the world points.
    H_undist_cornerpoints = estimate_homography(undist_cornerpoints', cornerpoints_W');
    
    % Use this homography to transform the grid points onto the undistorted
    % image space
    p_image_un = H_undist_cornerpoints * p_grid_hom';
    p_image_un = transpose(p_image_un(1:2,:) ./ (ones(2,1)*p_image_un(3,:)));
    p_image_un(:,3) = 1;
    
    % Apply projection so the images will again be on the viewable image
    % plane
    p_image_dist = camera2Image(p_image_un, 0.0, -cal.kappa, 0.0, cal.kappa, 0.0, cal.f, cal.offset);
    p_image = p_image_dist';
else
    % Tranform grid points to the image plane
    p_image = H_cornerpoints * p_grid_hom';
    p_image = p_image(1:2,:) ./ (ones(2,1)*p_image(3,:));
end;

% Plot the resulting set of points
plot(p_image(1,:), p_image(2,:), 'x', 'color', [ .0 1. .0 ], 'linewidth', 2, 'Markersize',10);

%
%
%
% Now refine the rough points using the localmarkerdetector, also mark any
% points that is outside of projection or too close at the edge for
% reasonable detection as NaN
refined_p_image = zeros(grid.p_along_x * grid.p_along_y, 2);
p_is_valid = ones(grid.p_along_x * grid.p_along_y,1);
for i=1:grid.p_along_x * grid.p_along_y
    % Check if the point is in the image
    if p_image(1,i) <= (size(A,2)-2) && p_image(1,i) > (0+2) && p_image(2,i) <= (size(A,1)-2) && p_image(2,i) > (0+2)
        % If yes, refine through localmarkerdetector
        [x y] = localmarkerdetector(p_image(1,i), p_image(2,i), A, 5);
    else
        % If no, set to NaN
        x = NaN;
        y = NaN;
    end;
    
    % Set the point, and mark if its valid
    refined_p_image(i, :) = [x y];
    if isnan(x) || isnan(y)
        p_is_valid(i) = 0;
    end;
end;

% Plot the refined set of points
plot(refined_p_image(:,1), refined_p_image(:,2), 'x', 'color', [ 1. .0 1. ], 'linewidth', 2, 'Markersize', 10);

%
%
%
% Now delete all grid points that would have a corresponding NaN image point
do_delete = true;
while(do_delete == true)
    found_non_valid = false;
    for i=1:size(p_is_valid)
        if (p_is_valid(i) == 0)
            p_is_valid(i) = [];
            refined_p_image(i,:) = [];
            p_grid(i,:) = [];
            found_non_valid = true;
            break;
        end;
    end;
    
    if found_non_valid == false
        break;
    end;
end;
% Assert equal number of points
assert(size(p_grid,1) == size(refined_p_image,1));

drawnow;
hold off;

% Finally, add a zero as third dimension for the grid points ...
p_grid(:,3) = 0;
% ... and the refined set of points as the one to return
p_image = refined_p_image;
end