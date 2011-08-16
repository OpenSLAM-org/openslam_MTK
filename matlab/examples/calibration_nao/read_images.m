% READ_IMAGES

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

global pattern;

pattern = [];
pattern.p_along_y = 9;
pattern.p_along_x = 7;
pattern.space_betw_p_x = 0.11125; % in m
pattern.space_betw_p_y = 0.11125; % in m
pattern.p_center_x = 2;
pattern.p_center_y = 3;

cal = [];
cal.f = [380 380];
cal.offset = [160 120];
cal.kappa = 0.00;

use_cal = true;

global both_cam_images top_only_images bottom_only_images num_images;

both_cam_images = input('Wie viele GEMEINSAME Bilder?', 's');
both_cam_images = str2num(both_cam_images);

for i_n=1:both_cam_images
    % Dateinamen holen und einlesen
    name = strcat(num2str(i_n), '.jpg');    
    A = double(imread(['o-' name]));

    [p_grid p_image] = calculate_image_points(A, pattern, cal, use_cal, 0);
    
    eval(['global p_grid_top' num2str(i_n) ' ;']);
    eval(['global p_image_top' num2str(i_n) ' ;']);  
    eval(['p_grid_top' num2str(i_n) ' = p_grid;']);
    eval(['p_image_top' num2str(i_n) ' = p_image;']);
    
    input('');
    
    % Dateinamen holen und einlesen
    name = strcat(num2str(i_n), '.jpg');
    A = double(imread(['u-' name]));
    
    [p_grid p_image] = calculate_image_points(A, pattern, cal, use_cal, 1);
    
    % Just the image data, assuming that the grid coordinates are the same
    eval(['global p_grid_bottom' num2str(i_n) ' ;']);
    eval(['global p_image_bottom' num2str(i_n) ' ;']);  
    eval(['p_image_bottom' num2str(i_n) ' = p_image;']);
    eval(['p_grid_bottom' num2str(i_n) ' = p_grid;']);
    
    input('');
end;

top_only_images = input('Wie viele Bilder von der OBEREN Kamera?', 's');
top_only_images = str2num(top_only_images);

for i_n=1:top_only_images
    % Dateinamen holen und einlesen
    name = strcat(num2str(i_n), '.jpg');    
    A = double(imread(['oben_only-' name]));

    [p_grid p_image] = calculate_image_points(A, pattern, cal, use_cal, 2);
    
    eval(['global p_grid_top_only' num2str(i_n) ' ;']);
    eval(['global p_image_top_only' num2str(i_n) ' ;']);  
    eval(['p_grid_top_only' num2str(i_n) ' = p_grid;']);
    eval(['p_image_top_only' num2str(i_n) ' = p_image;']);
    
    input('');   
end;

bottom_only_images = input('Wie viele Bilder von der UNTEREN Kamera?', 's');
bottom_only_images = str2num(bottom_only_images);

for i_n=1:bottom_only_images
    % Dateinamen holen und einlesen
    name = strcat(num2str(i_n), '.jpg');    
    A = double(imread(['unten_only-' name]));

    [p_grid p_image] = calculate_image_points(A, pattern, cal, use_cal, 2);
    
    eval(['global p_grid_bottom_only' num2str(i_n) ' ;']);
    eval(['global p_image_bottom_only' num2str(i_n) ' ;']);  
    eval(['p_grid_bottom_only' num2str(i_n) ' = p_grid;']);
    eval(['p_image_bottom_only' num2str(i_n) ' = p_image;']);
    
    input('');   
end;

num_images = both_cam_images;

global image_size_y image_size_x;
[image_size_y, image_size_x] = size(A);