% REMOVE_OUTLIERS

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

[r J r_orig] = o.fun(X);

res_xy = [r_orig(1:2:end) r_orig(2:2:end)];
dist = sqrt(sum(res_xy.^2, 2));

outlier_idx = find(dist > 1.5)

cnt = 0;
for i=1:num_images
    eval(['p_top = transpose(p_grid_top' num2str(i) ');']);
    eval(['p_bottom = transpose(p_grid_bottom' num2str(i) ');']);
    rmv_top = 0;
    for j=1:length(outlier_idx)       
        if outlier_idx(j) > cnt && outlier_idx(j) < (cnt+length(p_top(:))/3)+1
            eval(['p_grid_top' num2str(i) '(outlier_idx(' num2str(j) ')-' num2str(cnt+rmv_top) ', :) = [];']);
            eval(['p_image_top' num2str(i) '(outlier_idx(' num2str(j) ')-' num2str(cnt+rmv_top) ', :) = [];']);
            rmv_top = rmv_top + 1;
        end;
    end;
    cnt = cnt + length(p_top(:)) / 3;
    rmv_bottom = 0;
    for j=1:length(outlier_idx)
        if outlier_idx(j) > cnt && outlier_idx(j) < (cnt+length(p_bottom(:))/3)+1
            eval(['p_grid_bottom' num2str(i) '(outlier_idx(' num2str(j) ')-' num2str(cnt+rmv_bottom) ', :) = [];']);
            eval(['p_image_bottom' num2str(i) '(outlier_idx(' num2str(j) ')-' num2str(cnt+rmv_bottom) ', :) = [];']);
            rmv_bottom = rmv_bottom + 1;
        end;
    end;
    cnt = cnt + length(p_bottom(:)) / 3;       
end;

for i=1:top_only_images    
    eval(['p = transpose(p_grid_top_only' num2str(i) ');']);
    rmv = 0;
    for j=1:length(outlier_idx)
        if outlier_idx(j) > cnt && outlier_idx(j) < (cnt+length(p(:))/3)+1
            eval(['p_grid_top_only' num2str(i) '(outlier_idx(' num2str(j) ')-' num2str(cnt+rmv) ', :) = [];']);
            eval(['p_image_top_only' num2str(i) '(outlier_idx(' num2str(j) ')-' num2str(cnt+rmv) ', :) = [];']);
            rmv = rmv + 1;
        end;
    end;    
    cnt = cnt + length(p(:))/3;
end;

for i=1:bottom_only_images    
    eval(['p = transpose(p_grid_bottom_only' num2str(i) ');']);   
    rmv = 0;
    for j=1:length(outlier_idx)
        if outlier_idx(j) > cnt && outlier_idx(j) < (cnt+length(p(:))/3)+1                        
            eval(['p_grid_bottom_only' num2str(i) '(outlier_idx(' num2str(j) ')-' num2str(cnt+rmv) ', :) = [];']);
            eval(['p_image_bottom_only' num2str(i) '(outlier_idx(' num2str(j) ')-' num2str(cnt+rmv) ', :) = [];']);           
            rmv = rmv + 1;
        end;
    end;        
    cnt = cnt + length(p(:))/3;
end;
