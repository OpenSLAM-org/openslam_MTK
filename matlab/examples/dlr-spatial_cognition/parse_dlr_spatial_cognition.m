% PARSE_DLR_SPATIAL_COGNITION

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

function [odo lm] = parse_dlr_spatial_cognition(filename, n)
    fprintf('parsing data set %s...', filename);
    odo = [];
    lm = [];
    
    f = fopen(filename);
    
    i = 0;
    
    % read line-wise
    
    while ~feof(f)
        line = fgetl(f);
        if (line(1) == '#') || (length(line) < 4)
            % skip
            
        elseif all(line(1:4) == 'STEP')
            data = line(35:end);
            step = textscan(data, '%f%f%f%f%f%f%f%f%f');
            
            i = i + 1;
            if n > 0 && i > n
                return;
            end
        
            odo = [odo; i i+1 step{1} step{2} step{3} step{4} step{5} step{6} step{7} step{8} step{9}];
            
        elseif all(line(1:10) == 'LANDMARK_C')
            data = line(12:end);
            l = textscan(data, '%f%f%f%f%f%f%f');
            
            %disp(['landmark at ' num2str(i)]);
            lm = [lm; (i+1) l{7} l{1} l{2} l{4} l{5} l{6}];
            
        end
    end
    
    fprintf(' done.\n');
end
