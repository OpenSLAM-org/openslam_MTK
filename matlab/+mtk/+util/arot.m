% AROT determine scaled axis represenation from 3x3 rotation matrix Q.
% http://en.wikipedia.org/wiki/Logarithm_of_a_matrix

%  Copyright (c) 2010, 2011 DFKI GmbH
%  All rights reserved
%
%  Author: Udo Frese <udo.frese@dfki.de>
%          Rene Wagner <rene.wagner@dfki.de>
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

function axis = arot(Q)
    if size(Q) ~= [3 3]
        error('size mismatch')
    end
    
    c = (trace(Q) - 1) / 2.;
    if c < -1
        c = -1;
    elseif c > 1
        c = 1;
    end
       
    angle = acos(c);
    if (angle ~= 0)
        % normalize axis
        factor = angle / (2 * sin(angle));
    else
        % in the limit
        factor = 1/2.;
    end
    
    axis = factor * [Q(3,2) - Q(2,3);
                     Q(1,3) - Q(3,1);
                     Q(2,1) - Q(1,2)];
end
