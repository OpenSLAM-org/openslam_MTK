% ROT determine 3x3 rotation matrix Q from scaled axis representation.
% http://mathworld.wolfram.com/RodriguesRotationFormula.html

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

function Q = rot(axis)

    if any(size(axis) ~= [3 1])
        error('size mismatch')
    end
    
    x = axis(1,1);
    y = axis(2,1);
    z = axis(3,1);
    
    angle = sqrt(x * x + y * y + z * z);
    
    if (angle > 0)
        x = x / angle;
        y = y / angle;
        z = z / angle;
    end
    
    s = sin(angle);
    c = cos(angle);
    
    v = 1 - c;
    xyv = x * y * v;
    yzv = y * z * v;
    xzv = x * z * v;
    
    Q = [x * x * v + c, xyv - z * s, xzv + y * s;
         xyv + z * s, y * y * v + c, yzv - x * s;
         xzv - y * s, yzv + x * s, z * z * v + c];
end
