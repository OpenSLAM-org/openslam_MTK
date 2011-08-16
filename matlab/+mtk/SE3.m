% SE3 transforms in 3D
%   The vectorized rotation space is formed by relative rotations in world
%   coordinates in scaled axis representation

%  Copyright (c) 2010, 2011 DFKI GmbH
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

classdef SE3 < handle
    
    properties (Constant)
        DOF = 6;
    end
    properties
        pos = zeros(3,1);
        Q = diag(ones(3,1), 0);
    end
    
    methods
        function obj = SE3(p, Q)
            if nargin > 1
                obj.pos = p;
                obj.Q = Q;
            elseif nargin > 0 && isa(p, 'mtk.SE3')
                obj.copyfrom(p);
            elseif nargin > 0
                obj.pos = p(1:3,4);
                obj.Q = p(1:3, 1:3);
            end
        end
        
        function copyfrom(obj, other)
            obj.pos = other.pos;
            obj.Q = other.Q;
        end
        
        function A = transform(obj)
            A = [obj.Q, obj.pos;
                 0 0 0 1];
        end
        
        function v = tolocal(obj, v)
            v = inv(obj.Q) * (v - obj.pos);
        end
        
        function v = fromlocal(obj, v)
            v = obj.Q * v + obj.pos;
        end
        
        function r = plus(s, d)
            if ~(isa(s, 'mtk.SE3') && isa(d, 'numeric'))
                error('undefined')
            else
                % A\B is roughly the same as inv(A)*B
                r = mtk.SE3(s.pos + d(1:3), ...
                            s.Q * mtk.util.rot(s.Q \ d(4:6)));
            end
        end
        
        function r = minus(s2, s1)
            if ~(isa(s2, 'mtk.SE3') && isa(s1, 'mtk.SE3'))
                s2_ = class(s2)
                s1_ = class(s1)
                error('undefined')
            else
                % A\B is roughly the same as inv(A)*B
                r = [s2.pos - s1.pos;
                     s1.Q * mtk.util.arot(s1.Q \ s2.Q)];
            end
        end
        
        function s = size(obj, k)
            s = [obj.DOF 1];
            if nargin > 1
                s = s(k);
            end
        end
        
        function l = length(obj)
            l = obj.DOF;
        end
        
        function l = numel(obj)
            l = obj.DOF;
        end
    end
end