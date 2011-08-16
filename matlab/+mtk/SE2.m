% SE2 transforms in the 2D plane    

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

classdef SE2 < handle

    properties (Constant)
        DOF = 3;
    end
    properties
        pos = zeros(2,1);
        phi = 0;
    end
    
    methods
        function obj = SE2(p, phi)
            if nargin > 0
                if isa(p, 'mtk.SE2')
                    obj.copyfrom(p);
                elseif nargin > 1
                    obj.pos = p;
                    obj.phi = phi;
                end
            end
        end
        
        function copyfrom(obj, other)
            obj.pos = other.pos;
            obj.phi = other.phi;
        end
        
        function Q = rotation(obj)
            c = cos(obj.phi);
            s = sin(obj.phi);
            Q = [c -s; s c];
        end
        
        function A = transform(obj)
            A = [obj.rotation(), obj.pos;
                 0 0 1];
        end
        
        function v = fromlocal(obj, v)
            v = obj.rotation() * v + obj.pos;
        end
        
        function v = tolocal(obj, v)
            v = obj.rotation()' * (v - obj.pos);
        end
        
        function set.phi(obj, value)
            obj.phi = mtk.util.normalize_angle(value);
        end
        
        function r = plus(a, b)
            if ~(isa(a, 'mtk.SE2') && isa(b, 'numeric'))
                error('undefined')
            else
                r = mtk.SE2(a.pos + b(1:2), a.phi + b(3));
            end
        end
        
        function r = minus(a, b)
            if ~(isa(a, 'mtk.SE2') && isa(b, 'mtk.SE2'))
                error('undefined')
            else
                r = [a.pos - b.pos;
                     mtk.util.normalize_angle(a.phi - b.phi)];
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