% Rn R^n vectors

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

classdef Rn < handle
    properties (SetAccess = private)
        DOF
    end
    properties
        vec
    end
    
    methods
        function obj = Rn(p, dim)
            if isa(p, 'mtk.Rn')
                obj.vec = p.vec;
            elseif nargin > 1
                obj.vec = zeros(dim, 1);
            else
                obj.vec = p;
            end
            obj.DOF = size(obj.vec,1);
        end
        
        function copyfrom(obj, other)
            obj.vec = other.vec;
        end
        
        function set.vec(obj, value)
            if ~isempty(obj.vec) && any(size(value) ~= size(obj.vec))
                error('size mismatch')
            else
                obj.vec = value;
            end
        end
        
        function r = plus(a, b)
            if ~(isa(a, 'mtk.Rn') && isa(b, 'numeric'))
                error('undefined')
            else
                r = mtk.Rn(a.vec + b);
            end
        end
        
        function r = minus(a, b)
            if ~(isa(a, 'mtk.Rn') && isa(b, 'mtk.Rn'))
                error('undefined')
            else
                r = a.vec - b.vec;
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