% SO3 orientations in 3D
%   The vectorized space is formed by relative rotations in world
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

classdef SO3 < handle
    
    properties (Constant)
        DOF = 3;
    end
    properties
        Q = diag(ones(3,1), 0);
    end
    
    methods
        function obj = SO3(p)
            if nargin > 0
                if isa(p, 'mtk.SO3')
                    obj.copyfrom(p);
                else
                    obj.Q = p;
                end
            end
        end
        
        function copyfrom(obj, other)
            obj.Q = other.Q;
        end
        
        function r = plus(s, d)
            if ~(isa(s, 'mtk.SO3') && isa(d, 'numeric'))
                error('undefined')
            else
                % A\B is roughly the same as inv(A)*B
                r = mtk.SO3(s.Q * mtk.util.rot(s.Q \ d));
            end
        end
        
        function r = minus(s2, s1)
            if ~(isa(s2, 'mtk.SO3') && isa(s1, 'mtk.SO3'))
                error('undefined')
            else
                % A\B is roughly the same as inv(A)*B
                r = s1.Q * mtk.util.arot(s1.Q \ s2.Q);
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