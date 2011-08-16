% SO2 orientations in the 2D plane    

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

classdef SO2 < handle
    properties (Constant)
        DOF = 1;
    end
    properties
        phi = 0;
    end
    
    methods
        function obj = SO2(p)
            if nargin > 0
                if isa(p, 'mtk.SO2')
                    obj.copyfrom(p);
                else
                    obj.phi = p;
                end
            end
        end
        
        function copyfrom(obj, other)
            obj.phi = other.phi;
        end
        
        function set.phi(obj, value)
            obj.phi = mtk.util.normalize_angle(value);
        end
        
        function r = plus(a, b)
            if ~(isa(a, 'mtk.SO2') && isa(b, 'numeric'))
                error('undefined')
            else
                r = mtk.SO2(a.phi + b);
            end
        end
        
        function r = minus(a, b)
            if ~(isa(a, 'mtk.SO2') && isa(b, 'mtk.SO2'))
                error('undefined')
            else
                r = mtk.util.normalize_angle(a.phi - b.phi);
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
