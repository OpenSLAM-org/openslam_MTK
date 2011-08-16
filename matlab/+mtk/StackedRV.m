% STACKEDRV A random variable consisting of stacked manifolds
%   Intended for internal use only.

%  Copyright (c) 2011, DFKI GmbH
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

classdef StackedRV < handle
    
    properties (Access = public)
        DOF = 0
        rvs = {}
        rvs_size = 0
        rvs_capacity = 0
    end
    
    methods
        function obj = StackedRV()
            obj.DOF = 0;
            obj.rvs = {};
            obj.rvs_size = 0;
            obj.rvs_capacity = 0;
        end
        
        function reserve(obj, n)
            if obj.rvs_capacity < n
                old = obj.rvs;
                obj.rvs = cell(n,1);
                obj.rvs{n,1} = [];
                obj.rvs{1:obj.rvs_size} = old{1:obj.rvs_size};
                obj.rvs_capacity = n;
            end
        end
        
        function handle = append_random_var(obj, x)
            %d = obj.DOF
            %rs = obj.rvs_size
            obj.DOF = obj.DOF + length(x);
            obj.rvs_size = obj.rvs_size + 1;
            if obj.rvs_size > obj.rvs_capacity
                %siz = obj.rvs_size
                %cap = obj.rvs_capacity
                obj.rvs = [obj.rvs; {x}];
                obj.rvs_capacity = obj.rvs_capacity + 1;
            else
                obj.rvs{obj.rvs_size} = x;
            end
            handle = obj.rvs_size;
        end
        
        function handles = append_random_vars(obj, xs)
            start = obj.rvs_size + 1;
            for i=1:length(xs)
                obj.DOF = obj.DOF + length(xs{i});
            end
            % FIXME: allow for mixing append_random_var and _vars
            obj.rvs = [obj.rvs; xs];
            obj.rvs_size = size(obj.rvs, 1);
            obj.rvs_capacity = obj.rvs_size;
            handles = start:obj.rvs_size;
        end
        
        function rv = get_random_var(obj, handle)
            rv = obj.rvs{handle};
        end
    
        function r = plus(a, b)
            if ~(isa(a, 'mtk.StackedRV') && isa(b, 'numeric'))
                error('undefined')
            elseif length(b) ~= a.DOF
                error('size mismatch')
            else
                r = mtk.StackedRV();
                r.DOF = a.DOF;
                r_rvs = cell(size(a.rvs));
                r.rvs_size = a.rvs_size;
                r.rvs_capacity = a.rvs_capacity;
                a_rvs = a.rvs;
                i = 1;
                for k=1:a.rvs_size
                    xk = a_rvs{k};
                    xk_dof = length(xk);
                    r_rvs{k} = xk + b(i:i+xk_dof-1);
                    i = i + xk_dof;
                end
                r.rvs = r_rvs;
            end
        end
        
        function r = minus(a, b)
            if ~(isa(a, 'mtk.StackedRV') && isa(b, 'mtk.StackedRV'))
                error('undefined')
            elseif a.DOF ~= b.DOF
                error('size mismatch')
            else
                r = zeros(a.DOF,1);
                i = 1;
                for k=1:a.rvs_size
                    p = a.rvs{k};
                    q = b.rvs{k};
                    p_dof = length(p);
                    if p_dof ~= length(q)
                        error('size mismatch')
                    end
                    r(i:i+p_dof-1) = p - q;
                    i = i + p_dof;
                end
            end
        end
        
        function s = scale(obj)
            s = zeros(obj.DOF,1);
            i = 1;
            for k=1:obj.rvs_size
                p = obj.rvs{k};
                dof = length(p);
                s(i:i+dof-1) = mtk.scale(p);
                i = i + dof;
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