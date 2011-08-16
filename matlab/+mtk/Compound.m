% COMPOUND A compound manifold consisting of a sequence of sub-manifolds.
%   Use the mtk.make_compound helper function instead.

%  Copyright (c) 2010, DFKI GmbH
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

classdef Compound < hgsetget & dynamicprops

    properties (SetAccess = private)
        DOF
    end
    
    properties (Access = private)
        typespec
        props = {}
        propnames = {}
        idxdict
        epropnames = {}
    end
    
    methods
        function obj = Compound(p)
            copy = isa(p, 'mtk.Compound');
            if copy
                typespec = p.typespec;
            else
                typespec = p;
            end
            obj.typespec = typespec;
            nargs = size(typespec, 2);
            state = 0;
            k = 1;
            while k <= nargs
                switch state
                    case 0
                        if ~isa(typespec{k}, 'char')
                            %cp = copy
                            %spec = typespec
                            error('expected property name in typespec')
                        elseif strcmp(typespec{k}, 'extra')
                            state = 2;
                        else
                            prop = typespec{k};
                            obj.propnames = [obj.propnames {prop}];
                            state = 1;
                        end
                    case 1
                        cl = typespec{k};
                        addprop(obj, prop);
                        if copy
                            o = cl(get(p, prop));
                            if k < nargs && isa(typespec{k+1}, 'numeric')
                                k = k+1;
                            end
                        else
                            if k < nargs && isa(typespec{k+1}, 'numeric')
                                dof = typespec{k+1};
                                o = cl('dof', dof);
                                k = k+1;
                            else
                                o = cl();
                            end
                        end
                        set(obj, prop, o);
                        obj.props = [obj.props {o}];
                        state = 0;
                    case 2
                        if ~isa(typespec{k}, 'char')
                            error(['expected extra property name in ' ...
                                   'typespec'])
                        else
                            if  k >= nargs
                                error(['expected initial value of ' ...
                                       'extra property'])
                            end
                            eprop = typespec{k};
                            if isa(typespec{k+1}, 'numeric')
                                val = typespec{k+1};
                            elseif isa(typespec{k+1}, 'function_handle')
                                f = typespec{k+1};
                                val = f();
                            else
                                % deep copy
                                val = feval(class(typespec{k+1}));
                            end
                            addprop(obj, eprop);
                            set(obj, eprop, val);
                            obj.epropnames = [obj.epropnames {eprop}];
                        end
                        k = k+1;
                end
                k = k + 1;
            end
            nprops = length(obj.props);
            propidx = cell(nprops, 1);
            % start at idx 0 to make recusion work and fix it up in
            % non-recursive wrapper
            v = 0;
            for k=1:length(obj.props)
                propidx{k} = v;
                p = obj.props{k};
                v = v + p.DOF;
            end
            dict = reshape({obj.propnames{:};propidx{:}},2,[]);
            obj.idxdict = struct(dict{:});
            
            obj.DOF = v;
        end
        
        function copyfrom(obj, other)
        % copies property values from other object
            if ~isa(other, 'mtk.Compound') %|| other.typespec ~= obj.typespec
                error('type mismatch')
            end

            % copy sub-manifold properties
            for k=1:size(obj.props, 2)
                obj.props{k}.copyfrom(other.props{k});
            end
            
            % copy extra properties
            for k=1:length(other.epropnames)
                p = other.epropnames{k};
                val = get(other, p);
                
                if isa(val, 'numeric')
                    set(obj, p, val);
                else
                    % attempt deep copy
                    % FIXME: requiring copyfrom to be implemented
                    %        is not pretty...
                    val2 = get(obj, p);
                    val2.copyfrom(val);
                    set(obj, p, val2);
                end
            end
        end
        
        function r = plus(a, b)
            if ~(isa(a, 'mtk.Compound') && isa(b, 'numeric'))
                error('undefined')
            elseif length(b) ~= a.DOF
                error('size mismatch')
            else
                r = mtk.Compound(a);
                i = 1;
                for k=1:size(a.props, 2)
                    p = a.props{k};
                    r.props{k}.copyfrom(p + b(i:i+p.DOF-1));
                    i = i + p.DOF;
                end
            end
        end
        
        function r = minus(a, b)
            if ~(isa(a, 'mtk.Compound') && isa(b, 'mtk.Compound'))
                error('undefined')
            elseif a.DOF ~= b.DOF
                error('size mismatch')
            else
                r = zeros(a.DOF,1);
                i = 1;
                for k=1:size(a.props, 2)
                    p = a.props{k};
                    q = b.props{k};
                    dof = p.DOF;
                    if dof ~= q.DOF
                        error('size mismatch')
                    end
                    r(i:i+dof-1) = p - q;
                    i = i + dof;
                end
            end
        end
        
        function idx = idxOf(obj, varargin)
            idx = 1 + obj.idxOfRec(varargin);
        end
        
        function display(obj)
            disp('');
            disp([inputname(1),' = ']);
            disp('  ');
            disp('  Compound manifold handle');
            disp('  ');
            disp(['  DOF: ' num2str(obj.DOF)]);
            if ~isempty(obj.props)
                disp('  ');
                disp('  Sub-manifold properties:');
                for k=1:length(obj.props)
                    p = obj.props{k};
                    fprintf('    %s: %s (%d DOF)\n', obj.propnames{k}, class(p), p.DOF);
                end
            end
            if ~isempty(obj.epropnames)
                disp('  ');
                disp('  Extra properties:');
                for k=1:length(obj.epropnames)
                    p = obj.epropnames{k};
                    fprintf('    %s: ', p);
                    disp(get(obj, p));
                end
            end
        end
        
        function s = scale(obj)
            s = zeros(obj.DOF, 1);
            i = 1;
            for k=1:size(obj.props, 2)
                p = obj.props{k};
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

    methods (Access = private)
        function idx = idxOfRec(obj, proppath)
            propname = proppath{1};
            prop = get(obj, propname);
            if length(proppath) < 2
                idx = obj.idxdict.(propname);
            else
                idx = obj.idxdict.(propname) + prop.idxOfRec(proppath(2:length(proppath)));
            end
        end
    end
end