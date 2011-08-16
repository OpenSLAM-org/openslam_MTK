% OPTIMIZATIONPROBLEM A graph optimization problem.
    
%   You can collect random variables and measurements using this
%   class. A stacked random variable x0 and the "big" error function fun
%   required by the optimization algorithm is generated automatically
%   from these.
    
%  Copyright (c) 2011 DFKI GmbH
%  All rights reserved
%
%  Author: Rene Wagner <rene.wagner@dfki.de>
%          Oliver Birbach <oliver.birbach@dfki.de>
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

classdef OptimizationProblem < handle
    
    properties (Access = private)
        totaldof = 0
        x_stacked_idx = []
        x_stacked_idx_current = 1
        J_nonzero = 0
    end
    
    properties
        X0 = mtk.StackedRV()
        measurements = {}
    end
    
    methods
        function [r J r_orig] = fun(obj, x)
            % pre-allocate r based on total DOF
            r = zeros(obj.totaldof, 1);
            if nargout > 2
                r_orig = zeros(obj.totaldof, 1);
            end

            % sparse numerical approximation of J
            m = obj.totaldof;
            n = length(x);
            
            J_ = zeros(obj.J_nonzero, 3);
            
            eps_ = sqrt(eps);
            eps_2 = 0.5 * eps_;
            e = eye(n);

            q = 1;
            i = 1;
            % for each measurement
            for k=1:length(obj.measurements)
                m = obj.measurements{k};
                z = m{1};
                f = m{2};
                rvhandles = m{3};
                userdata = m{4};
                
                fixed_sigma = length(m) > 4;
                
                rvs = cellfun(@(h) x.rvs{h}, rvhandles, 'UniformOutput', false);
                if fixed_sigma
                    z_ = f(rvs{:}, userdata{:});
                else
                    [z_ sigma] = f(rvs{:}, userdata{:});
                end
                
                dof = length(z);
               
                rk_orig = z_ - z;
                
                if fixed_sigma
                    rk = mtk.normalize_innovation(rk_orig, m{5}, fixed_sigma);
                else
                    rk = mtk.normalize_innovation(rk_orig, sigma);
                end
                
                % for each dependent random variable
                for rv=1:length(rvhandles)
                    nstart = obj.x_stacked_idx(rvhandles{rv});
                    rv_dof = length(rvs{rv}); 
                    for j=1:rv_dof
                        ej = e(1:rv_dof,j);
                        %x1 = x + ((-eps_2) * ej);
                        %x2 = x + (eps_2 * ej);
                        %J(:,j) = (obj.f(x2) - obj.f(x1)) / eps_;
                        
                        rvs2 = {rvs{:}};                       
                        rvs2{rv} = rvs{rv} + eps_ * ej;

                        if fixed_sigma
                            z2_ = f(rvs2{:}, userdata{:});
                            rk2 = mtk.normalize_innovation(z2_ - z, m{5}, fixed_sigma);
                        else
                            [z2_ sigma2] = f(rvs2{:}, userdata{:});
                            rk2 = mtk.normalize_innovation(z2_ - z, sigma2);
                        end
                        
                        d = (rk2 - rk) / eps_;
                        
                        % NOTE: Given the spconvert call below, the following is the equivalent of
                        %       J(i:i+dof-1, nstart+j-1) = d;

                        J_(q:q+dof-1,:) = [(i:i+dof-1)', (nstart+j-1) * ones(dof,1), d];
                        
                        q = q + dof;
                    end
                end

                r(i:i+dof-1,1) = rk;
                if nargout > 2
                    r_orig(i:i+dof-1,1) = rk_orig;
                end
                
                i = i + dof;
            end
            
            % "sparsify"
            r = sparse(r);
            J = spconvert(J_);

            % Check if the number of non-zero elements we know we have
            % computed is always less or equal the number of non-zero
            % elements in the converted sparse Jacobian.
            if (nnz(J) > obj.J_nonzero)
                error(['The number of non-zero elements in the Jacobian is ' ...
                'larger than expected.']);
            end;
        end
        
        function reserve_random_vars(obj, n)
            obj.X0.reserve(n);
        end
        
        function handle = add_random_var(obj, x)
            obj.x_stacked_idx = [obj.x_stacked_idx obj.x_stacked_idx_current]; % FIXME: move to StackedRV ?
            obj.x_stacked_idx_current = obj.x_stacked_idx_current + length(x);
            handle = obj.X0.append_random_var(x);
        end
        
        function handles = add_random_vars(obj, xs)
            for i=1:length(xs)
                % FIXME: move to StackedRV ?
                obj.x_stacked_idx = [obj.x_stacked_idx obj.x_stacked_idx_current];
                obj.x_stacked_idx_current = obj.x_stacked_idx_current + length(xs{i});
            end
            handles = obj.X0.append_random_vars(xs);
        end
        
        function add_measurement(obj, z, f, rvhandles, userdata, sigma)
        % where [^z, sigma] = f(rvs{:}, userdata{:})
            if nargin > 5
                if size(sigma, 1) > 1
                    sigma = inv(chol(sigma, 'lower'));
                 end
                obj.measurements = [obj.measurements; {{z f rvhandles userdata sigma}}];
            else
                obj.measurements = [obj.measurements; {{z f rvhandles userdata}}];
            end
            % store DOF of measurement
            z_dof = length(z);
            obj.totaldof = obj.totaldof + z_dof;
        
            x = obj.X0;
            rvs = cellfun(@(h) x.rvs{h}, rvhandles, 'UniformOutput', false);
            rvs_totaldof = 0;
            for j=1:length(rvs)
                rvs_totaldof = rvs_totaldof + length(rvs{j});
            end
            obj.J_nonzero = obj.J_nonzero + z_dof * rvs_totaldof;
        end
        
        function add_measurements(obj, ms)
            x = obj.X0;
            added_dof = 0;
            added_J_nonzero = 0;
            for i=1:length(ms)
                % store DOF of measurement
                m = ms{i};
                z = m{1};
                z_dof = length(z);
                added_dof = added_dof + z_dof;
                
                % pre-compute Cholesky decomposition if needed
                if length(m) > 4
                    s = m{5};
                    if size(s, 1) > 1
                        ms{i}{5} = inv(chol(s, 'lower'));
                    end
                end

                % store number of non-zero Jacobian entries
                rvhandles = m{3};
                rvs = cellfun(@(h) x.rvs{h}, rvhandles, 'UniformOutput', false);
                rvs_totaldof = 0;
                for j=1:length(rvs)
                    rvs_totaldof = rvs_totaldof + length(rvs{j});
                end
                added_J_nonzero = added_J_nonzero + z_dof * rvs_totaldof;
            end
            obj.totaldof = obj.totaldof + added_dof;
            obj.J_nonzero = obj.J_nonzero + added_J_nonzero;
            obj.measurements = [obj.measurements; ms];
        end
    end
end


