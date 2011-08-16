%NRLM Levenberg-Marquardt algorithm
%   Levenberg Marquardt implementation as suggested by W.H. Press, S.A.
%   Teukolsky, W.T. Vetterling, B.P. Flannery; Numerical Recipes: The Art
%   of Scientific Computing, 3rd edition, Cambridge University Press, 2007
%   
%   This function finds xm = argmin{f(x)} where x is an nd-vector and
%   f(x) = sum(r(x)^2). The residual vector r and its Jacobian J have to
%   be provided by fun for any x. x0 is the initial guess.
%

%  Copyright (c) 2011 DFKI GmbH
%  All rights reserved
%
%  Author: Oliver Birbach <oliver.birbach@dfki.de>
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

function [ x chisq alpha ] = nrlm( fun, x0, params, plot )

done = 0;
alamda = .001;
x = x0;

if nargin > 2
    ITMAX = params(1);
    NDONE = params(2);
    tol = params(3);
else
    ITMAX = 25;
    NDONE = 4;
    tol = 1.e-3;
end

if nargin > 3
    plot(x);
end

tStart = tic;

% Initialization
[r J] = feval(fun,x);
chisq = sum(r.^2);
ochisq = chisq;

% Create matrix and vector
DOF = size(J,2);
alpha = sparse(DOF, DOF);

for iter=1:ITMAX
    if nargin > 3
        plot(x);
    end
    
    % Print status
    fprintf('%d       \t%12.6e\t', iter, full(ochisq(1,1)));
    toc(tStart);
    tStart = tic;
    
    % Last pass. Use zero alamda
    if (done == NDONE)
        alamda = 0.;
    end;
    
    % Alter linearized fitting matrix, by 
    % augmenting diagonal elements
    alpha = J' * J + alamda * sparse(eye(DOF));
    beta = J' * r;
    
    % Matrix solution (CHOLMOD)
    da = alpha \ -beta; 
    
    % alternative (CSparse)
    %da = cs_cholsol(alpha, full(-beta));
   
    % Converged. Return
    if (done == NDONE)
        return;
    end;
    
    % Did the trial succeed?
    x_new = x + da;
    
    [r_new J_new] = feval(fun, x_new);
    chisq = sum(r_new.^2);
    if (abs(chisq - ochisq) < max(tol, tol * chisq))
        done = done + 1;
        if nargout < 3 && done == NDONE
            % Converged. alpha not requested. Update x and return
            if chisq < ochisq
                x = x_new;
            end
            return;
        end
    end;
    
    if (chisq < ochisq) % Success, accept new solution
        alamda = alamda * .1;%sqrt(.1);
        ochisq = chisq;
        r = r_new;
        J = J_new;
        x = x_new;
    else % Failure, increase lambda
        alamda = alamda * 10.;%sqrt(10.);
        chisq = ochisq;
    end;
end;
end
