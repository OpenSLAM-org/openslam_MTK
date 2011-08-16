% GN Gauss-Newton algorithm
%   This function finds xm = argmin{f(x)} where x is an nd-vector and
%   f(x) = sum(r(x)^2). The residual vector r and its Jacobian J have to
%   be provided by fun for any x. x0 is the initial guess.

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

function [ x chisq alpha ] = gn( fun, x0, params, plot )

done = 0;
x = x0;

if nargin > 2
    ITMAX = params(1);
    gain_converged = params(2);
else
    ITMAX = 25;
    gain_converged = 1e-5;
end

if nargin > 3
    plot(x);
end

tStart = tic;

% Initialization
[r J] = feval(fun,x);
chisq = sum(r.^2);

% Create matrix and vector
DOF = size(J,2);
alpha = sparse(DOF, DOF);

for iter=1:ITMAX
    if nargin > 3
        plot(x);
    end
    
    % Print status
    fprintf('%d       \t%12.6e\t', iter, full(chisq(1,1)));
    toc(tStart);
    tStart = tic;

    % solve linear system
    alpha = J' * J;
    beta = J' * r;
    da = alpha \ -beta; 
    %da = cs_cholsol(alpha, full(-beta));
   
    x = x + da;
    
    [r J] = feval(fun, x);
    
    ochisq = chisq;
    chisq = sum(r.^2);
    
    gain = (ochisq - chisq)/chisq;
    if gain > 0 && gain < gain_converged
        % converged
        return;
    end
end;
end
