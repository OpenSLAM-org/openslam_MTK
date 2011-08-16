% SCALE Determines the scale of a manifold x.
%   For vectors this is abs(x). For other manifolds it gives an idea of
%   the magnitude of meaningful changes.

%  Copyright (c) 2011 DFKI GmbH
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

function [ s ] = scale( x )

if isa(x, 'mtk.StackedRV') || isa(x, 'mtk.Compound')
    s = x.scale();
elseif isa(x, 'mtk.Rn')
    s = abs(x.vec);
elseif isa(x, 'mtk.SO2') || isa(x, 'mtk.SO3')
    dof = size(x, 1);
    s = pi / (2. * sqrt(dof)) * ones(dof, 1);
elseif isa(x, 'mtk.SE2')
    s = [abs(x.pos); pi / 2.];
elseif isa(x, 'mtk.SE3')
    s = [abs(x.pos); pi / (2. * sqrt(3)) * ones(3, 1)];
elseif isa(x, 'numeric')
    s = abs(x);
else
    error('undefined');
end
