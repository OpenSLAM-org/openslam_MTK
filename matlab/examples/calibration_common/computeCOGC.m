% COMPUTECOGC

%  Copyright (c) 2011 DFKI GmbH
%  All rights reserved
%
%  Author: Oliver Birbach <oliver.birbach@dfki.de>
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

function [x, y, rms] = computeCOGC(x0, y0, threshold, radius, sigma, image, width, height)
    rX = radius * sigma;
    rY = rX;
    
    if x0 - rX < 1
        rX = x0 - 1;
    end;
    
    if x0 + rX > width - 2
        rX = width - 2 - x0;
    end;
    
    if y0 - rY < 1
        rY = y0 -1;
    end;
    if y0 + rY > height - 2
        rY = height - 2 - y0;
    end;
    
    sumggt00=0;
    sumggt01=0;
    sumggt11=0; % symmetric 2*2 matrix ggT
    
    sumggtx0=0;
    sumggtx1=0; % ggTx vector
  
    s2 =2*sigma*sigma;
    eps = 1;
    w=0;
    delta = 1;
  
    for xx=-rX:delta:rX
        for yy=-rY:delta:rY
            weight = exp(-(xx*xx+yy*yy)/s2);
            gradX = (interpolated(x0+xx+eps, y0+yy, image) - interpolated(x0+xx-eps, y0+yy, image))/(2*eps);
            gradY = (interpolated(x0+xx, y0+yy+eps, image) - interpolated(x0+xx, y0+yy-eps, image))/(2*eps);

            sumggt00 = sumggt00 + weight * gradX * gradX;
            sumggt01 = sumggt01 + weight * gradX * gradY;
            sumggt11 = sumggt11 + weight * gradY * gradY;

            sumggtx0 = sumggtx0 + weight*(xx*gradX+yy*gradY)*gradX;
            sumggtx1 = sumggtx1 + weight*(xx*gradX+yy*gradY)*gradY;
            
            w = w + weight;
        end;
    end;

    p = - sumggt00-sumggt11;
    q = sumggt00*sumggt11 - sumggt01*sumggt01;
    p = double(p);
    q = double(q);
    lambda0 = sqrt(-p/2-sqrt(p*p/4-q));
    rms = lambda0;

    % 'lambda0' is the sqrt of the smaller eigenvalue of Integral('ggt')
    % thus defining contrast in the worst direction
    if lambda0 < threshold
        x = x0;
        y = y0;
    end;

  % Now compute I(ggt)^-1*I(ggtx) by computing the inverse of sumggt??
  det  = sumggt00*sumggt11 - sumggt01*sumggt01;
  xSol = (+sumggt11*sumggtx0 - sumggt01*sumggtx1) / det;
  ySol = (-sumggt01*sumggtx0 + sumggt00*sumggtx1) / det;
  x = x0+xSol;
  y = y0+ySol;