function [err,ss,tt,w,sdim,eigval,info] = mjdgges(e,d,qz_criterium, fake)
%function [err,ss,tt,w,sdim,eigval,info] = mjdgges(e,d,qz_criterium)
% QZ decomposition, Sims' codes are used.
%
% INPUTS
%   e            [double] real square (n*n) matrix.
%   d            [double] real square (n*n) matrix.
%   qz_criterium [double] scalar (1+epsilon).
%
% OUTPUTS
%   err          [double]  scalar: 1 indicates failure, 0 indicates success
%   ss           [complex] (n*n) matrix.
%   tt           [complex] (n*n) matrix.
%   w            [complex] (n*n) matrix.
%   sdim         [integer] scalar.
%   eigval       [complex] (n*1) vector.
%   info         [integer] scalar.
%
% ALGORITHM
%   Sims's qzdiv routine is used.
%
% SPECIAL REQUIREMENTS
%   none.

% Copyright (C) 1996-2017 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

% Check number of inputs and outputs.
if nargin>4 || nargin<2 || nargout>7 || nargout==0
    error('MJDGGES: takes 2 or 3 input arguments and between 1 and 7 output arguments.')
end

% Check the first two inputs.
[me,ne] = size(e);
[md,nd] = size(d);
if ( ~isreal(e) || ~isreal(d) || me~=ne || md~=nd || me~=nd)
    % info should be negative in this case, see dgges.f.
    error('MJDGGES requires two square real matrices of the same dimension.')
end

% Set default value of qz_criterium.
if nargin <3
    qz_criterium = 1 + 1e-6;
end

info = 0;

% Initialization of the output arguments.
ss = zeros(ne,ne);
tt = zeros(ne,ne);
w  = zeros(ne,ne);
sdim   = 0;
eigval = zeros(ne,1);

% Computational part.
try
    if isoctave()
        % Force complex QZ factorization.
        e = complex(e);
        d = complex(d);
    end
    [ss,tt,qq,w,a,b,c] = qz(e, d);
    warning_old_state = warning;
    warning off;
    [tt,ss,qq,w] = qzdiv(qz_criterium, tt, ss, qq, w);
    eigval = diag(ss)./diag(tt);
    warning(warning_old_state);
    sdim = sum(abs(eigval) < qz_criterium);
catch
    info = 1;% Not as precise as lapack's info!
end

err = 0;