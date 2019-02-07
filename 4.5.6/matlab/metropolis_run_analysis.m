function metropolis_run_analysis(M,basetopt,j)
%function metropolis_run_analysis(M)
% analizes Metropolis runs
%
% INPUTS
%   M:         (struct)  Model structure
%   basetopt:  (struct)  Estimated parameters structure
%   j:         (int)     Index of estimated paramter
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2003-2017 Dynare Team
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

load([M.fname '/metropolis/' M.fname '_mh_history'])
nblck = record.Nblck;
ndraws = sum(record.MhDraws(:,1));

logPost = [];
params = [];
blck = 1;
for i=1:record.LastFileNumber
    fname = [M.fname '/metropolis/' M.fname '_mh' int2str(i) '_blck' ...
             int2str(blck) '.mat'];
    if exist(fname,'file')
        o=load(fname);
        logPost = [logPost; o.logpo2];
        params  = [params; o.x2];
    end
end

figure;
subplot(2,1,1)
plot(logPost)
subplot(2,1,2)
plot(params(:,j))
title(['parameter ' int2str(j)])