function [mhname,info] = get_name_of_the_last_mh_file(M_)
% This function returns the name of the last mh file and test if the metropolis was completed.
%
% INPUTS
%   M_       [structure]   Dynare structure specifying the model.
%
% OUTPUTS
%  mhname    [string]      Name of the last mh file (with complete path).
%  info      [integer]     Scalar. If info is equal to 1 then the predicted name of the last
%                          metropolis hastings matches the name of the name of the last mh
%                          file. Otherwise info is equal to zero (a likely reason for this is
%                          that the mcmc simulations were not completed).

% Copyright (C) 2008-2017 Dynare Team
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

mhname = [];
info = 1;

MetropolisFolder = CheckPath('metropolis',M_.dname);
ModelName = M_.fname;
BaseName = [MetropolisFolder filesep ModelName];

load_last_mh_history_file(MetropolisFolder, ModelName);

mh_number = record.LastFileNumber ;
bk_number = record.Nblck ;

predicted_mhname = [ BaseName '_mh' int2str(mh_number) '_blck' int2str(bk_number) '.mat' ] ;

all_mh_files = dir([BaseName '_mh*_blck*' ]);
[junk,idx] = sort([all_mh_files.datenum]);
mhname = all_mh_files(idx(end)).name;

if ~strcmpi(mhname,predicted_mhname)
    info = 0;
end