function CutSample(M_, options_, estim_params_)
% function CutSample(M_, options_, estim_params_)
% Takes a subset from metropolis draws by storing the required information
% like the first MH-file to be loaded and the first line in that file to be
% loaded into the record structure saved on harddisk into the
% _mh_history-file
%
% INPUTS
%   M_               [structure]    Dynare model structure
%   options_         [structure]    Dynare options structure
%   estim_params_    [structure]    Parameter structure
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2005-2017 Dynare Team
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

npar = estim_params_.np+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.nvx;

% Get the path to the metropolis files.
MetropolisFolder = CheckPath('metropolis',M_.dname);

% Get the (base) name of the mod file.
ModelName = M_.fname;

% Load the last mh-history file.
load_last_mh_history_file(MetropolisFolder, M_.fname);

% Get the list of files where the mcmc draw are saved.
mh_files = dir([ MetropolisFolder ,filesep, M_.fname '_mh*.mat' ]);

if ~length(mh_files)
    error('Estimation::mcmc: I can''t find MH file to load here!')
end

TotalNumberOfMhFiles = sum(record.MhDraws(:,2));
TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);
FirstDraw = max(1,floor(options_.mh_drop*TotalNumberOfMhDraws)+1);
FirstMhFile = ceil(FirstDraw/MAX_nruns);
FirstLine = FirstDraw-(FirstMhFile-1)*MAX_nruns;
record.KeepedDraws.FirstMhFile = FirstMhFile;
record.KeepedDraws.FirstLine = FirstLine;
if (TotalNumberOfMhFiles-1)-(FirstMhFile+1)+1 > 0
    record.KeepedDraws.Distribution = [ MAX_nruns-FirstLine+1 ; ...
                        ones((TotalNumberOfMhFiles-1)-(FirstMhFile+1)+1,1)*MAX_nruns ; ...
                        record.MhDraws(end,3) ];
elseif TotalNumberOfMhFiles == 1
    record.KeepedDraws.Distribution = [];
elseif TotalNumberOfMhFiles == 2 && FirstMhFile > 1
    record.KeepedDraws.Distribution = [MAX_nruns-FirstLine+1 ; record.MhDraws(end,3)];
end

% Save updated mh-history file.
update_last_mh_history_file(MetropolisFolder, ModelName, record);

fprintf('Estimation::mcmc: Total number of MH draws per chain: %d.\n',TotalNumberOfMhDraws);
fprintf('Estimation::mcmc: Total number of generated MH files: %d.\n',TotalNumberOfMhFiles);
fprintf('Estimation::mcmc: I''ll use mh-files %d to %d.\n',FirstMhFile,TotalNumberOfMhFiles);
fprintf('Estimation::mcmc: In MH-file number %d I''ll start at line %d.\n',FirstMhFile,FirstLine);
fprintf('Estimation::mcmc: Finally I keep %d draws per chain.\n',TotalNumberOfMhDraws-FirstDraw+1);
skipline()