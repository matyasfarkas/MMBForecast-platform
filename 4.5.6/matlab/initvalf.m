function initvalf(fname_)
% function initvalf(fname_)
%
% Reads an initial path from the 'fname_' file for exogenous and endogenous variables
%
% INPUTS
%    fname_:         name of the function or file containing the data
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    All variables local to this function have an underscore appended to
%    their name, to minimize clashes with model variables loaded by this function.

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

global M_ oo_ options_

series_ = 1;

[directory,basename,extension] = fileparts(fname_);

% Auto-detect extension if not provided
if isempty(extension)
    if exist([basename '.m'],'file')
        extension = '.m';
    elseif exist([basename '.mat'],'file')
        extension = '.mat';
    elseif exist([basename '.xls'],'file')
        extension = '.xls';
    elseif exist([basename '.xlsx'],'file')
        extension = '.xlsx';
    else
        error(['Can''t find datafile: ' basename '.{m,mat,xls,xlsx}']);
    end
end

fullname = [basename extension];

if ~exist(fullname)
    error(['Can''t find datafile: ' fullname ]);
end

switch (extension)
  case '.m'
    eval(basename);
  case '.mat'
    load(basename);
  case { '.xls', '.xlsx' }
    [data_,names_v_]=xlsread(fullname); % Octave needs the extension explicitly
    series_=0;
  otherwise
    error(['Unsupported extension for datafile: ' extension])
end

options_.initval_file = 1;
oo_.endo_simul = [];
oo_.exo_simul = [];

for i_=1:size(M_.endo_names,1)
    if series_ == 1
        x_ = eval(M_.endo_names(i_,:));
        if size(x_,2)>size(x_,1) %oo_.endo_simul must be collection of row vectors
            oo_.endo_simul = [oo_.endo_simul; x_];
        else %transpose if column vector
            oo_.endo_simul = [oo_.endo_simul; x_'];
        end
    else
        k_ = strmatch(deblank(M_.endo_names(i_,:)),names_v_,'exact');
        if isempty(k_)
            error(['INITVAL_FILE: ' deblank(M_.endo_names(i_,:)) ' not found'])
        end
        x_ = data_(:,k_);
        oo_.endo_simul = [oo_.endo_simul; x_'];
    end
end

for i_=1:size(M_.exo_names,1)
    if series_ == 1
        x_ = eval(M_.exo_names(i_,:) );
        if size(x_,2)>size(x_,1) %oo_.endo_simul must be collection of row vectors
            oo_.exo_simul = [oo_.exo_simul x_'];
        else %if column vector
            oo_.exo_simul = [oo_.exo_simul x_];
        end
    else
        k_ = strmatch(deblank(M_.exo_names(i_,:)),names_v_,'exact');
        if isempty(k_)
            error(['INITVAL_FILE: ' deblank(M_.exo_names(i_,:)) ' not found'])
        end
        x_ = data_(:,k_);
        oo_.exo_simul = [oo_.exo_simul x_];
    end
end