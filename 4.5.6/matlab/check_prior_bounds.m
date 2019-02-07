function check_prior_bounds(xparam1,bounds,M_,estim_params_,options_,bayestopt_)
% function check_prior_bounds(xparam1,bounds,M_,estim_params_,options_)
% checks the parameter vector of violations of the prior bounds
% Inputs:
%   -xparam1        [double]    vector of parameters to be estimated (initial values)
%   -bounds         [vector]    vector containing the lower and upper
%   bounds
%   -M_             [structure] characterizing the model.
%   -estim_params_  [structure] characterizing parameters to be estimated
%   -options_       [structure] characterizing the options
%   -bayestopt_     [structure] characterizing priors

% Copyright (C) 2013-2017 Dynare Team
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

outside_bound_pars=find(xparam1 < bounds.lb | xparam1 > bounds.ub);
if ~isempty(outside_bound_pars)
    for ii=1:length(outside_bound_pars)
        outside_bound_par_names{ii,1}=get_the_name(outside_bound_pars(ii),0,M_,estim_params_,options_);
    end
    disp_string=[outside_bound_par_names{1,:}];
    for ii=2:size(outside_bound_par_names,1)
        disp_string=[disp_string,', ',outside_bound_par_names{ii,:}];
    end
    error(['Initial value(s) of ', disp_string ,' are outside parameter bounds. Potentially, you should set prior_trunc=0. If you used the mode_file-option, check whether your mode-file is consistent with the priors.'])
end
inadmissible_inverse_gamma_values=find(bayestopt_.pshape==4 & xparam1 == 0);
if ~isempty(inadmissible_inverse_gamma_values)
    for ii=1:length(inadmissible_inverse_gamma_values)
        inadmissible_inverse_gamma_par_names{ii,1}=get_the_name(inadmissible_inverse_gamma_values(ii),0,M_,estim_params_,options_);
    end
    disp_string=[inadmissible_inverse_gamma_par_names{1,:}];
    for ii=2:size(inadmissible_inverse_gamma_par_names,1)
        disp_string=[disp_string,', ',inadmissible_inverse_gamma_par_names{ii,:}];
    end
    error(['Initial value(s) of ', disp_string ,' is zero. This is not allowed when using an inverse gamma prior.\n'])
end