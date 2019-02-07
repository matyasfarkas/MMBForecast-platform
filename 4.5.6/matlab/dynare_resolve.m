function [A,B,ys,info,Model,DynareOptions,DynareResults] = dynare_resolve(Model,DynareOptions,DynareResults,mode)
% Computes the linear approximation and the matrices A and B of the transition equation.

%@info:
%! @deftypefn {Function File} {[@var{A},@var{B},@var{ys},@var{info},@var{Model},@var{DynareOptions},@var{DynareResults}] =} resol (@var{Model},@var{DynareOptions},@var{DynareResults})
%! @anchor{dynare_resolve}
%! @sp 1
%! Computes the linear approximation and the matrices A and B of the transition equation.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item Model
%! Matlab's structure describing the model (initialized by dynare, see @ref{M_}).
%! @item DynareOptions
%! Matlab's structure describing the options (initialized by dynare, see @ref{options_}).
%! @item DynareResults
%! Matlab's structure gathering the results (initialized by dynare, see @ref{oo_}).
%! @item mode
%! Passed argument if restricted state-space is required, not passed otherwise
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item A
%! Matrix of doubles, transition matrix of the state equation.
%! @item B
%! Matrix of doubles, matrix relating the endogenous variables to the innovations in the state equation.
%! @item ys
%! Vector of doubles, steady state level of the endogenous variables in declaration order
%! @item info
%! Integer scalar, error code as given by @ref{resol}.
%! @item Model
%! Matlab's structure describing the model (initialized by dynare, see @ref{M_}).
%! @item DynareOptions
%! Matlab's structure describing the options (initialized by dynare, see @ref{options_}).
%! @item DynareResults
%! Matlab's structure gathering the results (initialized by dynare, see @ref{oo_}).
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{dsge_likelihood}, @ref{DsgeLikelihood_hh}, @ref{DsgeVarLikelihood}, @ref{dsge_posterior_kernel}, @ref{DsgeSmoother}, @ref{dynare_sensitivity}, @ref{gsa/thet2tau}, @ref{gsa/stab_map}, @ref{identification_analysis}, @ref{imcforecast}, @ref{thet2tau}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @ref{resol}, @ref{kalman_transition_matrix}
%! @end deftypefn
%@eod:

% Copyright (C) 2001-2017 Dynare Team
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

[dr,info,Model,DynareOptions,DynareResults] = resol(0,Model,DynareOptions,DynareResults);
DynareResults.dr = dr;

if info(1) > 0
    A = [];
    if nargout>1
        B = [];
        if nargout>2
            ys = [];
        end
    end
    return
end

switch nargin
  case 3
    endo_nbr = Model.endo_nbr;
    nstatic = Model.nstatic;
    nspred = Model.nspred;
    iv = (1:endo_nbr)';
    if DynareOptions.block == 0
        ic = [ nstatic+(1:nspred) endo_nbr+(1:size(DynareResults.dr.ghx,2)-nspred) ]';
    else
        ic = DynareResults.dr.restrict_columns;
    end
  case 4
    iv = DynareResults.dr.restrict_var_list;
    ic = DynareResults.dr.restrict_columns;
  otherwise
    error('dynare_resolve:: Error in the calling sequence!')
end

if nargout==1
    A = kalman_transition_matrix(DynareResults.dr,iv,ic,Model.exo_nbr);
    return
end

[A,B] = kalman_transition_matrix(DynareResults.dr,iv,ic,Model.exo_nbr);
ys = DynareResults.dr.ys;
