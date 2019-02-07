function [eigenvalues_,result,info] = check(M, options, oo)
% Checks determinacy conditions by computing the generalized eigenvalues.

%@info:
%! @deftypefn {Function File} {[result,info] =} check (@var{M},@var{options},@var{oo})
%! @anchor{check}
%! @sp 1
%! Checks determinacy conditions by computing the generalized eigenvalues.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item M
%! Matlab's structure describing the model (initialized by dynare).
%! @item options
%! Matlab's structure describing the options (initialized by dynare).
%! @item oo
%! Matlab's structure gathering the results (initialized by dynare).
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item eigenvalues_
%! Eigenvalues of the model.
%! @item result
%! Integer scalar equal to one (BK conditions are satisfied) or zero (otherwise).
%! @item info
%! Integer scalar, error code as returned by @ref{resol}.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{smm_objective}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @ref{resol}
%! None.
%! @end deftypefn
%@eod:

% Copyright (C) 2001-2014 Dynare Team
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


if ~options.initval_file && M.exo_nbr > 1
    oo.exo_simul = ones(M.maximum_lead+M.maximum_lag+1,1)*oo.exo_steady_state';
end

options.order = 1;

if isempty(options.qz_criterium)
    options.qz_criterium = 1+1e-6;
end

oo.dr=set_state_space(oo.dr,M,options);

[dr,info,M,options,oo] = resol(1,M,options,oo);

if info(1) ~= 0 && info(1) ~= 3 && info(1) ~= 4
    print_info(info, 0, options);
end

eigenvalues_ = dr.eigval;
[m_lambda,i]=sort(abs(eigenvalues_));
n_explod = nnz(abs(eigenvalues_) > options.qz_criterium);

% Count number of forward looking variables
if ~options.block
    nyf = M.nsfwrd;
else
    nyf = 0;
    for j = 1:length(M.block_structure.block)
        nyf = nyf + M.block_structure.block(j).n_forward + M.block_structure.block(j).n_mixed;
    end
end

result = 0;
if (nyf == n_explod) && (dr.full_rank)
    result = 1;
end

if options.noprint == 0
    skipline()
    disp('EIGENVALUES:')
    disp(sprintf('%16s %16s %16s\n','Modulus','Real','Imaginary'))
    z=[m_lambda real(eigenvalues_(i)) imag(eigenvalues_(i))]';
    disp(sprintf('%16.4g %16.4g %16.4g\n',z))
    disp(sprintf('\nThere are %d eigenvalue(s) larger than 1 in modulus ', n_explod));
    disp(sprintf('for %d forward-looking variable(s)',nyf));
    skipline()
    if result
        disp('The rank condition is verified.')
    else
        disp('The rank condition ISN''T verified!')
    end
    skipline()
end
