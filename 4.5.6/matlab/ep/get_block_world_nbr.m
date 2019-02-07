function [block_nbr,world_nbr] = get_block_world_nbr(algo,nnodes,order,periods)

% Copyright (C) 2014 Dynare Team
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


switch algo
  case 0
    world_nbr = nnodes^order;
    block_nbr = 1+(nnodes^(order+1)-nnodes)/(nnodes-1)+(periods-order)*world_nbr;
  case 1
    world_nbr = 1+(nnodes-1)*order;
    block_nbr = (order+(nnodes-1)*(order-1)*order/2+(periods-order)* ...
                 world_nbr);
  otherwise
    error('This case is not supposed to happen')
end