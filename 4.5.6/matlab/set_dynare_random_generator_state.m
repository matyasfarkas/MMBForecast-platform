function [state_u,state_n] = set_dynare_random_generator_state(state_u,state_n)
% Wet state of Matlab/Octave random generator depending on matlab
% (octave) version.
% In older versions, Matlab kept one generator for uniformly distributed numbers and
% one for normally distributed numbers.
% For backward compatibility, we return two vectors, but, in recent
% versions of Matlab and in Octave, we return two identical vectors.
%
% Copyright (C) 2010-2017 Dynare Team
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

matlab_random_streams = ~(isoctave || matlab_ver_less_than('7.7'));

if matlab_random_streams% Use new matlab interface.
    if matlab_ver_less_than('7.12')
        s = RandStream.getDefaultStream();
    else
        s = RandStream.getGlobalStream();
    end
    if isequal(s.Type,'legacy')
        rand('state',state_u);
        randn('state',state_n);
    else
        if ~isequal(state_u,state_n)
            error(['You are using the new Matlab RandStream mechanism ' ...
                   'with a single random generator, but the values ' ...
                   'of the state of the uniformly ' ...
                   'distributed numbers and of the state of the ' ...
                   'normally distributed numbers are different. Something must be ' ...
                   'wrong, such as reloading old Metropolis runs, ' ...
                   'computed on a different version of Matlab. If you ' ...
                   'don''t understand the origin of the problem, ' ...
                   'please, contact Dynare''s development team.'])
        end
        s.State = state_u;
        if matlab_ver_less_than('7.12')
            RandStream.setDefaultStream(s);
        else
            RandStream.setGlobalStream(s);
        end
    end
else% Use old matlab interface.
    rand('state',state_u);
    randn('state',state_n);
end