function [res,jac,domerr] = mcpath_function(func,z,jacflag,varargin)
domerr = 0;

if jacflag
    [res,jac] = func(z,varargin{:});
else
    res = func(z,varargin{:});
end
