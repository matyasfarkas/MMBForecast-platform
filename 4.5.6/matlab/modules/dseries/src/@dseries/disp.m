function disp(A)

%@info:
%! @deftypefn {Function File} disp (@var{A})
%! @anchor{@dseries/disp}
%! @sp 1
%! Overloads the disp method for the Dynare time series class (@ref{dseries}).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item A
%! Dynare time series object instantiated by @ref{dseries}.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! None
%! @end deftypefn
%@eod:

separator = repmat(' | ', nobs(A)+1,1);
vspace = ' ';
TABLE = ' ';
for t=1:nobs(A)
    TABLE = char(TABLE, date2string(A.dates(t)));
end
for i = 1:vobs(A)
    TABLE = horzcat(TABLE,separator);
    tmp = A.name{i};
    for t=1:nobs(A)
        tmp = char(tmp,num2str(A.data(t,i)));
    end
    TABLE = horzcat(TABLE, tmp);
end
disp(vspace)
disp([inputname(1) ' is a dseries object:'])
disp(vspace);
disp(TABLE);
disp(vspace);