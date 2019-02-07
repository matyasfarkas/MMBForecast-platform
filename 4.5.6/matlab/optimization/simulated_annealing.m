function [xopt, fopt,exitflag, n_accepted_draws, n_total_draws, n_out_of_bounds_draws, t, vm] = ...
    simulated_annealing(fcn,x,optim,lb,ub,varargin)
% function [xopt, fopt,exitflag, n_accepted_draws, n_total_draws, n_out_of_bounds_draws, t, vm] = ...
%     simulated_annealing(fcn,x,optim,lb,ub,varargin)
%
% Implements the continuous simulated annealing global optimization
% algorithm described in Corana et al. (1987)
%
%  A very quick (perhaps too quick) overview of SA:
%     SA tries to find the global optimum of an N dimensional function.
%  It moves both up and downhill and as the optimization process
%  proceeds, it focuses on the most promising area.
%     To start, it randomly chooses a trial point within the step length
%  VM (a vector of length N) of the user selected starting point. The
%  function is evaluated at this trial point and its value is compared
%  to its value at the initial point.
%     In a maximization problem, all uphill moves are accepted and the
%  algorithm continues from that trial point. Downhill moves may be
%  accepted  the decision is made by the Metropolis criteria. It uses T
%  (temperature) and the size of the downhill move in a probabilistic
%  manner. The smaller T and the size of the downhill move are, the more
%  likely that move will be accepted. If the trial is accepted, the
%  algorithm moves on from that point. If it is rejected, another point
%  is chosen instead for a trial evaluation.
%     Each element of VM periodically adjusted so that half of all
%  function evaluations in that direction are accepted.
%     A fall in T is imposed upon the system with the RT option by
%  T(i+1) = RT*T(i) where i is the ith iteration. Thus, as T declines,
%  downhill moves are less likely to be accepted and the percentage of
%  rejections rise. Given the scheme for the selection for VM, VM falls.
%  Thus, as T declines, VM falls and SA focuses upon the most promising
%  area for optimization.
%
%  The importance of the parameter T (initial_temperature):
%     The parameter T is crucial in using SA successfully. It influences
%  VM, the step length over which the algorithm searches for optima. For
%  a small intial T, the step length may be too small  thus not enough
%  of the function might be evaluated to find the global optima. The user
%  should carefully examine VM in the intermediate output (set verbosity =
%  1) to make sure that VM is appropriate. The relationship between the
%  initial temperature and the resulting step length is function
%  dependent.
%     To determine the starting temperature that is consistent with
%  optimizing a function, it is worthwhile to run a trial run first. Set
%  RT = 1.5 and T = 1.0. With RT > 1.0, the temperature increases and VM
%  rises as well. Then select the T that produces a large enough VM.
%
%  Input Parameters:
%    Note: The suggested values generally come from Corana et al. To
%          drastically reduce runtime, see Goffe et al., pp. 90-1 for
%          suggestions on choosing the appropriate RT and NT.
%
%    fcn - function to be optimized.
%    x - The starting values for the variables of the function to be
%        optimized. (N)
%    optim:     Options structure with fields
%
%       optim.maximizer_indicator - Denotes whether the function should be maximized or
%           minimized. A value =1 denotes maximization while a
%           value =0 denotes minimization. Intermediate output (see verbosity)
%           takes this into account.
%       optim.RT - The temperature reduction factor
%       optim.TolFun - Error tolerance for termination. If the final function
%           values from the last neps temperatures differ from the
%           corresponding value at the current temperature by less than
%           optim.TolFun and the final function value at the current temperature
%           differs from the current optimal function value by less than
%           optim.TolFun, execution terminates and exitflag = 0 is returned.
%       optim.ns - Number of cycles. After NS*N function evaluations, each
%           element of VM is adjusted so that approximately half of
%           all function evaluations are accepted. The suggested value
%           is 20.
%       optim.nt - Number of iterations before temperature reduction. After
%           NT*NS*N function evaluations, temperature (T) is changed
%           by the factor optim.RT. Value suggested by Corana et al. is
%           max(100, 5*n). See Goffe et al. for further advice.
%       optim.neps - Number of final function values used to decide upon termi-
%           nation. See optim.TolFun. Suggested value is 4.
%       optim.MaxIter - Maximum number of function evaluations. If it is
%           exceeded, exitflag = 1.
%       optim.step_length_c - Vector that controls the step length adjustment. The suggested
%           value for all elements is 2.0.
%       optim.verbosity - controls printing inside SA.
%           Values: 0 - Nothing printed.
%                     1 - Function value for the starting value and summary results before each temperature
%                         reduction. This includes the optimal function value found so far, the total
%                         number of moves (broken up into uphill, downhill, accepted and rejected), the
%                         number of out of bounds trials, the number of new optima found at this
%                         temperature, the current optimal X and the step length VM. Note that there are
%                         N*NS*NT function evalutations before each temperature reduction. Finally, notice is
%                         is also given upon achieveing the termination criteria.
%                     2 - Each new step length (VM), the current optimal X (XOPT) and the current trial X (X). This
%                         gives the user some idea about how far X strays from XOPT as well as how VM is adapting
%                         to the function.
%                     3 - Each function evaluation, its acceptance or rejection and new optima. For many problems,
%                         this option will likely require a small tree if hard copy is used. This option is best
%                         used to learn about the algorithm. A small value for optim.MaxIter is thus recommended when
%                         using optim.verbosity = 3.
%    optim.initial_temperature  initial temperature. See Goffe et al. for advice.
%    optim.initial_step_length  (VM) step length vector. On input it should encompass the
%                               region of interest given the starting value X. For point
%                               X(I), the next trial point is selected is from X(I) - VM(I)
%                               to  X(I) + VM(I). Since VM is adjusted so that about half
%                               of all points are accepted, the input value is not very
%                               important (i.e. is the value is off, SA adjusts VM to the
%                               correct value)
%
%    lb - The lower bound for the allowable solution variables.
%    ub - The upper bound for the allowable solution variables.
%         If the algorithm chooses X(I) < LB(I) or X(I) > UB(I),
%         I = 1, N, a point is from inside is randomly selected.
%         This focuses the algorithm on the region inside UB and LB.
%         Unless the user wishes to concentrate the search to a par-
%         ticular region, UB and LB should be set to very large positive
%         and negative values, respectively. Note that the starting
%         vector X should be inside this region. Also note that LB and
%         UB are fixed in position, while VM is centered on the last
%         accepted trial set of variables that optimizes the function.
%
%
% Input/Output Parameters:
%
%  Output Parameters:
%    xopt - The variables that optimize the function. (N)
%    fopt - The optimal value of the function.
%    exitflag - The error return number.
%          Values: 0 - Normal return  termination criteria achieved.
%                  1 - Number of function evaluations (NFCNEV) is
%                      greater than the maximum number (optim.MaxIter).
%                  2 - The starting value (X) is not inside the
%                      bounds (LB and UB).
%                  3 - The initial temperature is not positive.
%                  99 - Should not be seen  only used internally.
%    n_accepted_draws - The number of accepted function evaluations.
%    n_total_draws - The total number of function evaluations. In a minor
%             point, note that the first evaluation is not used in the
%             core of the algorithm  it simply initializes the
%             algorithm.
%    n_out_of_bounds_draws - The total number of trial function evaluations that
%            would have been out of bounds of LB and UB. Note that
%            a trial point is randomly selected between LB and UB.
%    t:     On output, the final temperature.
%    vm:    Final step length vector
%
% Algorithm:
%  This routine implements the continuous simulated annealing global
%  optimization algorithm described in Corana et al.'s article
%  "Minimizing Multimodal Functions of Continuous Variables with the
%  "Simulated Annealing" Algorithm" in the September 1987 (vol. 13,
%  no. 3, pp. 262-280) issue of the ACM Transactions on Mathematical
%  Software.
%
% For modifications to the algorithm and many details on its use,
%  (particularly for econometric applications) see Goffe, Ferrier
%  and Rogers, "Global Optimization of Statistical Functions with
%  Simulated Annealing," Journal of Econometrics, vol. 60, no. 1/2,
%  Jan./Feb. 1994, pp. 65-100.
%
%  Based on the Matlab code written by Thomas Werner (Bundesbank December
%  2002), which in turn is based on the GAUSS version of Bill Goffe's simulated annealing
%  program for global optimization, written by E.G.Tsionas (9/4/95).
%
% Copyright (C) 1995 E.G.Tsionas
% Copyright (C) 1995-2002 Thomas Werner
% Copyright (C) 2002-2015 Giovanni Lombardo
% Copyright (C) 2015-2017 Dynare Team
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

c=optim.step_length_c;
t=optim.initial_temperature;
vm=optim.initial_step_length;
n=size(x,1);
xp=zeros(n,1);
%*  Set initial values.*
n_accepted_draws=0;
n_out_of_bounds_draws=0;
n_total_draws=0;
exitflag=99;
xopt=x;
nacp=zeros(n,1);
fstar=1e20*ones(optim.neps,1);
%* If the initial temperature is not positive, notify the user and abort. *
if(t<=0.0)
    fprintf('\nThe initial temperature is not positive. Reset the variable t\n');
    exitflag=3;
    return
end
%*  If the initial value is out of bounds, notify the user and abort. *
if(sum(x>ub)+sum(x<lb)>0)
    fprintf('\nInitial condition out of bounds\n');
    exitflag=2;
    return
end
%*  Evaluate the function with input x and return value as f. *
f=feval(fcn,x,varargin{:});
%*
%  If the function is to be minimized, switch the sign of the function.
%  Note that all intermediate and final output switches the sign back
%  to eliminate any possible confusion for the user.
%*
if(optim.maximizer_indicator==0)
    f=-f;
end
n_total_draws=n_total_draws+1;
fopt=f;
fstar(1)=f;
if(optim.verbosity >1)
    disp '  ';
    disp(['initial x    ' num2str(x(:)')]);
    if(optim.maximizer_indicator)
        disp(['initial f    ' num2str(f)]);
    else
        disp(['initial f    ' num2str(-f)]);
    end
end
%  Start the main loop. Note that it terminates if (i) the algorithm
%  succesfully optimizes the function or (ii) there are too many
%  function evaluations (more than optim.MaxIter).

while (1>0)
    nup=0;
    nrej=0;
    nnew=0;
    ndown=0;
    lnobds=0;
    m=1;
    while m<=optim.nt
        j=1;
        while j<=optim.ns
            h=1;
            while h<=n
                %*  Generate xp, the trial value of x. Note use of vm to choose xp. *
                i=1;
                while i<=n
                    if(i==h)
                        xp(i)=x(i)+(rand(1,1)*2.0-1.0)*vm(i);
                    else
                        xp(i)=x(i);
                    end
                    %*  If xp is out of bounds, select a point in bounds for the trial. *
                    if((xp(i)<lb(i) || xp(i)>ub(i)))
                        xp(i)=lb(i)+(ub(i)-lb(i))*rand(1,1);
                        lnobds=lnobds+1;
                        n_out_of_bounds_draws=n_out_of_bounds_draws+1;
                        if(optim.verbosity >=3)
                            if exist('fp','var')
                                print_current_invalid_try(optim.maximizer_indicator,xp,x,fp,f);
                            end
                        end
                    end
                    i=i+1;
                end
                %*  Evaluate the function with the trial point xp and return as fp. *
                % fp=feval(fcn,xp,listarg);
                fp=feval(fcn,xp,varargin{:});
                if(optim.maximizer_indicator==0)
                    fp=-fp;
                end
                n_total_draws=n_total_draws+1;
                if(optim.verbosity >=3)
                    print_current_valid_try(optim.maximizer_indicator,xp,x,fp,f);
                end
                %*  If too many function evaluations occur, terminate the algorithm. *
                if(n_total_draws>=optim.MaxIter)
                    fprintf('Too many function evaluations; consider\n');
                    fprintf('increasing optim.MaxIter or optim.TolFun or decreasing\n');
                    fprintf('optim.nt or optim.rt. These results are likely to be poor\n');
                    if(optim.maximizer_indicator==0)
                        fopt=-fopt;
                    end
                    exitflag=1;
                    return
                end
                %*  Accept the new point if the function value increases. *
                if(fp>=f)
                    if(optim.verbosity >=3)
                        fprintf('point accepted\n');
                    end
                    x=xp;
                    f=fp;
                    n_accepted_draws=n_accepted_draws+1;
                    nacp(h)=nacp(h)+1;
                    nup=nup+1;
                    %*  If greater than any other point, record as new optimum. *
                    if(fp>fopt)
                        if(optim.verbosity >=3)
                            fprintf('new optimum\n');
                        end
                        xopt=xp;
                        fopt=fp;
                        nnew=nnew+1;
                    end
                    %*
                    % If the point is lower, use the Metropolis criteria to decide on
                    % acceptance or rejection.
                    %*
                else
                    p=exp((fp-f)/t);
                    pp=rand(1,1);
                    if(pp<p)
                        if(optim.verbosity >=3)
                            if(optim.maximizer_indicator)
                                fprintf('though lower, point accepted\n');
                            else
                                fprintf('though higher, point accepted\n');
                            end
                        end
                        x=xp;
                        f=fp;
                        n_accepted_draws=n_accepted_draws+1;
                        nacp(h)=nacp(h)+1;
                        ndown=ndown+1;
                    else
                        nrej=nrej+1;
                        if(optim.verbosity >=3)
                            if(optim.maximizer_indicator)
                                fprintf('lower point rejected\n');
                            else
                                fprintf('higher point rejected\n');
                            end
                        end
                    end
                end
                h=h+1;
            end
            j=j+1;
        end
        %*  Adjust vm so that approximately half of all evaluations are accepted. *
        i=1;
        while i<=n
            ratio=nacp(i)/optim.ns;
            if(ratio>.6)
                vm(i)=vm(i)*(1.+c(i)*(ratio-.6)/.4);
            elseif(ratio<.4)
                vm(i)=vm(i)/(1.+c(i)*((.4-ratio)/.4));
            end
            if(vm(i)>(ub(i)-lb(i)))
                vm(i)=ub(i)-lb(i);
            end
            i=i+1;
        end
        if(optim.verbosity >=2)
            fprintf('intermediate results after step length adjustment\n');
            fprintf('new step length(vm)  %4.3f', vm(:)');
            fprintf('current optimal x    %4.3f', xopt(:)');
            fprintf('current x            %4.3f', x(:)');
        end
        nacp=zeros(n,1);
        m=m+1;
    end
    if(optim.verbosity >=1)
        print_intermediate_statistics(optim.maximizer_indicator,t,xopt,vm,fopt,nup,ndown,nrej,lnobds,nnew);
    end
    %*  Check termination criteria. *
    quit=0;
    fstar(1)=f;
    if((fopt-fstar(1))<=optim.TolFun)
        quit=1;
    end
    if(sum(abs(f-fstar)>optim.TolFun)>0)
        quit=0;
    end
    %*  Terminate SA if appropriate. *
    if(quit)
        exitflag=0;
        if(optim.maximizer_indicator==0)
            fopt=-fopt;
        end
        if(optim.verbosity >=1)
            fprintf('SA achieved termination criteria.exitflag=0\n');
        end
        return
    end
    %*  If termination criteria are not met, prepare for another loop. *
    t=optim.rt*t;
    i=optim.neps;
    while i>=2
        fstar(i)=fstar(i-1);
        i=i-1;
    end
    f=fopt;
    x=xopt;
    %*  Loop again. *
end

end

function  print_current_invalid_try(max,xp,x,fp,f)
fprintf('\n');
disp(['Current x    ' num2str(x(:)')]);
if(max)
    disp(['Current f    ' num2str(f)]);
else
    disp(['Current f    ' num2str(-f)]);
end
disp(['Trial x      ' num2str(xp(:)')]);
disp 'Point rejected since out of bounds';
end

function print_current_valid_try(max,xp,x,fp,f)

disp(['Current x    ' num2str(x(:)')]);
if(max)
    disp(['Current f   ' num2str(f)]);
    disp(['Trial x     ' num2str(xp(:)')]);
    disp(['Resulting f ' num2str(fp)]);
else
    disp(['Current f   ' num2str(-f)]);
    disp(['Trial x     ' num2str(xp(:)')]);
    disp(['Resulting f ' num2str(-fp)]);
end
end


function  print_intermediate_statistics(max,t,xopt,vm,fopt,nup,ndown,nrej,lnobds,nnew)

totmov=nup+ndown+nrej;
fprintf('\nIntermediate results before next temperature reduction\n');
disp(['current temperature         ' num2str(t)]);

if(max)
    disp(['Max function value so far   ' num2str(fopt)]);
    disp(['Total moves                 ' num2str(totmov)]);
    disp(['Uphill                      ' num2str(nup)]);
    disp(['Accepted downhill           ' num2str(ndown)]);
    disp(['Rejected downhill           ' num2str(nrej)]);
    disp(['Out of bounds trials        ' num2str(lnobds)]);
    disp(['New maxima this temperature ' num2str(nnew)]);
else
    disp(['Min function value so far   ' num2str(-fopt)]);
    disp(['Total moves                 ' num2str(totmov)]);
    disp(['Downhill                    ' num2str(nup)]);
    disp(['Accepted uphill             ' num2str(ndown)]);
    disp(['Rejected uphill             ' num2str(nrej)]);
    disp(['Trials out of bounds        ' num2str(lnobds)]);
    disp(['New minima this temperature ' num2str(nnew)]);
end
xopt1=xopt(1:round(length(xopt)/2));
xopt2=xopt(round(length(xopt)/2)+1:end);

disp(['Current optimal x1           ' num2str(xopt1')]);
disp(['Current optimal x2           ' num2str(xopt2')]);
disp(['Strength(vm)                ' num2str(vm')]);
fprintf('\n');
end