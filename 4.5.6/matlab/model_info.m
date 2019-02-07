function model_info(varargin)
%function model_info;

% Copyright (C) 2008-2017 Dynare Team
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

global M_;
if sum(strcmp(varargin,'static')) > 0
    static_ = 1;
else
    static_ = 0;
end
if sum(strcmp(varargin,'incidence')) > 0
    incidence = 1;
else
    incidence = 0;
end
if static_
    fprintf('                                          Informations about %s (static model)\n',M_.fname);
    block_structre_str = 'block_structure_stat';
    nb_leadlag = 1;
else
    fprintf('                                          Informations about %s (dynamic model)\n',M_.fname);
    block_structre_str = 'block_structure';
    nb_leadlag = 3;
end
if(isfield(M_,block_structre_str))
    if static_
        block_structure = M_.block_structure_stat;
    else
        block_structure = M_.block_structure;
    end
    fprintf(strcat('                                          ===================',char(ones(1,length(M_.fname))*'='),'\n\n'));
    nb_blocks=length(block_structure.block);
    fprintf('The model has %d equations and is decomposed in %d blocks as follow:\n',M_.endo_nbr,nb_blocks);
    fprintf('===============================================================================================================\n');
    fprintf('| %10s | %10s | %30s | %14s | %31s |\n','Block no','Size','Block Type','   Equation','Dependent variable');
    fprintf('|============|============|================================|================|=================================|\n');
    for i=1:nb_blocks
        size_block=length(block_structure.block(i).equation);
        if(i>1)
            fprintf('|------------|------------|--------------------------------|----------------|---------------------------------|\n');
        end
        for j=1:size_block
            if(j==1)
                fprintf('| %10d | %10d | %30s | %14d | %-6d %24s |\n',i,size_block,Sym_type(block_structure.block(i).Simulation_Type),block_structure.block(i).equation(j),block_structure.block(i).variable(j),M_.endo_names(block_structure.block(i).variable(j),:));
            else
                fprintf('| %10s | %10s | %30s | %14d | %-6d %24s |\n','','','',block_structure.block(i).equation(j),block_structure.block(i).variable(j),M_.endo_names(block_structure.block(i).variable(j),:));
            end
        end
    end
    fprintf('===============================================================================================================\n');
    fprintf('\n');
    if static_
        fprintf('%-30s %s','the variable','is used in equations Contemporaneously');
        if(size(block_structure.incidence.sparse_IM,1)>0)
            IM=sortrows(block_structure.incidence.sparse_IM,2);
        else
            IM=[];
        end
        size_IM=size(IM,1);
        last=99999999;
        for i=1:size_IM
            if(last~=IM(i,2))
                fprintf('\n%-30s',M_.endo_names(IM(i,2),:));
            end
            fprintf(' %5d',IM(i,1));
            last=IM(i,2);
        end
        fprintf('\n\n');
    else
        for k=1:M_.maximum_endo_lag+M_.maximum_endo_lead+1
            if(k==M_.maximum_endo_lag+1)
                fprintf('%-30s %s','the variable','is used in equations Contemporaneously');
            elseif(k<M_.maximum_endo_lag+1)
                fprintf('%-30s %s %d','the variable','is used in equations with lag ',M_.maximum_endo_lag+1-k);
            else
                fprintf('%-30s %s %d','the variable','is used in equations with lead ',k-(M_.maximum_endo_lag+1));
            end
            if(size(block_structure.incidence(k).sparse_IM,1)>0)
                IM=sortrows(block_structure.incidence(k).sparse_IM,2);
            else
                IM=[];
            end
            size_IM=size(IM,1);
            last=99999999;
            for i=1:size_IM
                if(last~=IM(i,2))
                    fprintf('\n%-30s',M_.endo_names(IM(i,2),:));
                end
                fprintf(' %5d',IM(i,1));
                last=IM(i,2);
            end
            fprintf('\n\n');
        end
    end

    %printing the gross incidence matrix
    IM_star = char([kron(ones(M_.endo_nbr, M_.endo_nbr-1), double(blanks(3))) double(blanks(M_.endo_nbr)')]);
    for i = 1:nb_leadlag
        n = size(block_structure.incidence(i).sparse_IM,1);
        for j = 1:n
            if ismember(block_structure.incidence(i).sparse_IM(j,2), M_.state_var)
                IM_star(block_structure.incidence(i).sparse_IM(j,1), 3 * (block_structure.incidence(i).sparse_IM(j,2) - 1) + 1) = 'X';
            else
                IM_star(block_structure.incidence(i).sparse_IM(j,1), 3 * (block_structure.incidence(i).sparse_IM(j,2) - 1) + 1) = '1';
            end
        end
    end
    seq = 1: M_.endo_nbr;
    blank = [ blanks(size(M_.endo_names,2)); blanks(size(M_.endo_names,2))];
    for i = 1:M_.endo_nbr
        if i == 1
            var_names = [blank; M_.endo_names(i,:)];
        else
            var_names = [var_names; blank; M_.endo_names(i,:)];
        end
    end
    if incidence
        topp = [char(kron(double(blanks(ceil(log10(M_.endo_nbr)))),ones(size(M_.endo_names,2),1))) var_names' ];
        bott = [int2str(seq') blanks(M_.endo_nbr)' blanks(M_.endo_nbr)' IM_star];
        fprintf('\n                                          Gross incidence matrix\n');
        fprintf('                                          =======================\n');
        disp([topp; bott]);

        %printing the reordered incidence matrix
        IM_star_reordered = char([kron(ones(M_.endo_nbr, M_.endo_nbr-1), double(blanks(3))) double(blanks(M_.endo_nbr)')]);
        eq(block_structure.equation_reordered) = seq;
        va(block_structure.variable_reordered) = seq;
        barre_blank = [ barre(size(M_.endo_names,2)); blanks(size(M_.endo_names,2))];
        cur_block = 1;
        for i = 1:M_.endo_nbr
            past_block = cur_block;
            while ismember(block_structure.variable_reordered(i), block_structure.block(cur_block).variable) == 0
                cur_block = cur_block + 1;
            end
            if i == 1
                var_names = [blank; M_.endo_names(block_structure.variable_reordered(i),:)];
            else
                if past_block ~= cur_block
                    var_names = [var_names; barre_blank; M_.endo_names(block_structure.variable_reordered(i),:)];
                else
                    var_names = [var_names; blank; M_.endo_names(block_structure.variable_reordered(i),:)];
                end
            end
        end
        topp = [char(kron(double(blanks(ceil(log10(M_.endo_nbr)))),ones(size(M_.endo_names,2),1))) var_names' ];
        n_state_var = length(M_.state_var);
        IM_state_var = zeros(n_state_var, n_state_var);
        inv_variable_reordered(block_structure.variable_reordered) = 1:M_.endo_nbr;
        state_equation = block_structure.equation_reordered(inv_variable_reordered(M_.state_var));
        for i = 1:nb_leadlag
            n = size(block_structure.incidence(i).sparse_IM,1);
            for j = 1:n
                [tf, loc] = ismember(block_structure.incidence(i).sparse_IM(j,2), M_.state_var);
                if tf
                    IM_star_reordered(eq(block_structure.incidence(i).sparse_IM(j,1)), 3 * (va(block_structure.incidence(i).sparse_IM(j,2)) - 1) + 1) = 'X';
                    [tfi, loci] = ismember(block_structure.incidence(i).sparse_IM(j,1), state_equation);
                    if tfi
                        IM_state_var(loci, loc) = 1;
                    end
                else
                    IM_star_reordered(eq(block_structure.incidence(i).sparse_IM(j,1)), 3 * (va(block_structure.incidence(i).sparse_IM(j,2)) - 1) + 1) = '1';
                end
            end
        end
        fprintf('1: non nul element, X: non nul element related to a state variable\n');

        cur_block = 1;
        i_last = 0;
        block = {};
        for i = 1:n_state_var
            past_block = cur_block;
            while ismember(M_.state_var(i), block_structure.block(cur_block).variable) == 0
                cur_block = cur_block + 1;
            end
            if (past_block ~= cur_block) || (past_block == cur_block && i == n_state_var)
                block(past_block).IM_state_var(1:(i - 1 - i_last), 1:i - 1) = IM_state_var(i_last+1:i - 1, 1:i - 1);
                i_last = i - 1;
            end
        end
        cur_block = 1;
        for i = 1:M_.endo_nbr
            past_block = cur_block;
            while ismember(block_structure.variable_reordered(i), block_structure.block(cur_block).variable) == 0
                cur_block = cur_block + 1;
            end
            if past_block ~= cur_block
                for j = 1:i-1
                    IM_star_reordered(j, 3 * (i - 1) - 1) = '|';
                end
            end
        end

        bott = [int2str(block_structure.equation_reordered') blanks(M_.endo_nbr)' blanks(M_.endo_nbr)' IM_star_reordered];
        fprintf('\n                                          Reordered incidence matrix\n');
        fprintf('                                          ==========================\n');
        disp([topp; bott]);
        fprintf('1: non nul element, X: non nul element related to a state variable\n');
    end
else
    fprintf('There is no block decomposition of the model.\nUse ''block'' model''s option.\n');
end

function ret=Sym_type(type)
UNKNOWN=0;
EVALUATE_FORWARD=1;
EVALUATE_BACKWARD=2;
SOLVE_FORWARD_SIMPLE=3;
SOLVE_BACKWARD_SIMPLE=4;
SOLVE_TWO_BOUNDARIES_SIMPLE=5;
SOLVE_FORWARD_COMPLETE=6;
SOLVE_BACKWARD_COMPLETE=7;
SOLVE_TWO_BOUNDARIES_COMPLETE=8;
EVALUATE_FORWARD_R=9;
EVALUATE_BACKWARD_R=10;
switch (type)
  case (UNKNOWN)
    ret='UNKNOWN                     ';
  case {EVALUATE_FORWARD,EVALUATE_FORWARD_R}
    ret='EVALUATE FORWARD            ';
  case {EVALUATE_BACKWARD,EVALUATE_BACKWARD_R}
    ret='EVALUATE BACKWARD            ';
  case SOLVE_FORWARD_SIMPLE
    ret='SOLVE FORWARD SIMPLE        ';
  case SOLVE_BACKWARD_SIMPLE
    ret='SOLVE BACKWARD SIMPLE        ';
  case SOLVE_TWO_BOUNDARIES_SIMPLE
    ret='SOLVE TWO BOUNDARIES SIMPLE  ';
  case SOLVE_FORWARD_COMPLETE
    ret='SOLVE FORWARD COMPLETE      ';
  case SOLVE_BACKWARD_COMPLETE
    ret='SOLVE BACKWARD COMPLETE      ';
  case SOLVE_TWO_BOUNDARIES_COMPLETE
    ret='SOLVE TWO BOUNDARIES COMPLETE';
end


function ret = barre(n)
s = [];
for i=1:n
    s = [s '|'];
end
ret = s;
