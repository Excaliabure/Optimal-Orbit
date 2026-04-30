
%{
@file    simplexmethod.m
@author  Alexzia Gutierrez / UCMerced (email)
@date    2026-04-23
@brief   geenral simplex method (Alexzia)

@section This is an application of the simplex method to a blackbox
function.
%}

%{ 
Bust, faulty implimentation and miscommunication (my bad :3) 
Due to the highly non-linear nature and iterative process of the
objective function and the "black box" nature of it, it is necaessary to
expand upon the simplex algorithm and use something called "Nelder–Mead
method" also known as "downhill simplex" method.

Nedler-Mead -> https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method

TODO:
    Impliment Nedler-Mead method

%} 



% HWEloooo


function [x_opt, z_opt] = simplexmethod(A, b, c)
    %{
    Params
        A - 
    
    Return
        x_optim - optimal variables
        z_optim - optimal slack variables? (Needs verification)

    %}
    
    % Needs comments
    [m, n] = size(A);
    
    % Slack variables?
    A = [A eye(m)];
    c = [c; zeros(m,1)];
    
    % tableau
    tableau = [A b];
    tableau = [tableau; -c' 0];
    
    
    basic_vars = n+1:n+m;
    
    while true
        [min_val, pivot_col] = min(tableau(end,1:end-1));
        
        if min_val >= 0
            break;
        end
        
        % test for ratios
        ratios = inf(m,1);
        for i = 1:m
            if tableau(i,pivot_col) > 0
                ratios(i) = tableau(i,end) / tableau(i,pivot_col);
            end
        end
        
        [~, pivot_row] = min(ratios);
        
        if isinf(ratios(pivot_row))
            error('Unbounded solution');
        end
        
        % Pivot
        pivot = tableau(pivot_row, pivot_col);
        tableau(pivot_row,:) = tableau(pivot_row,:) / pivot;
        
        for i = 1:m+1
            if i ~= pivot_row
                tableau(i,:) = tableau(i,:) - ...
                    tableau(i,pivot_col)*tableau(pivot_row,:);
            end
        end
        
        basic_vars(pivot_row) = pivot_col;
    end
    
    
    x_opt = zeros(n,1);
    
    for i = 1:m
        if basic_vars(i) <= n
            x_opt(basic_vars(i)) = tableau(i,end);
        end
    end
    
    
    x_optim = x_opt(1:n-m);
    z_opt = tableau(end,end);
    
    % End of simplex
    
end